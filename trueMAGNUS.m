function [nodes, weights] = NC(k, h)
    % k: order
    % h: length of the interval
    % output: nodes and weights of Newton-Cotes quadrature

    % Equally spaced nodes
    nodes = sym(linspace(0, h, k));

    if k == 2
        % Trapezoidal rule
        weights = [sym(0.5), sym(0.5)] * h;
    elseif k == 3
        % Simpson's rule
        weights = [sym(1), sym(4), sym(1)] * h / 6;
    elseif k == 4
        % Simpson's 3/8 rule
        weights = [sym(1), sym(3), sym(3), sym(1)] * h / 8;
    elseif k == 5
        % Boole's rule
        weights = [sym(7), sym(32), sym(12), sym(32), sym(7)] * h /90;
    else
        error('Higher-order Newton-Cotes weights not implemented.');
    end
end

function [nodes, weights] = GL(k, h)
    % k: order
    % h: length of the interval
    % output: nodes and weights of Gauss-Legendre quadrature
    % Use symbolic arithmetic for high precision
    syms x;
    % Legendre polynomial of order k
    Pk = legendreP(k, x);
    
    % Compute nodes (roots of Legendre polynomial)
    nodes = solve(Pk, x); % Preserve symbolic form
    nodes = sort(nodes); % Ensure correct ordering
    
    % Compute weights
    weights = sym(zeros(k, 1));
    for i = 1:k
        % Derivative of Legendre polynomial
        dPk = diff(Pk, x);
        weights(i) = simplify((2 / ((1 - nodes(i)^2) * (subs(dPk, nodes(i))^2))) * h / 2);
    end
    
    % Transform nodes to the interval [0, h]
    nodes = simplify(h / 2 * nodes + h / 2);
end

function [weights, nodes, RQ] = alpha(s, k, h, quad)
    % s: order of the method
    % k: order of the quadrature rule
    % h: length of the interval
    % quad: GL for Gauss-Legendre or NC for Newton-Cotes)
    % 
    % output: matrix RQ with coefficients to express alpha's through A(t_i)

    h = sym(h);

    % Compute weights and nodes based on the quadrature type
    if strcmp(quad, 'GL')
        [nodes, weights] = GL(k, h);
    elseif strcmp(quad, 'NC')
        [nodes, weights] = NC(k, h);
    else
        error('Unsupported quadrature type. Choose "GL" or "NC".');
    end

    % Transition matrix Q for [0, h]
    Q = sym(zeros(s, k));
    for i = 1:s
        for j = 1:k
            Q(i, j) = weights(j) * (nodes(j) - sym(0.5))^(i-1);
        end
    end

    % matrix T with symbolic arithmetic
    T = sym(zeros(s, s));
    for i = 1:s
        for j = 1:s
            T(i, j) = (1 - (-1)^(i-1+j)) / ((i-1 + j) * 2^(i-1 + j));
        end
    end

    % Compute R as the inverse of T using symbolic arithmetic
    %R = inv(T);
    RQ = T\Q;
end

function C = Commutator(A, B)
    C = A*B - B*A;
end

function error = Magnus2err(Y, A, y0, T, h, quad)
    %draws a plot of the error between the exact solution y0 and an
    %approximate soloution obtained through exponential midpoint rule and MG2GL
    tau = 0:h:T;
    n = T/h+1;
    error = zeros(n, 1);
    syms t;
    if strcmp(quad, 'NC')
        for i = 1:n
            error(i) = subs(y0, t, tau(i)) - Y(1, 1);
            Y = expm(subs(A, t, (tau(i) + h/2)/2))*Y;
        end
    elseif strcmp(quad, 'GL')
        A1 = subs(A, t, tau(1));
        A2 = subs(A, t, tau(2));
        error(1) = subs(y0, t, tau(1)) - Y(1, 1);
        Y = vpa(expm((A1 + A2)/2)*Y);
        for i = 2:(n-1)
            A1 = A2;
            A2 = subs(A, tau(i+1));
            error(i) = vpa(subs(y0, t, tau(i)) - Y(1, 1));
            Y = vpa(expm((A1 + A2)/2)*Y);
        end
    else
        error('Unsupported quadrature type. Choose "GL" or "NC".');
    end
    plot(tau, error,'b');
    grid;
    xlabel('t');
    ylabel('y');
end

function error = Magnus4err(Y, A, ep, y0, T, h, k, quad)
    %draws a plot of the error between the exact solution y0 and an
    %approximate soloution obtained through Magnus series method 
    %of the 4th order with stepsize h
    tau = 0:h:T;
    n = T/h+1;
    [b0, c0, RQ0] = alpha(2, k, 1, quad);
    error = sym(zeros(n, 1));
    syms t;
    for i = 1:n
        error(i) = subs(y0, t, tau(i)) - vpa(Y(1, 1)); 
        alpha1 = sym(zeros(size(A)));
        alpha2 = sym(zeros(size(A)));
        for j = 1:k
            Aa_j = [0 1;-(1+ep*Y(1, 1)^2) 0];
            alpha1 = alpha1 + RQ0(1, j) * Aa_j;
            alpha2 = simplify(alpha2 + (RQ0(2, j) )* Aa_j); %+ b0(j)*tau(j)
        end
        Omega = simplify(h*alpha1 - h^2*(alpha1 * alpha2 - alpha2 * alpha1) / 12);
        Y = vpa(expm(Omega)*Y);
    end
    plot(tau, error,'b');
    syms p;
    p = sym(h);
    title(strcat('Error of the 4th order Magnus method with h=', char(p)));
    grid;
    xlabel('t');
    ylabel('y');
end


function error = Magnus6err(Y, A, y0, T, h, k, quad)
    %draws a plot of the error between the exact solution y0 and an
    %approximate soloution obtained through Magnus series method 
    %of the 6th order with stepsize h
    tau = 0:h:T;
    n = T/h+1;
    [b0, c0, RQ0] = alpha(3, k, 1, quad);
    error = sym(zeros(n, 1));
    syms t;
    for i = 1:n
        error(i) = subs(y0, t, tau(i)) - vpa(Y(1, 1)); 
        alpha1 = sym(zeros(size(A)));
        alpha2 = sym(zeros(size(A)));
        alpha3 = sym(zeros(size(A)));
        for j = 1:k 
            Aa_j = subs(A, t, tau(i) + h*c0(j));
            alpha1 = alpha1 + h*RQ0(1, j)*Aa_j;
            alpha2 = alpha2 + h*RQ0(2, j)*Aa_j;
            alpha3 = alpha3 + h*RQ0(3, j)*Aa_j;
        end
        beta1 = Commutator(alpha1, alpha2);
        beta2 = - Commutator(alpha1, 2*alpha3 + beta1)/60;
        Omega = (alpha1 + alpha3/12 + Commutator(beta1 - 20*alpha1 - alpha3, alpha2 + beta2)/240);
        Y = vpa(expm(Omega)*Y);
    end
    plot(tau, error,'b');
    syms p;
    p = sym(h);
    title(strcat('Error of the 6th order Magnus method with h=', char(p)));
    grid;
    xlabel('t');
    ylabel('y');
    hAx=gca;
    hYAx = hAx.YAxis;
    hYExpText = hYAx.NodeChildren(1);
    hYExpText.FontSize = 20; 
end



syms t;
Y = [1; 0];
A = [0, 1; -t, 0];
ep = 0.00001;
y0 = cos(t);
y1 = (airy(3,0)*airy(-t)-airy(1,0)*airy(2,-t))/(airy(0)*airy(3,0)-airy(1,0)*airy(2,0));
y2 = 1/2*pi/gamma(3/4)*t^(1/2)*besselj(1/4,1/2*t^2);
y3 = cos(t)+ep*(cos(3*t)/32-cos(t)/32-(3*t*sin(t))/8);
T = 1000;
h = 1/16;
tau = 0:h:T;
s = 3;
k = 2;
quad='GL';
%[nodes, weights, RQ] = alpha(3, k, 1, quad)
%[nodes, weights] = GL(k, h);

error = Magnus4err(Y, A, ep, y3, T, h, k, quad);
