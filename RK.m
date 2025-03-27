function error = gauss_legendre_rk(Y0, A_sym, y0, T, h, epsilon)
    syms y;
    A_func = matlabFunction(A_sym, 'Vars', {y}); % Ensure A_func returns a matrix
    
    m = length(Y0); % Dimension of the system
    N = round(T/h); % Number of time steps
    Y = zeros(m, N+1); % Store solution
    Y(:,1) = Y0;
    error = zeros(1, N+1); % Store first component for plotting
    error(1) = subs(y0, 0) - Y0(1);
    
    % Butcher Tableau values for Gauss-Legendre 2-stage
    c1 = (3 - sqrt(3)) / 6;
    c2 = (3 + sqrt(3)) / 6;
    a11 = 1/4; a12 = (3 - 2*sqrt(3)) / 12;
    a21 = (3 + 2*sqrt(3)) / 12; a22 = 1/4;
    b1 = 1/2; b2 = 1/2;
    
    for n = 1:N
        Yn = Y(:,n);
        
        % Define nonlinear system for K1, K2
        F = @(K) [K(1:m,1) - A_func(Yn(1) + h*(a11*K(1:m,1) + a12*K(1:m,2))) * (Yn + h*(a11*K(1:m,1) + a12*K(1:m,2)));
                  K(1:m,2) - A_func(Yn(1) + h*(a21*K(1:m,1) + a22*K(1:m,2))) * (Yn + h*(a21*K(1:m,1) + a22*K(1:m,2)))];
        
        % Initial guess for K (set to zero matrix)
        K_init = zeros(m, 2);
        
        % Solve for K1 and K2
        K_sol = fsolve(@(K) reshape(F(reshape(K, m, 2)), [], 1), reshape(K_init, [], 1));
        K_sol = reshape(K_sol, m, 2);
        
        % Compute Y_{n+1}
        Y(:,n+1) = Y(:,n) + h * (b1 * K_sol(:,1) + b2 * K_sol(:,2));
        
        % Store first component
        error(n+1) = subs(y0, n*h) - Y(1, n+1);
    end
    
    % Plot first component
    t_vals = 0:h:T;
    plot(t_vals, error, 'r');
    grid;
    xlabel('t'); 
    ylabel('Error');
    title(['Error of 2-stage Gauss-Legendre RK method with h=', num2str(h)]);
end

% Define parameters for Duffing equation
epsilon = 0.0001; % Nonlinearity parameter
syms y;
A_sym = [0, 1; -(1 + epsilon*y^2), 0];

% Initial conditions
Y0 = [1; 0];
T = 100;
h = 1/16;
y0 = cos(0); % Initial true solution approximation

% Run solver
error = gauss_legendre_rk(Y0, A_sym, y0, T, h, epsilon);