from manim import *
import numpy as np


class MatrixODEScheme(MovingCameraScene):
    def construct(self):
        # Setup axes
        axes = Axes(
            x_range=[0, 2 * PI, PI / 2],
            y_range=[-1.5, 1.5, 0.5],
            axis_config={"color": BLUE},
        )
        axes_labels = axes.get_axis_labels(x_label="t", y_label="Y(t)")
        self.play(Create(axes), Write(axes_labels))
        #self.wait(1)

        h = (2 * PI / 20)
        # Discretize time interval
        time_points = [n * (2 * PI / 20) for n in range(21)]
        labels = [MathTex(f"t_{{{n}}}").next_to(axes.c2p(t, 0), DOWN) for n, t in enumerate(time_points)]

        # Add notches and labels
        for label, t in zip(labels, time_points):
            notch = Line(axes.c2p(t, -0.1), axes.c2p(t, 0.1), color=BLUE)
            self.play(Create(notch), Write(label), run_time=0.1)  # Faster appearance
            self.wait(0.1)

        # Plot exact solution
        exact_solution = axes.plot(lambda t: np.cos(t), color=BLUE, stroke_opacity=0.5, stroke_width=2)
        exact_solution.set_dash([0.1, 0.1])
        self.play(Create(exact_solution))
        #self.wait(1)

        # Initial point y_0
        y_points = [axes.c2p(t, np.cos(t)) for t in time_points]
        y_labels = [MathTex(f"y_{{{n}}}").next_to(y_points[n], UP) for n in range(21)]

        y0_dot = Dot(point=y_points[0], color=RED)
        y0_label = y_labels[0]
        self.play(Create(y0_dot), Write(y0_label))
        #self.wait(1)

        # Group axes and graph for zooming
        #graph_group = VGroup(axes, exact_solution, y0_dot, y0_label, *[notch for notch in notches], *labels)

        self.camera.frame.save_state()
        # Zoom in to t_0 to t_1
        self.play(self.camera.frame.animate.set(width=7*h).move_to(axes.coords_to_point(time_points[0] + h/2, np.cos(time_points[0]))))
        #self.wait(1)

        # Display formula
        formula1 = MathTex(r"y_1 = \exp(\Omega(t_1, t_0)) y_0").move_to(axes.coords_to_point(time_points[0] + h/2, np.cos(time_points[0]) - h/2))
        formula1.scale(0.2)
        self.play(Write(formula1))
        #self.wait(1)

        # Pointer and y_1
        pointer = Arrow(start=y_points[0], end=y_points[1], color=RED, buff=0.1, stroke_opacity=0.5)
        self.play(Create(pointer))
        #self.wait(0.5)
        y1_dot = Dot(point=y_points[1], color=RED)
        y1_label = y_labels[1]
        self.play(Create(y1_dot), Write(y1_label))
        #self.wait(1)

        # Fade out formula
        self.play(FadeOut(formula1))
        #self.wait(0.5)

        # Zoom out
        #self.play(Restore(self.camera.frame))
        #self.wait(1)

        # Repeat for t_1 to t_2 and t_2 to t_3
        for n in range(2, 4):
            # Zoom in
            self.play(self.camera.frame.animate.set(width=7 * h).move_to(axes.coords_to_point(time_points[n-1] + h/2, np.cos(time_points[n-1]))))
            #self.wait(1)

            # Display formula
            formula2 = MathTex(f"y_{{{n}}} = \exp(\Omega(t_{{{n}}}, t_{{{n-1}}})) y_{{{n-1}}}").move_to(axes.coords_to_point(time_points[n-1] - h/5, np.cos(time_points[n-1]) - 3*h/4))
            formula2.scale(0.2)
            self.play(Write(formula2))
            #self.wait(1)

            # Pointer and y_n
            pointer = Arrow(start=y_points[n-1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer))
            self.wait(0.2)
            y_dot = Dot(point=y_points[n], color=RED)
            y_label = y_labels[n]
            self.play(Create(y_dot), Write(y_label), run_time=0.4)
            #self.wait(0.5)

            # Fade out formula
            self.play(FadeOut(formula2))
            #self.wait(0.3)

        # Zoom out
        self.play(Restore(self.camera.frame))
        #self.wait(1)

        # Display formula in the corner
        formula = MathTex(r"Y(t_N) = \prod_{n=1}^N \exp(\Omega(t_n, t_{n-1})) Y_0").to_corner(UP + RIGHT)
        self.play(Write(formula.to_corner(UP + RIGHT)))
        #self.wait(1)

        # Animate remaining points faster
        for n in range(4, 6):
            pointer = Arrow(start=y_points[n-1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer), run_time=0.3)
            y_dot = Dot(point=y_points[n], color=RED)
            y_label = y_labels[n]
            self.play(Create(y_dot), Write(y_label), run_time=0.2)

        for n in range(6, 8):
            pointer = Arrow(start=y_points[n - 1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer), run_time=0.3)
            y_dot = Dot(point=y_points[n], color=RED)
            self.play(Create(y_dot), run_time=0.2)

        for n in range(8, 13):
            pointer = Arrow(start=y_points[n-1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer), run_time=0.3)
            y_dot = Dot(point=y_points[n], color=RED)
            y_label = y_labels[n]
            self.play(Create(y_dot), Write(y_label), run_time=0.2)

        for n in range(13, 15):
            pointer = Arrow(start=y_points[n - 1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer), run_time=0.3)
            y_dot = Dot(point=y_points[n], color=RED)
            self.play(Create(y_dot), run_time=0.2)

        for n in range(15, 21):
            pointer = Arrow(start=y_points[n-1], end=y_points[n], color=RED, buff=0.1, stroke_opacity=0.5)
            self.play(Create(pointer), run_time=0.3)
            y_dot = Dot(point=y_points[n], color=RED)
            y_label = y_labels[n]
            self.play(Create(y_dot), Write(y_label), run_time=0.2)

        self.wait(2)