""" Poincare disk and Klein disk """
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from matplotlib.patches import Circle
from scipy.spatial.transform import Rotation
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d

""" Global variables """
x_klein = 0.
y_klein = 0.

""" Animation control """
is_play = False
is_tilt = False
is_precession = False

""" Axis vectors """
vector_x_axis = np.array([1., 0., 0.])
vector_y_axis = np.array([0., 1., 0.])
vector_z_axis = np.array([0., 0., 1.])

""" Create figure and axes """
title_ax0 = "Poincare disk and Klein disk"
title_tk = title_ax0

x_min = -4.
x_max = 4.
y_min = -4.
y_max = 4.
z_min = -4.
z_max = 4.

fig = Figure()
ax0 = fig.add_subplot(111, projection="3d")
ax0.set_box_aspect((4, 4, 4))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel("x")
ax0.set_ylabel("y")
ax0.set_zlabel("z")
ax0.set_xlim(x_min, x_max)
ax0.set_ylim(y_min, y_max)
ax0.set_zlim(z_min, z_max)

""" Embed in Tkinter """
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill="both")

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

""" Global objects of Tkinter """
var_x_klein = tk.StringVar(root)
var_y_klein = tk.StringVar(root)

""" Classes and functions """


class Counter:
    def __init__(self, is3d=None, ax=None, xy=None, z=None, label=""):
        self.is3d = is3d if is3d is not None else False
        self.ax = ax
        self.x, self.y = xy[0], xy[1]
        self.z = z if z is not None else 0
        self.label = label

        self.count = 0

        if not is3d:
            self.txt_step = self.ax.text(self.x, self.y, self.label + str(self.count))
        else:
            self.txt_step = self.ax.text2D(self.x, self.y, self.label + str(self.count))
            self.xz, self.yz, _ = proj3d.proj_transform(self.x, self.y, self.z, self.ax.get_proj())
            self.txt_step.set_position((self.xz, self.yz))

    def count_up(self):
        self.count += 1
        self.txt_step.set_text(self.label + str(self.count))

    def reset(self):
        self.count = 0
        self.txt_step.set_text(self.label + str(self.count))

    def get(self):
        return self.count


class Hyperboloid:
    def __init__(self, ax, cmap, edge_color, alpha):
        self.ax = ax
        self.cmap = cmap
        self.edge_color = edge_color
        self.alpha = alpha

        # Define the range of parameters
        u = np.linspace(-2, 2, 100)  # Controls the vertical expansion
        v = np.linspace(0, 2 * np.pi, 100)  # Circular rotation
        U, V = np.meshgrid(u, v)

        # Hyperboloid equation using sinh and cosh
        a, b = 1, 1
        X = a * np.sinh(U) * np.cos(V)  # X-coordinate
        Y = a * np.sinh(U) * np.sin(V)  # Y-coordinate
        Z = b * np.cosh(U)  # Z-coordinate

        # Create plot
        surface = self.ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor=self.edge_color, alpha=self.alpha)


class UnitSphere:
    def __init__(self, ax, cmap, edge_color, alpha):
        self.ax = ax
        self.cmap = cmap
        self.edge_color = edge_color
        self.alpha = alpha

        # Define the range of parameters
        phi = np.linspace(0, np.pi, 50)  # Elevation angle
        theta = np.linspace(0, 2 * np.pi, 50)  # Azimuth angle
        PHI, THETA = np.meshgrid(phi, theta)

        # Unit sphere
        X = np.sin(PHI) * np.cos(THETA)  # X-coordinate
        Y = np.sin(PHI) * np.sin(THETA)  # Y-coordinate
        Z = np.cos(PHI)  # T-coordinate

        # Create plot
        surface = self.ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor=self.edge_color, alpha=self.alpha)


class CircularDisk:
    def __init__(self, ax, offset, color, edge_color, alpha):
        self.ax = ax
        self.offset = offset
        self.color = color
        self.edge_color = edge_color
        self.alpha = alpha

        self .radius = 1.

        # Define the range of parameters
        disk_theta = np.linspace(0, 2 * np.pi, 50)
        disk_r = np.linspace(0, self.radius, 50)
        Disk_R, Disk_Theta = np.meshgrid(disk_r, disk_theta)

        # Disk at z = 0
        X = Disk_R * np.cos(Disk_Theta)
        Y = Disk_R * np.sin(Disk_Theta)
        Z = np.zeros_like(X) + self.offset

        # Create plot
        surface = self.ax.plot_surface(X, Y, Z, color=self.color, edgecolor="none", alpha=self.alpha)

        # Circle edge
        self.angle_space = np.arange(0, 360)

        self.x_circle = self.radius * np.cos(self.angle_space * np.pi / 180.)
        self.y_circle = self.radius * np.sin(self.angle_space * np.pi / 180.)
        self.z_circle = self.angle_space * 0.  + self.offset

        self.plt_circle, = self.ax.plot(self.x_circle, self.y_circle, self.z_circle,
                                        linewidth=1, linestyle="-", color=self.edge_color)


class LineKlein:
    def __init__(self, ax, color, alpha):
        self.ax = ax
        self.color = color
        self.alpha = alpha

        self.r = 1.
        self.x = 0.
        self.discriminant = self.r ** 2 - self.x ** 2

        self.x, self.y, self.z0 = 0., -1., 0.
        self.u, self.v, self.w0 = 0., 1., 0.
        self.line_klein_z0, = self.ax.plot([self.x, self.u], [self.y, self.v], [self.z0, self.w0],
                                           linewidth=0.5, color=self.color, linestyle="--")

        self.x, self.y, self.z1 = 0., -1., 1.
        self.u, self.v, self.w1 = 0., 1., 1.
        self.line_klein_z1, = self.ax.plot([self.x, self.u], [self.y, self.v], [self.z1, self.w1],
                                           linewidth=1, color=self.color, linestyle="-")

        self.y_marker = 0.
        self.marker, = self.ax.plot(self.x, self.y_marker, self.z1, marker="o", markersize=3, color=self.color)

    def set_x(self, value):
        if -1. <= value <= 1.:
            self.x = value
            self.u = self.x

            self.discriminant = self.r ** 2 - self.x ** 2

            if self.discriminant >= 0.:
                y_c = np.sqrt(self.discriminant)
                v_c = -np.sqrt(self.discriminant)

                self.line_klein_z0.set_data_3d([self.x, self.u], [y_c, v_c], [self.z0, self.w0])
                self.line_klein_z1.set_data_3d([self.x, self.u], [y_c, v_c], [self.z1, self.w1])

                if self.y_marker > y_c:
                    self.y_marker = y_c

                if self.y_marker < v_c:
                    self.y_marker = v_c

                self.marker.set_data_3d([self.x], [self.y_marker], [self.z1])

    def set_y(self, value):
        if -1. <= value <= 1.:
            self.y_marker = value

            self.discriminant = self.r ** 2 - self.x ** 2

            if self.discriminant >= 0.:
                y_c = np.sqrt(self.discriminant)
                v_c = -np.sqrt(self.discriminant)

                if self.y_marker > y_c:
                    self.y_marker = y_c

                if self.y_marker < v_c:
                    self.y_marker = v_c

                self.marker.set_data_3d([self.x], [self.y_marker], [self.z1])


class CircleSphere:
    def __init__(self, ax, color, alpha):
        self.ax = ax
        self.color = color
        self.alpha = alpha

        self.radius = 1.
        self.x = 0.

        self.angle_space = np.arange(0, 360)

        self.x_circle = self.angle_space * 0.
        self.y_circle = self.radius * np.cos(self.angle_space * np.pi / 180.)
        self.z_circle = self.radius * np.sin(self.angle_space * np.pi / 180.)

        self.plt_circle, = self.ax.plot(self.x_circle, self.y_circle, self.z_circle,
                                        linewidth=1, linestyle="-", color=self.color)

    def set_x(self, value):
        if -1. <= value <= 1.:
            self.x = value

            r = 1.

            discriminant = r ** 2 - self.x ** 2
            if discriminant >= 0.:
                self.radius = np.sqrt(discriminant)

                self.x_circle = self.angle_space * 0. + self.x
                self.y_circle = self.radius * np.cos(self.angle_space * np.pi / 180.)
                self.z_circle = self.radius * np.sin(self.angle_space * np.pi / 180.)

                self.plt_circle.set_data_3d(np.array(self.x_circle),
                                            np.array(self.y_circle),
                                            np.array(self.z_circle))


class DotKlein:
    def __init__(self, ax, color, alpha):
        self.ax = ax
        self.color = color
        self.alpha = alpha

        self.x, self.y, self.z1 = 0., 0., 1.
        self.marker, = self.ax.plot(self.x, self.y, self.z1, marker="o", markersize=3, color=self.color)

    def set_x(self, value):
        if -1. <= value <= 1.:
            self.x = value


def set_x(value):
    line_klein.set_x(value)
    circle_sphere.set_x(value)


def set_y(value):
    line_klein.set_y(value)
    # circle_sphere.set_y(value)


def create_parameter_setter():
    # Position of a line on Klein disk
    frm_klein = ttk.Labelframe(root, relief="ridge", text="x, y on Klein disk)", labelanchor='n')
    frm_klein.pack(side="left", fill=tk.Y)

    lbl_x_klein = tk.Label(frm_klein, text="x (line)")
    lbl_x_klein.pack(side='left')

    # var_x_klein = tk.StringVar(root)
    var_x_klein.set(str(x_klein))
    spn_x_klein = tk.Spinbox(
        frm_klein, textvariable=var_x_klein, format="%.2f", from_=-1, to=1, increment=0.02,
        command=lambda: set_x(float(var_x_klein.get())), width=5
    )
    spn_x_klein.pack(side="left")

    lbl_y_klein = tk.Label(frm_klein, text="y (dot)")
    lbl_y_klein.pack(side='left')

    # var_y_klein = tk.StringVar(root)
    var_y_klein.set(str(y_klein))
    spn_y_klein = tk.Spinbox(
        frm_klein, textvariable=var_y_klein, format="%.2f", from_=-1, to=1, increment=0.02,
        command=lambda: set_y(float(var_y_klein.get())), width=5
    )
    spn_y_klein.pack(side="left")


def create_animation_control():
    frm_anim = ttk.Labelframe(root, relief="ridge", text="Animation", labelanchor="n")
    frm_anim.pack(side="left", fill=tk.Y)
    btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
    btn_play.pack(side="left")
    btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
    btn_reset.pack(side="left")
    # btn_clear = tk.Button(frm_anim, text="Clear path", command=lambda: aaa())
    # btn_clear.pack(side="left")


def create_center_lines():
    ln_axis_x = art3d.Line3D([x_min, x_max], [0., 0.], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_x)
    ln_axis_y = art3d.Line3D([0., 0.], [y_min, y_max], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_y)
    ln_axis_z = art3d.Line3D([0., 0.], [0., 0.], [z_min, z_max], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_z)


def draw_static_diagrams():
    create_center_lines()


def update_diagrams():
    pass


def reset():
    global is_play
    cnt.reset()
    if is_play:
        is_play = not is_play


def switch():
    global is_play
    is_play = not is_play


def update(f):
    if is_play:
        cnt.count_up()
        update_diagrams()


""" main loop """
if __name__ == "__main__":
    cnt = Counter(ax=ax0, is3d=True, xy=np.array([x_min, y_max]), z=z_max, label="Step=")
    draw_static_diagrams()
    create_animation_control()
    create_parameter_setter()

    hyperboloid = Hyperboloid(ax0, "viridis", "none", 0.2)
    unit_sphere = UnitSphere(ax0, "plasma", "none", 0.2)
    poincare_disk = CircularDisk(ax0, 0., "pink", "pink", 0.4)
    klein_disk = CircularDisk(ax0, 1., "skyblue", "skyblue", 0.4)

    line_klein = LineKlein(ax0, "blue", 1)

    circle_sphere = CircleSphere(ax0, "magenta", 1)

    anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
    root.mainloop()