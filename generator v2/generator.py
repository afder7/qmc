import random
from math import *
from scipy.stats import qmc
import matplotlib.pyplot as plt
import signal


number_of_points = 1024
dimension = 20
iter_count = 1000
delta_t = 0.01
kappa = 200

it = 0

G = 1
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_U = 0

counter = 0
x_prev = []
x_0 = []


potential_energy_list = []
l2_discrepancy_list = []


def dist(x_i, x_j):
    dim = len(x_i)
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


# points is a set of num points of class Point in a given dim
# U_{1,1}
def potential_energy(points):
    num = number_of_points
    dim = dimension

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


def force(points):
    num = number_of_points
    dim = dimension
    f = [[0 for __ in range(dim)] for _ in range(num)]

    for i in range(num):
        for j in range(num):
            if j != i:
                for k in range(dim):
                    sign = -1 if points[i][k] < points[j][k] else 1
                    delta = abs(points[i][k] - points[j][k])
                    a_ijk = sign * delta * (1 - delta) * (1 - 2 * delta)
                    f[i][k] += -G * a_ijk / dist(points[i], points[j]) ** 3

    return f


class MovingPoints:
    def __init__(self, point_count=number_of_points, dim=dimension):
        self.num = point_count
        self.dim = dim
        self.points = []
        self.velocities = []

        self.create_points()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_U

        for i in range(self.num):
            point = [random.random() for _ in range(self.dim)]
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        U = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(U) * delta_t ** 2 / scalar
        c = -abs(U) * sqrt(kappa) * delta_t / scalar
        cur_U = U

        # calculating x_1
        f = force(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                self.points[i][k] = self.points[i][k] + 0.5 * (1 / m) * f[i][k] * delta_t ** 2

        save_points(self, "random")

    def update(self):
        global x_prev
        f = force(self.points)

        new_x_prev = self.points
        for i in range(self.num):
            for k in range(self.dim):
                x_j = self.points[i][k]
                self.points[i][k] = (-f[i][k] + A1 * x_j - A2 * x_prev[i][k]) / A0
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1

        x_prev = new_x_prev.copy()

        global cur_U
        cur_U = potential_energy(self.points)
        # print("p1               ", f"({point_array[0][0]}, {point_array[0][1]})")


def save_points(app, type):
    with open(f"generation/{type}_num_{app.num}_dim_{app.dim}.txt", "w") as fi:
        for point in app.points:
            for coord in point:
                if coord > 1:
                    fi.write(str(coord - int(coord)) + " ")
                elif coord < 0:
                    fi.write(str(coord - int(coord) + 1) + " ")
                else:
                    fi.write(str(coord) + " ")
            fi.write("\n")


def build_graph():
    plt.figure(figsize=(8, 6))

    plt.subplot(2, 1, 1)  # 2 строки, 1 столбец, позиция 1
    plt.plot(list(range(it)), potential_energy_list, 'g-', label='U')
    plt.legend()

    plt.subplot(2, 1, 2)  # 2 строки, 1 столбец, позиция 2
    plt.plot(list(range(it)), l2_discrepancy_list, 'm--', label='L2')
    plt.legend()

    plt.tight_layout()
    plt.show()


app = MovingPoints(dim=dimension, point_count=number_of_points)


def save_on_exit(signum, frame):
    print(f"\nSaving results of current iteration #{it}")
    save_points(app, "modeled")
    build_graph()
    exit(0)


signal.signal(signal.SIGINT, save_on_exit)
signal.signal(signal.SIGTERM, save_on_exit)


for _ in range(iter_count):
    app.update()

    l2_disc = qmc.discrepancy(app.points, method='L2-star')
    potential_energy_list.append(cur_U)
    l2_discrepancy_list.append(l2_disc)

    if it % 10 == 0:
        print("Iteration        ", it)
        print(f"Potential Energy", cur_U)
        print(f"L2 discrepancy  ", l2_disc)
        print("")
    it += 1

save_points(app, "modeled")
build_graph()
