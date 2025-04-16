import random
from math import *
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt


number_of_points = 128
dimension = 5
iter_count = 500
delta_t = 0.01
kappa = 200

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


class Point:
    def __init__(self, coords):
        self.coords = coords


def dist(x_i, x_j):
    dim = len(x_i)
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


# points is a set of num points of class Point in a given dim
# U_{1,1}
def potential_energy(points):
    num = len(points)
    dim = len(points[0].coords)

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i].coords, points[j].coords)

    return U


def force(points):
    num = len(points)
    dim = len(points[0].coords)
    f = [[0 for __ in range(dim)] for _ in range(num)]

    for i in range(num):
        for j in range(num):
            d = dist(points[i].coords, points[j].coords)
            if j != i:
                for k in range(dim):
                    sign = -1 if points[i].coords[k] < points[j].coords[k] else 1
                    delta = abs(points[i].coords[k] - points[j].coords[k])
                    a_ijk = sign * delta * (1 - delta) * (1 - 2 * delta)
                    f[i][k] += -G * a_ijk / (d ** 3)

    return f


class MovingPoints:
    def __init__(self, point_count=100, dim=2):
        self.num = point_count
        self.dim = dim
        self.points = []
        self.velocities = []

        self.create_points()

        self.update()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_U

        for i in range(self.num):
            point = Point(coords=[random.random() for _ in range(self.dim)])
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        U = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i.coords:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(U) * delta_t ** 2 / scalar
        c = -abs(U) * sqrt(kappa) * delta_t / scalar
        cur_U = U

        # calculating x_1
        f = force(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                self.points[i].coords[k] = self.points[i].coords[k] + 0.5 * (1 / m) * f[i][k] * delta_t ** 2

        save_points(self, "random")

    def update(self):
        global x_prev
        f = force(self.points)

        new_x_prev = self.points
        for i in range(self.num):
            for k in range(self.dim):
                x_j = self.points[i].coords[k]
                self.points[i].coords[k] = (-f[i][k] + A1 * x_j - A2 * x_prev[i].coords[k]) / A0
                if self.points[i].coords[k] > 1:
                    self.points[i].coords[k] = self.points[i].coords[k] - int(self.points[i].coords[k])
                elif self.points[i].coords[k] < 0:
                    self.points[i].coords[k] = self.points[i].coords[k] - int(self.points[i].coords[k]) + 1

        x_prev = new_x_prev.copy()

        global cur_U
        cur_U = potential_energy(self.points)
        point_array = []
        for pt in self.points:
            point_array.append(pt.coords)
        point_array = np.array(point_array)
        l2_disc = qmc.discrepancy(point_array, method='L2-star')
        potential_energy_list.append(cur_U)
        l2_discrepancy_list.append(l2_disc)

        print(f"Potential Energy", cur_U)
        print(f"L2 discrepancy  ", l2_disc)
        # print("p1               ", f"({point_array[0][0]}, {point_array[0][1]})")


def save_points(app, type):
    with open(f"{type}_num_{app.num}_dim_{app.dim}.txt", "w") as fi:
        for point in app.points:
            for coord in point.coords:
                if coord > 1:
                    fi.write(str(coord - int(coord)) + " ")
                elif coord < 0:
                    fi.write(str(coord - int(coord) + 1) + " ")
                else:
                    fi.write(str(coord) + " ")
            fi.write("\n")


app = MovingPoints(dim=dimension, point_count=number_of_points)

for i in range(iter_count):
    app.update()

save_points(app, "modeled")

plt.figure(figsize=(8, 6))

plt.subplot(2, 1, 1)  # 2 строки, 1 столбец, позиция 1
plt.plot(list(range(iter_count + 1)), potential_energy_list, 'g-', label='U')
plt.legend()

plt.subplot(2, 1, 2)  # 2 строки, 1 столбец, позиция 2
plt.plot(list(range(iter_count + 1)), l2_discrepancy_list, 'm--', label='L2')
plt.legend()

plt.tight_layout()
plt.show()
