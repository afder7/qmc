import numpy as np
from scipy.stats import qmc
from math import *

num = 128
dim = 5


def dist(x_i, x_j):
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


def potential_energy(points):
    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


with open(f"random_num_{num}_dim_{dim}.txt", "r") as f:
    points = []
    k = 0
    for i in f.readlines():
        points.append(list(map(float, i.strip().split())))
        k += 1
    points = np.array(points)
    r_l2_disc = qmc.discrepancy(points, method='L2-star')
    r_u = potential_energy(points)

with open(f"modeled_num_{num}_dim_{dim}.txt", "r") as f:
    points = []
    for i in f.readlines():
        points.append(list(map(float, i.strip().split())))
    points = np.array(points)
    m_l2_disc = qmc.discrepancy(points, method='L2-star')
    m_u = potential_energy(points)


print("Random L2-star discrepancy: ", r_l2_disc)
print("Modeled L2-star discrepancy:", m_l2_disc)
print("")

print("Random Potential Energy: ", r_u)
print("Modeled Potential Energy:", m_u)
print("")
