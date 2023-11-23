#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:18:13 2023

@author: fredriklorensson
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import linprog

# Definiera beslutsvariabler och koefficienter för objektivfunktion
c = np.array([-700, -1000])

# Definiera koefficienter för olikhetsbegränsningarna Ax <= b
A = np.array([[3, 5], [1, 3], [2, 2]])

# Definiera högerledet för olikhetsbegränsningarna Ax <= b
b = np.array([3900, 2100, 2200])

# 2) Rita linjerna för olikhetsbegränsningarna
x = np.linspace(0, 1200, 100)

# Rita linjerna för olikhetsbegränsningarna
#(b[0]−A[0,0]⋅x)/A[0,1]
plt.plot(x, (b[0] - A[0, 0] * x) / A[0, 1], label=r'$3x_1 + 5x_2 \leq 3900$')
plt.plot(x, (b[1] - A[1, 0] * x) / A[1, 1], label=r'$x_1 + 3x_2 \leq 2100$')
plt.plot(x, (b[2] - A[2, 0] * x) / A[2, 1], label=r'$2x_1 + 2x_2 \leq 2200$')

# Ställ in förklaring och y-gränser
plt.legend()
plt.ylim(0, 1200)

# Visa plotten
plt.show()

# 3) Hitta hörn på den möjliga regionen
vertices = np.array([
    [600, 600],
    [800, 400],
    [1000, 0],
    [0, 800],
    [0, 400]
])

# 4) Skapa och plotta konvext hölje
hull = ConvexHull(vertices)
plt.plot(vertices[:, 0], vertices[:, 1], 'o')  # Plotta hörn
for simplex in hull.simplices:
    plt.plot(vertices[simplex, 0], vertices[simplex, 1], 'k-')  # Plotta konvext hölje

plt.ylim(0, 800)
plt.show()

# 5) Använd linprog för att hitta optimala lösningen
result = linprog(c, A_ub=A, b_ub=b)
optimal_x = result.x
max_profit = -result.fun  # Negation eftersom linprog löser minimiseringsproblem

print("Optimala värden för beslutsvariablerna (x1, x2):", optimal_x)
print("Maximal vinst:", max_profit)

# 6) Use linprog to find the maximum
result = linprog(c, A_ub=A, b_ub=b)
optimal_x = result.x
max_profit = -result.fun  # Negation eftersom linprog löser minimiseringsproblem, max värde för funktionen

print("Optimal solution using linprog:")
print("Optimal values of x:", optimal_x)
print("Maximum profit:", max_profit)