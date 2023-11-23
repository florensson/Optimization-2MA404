#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:18:13 2023

@author: fredriklorensson
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d

# i) Define the model with decision variables x1, x2 and create A, b, and c using numpy in Python.

# Coefficients of the objective function c^T x
c = np.array([-700, -1000])

# Coefficients of the inequality constraints Ax <= b
A = np.array([[3, 5], [1, 3], [2, 2]])

# Right-hand side of the inequality constraints Ax <= b
b = np.array([3900, 2100, 2200])

# ii) Plot the constraint lines

x = np.linspace(0, 1200, 100)

# Line equations for each constraint
y1 = (b[0] - A[0, 0] * x) / A[0, 1]
y2 = (b[1] - A[1, 0] * x) / A[1, 1]
y3 = (b[2] - A[2, 0] * x) / A[2, 1]

# Plot the constraint lines
plt.plot(x, y1, label=r'$3x_1 + 5x_2 \leq 3900$')
plt.plot(x, y2, label=r'$x_1 + 3x_2 \leq 2100$')
plt.plot(x, y3, label=r'$2x_1 + 2x_2 \leq 2200$')

# Set legend and y-limits
plt.legend()
plt.ylim(0, 1200)

# Show the plot
plt.show()

# iii) Find vertices of the feasible region

# Vertices (intersection points of constraint lines)
vertex1 = np.array([600, 600])  # Example, adjust based on the plot
vertex2 = np.array([800, 400])  # Example, adjust based on the plot
vertex3 = np.array([1000, 0])   # Example, adjust based on the plot
vertex4 = np.array([0, 800])    # Example, adjust based on the plot
vertex5 = np.array([0, 400])    # Example, adjust based on the plot

vertices = np.array([vertex1, vertex2, vertex3, vertex4, vertex5])

# iv) Plot the convex hull and level curves

# Create convex hull
hull = ConvexHull(vertices)

# Plot the convex hull
convex_hull_plot_2d(hull)

# Plot level curves for z = c^T x
C_values = np.arange(1, 9) * 1e5
for k, C in enumerate(C_values, start=1):
    z_values = (C - c[0] * x) / c[1]
    plt.plot(x, z_values, 'r', alpha=k / 8)

# Set y-limits and show the plot
plt.ylim(0, 800)
plt.show()
