#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:24:15 2023

@author: fredriklorensson
"""
import numpy as np
from scipy.optimize import linprog


# Primal problem coefficients
c_primal = np.array([-700, -1000])
A_primal = np.array([[3, 5], [1, 3], [2, 2]])
b_primal = np.array([3900, 2100, 2200])

# Dual problem coefficients
c_dual = np.array([3900, 2100, 2200])
A_dual = -A_primal.T
b_dual = np.array([-700, -1000])

# Solve primal problem
result_primal = linprog(c_primal, A_ub=A_primal, b_ub=b_primal)
x_optimal_primal = result_primal.x
z_optimal_primal = result_primal.fun

# Solve dual problem
result_dual = linprog(c_dual, A_ub=A_dual, b_ub=b_dual, method='highs')  # Använd method='highs'
y_optimal_dual = result_dual.x
w_optimal_dual = result_dual.fun

# Verify that the solutions are close enough
tolerance = 1e-6

if np.isclose(z_optimal_primal, w_optimal_dual, rtol=tolerance):
    print("The optimal solutions of the primal and dual problems are close enough.")
    print(f"Optimal primal solution (x): {x_optimal_primal}")
    print(f"Optimal dual solution (y): {y_optimal_dual}")
else:
    print("The optimal solutions of the primal and dual problems are not close enough.")
    print(f"Optimal primal solution (x): {x_optimal_primal}")
    print(f"Optimal dual solution (y): {y_optimal_dual}")

# Debugging print of the actual optimal values
print(f"Optimal primal objective value: {z_optimal_primal}")
print(f"Optimal dual objective value: {w_optimal_dual}")


c_primal = np.array([-700, -1000])
A_primal = np.array([[3, 5], [1, 3], [2, 2]])
b_primal = np.array([3900, 2100, 2200])

# Lös det primala problemet
result_primal = linprog(c_primal, A_ub=A_primal, b_ub=b_primal)
x_optimal_primal = result_primal.x

# Hämta skuggpriserna (dualvärden)
skuggpriser = result_primal.fun

# Visa resultaten
print("Optimal primal lösning (x):", x_optimal_primal)
print("Skuggpriser (dualvärden):", skuggpriser)

# Assuming shadow_prices contains the shadow prices obtained from the linear programming solution

# Find the index of the stage with the highest positive shadow price
max_positive_index = np.argmax(shadow_prices)

# Find the index of the stage with the most negative shadow price
min_negative_index = np.argmin(shadow_prices)

# Display recommendations
print(f"Recommendation: Invest the additional 100 hours in Stage {max_positive_index + 1}.")
print(f"Avoidance Recommendation: Avoid investing the additional hours in Stage {min_negative_index + 1}.")
Replace shadow_prices with the actual array of shadow prices obtained from your linear 



