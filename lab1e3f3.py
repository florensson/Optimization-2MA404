#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:52:38 2023

@author: fredriklorensson
"""

import numpy as np
from scipy.optimize import linprog

# Primal problem coefficients
c_primal = np.array([-700, -1000])
A_primal = np.array([[3, 5], [1, 3], [2, 2]])
b_primal = np.array([3900, 2100, 2200])

# Solve primal problem
result_primal = linprog(c_primal, A_ub=A_primal, b_ub=b_primal, method='highs')
x_optimal_primal = result_primal.x
z_optimal_primal = result_primal.fun
shadow_prices = result_primal.fun  # Get the shadow prices

# Find the index of the stage with the highest positive shadow price
max_positive_index = np.argmax(shadow_prices)

# Find the index of the stage with the most negative shadow price
min_negative_index = np.argmin(shadow_prices)

# ...
# ...

# Display recommendations
print("Optimal primal solution (x):", x_optimal_primal)
print("Optimal primal objective value:", z_optimal_primal)
print("Shadow prices (dual values):", shadow_prices)

# Find the index of the stage with the highest positive shadow price
max_positive_index = np.argmax(shadow_prices)

# Find the index of the stage with the most negative shadow price
min_negative_index = np.argmin(shadow_prices)

# Display recommendations
print(f"Recommendation: Invest the additional 100 hours in Stage {max_positive_index + 1}.")
print(f"Avoidance Recommendation: Avoid investing the additional hours in Stage {min_negative_index + 1}.")
