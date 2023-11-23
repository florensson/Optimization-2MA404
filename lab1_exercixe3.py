# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:48:49 2023

@author: Fredr
"""

import numpy as np
from scipy.optimize import linprog

# Function to solve primal and dual problems and make recommendations
def solve_and_recommend(c_primal, A_primal, b_primal, method_primal, method_dual):
    # Solve the primal problem
    result_primal = linprog(c_primal, A_ub=A_primal, b_ub=b_primal, method=method_primal)
    x_optimal_primal = result_primal.x
    z_optimal_primal = result_primal.fun
    shadow_prices = result_primal.fun  # Assume shadow prices are the same as primal

    # Solve the dual problem
    c_dual = np.array([3900, 2100, 2200])
    A_dual = -A_primal.T
    b_dual = np.array([-700, -1000])
    result_dual = linprog(c_dual, A_ub=A_dual, b_ub=b_dual, method=method_dual)
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

    # Make recommendations based on shadow prices
    max_positive_index = np.argmax(shadow_prices)
    min_negative_index = np.argmin(shadow_prices)
    print(f"Recommendation: Invest the additional 100 hours in Stage {max_positive_index + 1}.")
    print(f"Avoidance Recommendation: Avoid investing the additional hours in Stage {min_negative_index + 1}.")

# Example data for the primal problem
c_primal_example = np.array([-700, -1000])
A_primal_example = np.array([[3, 5], [1, 3], [2, 2]])
b_primal_example = np.array([3900, 2100, 2200])

# Call the function with example data
solve_and_recommend(c_primal_example, A_primal_example, b_primal_example, 'simplex', 'simplex')

# Code to handle the new type of TV (point v)
c_new_tv = np.array([7, 4, 2])
A_primal_new_tv = np.column_stack([A_primal_example, np.zeros_like(c_new_tv), c_new_tv])  # Corrected this line
b_primal_new_tv = np.append(b_primal_example, 100)  # Add 100 extra work hours for the new TV

# Call the function with the new TV type
solve_and_recommend(c_primal_example, A_primal_new_tv, b_primal_new_tv, 'simplex', 'simplex')

# Code to handle the new production line (point vi)
# Assume that inspection_times contains inspection time for each type of TV (A, B, C)
inspection_times = np.array([0.5, 0.75, 0.1 / 60])  # 0.1 minutes for TV of type C
A_primal_inspection = np.column_stack([A_primal_example, np.zeros_like(inspection_times), inspection_times])
b_primal_inspection = b_primal_example

# Call the function with the new production line
solve_and_recommend(c_primal_example, A_primal_inspection, b_primal_inspection, 'simplex', 'simplex')
