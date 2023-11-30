#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:18:13 2023

@author: fredriklorensson

Exercise 1: Warm-up LP-problem

A company produces two types of television sets, a cheap type (A) and an expensive type (B). The
company has a profit of 700 per unit for type A and 1000 per unit for type B. There are three stages
in the production process. In stage I, 3 hours are needed per unit of type A and 5 hours per unit of
type B. The total number of available working hours in stage I is 3900. In step II, 1 hour is worked
on each device of type A and for 3 hours on each device of type B. The total working time available
is 2100 hours. In step III, 2 working hours are needed for both types, and 2200 working hours are
available. How many TV sets of each type should the company produce to maximize its profit?
As a first step, we want to create matrices A; b and c such that the problem is formulated in terms
of
max c^T x = z
s.t. Ax =< b
x >= 0:
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import linprog

"""
i) Define the model with decision variables x1, x2 and create A; b and c using numpy in Python
"""

# Definiera beslutsvariabler och koefficienter för objektivfunktion
c = np.array([-700, -1000])

# Definiera koefficienter för olikhetsbegränsningarna Ax <= b
A = np.array([[3, 5], [1, 3], [2, 2]])

# Definiera högerledet för olikhetsbegränsningarna Ax <= b
b = np.array([3900, 2100, 2200])

"""
ii) The feasible region is a convex set restricted by the x1- and x2-axis and the lines
ai;1x1 + ai;2x2 = bi
;
for i = 1; 2; 3. We want to plot these lines to get a good overview of the problem. Create a
numpy-array using linspace that goes between 0 and 1200 with a suitable number of points, e.g.,
x = np.linspace(0, 1200, 100):
This now represents some x1-values.
Reformulate each line equation as yi:=1

ai;2

(bi  ai;1x) for i = 1; 2; 3 and compute the corresponding points of x2-values using x. Plot all lines in the same figure by importing matplotlib.pyplot
as plt in your script and then use plt.plot for each pair (x; yi). Be sure to set a legend by supplying plt.legend with a list of three strings to clearly see which constraint each line represents
(the order must follow the order of calls to plt.plot). Also set appropriate y-limits using
plt.ylim(0, 1200)
and call plt.show() in the end to visualize your figure
"""

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

"""

iii) Study the figure produced in (ii). There are five vertices that enclose the feasible region. Each
of these will be the intersection of two lines (including axis lines). Graphically determine their
values, i.e., see which two constraints are satisfied with equality for each vertex. Store the vertices
in a 5 by 2 numpy-array.
The feasible region X in our case is precisely the convex hull of the five vertices (smallest convex set containing them). So add an import of ConvexHull and convex_hull_plot_2d from
scipy.spatial in your script. Create a hull and plot it all in one go by
convex_hull_plot_2d(ConvexHull(vertices));
where vertices are your numpy-array of points. Adjust y-values with plt.ylim(0, 800). Note
that this will always begin a new figure and we will use this in the next step so delay the call to
plt.show() a little bit (or move it down later

"""


# 3) Hitta hörn på den möjliga regionen
vertices = np.array([
    [600, 600],
    [800, 400],
    [1000, 0],
    [0, 800],
    [0, 400]
])


'''

iv) Let’s also plot some level curves for z = c
T x in the convex hull plot. Using the previous x-array
for x1-values, figure out what the corresponding x2-values must be if z = C for
C = 1 * 105, 2 * 105,...,8 * 105
:

If the x2-values for a specific C := k * 105
is y(k), then plot all curves with a suitable for-loop
indexed by k and using
plt.plot(x, y(k), ’r’, alpha=k/8):
This will now plot eight level curves of z on top of your hull with an increasing red color gradient
for larger values of z. Determine from your plot in which vertex the maximum is obtained.
'''

# 4) Skapa och plotta konvext hölje
hull = ConvexHull(vertices)
plt.plot(vertices[:, 0], vertices[:, 1], 'o')  # Plotta hörn
for simplex in hull.simplices:
    plt.plot(vertices[simplex, 0], vertices[simplex, 1], 'k-')  # Plotta konvext hölje

plt.ylim(0, 800)
plt.show()

'''
v) Verify by checking the value of z in all extreme points. This can be done by looping over each
vertex v in your collection and using c @ v (numpy matrix multiplication).
'''

# 5) Använd linprog för att hitta optimala lösningen
result = linprog(c, A_ub=A, b_ub=b)
optimal_x = result.x
max_profit = -result.fun  # Negation eftersom linprog löser minimiseringsproblem

print("Optimala värden för beslutsvariablerna (x1, x2):", optimal_x)
print("Maximal vinst:", max_profit)

'''
vi) Now that we know the ground truth we will instead use scipy to find the maximum for us. Import
linprog and OptimizeResult from scipy.optimize. Have a good look at the documentation
of both (1 2)
(can also be found in their respective definitions). You will need to provide arguments
to linprog in precisely the right way and then know how to interpret the result in the returned
OptimizeResult (which is a dictionary). Recall that max z is equivalent to min(z). Find and
extract the maximum using your defined matrices A; b and 
'''

# 6) Use linprog to find the maximum
result = linprog(c, A_ub=A, b_ub=b)
optimal_x = result.x
max_profit = -result.fun  # Negation eftersom linprog löser minimiseringsproblem, max värde för funktionen

print("Optimal solution using linprog:")
print("Optimal values of x:", optimal_x)
print("Maximum profit:", max_profit)


'''
vii) Setup the problem in standard form (using equality constraints with slack variables). Again
make sure that you have a good grasp on the usage of linprog. Solve it and verify that the
optimal solution stays the samE
'''

# Objective coefficients for the minimization problem
c_minimization = -c

# Solve the linear programming problem
result = linprog(c_minimization, A_ub=A, b_ub=b)

# Extract the maximum value
#
max_value = -result.fun if result.fun < 0 else 0.0
tolerance = 1e-8  # Just a small tolerance value
if result.fun < -tolerance:
    print(f"Maximum value using linprog: {max_value}")
else:
    print("Linprog reported an optimal value very close to zero, possibly due to numerical precision.")
    print("Consider it as a zero or near-zero value.")

# Display the maximum value
print(f"Maximum value using linprog: {max_value}")


'''
Exercise 3: Sensitivity analysis
In this exercise we will perform a sensitivity analysis of the LP-problem in Exercise 1.
'''


'''
i) Formulate and solve the dual problem. Verify that the solutions are the same.
'''


# Primal problem coefficients
# we get the arrays from Exercise 1, i)
c_primal = c
A_primal = A
b_primal = b

# Dual problem coefficients
c_dual = np.array([3900, 2100, 2200])
A_dual = -A_primal.T
b_dual = np.array([-700, -1000])

# Solve primal problem
result_primal = linprog(c_primal, A_ub=A_primal, b_ub=b_primal)
x_optimal_primal = result_primal.x
z_optimal_primal = result_primal.fun

# Solve dual problem
result_dual = linprog(c_dual, A_ub=A_dual, b_ub=b_dual, method='highs')  # Use method='highs'
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


'''
ii) Compute the shadow price of the various constraints.
'''
# Assuming shadow_prices contains the shadow prices obtained from the linear programming solution
shadow_prices = result_primal.fun

# Find the index of the stage with the highest positive shadow price
max_positive_index = np.argmax(shadow_prices)

# Find the index of the stage with the most negative shadow price
min_negative_index = np.argmin(shadow_prices)

# Display recommendations
print(f"Recommendation: Invest the additional 100 hours in Stage {max_positive_index + 1}.")
print(f"Avoidance Recommendation: Avoid investing the additional hours in Stage {min_negative_index + 1}.")

'''
iii) Suppose that you are granted 100 extra working hours in either stage I, II or III. Where would
you invest this time, and where should you not invest it?
'''

'''
iv) Suppose that the company investigates the possibility to increase the price of the tv of type (B).
How much must the price increase for the optimal solution to change?
'''

'''
v) A team of marketing people and some technicians have considered a new kind of tv to produce.
Suppose that this new type of tv is called type (C). This is a premium tv and will return a
profit of 1350. However, it takes some time to construct and the times to produce it are (7; 4; 2)
for each separate stage. Add this type of tv to the problem, and find the reduced cost for the
type (C) tv. Consider the value of the reduced cost. Should the tv be produced or not? If so,
re-run the optimization to get a new optimum.
'''


'''
vi) The company’s new quality director have decided to add a new stage in the production for
quality inspection. She claims that it’ll take 0.5 hours to inspect the type (A) tv, and 0.75
hours to inspect the more advanced tv of type (B). One upside with producing tv of type (C)
is the advanced production which involves quality inspection already in the production, so the
inspection time is only 6 minutes. The CEO do not want to decrease production and is willing
to pay for the quality inspections. How many hours must the company add to this working line,
not to interfer with the current production amounts?
'''