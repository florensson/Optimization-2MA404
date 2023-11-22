import numpy as np
from scipy.optimize import linprog
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

'''
Exercise 1 from 2Ma404
A company is producing 2 types of television set

max c^T x = z
s.t. Ax =< b
x >= 0
        A    B    total:
Profit: 700  1000
Steg 1: -3   -5   3900
Steg 2: -1   -3   2100
Steg 3: -2   -2   2200           
'''

# Objective function
c = np.array([700, 1000]) # koefficienter för objektiva funktionen
A = np.array([
        [-3, -5], #Koefficienter för arbetstimmar i steg 1
        [-1, -3], #Koefficienter för arbetstimmar i steg 2
        [-2, -2]  #Koefficienter för arbetstimmar i steg 3
])

b = np.array([-3900, -2100, -2200]) #högst antal timmar tillgängliga

# Gränser för beslutsvariabler
x_bounds = (0, None)

# Lösning med LP-problem med hjälp av linprog
result = linprog(c, A_ub=A, b_ub=b, bounds=[x_bounds, x_bounds])

# Resultatet för x
max_profit = result.fun #minimizing problem 

print("Optimala värden för x_1, x_2 är:", result.x)
print("Maximal vinst:", max_profit)

### II 

x = np.linspace(0, 1200, 100)

# Reformulate each line equation as yi := 1 ai;2 (bi ai;1x) for i = 1; 2; 3
y1 = (b[0] - A[0, 1] * np.clip(x, 0, 1200)) / A[0, 1]
y2 = (b[1] - A[1, 1] * np.clip(x, 0, 1200)) / A[1, 1]
y3 = (b[2] - A[2, 1] * np.clip(x, 0, 1200)) / A[2, 1]

# Plot all lines in the same figure
plt.plot(x, y1, label='Constraint 1')
plt.plot(x, y2, label='Constraint 2')
plt.plot(x, y3, label='Constraint 3')

# Set a legend
plt.legend()

# Set appropriate y-limits
plt.ylim(0, 1200)

# Determine the vertices of the feasible region
vertices = []
for i in range(0, 3):
    for j in range(i + 1, 3):
        if (y1[i] == y1[j]) and (y2[i] == y2[j]):
            vertices.append([x[i], y1[i]])
        if (y1[i] == y2[j]) and (y2[i] == y3[j]):
            vertices.append([x[i], y2[j]])
        if (y1[i] == y3[j]) and (y2[j] == y3[j]):
            vertices.append([x[i], y3[j]])

# Create a numpy-array of the vertices
vertices = np.array(vertices)

# Plot the feasible region
plt.plot(vertices[:, 0], vertices[:, 1], 'r-', label='Feasible region')

# Show the plot
plt.show()
