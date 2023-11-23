#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 09:12:57 2023

@author: fredriklorensson
"""




import numpy as np
from numpy.random import rand, randn
from scipy.optimize import linprog, OptimizeResult
from time import perf_counter

n = 10
m = 10

#size 200
A = np.concatenate ([rand(m,n) + 1, np.eye(m)], axis = -1) 

#Ger en array med 10 ettor 1x10
b = np.ones(m)


c = np.concatenate ([randn (n), np.zeros(m)])


# Lös instansen med simplex-metoden och mät tiden
start_time = perf_counter()

#simplex använder sig av vertex och jämför med målfunktionen för att hitta en optimal lösning
result_simplex = linprog(c, A_ub=A, b_ub=b, method="simplex")

end_time = perf_counter()
elapsed_time = (end_time - start_time) * 1000  # Omvandla till millisekunder

# Visa resultat och tid
print("Optimal solution:", result_simplex.x)

#Får ett negativet värde
print("Optimal value:", result_simplex.fun)

# ger olika medelande om det gick bra, om det är unbound osv om det var några problem
print("Status:", result_simplex.message)

#Tiden
print(f"Elapsed time: {elapsed_time:.5f} milliseconds")


'''
Testkörning för att se näör vi når 1000ms
'''


def run_experiment(n_variables):
    elapsed_times = []
    n_trials = 5

    for _ in range(n_trials):
        # Samma som A,b och c som ovan (lägg in så de impoteras i funktionen istället)
        n, m = n_variables, n_variables
        A = np.concatenate([rand(m, n) + 1, np.eye(m)], axis=-1)
        b = np.ones(m)
        c = np.concatenate([randn(n), np.zeros(m)])

        # Lös instansen med simplex-metoden och mät tiden
        start_time = perf_counter()
        # Tydligen används HiGHS automatiskt om ingen metod anges
        result_simplex = linprog(c, A_ub=A, b_ub=b)
        end_time = perf_counter()
        elapsed_time = (end_time - start_time) * 1000  # Millisekunder
        elapsed_times.append(elapsed_time)

    # Beräkna genomsnittlig tid
    average_time = np.mean(elapsed_times)
    return average_time

# Målet är att hitta det lägsta n_variables där genomsnittlig tid överstiger 1000 ms (1 sekund)
n_variables = 1
target_time = 1000  # 1 sekund

while True:
    average_time = run_experiment(n_variables)
    print(f"Number of variables: {n_variables}, Average time: {average_time:.2f} ms")

    if average_time > target_time:
        break

    n_variables += 100

print(f"\nAverage time exceeds 1 second at m = n = {n_variables}.")


























