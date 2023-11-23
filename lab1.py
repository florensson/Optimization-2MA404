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


#eyes skapar en array med 1 på diagonalen och 0 överallt annars
a2 = np.eye(m)

#En matris med 10x10 i storlek som har randoms nummer
a1 = rand(m,n) + 1

#size 200, axis = -1 tar den sista dim
A = np.concatenate ([a1, a2], axis = 1) 

#Ger en array med 10 ettor 1x10 (vektor)
b = np.ones(m)

# column vector med först 10 random sen 10 nollor 20x1
c = np.concatenate ([randn (n), np.zeros(m)])


# Lös instansen med simplex-metoden och mät tiden
start_time = perf_counter()

#simplex använder sig av vertex och jämför med målfunktionen för att hitta en optimal lösning
result_simplex = linprog(c, A_ub=A, b_ub=b, method="simplex")

end_time = perf_counter()
elapsed_time = (end_time - start_time) * 1000  # Omvandla till millisekunder

# Visa resultat och tid
print("Bästa lösningen:", result_simplex.x)

#Får ett negativet värde
print("Bästa värdet:", result_simplex.fun)

# ger olika medelande om det gick bra, om det är unbound osv om det var några problem
print("Status:", result_simplex.message)

#Tiden
print("Tid:",elapsed_time,  "millisec")


'''
Testkörning för att se när vi når 1000ms
Normalt hade man nog lyft ut första delen i loopen och skapat en egen funktion, men jag är lat ibland
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


























