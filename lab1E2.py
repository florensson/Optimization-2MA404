#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 09:12:57 2023

@author: fredriklorensson
"""


import numpy as np
from numpy . random import rand , randn
from scipy . optimize import linprog , OptimizeResult
from time import perf_counter
# Example of a problem instance

'''
i) By formulating the matrix A and the vectors b, c in the above manner we guarantee a solution.
Why?

Svar: Jag vet inte 
'''

n,m = 10,10 #variabler som går att ändra
A = np.concatenate ([rand(m,n) + 1, np.eye(m)], axis =-1)
b = np.ones(m)
c = np.concatenate ([randn(n), np.zeros(m)])

'''

ii) We warm up by solving the example instance using the simplex method by supplying method =
"simplex" to linprog. Further, we can estimate the elapsed time for the procedure in milliseconds using perf_counter().
Running a few times shows that it takes roughly 3 ms on my machine. How about your machine?

Svar: ca 5-7ms på min maskin
'''


start_time = perf_counter()
result = linprog (-c, A_eq = A, b_eq = b, method = 'simplex')

elapsed_time = 1000*(perf_counter() - start_time)


'''
iii) Keep m = n and steadily increase the number of variables using the ’simplex’ method. For each
value, perform the experiment five times by re-sampling and then present the average elapsed
time. On my current machine the average time exceeds 1 second at m = n = 170. How high
must you go for the average time to exceed 1 second (1000 ms)?
'''

tid = []
test = 1
var = 1

for _ in range(test):
    n,m = var,var #variabler som går att ändra
    A = np.concatenate ([rand(m,n) + 1, np.eye(m)], axis =-1)
    b = np.ones(m)
    c = np.concatenate ([randn(n), np.zeros(m)])

    start_time = perf_counter()
    result = linprog (-c, A_eq = A, b_eq = b, method = 'simplex')

    elapsed_time = 1000*(perf_counter() - start_time)

while True:
    
    if tid < 1000:
        print(tid)
    else:
        break
    
    var += 1

'''
iv) Now we will change the method to a sophisticated and state-of-the-art algorithm (with its roots
in simplex). Supply the default LP-method HiGHS to linprog, either by removing the method
keyword or passing method = ’highs’. On my current machine the average time exceeds 1
second at m = n = 2800. How high must you now go for the 1 second average?


Här kommer vi köra samma kod som ovan med ta bort method = 'simplex' så vi kör Highs
'''


tid = []
test = 1
var = 1

for _ in range(test):
    n,m = var,var #variabler som går att ändra
    A = np.concatenate ([rand(m,n) + 1, np.eye(m)], axis =-1)
    b = np.ones(m)
    c = np.concatenate ([randn(n), np.zeros(m)])

    start_time = perf_counter()
    result = linprog (-c, A_eq = A, b_eq = b, method = 'simplex')

    elapsed_time = 1000*(perf_counter() - start_time)

while True:
    
    if tid < 1000:
        print(tid)
    else:
        break
    
    var += 1
























'''
# Extract fractional seconds ( high resolution )
start_time = perf_counter ()
# Optimize with ’simplex ’
result = linprog ( -c ,
A_eq =A ,
b_eq =b ,
method =’simplex’,
options ={ ’maxiter ’: 5000})
# Est. elapsed time in milliseconds
elapsed_time = 1000*( perf_counter () - start_time )
'''




'''
n = 10
m = 10

#def calc():
#    m,n = 10,10
#    A = np.concatenate ([a1, a2], axis = 1)
    
#eyes skapar en array med 1 på diagonalen och 0 överallt annars
a2 = np.eye(m)

#En matris med 10x10 i storlek som har randoms nummer, random nummer 0 < x < 1
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

'''
Testkörning för att se när vi når 1000ms
Normalt hade man nog lyft ut första delen i loopen och skapat en egen funktion, men jag är lat ibland
'''



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




'''





















