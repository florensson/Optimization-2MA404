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
(concateneringen av identitestmatrisen garanterar lösningar då det är ett krav för simplex lösningar.)

Linjärt oberoende kolumner i 
�
A säkerställer att varje hörnlösning är en punkt där precis 
�
n variabler (motsvarande antalet variabler i problemet) är bundna till icke-nollvärden.
Detta möjliggör en entydig uppsättning värden för variablerna och därmed en väldefinierad lösning.
Om kolumnerna inte är linjärt oberoende, kan det finnas oändligt många lösningar, och simplexmetoden kan misslyckas med att hitta den optimala lösningen.
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

tid = [] #Delclear array
test = 1 
var = 1
tick = 1 #keeps track of how many tries

while True:
    n, m = var, var  # Variabler som går att ändra
    A = np.concatenate([rand(m, n) + 1, np.eye(m)], axis=-1)
    b = np.ones(m)
    c = np.concatenate([randn(n), np.zeros(m)])

    start_time = perf_counter()
    result = linprog(-c, A_eq=A, b_eq=b, method='simplex')
    elapsed_time = 1000 * (perf_counter() - start_time)

    print(f"Attempt {tick}: The simplex method elapsed time: {elapsed_time} ms for a {var}x{var} matrix") # print out the nxm dim, the number of attemps and the time that i took for tht try

    tid.append(elapsed_time)

    if elapsed_time > 1000:
        break

    var += 10
    tick += 1

print(tid)

'''
iv) Now we will change the method to a sophisticated and state-of-the-art algorithm (with its roots
in simplex). Supply the default LP-method HiGHS to linprog, either by removing the method
keyword or passing method = ’highs’. On my current machine the average time exceeds 1
second at m = n = 2800. How high must you now go for the 1 second average?


Här kommer vi köra samma kod som ovan med ta bort method = 'simplex' så vi kör Highs
'''


tid = []
test = 1
var = 700
tick = 1 #keeps track of how many tries

while True:
    n, m = var, var  # Variabler som går att ändra
    A = np.concatenate([rand(m, n) + 1, np.eye(m)], axis=-1)
    b = np.ones(m)
    c = np.concatenate([randn(n), np.zeros(m)])

    start_time = perf_counter()
    result = linprog(-c, A_eq=A, b_eq=b,)
    elapsed_time = 1000 * (perf_counter() - start_time)

    print(f"Attempt {tick}: The simplex method elapsed time: {elapsed_time} ms for a {var}x{var} matrix") # print out the nxm dim, the number of attemps and the time that i took for tht try

    tid.append(elapsed_time)

    if elapsed_time > 1000:
        break

    var += 10
    tick += 1

print(tid)
