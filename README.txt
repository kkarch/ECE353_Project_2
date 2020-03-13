############################################

Authors:
Kevin Karch (KSK62) 
Madeline Cook (MMC378)

Drexel University

ECEC 353 - Assignment 2

GitHub Repository: https://github.com/kkarch/ECE353_Project_2

############################################

To Run:

----------------------------------------------------------------

Integration/trap.c

Compile as: gcc -o trap trap.c -O3 -std=c99 -Wall -lpthread -lm

Usage: ./trap a b n num_threads


Purpose: Calculate definite integral using trapezoidal rule.
Input:   a, b, n, num_threads
Output:  Estimate of integral from a to b of f(x) using n trapezoids, with num_threads.

Argument Definitions-
a=lower limit of integral
b=upper limit of integral
n=number of trapezoids to be used in calculation
num_threads=number of threads used to calculate integral

------------------------------------------------------------------

Jacobi/solver.c

Compile as: gcc -o solver solver.c solver_gold.c -O3 -Wall -std=c99 -lm -lpthread -D_GNU_SOURCE

Usage: ./solver grid_size number_of_threads [min_temp max_temp]


Purpose: Calculate heat diffusion using Jacobi method.
Input:   grid_size, number_of_threads, [min_temp max_temp]
Output:  Estimate of integral from a to b of f(x) using n trapezoids, with num_threads.

Argument Definitions-
grid_size = single dimension (n) of an n x n matrix
num_threads = number of threads used to calculate Jacobian Matrix
min_temp = lowest possible temperature of the applied heat row
max_temp = highest possible temperature of the applied heat row
