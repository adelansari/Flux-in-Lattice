# Flux-in-Lattice
Solving non-homogenous linear system using the following methods:
1. Triangular Factorization Method
2. Gaussian Elimination Method
3. Jacobi Iterative method
4. Gauss-Seidel Iterative Method
5. Successive Relaxation Method

The code "project.m" will perform the following tasks:
1. Ask the user to specify the length of the region or mesh size.
2. Specify the number of divisions or number of meshes.
3. Enter the macroscopic absorption cross section.
4. Enter the Diffusion coefficient.
5. Entering the flux values.
6. Displaying the flux coefficients in Matrix form.
7. Entering the source values (asks the user to choose between entering it manually at each division or specifying the source-function in terms of the variable rad).
8. Asking the user to choose the method to use for solving the matrix.
9. Check and show errors if any (maximum element in column X is zero and no exact solution for singular matrix)
10. Solve according to the selected method.
10. Display the number of iterations in each method.
11. Print the Flux and Norm in the last iteration.
12. Plot the flux & the norm for the specified method.
13. Ask the user whether to continue and try another method or exit the program.
