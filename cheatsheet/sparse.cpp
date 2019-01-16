// useful for sparse modification
// see ex01

/**
 * 
 * DIAGONAL MATRIX x vector O(n) multiplication
 * 
 * */
d;                          // diagonal matrix
x;                          // vector
Dx = d.cwiseProduct(x);     // before: O(n^2)
Dx = d.array() * x.array(); // after: O(n)

/**
 * 
 * SPARSE MATRICES page 54
 * 
 * */

/**
 * 
 * SPARSE LSE SOLVER: page 56 
 * 
 * */
