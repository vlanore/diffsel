/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	NOTE: if ANSI function prototypes are not supported, define NO_PROTOTYPES
|		  before including this file.
*/

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */



extern int  EigenRealGeneral (int n, double **a, double *v, double *vi, double **u, int *iwork, double *work);


/*--------------------------------------------------------------------------------------------------
|
|	EigenRealGeneral
|
|	Calculate eigenvalues and eigenvectors of a general real matrix assuming that all eigenvalues
|	are real, using routines from the public domain EISPACK package.
*/

	/*      n = order of a                                                      */
	/*    **a = input matrix in row-ptr representation; will be destroyed       */
	/*     *v = array of size 'n' to receive eigenvalues                        */
	/*    *vi = work vector of size 'n' for imaginary components of eigenvalues */
	/*    **u = matrix in row-ptr representation to receive eigenvectors        */
	/* *iwork = work vector of size 'n'                                         */
	/*  *work = work vector of size 'n'                                         */

	


extern int  InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv);
/*--------------------------------------------------------------------------------------------------
|
|	InvertMatrix
|
|	Invert matrix 'a' using LU-decomposition technique, storing inverse in 'a_inv'.  Matrix 'a'
|	is destroyed.  Returns ERROR if matrix is singular, NO_ERROR otherwise.
*/

	/*     **a = matrix represented as vector of row pointers      */
	/*       n = order of matrix                                   */
	/*    *col = work vector of size n                             */
	/*   *indx = work vector of size n                             */
	/* **a_inv = inverse of input matrix a (matrix a is destroyed) */



extern int  LUDecompose (double **a, int n, double *vv, int *indx, double *pd);



/*--------------------------------------------------------------------------------------------------
|
|	LUDecompose
|
|	Replace matrix 'a' with its LU-decomposition.  Returns ERROR if matrix is singular, NO_ERROR
|	otherwise.
*/

	/*   **a = the matrix whose LU-decomposition is wanted                    */
	/*     n = order of a                                                     */
	/*   *vv = work vector of size n (stores implicit scaling of each row)    */
	/* *indx => row permutation according to partial pivoting sequence        */
	/*   *pd => 1 if number of row interchanges was even, -1 if odd (NULL OK) */


