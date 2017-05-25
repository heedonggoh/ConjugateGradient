#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "macros.h"

/* i/o structure of the operator */
typedef ecode (*operator)(Vec Ax, Vec x, void* paramp);

/* initializes CG solver */
ecode ConjugateGradientInitialize(integer _size);

/* finalizes CG solver */
ecode ConjugateGradientFinalize();

/* sets the operator A */
ecode ConjugateGradientSetOperator(operator _inA);

/* solves Ax=b */
ecode ConjugateGradientSolve(Vec x, Vec b, void* paramp);

/* gets convergence flag, iteration count, relative error
 * cflag == 0 - not converged
 * cflag == 1 - converged
 * cflag == 2 - not Hermitian 
 * iter - number of iterations
 * err - relative error in L2, ||r||/||b|| */
ecode CojugateGradientGetConvInfo(integer* _cflag, integer* _iter, real* _err);

#endif
