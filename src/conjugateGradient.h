#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "macros.h"

typedef ecode (*operator)(Vec Ax, Vec x, void* paramp);

ecode ConjugateGradientInitialize(integer _size);
ecode ConjugateGradientFinalize();
ecode ConjugateGradientSetOperator(operator _inA);
ecode ConjugateGradientSolve(Vec x, Vec b, void* paramp);

/* cflag == 0 not converged
 * cflag == 1 converged
 * cflag == 2 not Hermitian */
ecode CojugateGradientGetConvInfo(integer* _cflag, integer* _iter, real* _err);

#endif
