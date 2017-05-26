#include "conjugateGradient.h"

operator cgAx;
ecode ierr;
integer size, maxiter, iter, cflag;
scalar rtr, nrtr, ptAp, alpha, beta;
Vec r,p,d;
real errtol,err;

ecode ConjugateGradientInitialize(integer _size)
{
  size = _size;
  maxiter = _size;
  errtol = 1.0e-10;
  ierr = VecCreate(PETSC_COMM_WORLD,&r);   CHKERRQ(ierr);
  ierr = VecSetSizes(r,PETSC_DECIDE,size); CHKERRQ(ierr);
  ierr = VecSetFromOptions(r);             CHKERRQ(ierr);
  ierr = VecDuplicate(r,&p);               CHKERRQ(ierr);
  ierr = VecDuplicate(r,&d);               CHKERRQ(ierr);
  return 0;
}

ecode ConjugateGradientFinalize()
{
  ierr = VecDestroy(&r); CHKERRQ(ierr);
  ierr = VecDestroy(&p); CHKERRQ(ierr);
  ierr = VecDestroy(&d); CHKERRQ(ierr);
  return 0;
}

ecode ConjugateGradientSolve(Vec x, Vec b, void* paramp)
{
  scalar refnorm;
  cflag = 0;
  ierr = VecCopy(b,r);               CHKERRQ(ierr);
  ierr = cgAx(d,x,paramp);            CHKERRQ(ierr);
  ierr = VecAXPY(r,-1.0,d);          CHKERRQ(ierr);
  ierr = VecCopy(r,p);               CHKERRQ(ierr);
  ierr = VecDot(r,r,&rtr);           CHKERRQ(ierr);
  ierr = VecNorm(b,NORM_2,&refnorm); CHKERRQ(ierr);
  for(iter=0;iter<maxiter;++iter){
    ierr = cgAx(d,p,paramp);          CHKERRQ(ierr);
    ierr = VecDot(p,d,&ptAp);        CHKERRQ(ierr);
    if(ptAp<0){ cflag = 2; break; }
    alpha = rtr/ptAp;
    ierr = VecAXPY(x,alpha,p);       CHKERRQ(ierr);
    ierr = VecAXPY(r,-alpha,d);      CHKERRQ(ierr);
    ierr = VecDot(r,r,&nrtr);        CHKERRQ(ierr);
    ierr = VecNorm(r,NORM_2,&err);   CHKERRQ(ierr);
    err = err/refnorm;
    if(err<errtol){ cflag = 1; break; }
    beta = nrtr/rtr;
    ierr= VecAYPX(p,beta,r);         CHKERRQ(ierr);
    rtr = nrtr;
  }
  return 0;
}

ecode CojugateGradientGetConvInfo(integer* _cflag, integer* _iter, real* _err)
{
  *_cflag = cflag;
  *_iter = iter;
  *_err = err;
  return 0;
}

ecode ConjugateGradientSetOperator(operator _cgAx)
{
  cgAx = _cgAx;
  return 0;
}

ecode CojugateGradientSetTolerance(real _errtol)
{
  errtol = _errtol;
  return 0;
}

ecode CojugateGradientSetMaxIter(real _maxiter)
{
  maxiter = _maxiter;
  return 0;
}
