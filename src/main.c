#include <unistd.h>
#include <sys/types.h>
#include "macros.h"
#include "conjugateGradient.h"

#define charSize 30

static char help[] = 
"\n--------------------------------------------------------------------------\n\
 Conjugate gradient solver by Heedong Goh <wellposed@gmail.com> \n\
--------------------------------------------------------------------------\n";

ecode ierr;
boolean flg;
char title[charSize];
integer cpuSize, cpuRank;

ecode GetProcInfo();
ecode A(Vec Ax, Vec x, integer nargs, ...);
ecode vA(Vec Ax, Vec x, integer nvargs, va_list vargs);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  integer size = 100;
  integer cflag, iter;
  Mat M;
  Vec x, b, refx;
  PetscRandom rctx;
  real norm,refnorm,err;
  integer Istart,Iend,Ii,loc1,loc2;
  scalar v;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);             CHKERRQ(ierr);
  ierr = GetProcInfo();                                          CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-size",&size,NULL);       CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);   CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,size); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);             CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);               CHKERRQ(ierr);
  ierr = VecDuplicate(x,&refx);            CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&M);                     CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,size,size); CHKERRQ(ierr);
  ierr = MatSetType(M,MATAIJ);                               CHKERRQ(ierr);
  ierr = MatSetUp(M);                                        CHKERRQ(ierr);
  
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);      CHKERRQ(ierr);
  ierr = VecSetRandom(refx,rctx);                        CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);          CHKERRQ(ierr);
  for (Ii=Istart; Ii<Iend; Ii++) {
    loc1 = Ii; loc2 = Ii+1;
    v = 2.0; ierr = MatSetValues(M,1,&loc1,1,&loc1,&v,INSERT_VALUES); CHKERRQ(ierr);
    v = -1.0; 
    if (loc2<size){
      ierr = MatSetValues(M,1,&loc1,1,&loc2,&v,INSERT_VALUES);        CHKERRQ(ierr);
      ierr = MatSetValues(M,1,&loc2,1,&loc1,&v,INSERT_VALUES);        CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
  ierr = A(b,refx,1,M);                          CHKERRQ(ierr);

  ierr = ConjugateGradientInitialize(size);              CHKERRQ(ierr);
  ierr = ConjugateGradientSetOperator(&vA);              CHKERRQ(ierr);
  ierr = ConjugateGradientSolve(x,b,1,M);                CHKERRQ(ierr);
  ierr = CojugateGradientGetConvInfo(&cflag,&iter,&err); CHKERRQ(ierr);
  ierr = ConjugateGradientFinalize();                    CHKERRQ(ierr);

  ierr = VecNorm(refx,NORM_2,&refnorm);                                    CHKERRQ(ierr);
  ierr = VecAXPY(x,-1.0,refx);                                             CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);                                          CHKERRQ(ierr);
  ierr = wprintf("matrix size = %d\n",size);                               CHKERRQ(ierr);
  ierr = wprintf("cflag = %d (== 0 converged; == 1 not converged; == 2 not Hermitian positive definite)\n",cflag);
  ierr = wprintf("iteration count = %d\ncg error = %e (||r||/||b||)\n",iter,err); CHKERRQ(ierr);
  ierr = wprintf("x error = %e (||x-xr||/||xr||)\n",norm/refnorm);                           CHKERRQ(ierr);
  ierr = wprintf("\n" BOLD "End %s" RESET "\n",title);                     CHKERRQ(ierr);

  ierr = VecDestroy(&refx);         CHKERRQ(ierr);
  ierr = VecDestroy(&x);            CHKERRQ(ierr);
  ierr = VecDestroy(&b);            CHKERRQ(ierr);
  ierr = MatDestroy(&M);            CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
  ierr = PetscFinalize();           CHKERRQ(ierr);
  return 0;
}

ecode GetProcInfo()
{
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&cpuSize);                           CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&cpuRank);                           CHKERRQ(ierr);
  ierr = wprintf("\n" BLUE "cpu size = %ld" RESET "\n",cpuSize);             CHKERRQ(ierr);
  ierr = syprintf(BLUE "cpu %d:\tpid = %ld, ppid = %ld" RESET "\n",
		  cpuRank,(long)getpid(),(long)getppid());                   CHKERRQ(ierr); 
  ierr = sypflush();                                                         CHKERRQ(ierr); 
  ierr = PetscOptionsGetString(NULL,NULL,"-title",title,sizeof(title),&flg); CHKERRQ(ierr); 
  if(!flg) strcat(title,"noname");
  ierr = wprintf("\n" BOLD "title:\t%s" RESET "\n",title);                   CHKERRQ(ierr);
  return 0;
}

ecode A(Vec Ax, Vec x, integer nvargs, ...)
{
  va_list vargs;
  va_start(vargs,nvargs);
  ierr = vA(Ax,x,nvargs,vargs); CHKERRQ(ierr);
  va_end(vargs);
  return 0;
}

ecode vA(Vec Ax, Vec x, integer nvargs, va_list vargs)
{
  Mat M = va_arg(vargs,Mat);
  return MatMult(M,x,Ax);
}
