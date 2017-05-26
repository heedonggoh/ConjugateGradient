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
ecode A(Vec Ax, Vec x, void* paramp);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{

  integer size = 500;
  integer cflag, iter;
  Mat M;
  Vec x, b, refx;
  PetscRandom rctx;
  real norm,refnorm,err;
  integer Istart,Iend,Ii;
  scalar v;

  ierr = PetscInitialize(&argc,&args,(char*)0,help); CHKERRQ(ierr);
  ierr = GetProcInfo();                              CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);   CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,size); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);             CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);               CHKERRQ(ierr);
  ierr = VecDuplicate(x,&refx);            CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&M);                     CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,size,size); CHKERRQ(ierr);
  ierr = MatSetType(M,MATAIJ);                               CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(M,size,NULL,size,NULL);   CHKERRQ(ierr);
  
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx); CHKERRQ(ierr);
  ierr = VecSetRandom(refx,rctx);                   CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);     CHKERRQ(ierr);
  for (Ii=Istart; Ii<Iend; Ii++) {
    v = 4.0*(1+Ii); ierr = MatSetValues(M,1,&Ii,1,&Ii,&v,INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
  ierr = A(b,refx,M);                            CHKERRQ(ierr);
  
  ierr = ConjugateGradientInitialize(size);              CHKERRQ(ierr);
  ierr = ConjugateGradientSetOperator(&A);               CHKERRQ(ierr);
  ierr = ConjugateGradientSolve(x,b,M);                  CHKERRQ(ierr);
  ierr = CojugateGradientGetConvInfo(&cflag,&iter,&err); CHKERRQ(ierr);
  ierr = ConjugateGradientFinalize();                    CHKERRQ(ierr);

  ierr = VecNorm(refx,NORM_2,&refnorm);                                    CHKERRQ(ierr);
  ierr = VecAXPY(x,-1.0,refx);                                             CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);                                          CHKERRQ(ierr);
  ierr = wprintf("cflag = %d\titer = %d\tcg error = %e\t",cflag,iter,err); CHKERRQ(ierr);
  ierr = wprintf("x error = %e\n",norm/refnorm);                           CHKERRQ(ierr);

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

ecode A(Vec Ax, Vec x, void* paramp)
{
  Mat M = paramp;
  return MatMult(M,x,Ax);
}
