#ifndef MACROS_H
#define MACROS_H

#include <petsc.h>

#define boolean PetscBool
#define integer PetscInt
#define real PetscReal
#define scalar PetscScalar
#define True PETSC_TRUE
#define False PETSC_FALSE
#define ecode PetscErrorCode

#define BOLD    "\033[1m"
#define BLACK   "\x1b[30m"
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define WHITE   "\x1b[37m"
#define BBLACK  "\x1b[40m"
#define BRED    "\x1b[41m"
#define BGREEN  "\x1b[42m"
#define BYELLOW "\x1b[43m"
#define BBLUE   "\x1b[44m"
#define BWHITE  "\x1b[47m"
#define RESET   "\x1b[0m"

#define wprintf(f_, ...) PetscPrintf(PETSC_COMM_WORLD,(f_),__VA_ARGS__)
#define sprintf(f_, ...) PetscPrintf(PETSC_COMM_SELF,(f_),__VA_ARGS__)
#define syprintf(f_, ...) PetscSynchronizedPrintf(PETSC_COMM_WORLD,(f_),__VA_ARGS__)
#define sypflush(f_) PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT)

#endif
