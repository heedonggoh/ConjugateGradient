## Generic Makefile for PETSc by Heedong Goh <wellposed@gmail.com> ##
## Name of the executable
NAME = run
## Source extension
SEXT = c
## Source path
SDIR = ./src
## Output path
ODIR = ./bin
## PETSc path
PETSC_DIR = /terra3/heedong/software/petsc-3.7.5
PETSC_ARCH = arch-icc-opt
## Additional compile option
CFLAGS = -Wall -Wextra -pedantic -std=c11
## DO NOT CHANGE ###################################################
TARGET = $(ODIR)/$(NAME)
.PHONY: all test
all: $(TARGET)

SRCS = $(wildcard $(SDIR)/*.$(SEXT))
OBJECTS = $(SRCS:.$(SEXT)=.o)
CLEANFILES = $(TARGET)

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

$(TARGET): $(OBJECTS) chkopts
	@mkdir -p $(ODIR)
	@$(CLINKER) -o $(TARGET) $(OBJECTS) $(PETSC_LIB)
	@$(RM) $(OBJECTS)
test:
	@echo "Source files:" $(SRCS)
	@echo "Output file:" $(TARGET)
	@echo "Compile options:" $(CLINKER)
	@echo "PETSc libraries:" $(PETSC_LIB)
####################################################################
