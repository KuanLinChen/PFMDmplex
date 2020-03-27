
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

OBJ =  geometry_structure.o linear_solvers_structure.o solver_poisson.o variable_structure.o solver_sepcies.o numerics_convection.o

main: main.o  ${OBJ}
	-${CLINKER} -o main main.o ${OBJ} ${PETSC_LIB}
	${RM} -f *.o 


