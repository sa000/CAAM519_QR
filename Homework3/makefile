
bin/qr: src/qr.o 
	-@${MKDIR} bin
	${CLINKER} -o $@ $^ ${PETSC_LIB}
	${RM} qr.o 

clean::
	${RM} bin/qr
	${RM} -r bin/qr.dSYM


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
