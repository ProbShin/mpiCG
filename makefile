MPICC    = mpiicpc
COMPILE  = mpiicc #mpiicc #mpicxx

MKL_MIC_ENABLE=1

MKL =    ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm


FLAGS= -DMKL_ILP64 -qopenmp -std=c++11 -I${MKLROOT}/include 

all: mycg  autorecv

mycg: basic_cg.o mmio.o 
		$(COMPILE) $(FLAGS) basic_cg.o mmio.o -o mycg $(MKL)
			
basic_cg.o:
	$(COMPILE) $(FLAGS) -c basic_cg.cpp   

mmio.o:
	$(COMPILE) $(FLAGS) -c mmio.cpp  -Wno-write-strings

autorecv: kkproj.o mmio.o
		$(COMPILE) $(FLAGS) kkproj.o mmio.o -o autorecv $(MKL)

kkproj.o:
	$(COMPILE) $(FLAGS) -c kkproj.cpp  

clean:
		rm -f mycg
		rm -f autorecv
		rm -f *.o
