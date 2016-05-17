MPICC    = mpiicpc
COMPILE  = mpiicc #mpiicc #mpicxx

MKL_MIC_ENABLE=1

MKL =    ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm


FLAGS= -DMKL_ILP64 -qopenmp -std=c++11 -I${MKLROOT}/include 

all: mycg

mycg: main.o mmio.o 
		$(COMPILE) $(FLAGS) main.o mmio.o -o mycg $(MKL)
			
main.o:
	$(COMPILE) $(FLAGS) -c main.cpp  

mmio.o:
	$(COMPILE) $(FLAGS) -c mmio.cpp  -Wno-write-strings

clean:
		rm -f mycg
		rm -f *.o
