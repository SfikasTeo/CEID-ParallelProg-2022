UNAME_S := $(shell uname -s)
CFLAGS = -O3
#specify weather aggresive optimization should be enabled

# OK... this might be incorrect, and we understand that putting more and more compiler options wont make the code magically faster.
# Having said that in my understanding using -march=native first of all makes the compiled code almost bound to my CPU but also enables
# -mavx2 -mfma , may also make changes to caching etc so there is no need to include those into the CFLAGS. 
# In our case -flto doen not seem to not hinder the performance from our (small) understanding.
# Lastly -ffast-math is a debatalbe flag risks explained at: https://stackoverflow.com/questions/7420665/what-does-gccs-ffast-math-actually-do.
# We opted to include it in the most aggresive of the 3 optimization scales because the performance boost is UNCANNY. 
# March=native also boosts the performance by a significant amount but the executables become architecturally bound. Lastly -flto for most cases is a nice addition.

aggresive_optimization=1
ifeq ($(UNAME_S),Linux)
	CC=gcc-11
	ifeq ($(aggresive_optimization),2)
		CFLAGS=-O3 -mfma -flto
	endif
	ifeq ($(aggresive_optimization),3)
		CFLAGS=-O3 -flto -ffast-math -march=native
	endif
endif
ifeq ($(UNAME_S),Darwin)
	CC=gcc-8
endif

#Set mpicc to mpicc if not already setted.
MPICC := mpicc
CFLAGS_THREADS = $(CFLAGS) -fopenmp
LIBS = -lm

#Run Mpi+openmp by mpiexec -n # executable #
all: multistart_serial multistart_omp multistart_hybrid multistart_mpi

multistart_serial:
	@$(CC) -O3 -o mult_hooke_seq multistart_hooke_seq.c $(LIBS)
	
multistart_omp:
	@$(CC) $(CFLAGS_THREADS) -o mult_hooke_omp multistart_hooke_omp.c $(LIBS)
	@$(CC) $(CFLAGS_THREADS) -o mult_hooke_tasks multistart_hooke_tasks.c $(LIBS)
	
multistart_mpi:
	@$(MPICC) $(CFLAGS) -o mult_hooke_mpi multistart_hooke_mpi.c $(LIBS)

multistart_hybrid:
	@$(MPICC) $(CFLAGS_THREADS) -o mult_hooke_mpi_omp multistart_hooke_mpi_omp.c $(LIBS)

cleanall: clean cleantxt

clean:
	@echo "Cleaning the executables..."
	@rm -f mult_hooke_seq
	@rm -f mult_hooke_omp
	@rm -f mult_hooke_tasks
	@rm -f mult_hooke_tasks2
	@rm -f mult_hooke_mpi
	@rm -f mult_hooke_mpi2
	@rm -f mult_hooke_mpi_omp
	@rm -f mult_hooke_mpi_omp2

cleantxt:
	@echo "Cleaning the outputs..."
	@rm -f *.txt

	
