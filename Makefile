CC = mpicc
CXX = mpicxx
LDFLAGS = -lpthread
CFLAGS = -O3 -march=native -std=gnu99
CXXFLAGS = -O3 -march=native -std=gnu++11
TARGETS = APSP_Pthread APSP_MPI_sync APSP_MPI_async APSP_Hybrid

.PHONY: all
all: $(TARGETS)

.PHONY: all
clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
