CCFLAGS = -g -Wall -std=c++11
CC = g++
Pthread = -lpthread
OMP = -fopenmp

all: pthread omp
	
pthread: LUpthreads.cpp
	mkdir -p ./bin
	$(CC) $(CCFLAGS) -o ./bin/LUpthreads LUpthreads.cpp $(Pthread)

omp: LUomp.cpp
	mkdir -p ./bin
	$(CC) $(CCFLAGS) -o ./bin/LUomp LUomp.cpp $(OMP)

clean:
	rm -rf ./bin
	