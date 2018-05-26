CC=gcc
CFLAGS=-lgmp
MATH=-lm

matrix : simpleMatrix.o matrix.o SargeNum.o MargeNum.o PargeNum.o Testing.o Reconstruction.o
	$(CC) -o matrix simpleMatrix.o matrix.o SargeNum.o MargeNum.o PargeNum.o Testing.o Reconstruction.o $(CFLAGS) $(MATH)

simpleMatrix.o : simpleMatrix.c matrix.h
	$(CC) -c simpleMatrix.c
matrix.o : matrix.c matrix.h
	$(CC) -c matrix.c
SargeNum.o : SargeNum.c matrix.h
	$(CC) -c SargeNum.c
MargeNum.o : MargeNum.c matrix.h
	$(CC) -c MargeNum.c
PargeNum.o : PargeNum.c matrix.h
	$(CC) -c PargeNum.c
Testing.o : Testing.c matrix.h
	$(CC) -c Testing.c
Reconstruction.o : Reconstruction.c matrix.h
	$(CC) -c Reconstruction.c
clean :
	rm matrix simpleMatrix.o matrix.o SargeNum.o MargeNum.o PargeNum.o Testing.o
