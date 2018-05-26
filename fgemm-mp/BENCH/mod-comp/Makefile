SRCS = 
OBJS = $(patsubst %.cpp, %.o, $(SRCS))

GXX = g++
GXXFLAGS = -Wall -O0 -g -march=native -D__FFLASFFPACK_HAVE_CBLAS -msse4.1 -fabi-version=6 -std=gnu++11 -fopenmp
LIBS = NTL-interface.o -lopenblas -lgivaro -lflint -lntl -lgmpxx -lgmp

.PHONY = clean execs

all:  $(OBJS) execs

clean:
	rm -f $(OBJS)
	rm -f bench
.cpp.o:
	$(GXX) -c $< -o $@ $(GXXFLAGS)

execs:  $(OBJS)
	$(GXX) bench.cpp -o bench $(OBJS) $(GXXFLAGS) $(LIBS)

