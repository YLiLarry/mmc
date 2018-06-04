SRCS = BlasCRT.cpp
OBJS = $(patsubst %.cpp, %.o, $(SRCS))

GXX = g++
GXXFLAGS = -std=c++11 -Wall -O3 -g -march=native -L/usr/local/atlas/lib
LIBS = -lntl -lgmp  -lcblas -latlas -lm

.PHONY = clean 

all:	 $(OBJS) mul

clean:
	rm -f $(OBJS)
	rm -f mul
.cpp.o:	
	$(GXX) -c $< -o $@ $(GXXFLAGS)

mul:	$(OBJS) mul.cpp
	$(GXX) -o mul mul.cpp $(OBJS) $(GXXFLAGS) $(LIBS)
