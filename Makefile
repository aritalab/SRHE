PREFIX=$(HOME)
CC=g++
AR=ar
# CFLAGS= -Wall -O3 
CFLAGS= -std=c++11 -g -Wall -O2 -pthread
LDFLAGS= -L. -lsbr2 -lntl -lgmp #-lrt

LOGSHE = ./log/log-she
PARAMS = ./params/params

HEADER = params.h nsgen.h ring.h she.h util.h

OBJ = nsgen.o ring_init.o ring.o she.o util.o

TESTPROGS = nsgenTest ringTest sheTest

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) -c $<

all: $(OBJ) libsbr2.a $(TESTPROGS)

clean:
	rm -f $(OBJ) *~ libsbr2.a $(TESTPROGS)

libsbr2.a: $(OBJ)
	$(AR) -rv libsbr2.a $(OBJ)

nsgenTest: nsgenTest.cpp libsbr2.a
	$(CC) $(CFLAGS) -o nsgenTest nsgenTest.cpp $(LDFLAGS)

ringTest: ringTest.cpp libsbr2.a
	$(CC) $(CFLAGS) -o ringTest ringTest.cpp $(LDFLAGS)

sheTest: sheTest.cpp libsbr2.a
	$(CC) $(CFLAGS) -o sheTest sheTest.cpp $(LDFLAGS)


