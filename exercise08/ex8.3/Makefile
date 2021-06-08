CC = g++
CFLAGS = -Wall -O3 --std=c++17


LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: qmc1d qmc1d_const

qmc1d.o: qmc1d.cpp
	$(CC) -Wall -c $< ${INCS}

qmc1d: qmc1d.o
	$(CC) ${CFLAGS} -o $@ $^ ${LIBS}

qmc1d_const.o: qmc1d_const.cpp
	$(CC) -Wall -c $< ${INCS}

qmc1d_const: qmc1d_const.o
	$(CC) ${CFLAGS} -o $@ $^ ${LIBS}



clear:
	rm *.o *~
