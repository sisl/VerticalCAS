CC=gcc
CFLAGS=-I /opt/OpenBLAS/include/ 
LDFLAGS= -L/opt/OpenBLAS/lib -lopenblas -lpthread -lm
DEPS=matrix.h nnet_VC.h split_VC.h

all: network_test_VC 
all: CFLAGS += -O3
all: CFLAGS += -std=c99
all: LDFLAGS += -O3

debug: network_test_VC
debug: CFLAGS += -DDEBUG -g
debug: CFLAGS += -std=c99
debug: LDFLAGS += -DDEBUG -g

network_test_VC: matrix.o nnet_VC.o network_test_VC.o split_VC.o
	$(CC) -o $@ $^ $(LDFLAGS)

c.o: 
	$(CC) $(CFLAGS) $<  -o $@

clean:
	rm -f *.o network_test_VC

