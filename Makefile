MPICC ?= mpicc 
CFLAGS += -Wall -I .
LDFLAGS += -lgsl -lgslcblas
NAME = exhuastive_search
OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))

NUM_NODES ?= 8
NODE_DEGREE ?= 3

CFLAGS += -DN=$(N) -DK=$(K) -g

all: $(NAME)

increment_test: graph.c graph.h
	gcc $(CFLAGS) -o $(NAME) -DINCREMENT_TEST graph.c $(LDFLAGS)

matrix_test: isomorphism.c isomorphism.h graph.h
	gcc $(CFLAGS) -o $(NAME) -DMATRIX_TEST isomorphism.c $(LDFLAGS)

clean:
	rm -f $(NAME) distance_test swap_test $(OBJECTS)

$(NAME): $(OBJECTS)
	$(MPICC) -g -o $@ $^ $(LDFLAGS)

$(OBJECTS): %.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

swap_test: node_swap.c distance.c priority_queue.c
	$(MPICC) $(CFLAGS) -DSWAP_TEST -o $@ $^

.PHONY: all clean
