MPICC ?= mpicc 
CFLAGS += -Wall -I .
LDFLAGS += -lgsl -lgslcblas
NAME = exhuastive_search
OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))

N ?= 6
K ?= 3

CFLAGS += -DNUM_NODES=$(N) -DNODE_DEGREE=$(K) -g

all: $(NAME)

increment_test: graph.c graph.h
	gcc $(CFLAGS) -o $(NAME) -DINCREMENT_TEST graph.c $(LDFLAGS)

matrix_test: isomorphism.c isomorphism.h graph.h graph.c
	gcc $(CFLAGS) -o $(NAME) -DMATRIX_TEST isomorphism.c graph.c $(LDFLAGS)

iso_test: isomorphism.c isomorphism.h graph.c graph.h
	gcc $(CFLAGS) -o $(NAME) -DISO_TEST isomorphism.c graph.c graph_set.c $(LDFLAGS)

graph_set_test: graph.c graph.h graph_set.c graph_set.h isomorphism.c isomorphism.h
	gcc $(CFLAGS) -o $(NAME) -DGRAPH_SET_TEST graph_set.c graph.c isomorphism.c $(LDFLAGS)

clean:
	rm -f $(NAME) distance_test swap_test $(OBJECTS)

$(NAME): $(OBJECTS)
	$(MPICC) -g -o $@ $^ $(LDFLAGS)

$(OBJECTS): %.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

swap_test: node_swap.c distance.c priority_queue.c
	$(MPICC) $(CFLAGS) -DSWAP_TEST -o $@ $^

.PHONY: all clean
