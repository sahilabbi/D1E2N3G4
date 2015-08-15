MPICC ?= mpicc 
CFLAGS += -Wall -O3 -I .
LDFLAGS += -lgsl -lgslcblas -lm
NAME = exhuastive_search
OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))

N ?= 8
K ?= 3
L ?= 3

CFLAGS += -DNUM_NODES=$(N) -DNODE_DEGREE=$(K) -DL=$(L) -g

all: $(NAME)

comp_graph: comp_graph.c 
	gcc $(CFLAGS) -o $(NAME) -DCOMP_GRAPH comp_graph.c

graph_test: graph_test.c  graph.c distance.c
	$(MPICC) $(CFLAGS) -o $@ $^

graph_master: graph_master.c  graph.c distance.c comp_graph.c
	$(MPICC) $(CFLAGS) -o $@ $^

ISO_check: ISO_checker.c ISO_main.c graph.c distance.c comp_graph.c
	$(MPICC) $(CFLAGS) -o $@ $^

increment_test: graph.c graph.h
	gcc $(CFLAGS) -o $(NAME) -DINCREMENT_TEST graph.c $(LDFLAGS)

increment_test1: graph.c graph.h
	mpicc $(CFLAGS) -o $(NAME) -DINCREMENT_TEST graph.c

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
