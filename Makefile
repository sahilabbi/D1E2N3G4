MPICC ?= mpicc 
CFLAGS += -Wall -O3 -I .
NAME = exhuastive_search
OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))

N ?= 64
K ?= 6

CFLAGS += -DN=$(N) -DK=$(K) -g

all: $(NAME)

clean:
	rm -f $(NAME) distance_test swap_test $(OBJECTS)

$(NAME): $(OBJECTS)
	$(MPICC) -g -o $@ $^ -lm

$(OBJECTS): %.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

swap_test: node_swap.c distance.c priority_queue.c
	$(MPICC) $(CFLAGS) -DSWAP_TEST -o $@ $^

.PHONY: all clean
