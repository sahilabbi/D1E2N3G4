#ifndef __DISTANCE__

#define __DISTANCE__

#include <stdlib.h>

#define GRAPH_INFINITY ((int)1000000)

typedef struct {
  int n;
  int *distances;
  int sum_of_distances;
  int diameter;
  int *k;
  int max_k;
} graph_info;

graph_info * getinfo(int * g, int n);
int * calc_dists(int * g, int n);
int sum_dists(int * d, int n);
int get_diameter(int * d, int n);
int * get_degrees(int * g, int n);
int get_max(int *list, int size);
float get_average(graph_info *info);

void copy_graph_info(graph_info *dst, graph_info * src);
void delete_graph_info(graph_info *g);

#endif
