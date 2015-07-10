#include "distance.h"
#include <string.h>
#include <stdio.h>

static int add_neighbors(int * list, int * g, int n, int * visited, int v); 

graph_info * getinfo(int * g, int n){
	graph_info *ret = malloc(sizeof(graph_info));
	ret->n = n;
	ret->distances = calc_dists(g, n);
	ret->sum_of_distances = sum_dists(ret->distances, n);
	ret->diameter = get_diameter(ret->distances, n);
	ret->k = get_degrees(g, n);
	ret->max_k = get_max(ret->k, n);
	return ret;
}

int * calc_dists(int *g, int n){
	int *ret = malloc(n * n * sizeof(int));
	int *current_nodes = malloc(n * sizeof(int));
	int current_nodes_size, next_nodes_size;
	int * next_nodes;
	int * visited = malloc(n * sizeof(int));
	int current_distance = 0;
	int i, j;
	int r = 0;
	for(i = 0;i < n; i++) for(j = 0; j < n; j++) ret[i * n + j] = GRAPH_INFINITY;
	for(;r < n; r++){
		for(i = 0; i < n; i++){
		   	current_nodes[i] = 0;
			visited[i] = 0;
		}
		current_distance = 0;
		current_nodes[0] = r;
		visited[r] = 1;
		current_nodes_size = 1;
		while(current_nodes_size > 0){
			for(i = 0; i < current_nodes_size; i++) ret[r*n + current_nodes[i]] = current_distance;
			next_nodes = malloc(n * sizeof(int));
			next_nodes_size = 0;
			for(i = 0; i < current_nodes_size; i++)
				next_nodes_size += add_neighbors(next_nodes + next_nodes_size, 
												  g, n, visited, 
												  current_nodes[i]);
			free(current_nodes);
			current_nodes = next_nodes;
			current_nodes_size = next_nodes_size;
			current_distance++;
		}
		
	}
	free(current_nodes);
	free(visited);
	return ret;
}

//adds the unvisited neighbors of node v in graph g to list
//returns the number of neighbors added
static int add_neighbors(int *list, int * g, int n, int * visited, int v){
	int i = 0;
	int k = 0;
	for(;i < n; i++){
		if(g[v*n + i] && !visited[i]){
			list[k++] = i;
			visited[i] = 1;
		}
	}
	return k;
}


int sum_dists(int *d, int n){
	int sum = 0;
	int i = 0;
	for(;i < n*n; i++) sum += d[i];
	return sum;
}

float get_average(graph_info *info){
	return (float) info->sum_of_distances /
		   (float) info->n /
		   (float)(info->n - 1);
}

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b));
#endif

int get_diameter(int * d, int n){
	int ret = 0;
	int i = 0, j = 0;
	for(; i < n; i++) for(j = 0;j < n; j++) ret = MAX(d[i*n + j], ret);
	return ret;
}

int *get_degrees(int *g, int n){
	int *ret = malloc(n * sizeof(int));
	int i,j,sum;
	for(i = 0; i < n; i++){
		sum = 0;
		for(j = 0; j < n; j++) sum += g[i*n + j];
		ret[i] = sum;
	}
	return ret;
}

int get_max(int * list, int size){
	int i = 0;
	int ret = 0;
	for(;i < size; i++) ret = MAX(ret, list[i]);
	return ret;
}

void delete_graph_info(graph_info *g){
	free(g->distances);
	free(g->k);
	free(g);
}

void copy_graph_info(graph_info *dst, graph_info *src){
	dst->n = src->n;
	memcpy(dst->distances, src->distances, src->n * src->n * sizeof(int));
	dst->sum_of_distances = src->sum_of_distances;
	dst->diameter = src->diameter;
	memcpy(dst->k, src->k, src->n * sizeof(int));
	dst->max_k = src->max_k;
}


/*
int main(){
	int * g = hypercube(3);
	graph_info *info = getinfo(g, 8);
	printf("Diameter: %d\nAverage: %f", info->diameter, (float)info->sum_of_distances / 56.0f);
	return 0;
}*/
