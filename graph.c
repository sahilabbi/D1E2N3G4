#include "graph.h"
#include <stdlib.h>
#include <string.h>

void copy_graph(graph_t * dst, graph_t * src){
    int i,j;
    for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			dst->adj[i][j] = src->adj[i][j];
		}
    }
}

//gives whether a node is one of the initial edges.
//these edges are not changed to avoid isomorphisms,
//as any graph obtained by altering them would be isomporphic
//to some other graph  obtained without altering them
static bool is_protected(int row, int column){
    if(row == 0) return true;
    if(row == column) return true;
    if(row - column != 1 && column - row != 1) return false;
    int x = row > column ? row : column;
    return x % 2 == K % 2;
}

//moves the cursor right to left on the row
//if the middle of the graph has been reached,
//then the cursor jumps up the the right edge of
//the row above
static void increment_cursor(int * i, int * j){
    *j -= 1;
    if(is_protected(*i,*j)){
	*i -= 1;
	*j = N-1;
	if(i == 0) return;
    }
}

bool is_regular(int * degree){
	int i;
	for(i = 0; i < N; i++){
		if(degree[i] != K) return false;
	}
	return true;
}

#include <stdio.h>

bool increment_graph(graph_t * graph){
    int degree[N];
    int i,j;
    for(i = 0; i < N; i++){
	int sum = 0;
	for(j = 0; j < N; j++){
	    sum += graph->adj[i][j];
	}
	degree[i] = sum;
    }
    i = N-2;
    j = N-1;
    //incase the first edge selected is one of the initial edges
    if(is_protected(i,j)) increment_cursor(&i,&j);
    //while(!is_regular(&degree[0])){
    while(graph->adj[i][j] == 1){
	graph->adj[i][j] = 0;
	graph->adj[j][i] = 0;
	degree[i]--;
	degree[j]--;
	increment_cursor(&i,&j);
	if(i == 0) return true;
    }
    graph->adj[i][j] = 1;
    graph->adj[j][i] = 1;
    degree[i]++;
    degree[j]++;
    printf("i: %d\nj: %d\n", i, j);
    //}
    return false;
}

graph_t * initial_graph(){
	graph_t * ret = (graph_t *) malloc(sizeof(graph_t));
	memset(&ret->adj[0][0], 0, sizeof(ret->adj));
	int i;
	for(i = 1; i <= K; i++){
		ret->adj[0][i] = 1;
		ret->adj[i][0] = 1;
	}

	for(i = K+1; i < N; i += 2){
		ret->adj[i][i+1] = 1;
		ret->adj[i+1][i] = 1;
	}
	
	return ret;
}

#include <stdio.h>
void print_graph(graph_t * g){
	int i,j;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			printf("%d ",g->adj[i][j]);
		}
		printf("\n");
	}
}

#ifdef INCREMENT_TEST

int main(){
	graph_t * g = initial_graph();
	print_graph(g);
	printf("\n");
	int i;
	for(i = 0; i < 5; i++){
		increment_graph(g);
		print_graph(g);
		printf("\n");
	}

	printf("is_protected(%d,%d)=%d\n",K+1,K+2,is_protected(K+1,K+2));
	return 0;
}


#endif
