#include "graph.h"

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
    int x,y;
    if(row < column) {
	x = row;
	y = column;
    }
    else{
	x = column;
	y = row;
    }
    return x - y == 1 && (x - K) % 2 == 1;
}

//moves the cursor right to left on the row
//if the middle of the graph has been reached,
//then the cursor jumps up the the right edge of
//the row above
static void increment_cursor(int * i, int * j){
    *j -= 1;
    if(is_protected(*i,*j)){
	i--;
	*j = N-1;
	if(i == 0) return;
    }
}

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
    while(graph->adj[i][j] == 1){
	graph->adj[i][j] = 0;
	graph->adj[j][i] = 0;
	increment_cursor(&i,&j);
	if(i == 0) return true;
    }
    return false;
}
