/*
 * comp_graph.cpp
 *
 *  Created on: July 24, 2015
 *      Author: devikakedia
 */

#include "graph.h"
#include "comp_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>




/*
int main()
{
	compress_graph c;
	graph_t mat;
	graph_t d_mat;

	int i,j;

	for(i=0; i<N; i++){
		for(j=0; j<N; j++){

			if(i==j || i>j)
				mat.adj[i][j] = 1;
			else
				mat.adj[i][j] = 0;
			printf("%d ",mat.adj[i][j]);
		}
		printf("\n");

	}
	Compress_graph(&mat, (comp_graph*)&c, 1);
	for(i=0; i<C; i++){
		printf("%d\n", c.comp[i]);
	}

	Compress_graph(&d_mat,(comp_graph*)&c, 0);
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%d ",d_mat.adj[i][j]);
		}
		printf("\n");
	}

	return 0;

}
*/




void
Compress_graph(graph_t* p, comp_graph *Compress, int is_compress)
{	

	int offset, n = 0;
	char flag[8] = {1,2,4,8,16,32,64,128};

	if(is_compress)
		memset(Compress->comp, 0, C);

		for(int i=0; i<N; i++){
		   for(int j=0; j<N; j++){
			offset = n/8;
			if(is_compress){
				if(p->adj[i][j])
					Compress->comp[offset]+=flag[n%8];
			} else {
				if(Compress->comp[offset] & flag[n%8])
					p->adj[i][j]=1;
				else p->adj[i][j]=0;
			}
			n++;
		  }
		}
	// printf("offset, n: %d %d\n", offset, n);
}
