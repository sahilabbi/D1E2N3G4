#ifndef __COMP_GRAPH__

#define __COMP_GRAPH__

#include "graph.h"

#define COMPRESS_SIZE ((int)(NUM_NODES * NUM_NODES / sizeof(char)))

// graphs kept in memory or saved to disk in comp_graph format. since all will have same diameter and avg_distance, no need to save per graph 
// this structure also used to send graph pairs across MPI channels for isomorphism checks
typedef struct
{
	unsigned char comp[COMPRESS_SIZE];
} comp_graph;


// this structure used to send graphs across MPI channel for distance calculation
typedef struct
{
	unsigned char comp[COMPRESS_SIZE];
	unsigned char	diameter;
        unsigned int 	sum_dist;
} compress_graph;

void Compress_graph(graph_t *p, comp_graph *Compress, int is_compress);

#define MAX_GRAPH_ARR		10000

#define COMP_DIST_TAG 	0 
#define CHECK_ISO_TAG	0

#define FPATH   "/tmp/fgraph/"

#define	_min(a,b)  ((a) < (b) ? (a) : (b))


#endif
