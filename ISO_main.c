//
//  ISO_main.c
//  
//

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <mpi.h>
// #include <assert.h>

#include "graph.h"
#include "distance.h"
#include "comp_graph.h"


comp_graph	cgraph[MAX_GRAPH_ARR];
comp_graph  	reduced_graphs[MAX_GRAPH_ARR];

static int file_size(FILE*);
static int read_arr_disk(comp_graph* pg, int ndx, char* suffix);
static int write_arr_disk(comp_graph *arr, int ngraphs, char* suffix, int ndx);
static int reduce_compressed_pairs(comp_graph *arr, int maxgraphs, int arr_size, int arr_offset,
					comp_graph* preduced, int* reduced_array_elements, int reduced_dim);

extern int num_ISOprocs;

int
process_ori_file(comp_graph* pgraph, int num_graphs, int* reduced_elements, int reduced_elements_count)
{
	int		nloop=0;
	static int *num_reduced_graphs = NULL;
	static comp_graph *cgr_reduced = NULL;

	if(!num_reduced_graphs){
		num_reduced_graphs = (int*)malloc(sizeof(int) * num_ISOprocs);
	}

	if(!num_reduced_graphs){
		printf("process-ori_file: num-reduced_graphs is NULL\n");
		return -1;
	}
	if(NULL == reduced_elements){
		reduced_elements = num_reduced_graphs; 
		reduced_elements_count=0;
	}

	if(!cgr_reduced)
		cgr_reduced = (comp_graph*)malloc(sizeof(comp_graph) * num_ISOprocs);

	if(!cgr_reduced){
		printf("process-ori_file: cgr_reduced is NULL\n");
		return -1;
	}

	int		total_reduced_graphs=0;

	int 		remaining_graphs, processed_graphs;

	remaining_graphs = num_graphs;
	processed_graphs = 0;

	int qty;
	while(remaining_graphs > 0){
		qty = _min(remaining_graphs, num_ISOprocs);
        
        	int n_reduced = reduce_compressed_pairs(pgraph, qty, num_graphs, processed_graphs,
					cgr_reduced, reduced_elements, reduced_elements_count);

		remaining_graphs -= qty;
		processed_graphs += qty;
		num_reduced_graphs[nloop++] = n_reduced;
		memcpy(pgraph + total_reduced_graphs, cgr_reduced, n_reduced);
		total_reduced_graphs += n_reduced;
	}
	if(total_reduced_graphs <= 1 || total_reduced_graphs == num_graphs)
		return total_reduced_graphs;

	process_ori_file(pgraph, total_reduced_graphs, (int*)num_reduced_graphs, nloop);

	return -1;
}
	

static int file_size(FILE* fp)
{
    struct stat fileStat;
    fstat(fileno(fp), &fileStat);
    return fileStat.st_size;
}
	

void ISOmain()
{
	int file_ndx=0;
    
    	for( ; ; ){
        
        	char suffix[3] = "ori";
        
        	int num_graphs = read_arr_disk(cgraph, file_ndx, suffix);

		if(num_graphs <= 0){
			printf("No more .ori files to read\n");
			break;
		}
		int niso_graphs = process_ori_file(cgraph, num_graphs, (int*)NULL, 0);
		if(niso_graphs >= 0)
				write_arr_disk(cgraph, niso_graphs, "iso", file_ndx);
        	file_ndx++;
    	}
}
    

// assumes that pg points to enough memory to accomodate the file content
static int
read_arr_disk(comp_graph* pg, int ndx, char* suffix)
{
    char    fpath[100];
    sprintf(fpath, "%s%d.%s", FPATH, ndx, suffix);
    FILE* fp = fopen(fpath, "r");

    if(fp){
	int fsiz = file_size(fp);

        fread((char* )pg, sizeof(comp_graph),  fsiz, fp);
        fclose(fp);
        return (fsiz / sizeof(comp_graph));
    }
    return 0;
}
    
static int
write_arr_disk(comp_graph *arr, int ngraphs, char* suffix, int ndx)
{
    char    fpath[100];
    sprintf(fpath, "%s%d.%s", FPATH, ndx, suffix);
    FILE* fp = fopen(fpath, "w");
    if(fp){
        fwrite((char* )arr, sizeof(comp_graph),  ngraphs, fp);
        fclose(fp);
        return 1;
    }
    return 0;
}
    
static int reduce_compressed_pairs(comp_graph *arr, int maxgraphs, int arr_size, int arr_offset,
					comp_graph* preduced, int* reduced_array_elements, int reduced_dim)
{
    MPI_Status status;
    int ierr;
    int num_non_iso=0;

    	if(maxgraphs > num_ISOprocs){
		printf("Error: reduce-pairs: too many ngraphs: %d\n", maxgraphs);
		return 0;
    	}

	int n=0, k=0, max_ndx=0, iso_node=1;
	int exclude_bound=0;

	static comp_graph* pairs = NULL;


	if(!pairs)
		pairs = (comp_graph*)malloc(num_ISOprocs*2*sizeof(comp_graph));
	if(!pairs){
		printf("Error: reduce-compressed_pairs: pairs is NULL\n");
		return 0;
	}
    

    	for(int i = arr_offset; i < arr_size-1  ; i++){
		exclude_bound=0;
		for(int j=0; j < reduced_dim; j++){ 
			exclude_bound += reduced_array_elements[j];
			if(i < exclude_bound){
				break;
			}
		}
        
		n=0;
        	for(k = i+1; k < arr_size; k++){
			if(exclude_bound && k < exclude_bound)	// do not pair elements from same set of reduced graphs
				continue;
			if(i > arr_offset && k > max_ndx)
				break;
            
            		pairs[2*n] = arr[i];
            		pairs[2*n+1] = arr[k];
	    		n++;
	    		if(i == arr_offset)
				max_ndx = k;
	    		if(n >= maxgraphs){
				break;
	    		}
        	}
		
        	ierr = MPI_Send(pairs, 2 * n * sizeof(comp_graph), MPI_CHAR,
                        iso_node++, COMP_DIST_TAG, MPI_COMM_WORLD);
	}

    	for(int j = 1; j < iso_node; j++){
                
                int index;
                ierr = MPI_Recv (&index, 1, MPI_INT,
                                 MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                if(index < 0){
                    //do nothing
                }
                else{
			printf("Recvd from %d\n", status.MPI_SOURCE);
                        preduced[num_non_iso++] = arr[status.MPI_SOURCE-1];
                }
                //if isomorphic, index = -1
                
    	}

    	return num_non_iso;
}
