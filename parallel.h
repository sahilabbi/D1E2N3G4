#ifndef __PARALLEL__

#define __PARALLEL__

//The generator node generates the graphs
//and sends them to other nodes to be checked.
//The generator sends out every regular graph,
//possibly repeats since isomorphism checking is hard
#define GENERATOR_RANK 0
//the manager node manages which threads are doing what 
//and also tells them where to send their data along
#define MANAGER_RANK 1

//each of these will perform the actions of their respective node(s)
void generator();
void manager();
void distance_calculator();
void iso_checker();

#endif
