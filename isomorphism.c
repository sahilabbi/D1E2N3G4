#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.00001

bool float_equals(float a, float b){
    float diff = a - b;
    if(diff < 0) diff = -diff;
    return diff < EPSILON;
}

typedef struct
{
    int index;
    double value;
} indexValuePair;





#ifdef MATRIX_TEST

int main(){
    
    return 0;
}
    
#endif

