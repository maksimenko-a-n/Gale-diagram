#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>  // uint32_t
#include <string.h> // memcmp()
#include <math.h>
#include <time.h>

typedef uint16_t Tpoint;
#define MAX_Tpoint 1023

#define DIM 5 // The maximum dimension of the Gale diagram (GD)
#define MAX_VERT 20 // The maximum number of vertices (cann't be greater than 32)
#define MAX_FACET 20 // For the acceleration of the processing (cann't be greater than 32)

const int Nall = pow(3, DIM); // The number of all points that can be added

/*
#ifndef GALE_PARAMS_H
#define GALE_PARAMS_H

// Extract the dimension and the number of vertices from the filename
int get_dv(char *filename, int &nv){
    char *pch = filename, *newpch;
    //pch = strchr(filename,'s');
    int dim = strtol(pch, &newpch, 10);
    if (pch == newpch) return 1;
	if (dim != DIM) return 2;
	if (*newpch != 'd') return 3;
    pch = newpch + 1;
    nv = strtol(pch, &newpch, 10);
    if (pch == newpch) return 4;
	if (nv < 1 || nv > MAX_VERT) return 5;
    return 0;
}
#endif // GALE_PARAMS_H 
*/