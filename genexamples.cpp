#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>  // int64_t

#define DIM 32
#define VERT 64

int nvert, dimension;
int vertices[VERT][DIM];

// Generate vertices of {-1,1}-qube
void gen_cube_vertices(char *name, int d){
    dimension = d;
    nvert = (1 << dimension);
	sprintf (name, "%dd%dvc.g", dimension, nvert);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++, x >>= 1)
            vertices[i][j] = 2 * (x&1) - 1;
	}	
}

// Generate all vectors with coordinates from {-2,-1,1,2}
void gen_12vectors(char *name, int d){
    dimension = d;
    nvert = 1 << (2*dimension);
	sprintf (name, "%dd%dv.g", dimension, nvert);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++, x >>= 2)
            vertices[i][j] = ((x&2) - 1) * (1 + (x&1));
	}	
}

// Generate vertices of the simplex with multiplicity mult
void gen_simplex_vertices(char *name, int d, int mult){
    dimension = d;
    nvert = (d+1) * mult;
	sprintf (name, "%dd%dvs.g", dimension, nvert);
	int i, j, v, m;
	for (i = 0, v = 0; i < dimension; i++)
        for (m = 0; m < mult; m++, v++){
            for (j = 0; j < dimension; j++)
                vertices[v][j] = 0;
            vertices[v][i] = 1;
        }
    for (m = 0; m < mult; m++, v++){
        for (j = 0; j < dimension; j++)
            vertices[v][j] = -1;
    }
}

// Gen vertices of the cross polytope with multiplicity mult
void gen_cross_vertices(char *name, int d, int mult){
    dimension = d;
    nvert = 2 * d * mult;
	sprintf (name, "%dd%dvcr.g", dimension, nvert);
	int i, j, v, m, s;
	for (i = 0, v = 0; i < dimension; i++){
        for (s = 1; s > -2; s -= 2){
            for (m = 0; m < mult; m++, v++){
                for (j = 0; j < dimension; j++)
                    vertices[v][j] = 0;
                vertices[v][i] = s;
            }
        }    
    }    
}

// Gen middles of edges of the cross polytope
void gen_ortahedr_edges(char *name, int d){
    dimension = d;
    nvert = 2 * d * (d-1);
	sprintf (name, "%dd%dvo.g", dimension, nvert);
	int i, j;
	for (i = 0; i < nvert; i++)
		for (j = 0; j < dimension; j++)
            vertices[i][j] = 0;
    int v = 0;
	for (i = 0; i < dimension-1; i++){
		for (j = i+1; j < dimension; j++){
            vertices[v][i] = 1;
            vertices[v][j] = 1;
            v++;
            vertices[v][i] = -1;
            vertices[v][j] = 1;
            v++;
            vertices[v][i] = 1;
            vertices[v][j] = -1;
            v++;
            vertices[v][i] = -1;
            vertices[v][j] = -1;
            v++;
        }    
	}	
}

int write_gale(char *fname){
    FILE *outf = fopen(fname, "w");
    if (outf == NULL){
        printf ("Cann't open file %s\n", fname);
        return 1;
    }    
    printf ("Open file %s for writing\n", fname);
	fprintf (outf, "%d %d\n", dimension, nvert);
	int i, j;
	for (i = 0; i < nvert; i++){
		for (j = 0; j < dimension; j++){
			fprintf (outf, " %2d", vertices[i][j]);
		}
		fprintf (outf, "\n");
	}
    fclose (outf);
    return 0;
}

	
//////////////
//
//    MAIN
//
/////////////


// вход программы - имя файла с описанием диаграммы Гейла
// выход программы - результаты тестирования

int main(int argc, char *argv[])
{
	char outfname[64];
    gen_12vectors(outfname, 2);
    write_gale(outfname);
    
    gen_cross_vertices(outfname, 7, 2);
    write_gale(outfname);
    
    gen_cube_vertices(outfname, 5);
    write_gale(outfname);
    
    gen_simplex_vertices(outfname, 10, 2);
    write_gale(outfname);
    
    gen_ortahedr_edges(outfname, 5);
    write_gale(outfname);
	return 0;
}
