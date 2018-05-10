#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define DIM 8
#define VERT 64

int nvert, dimension;
int vertices[VERT][DIM];


void gen_qube_vertices(char *name, int d){
    dimension = d;
    nvert = (1 << dimension);
	sprintf (name, "%dd%dvq.g", dimension, nvert);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++, x >>= 1)
            vertices[i][j] = 2 * (x&1) - 1;
	}	
}

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

void write_gale(FILE *outf){
	fprintf (outf, "%d %d\n", dimension, nvert);
	int i, j;
	for (i = 0; i < nvert; i++){
		for (j = 0; j < dimension; j++){
			fprintf (outf, " %2d", vertices[i][j]);
		}
		fprintf (outf, "\n");
	}
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
	//sprintf (outfname, "gen.g");
    for (int d = 2; d <= 6; d++){
        gen_qube_vertices(outfname, d);
        //gen_ortahedr_edges(outfname, d);
        FILE *outf = fopen(outfname, "w");
        if (outf == NULL){
            printf ("ERROR: Cann't open file %s\n", outfname);
            return 1;
        }
        printf ("Open file %s for writing\n", outfname);

        write_gale(outf);
        fclose (outf);
	}	
	return 0;
}
