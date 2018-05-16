#include <stdio.h>
#include <math.h>

#define MAX_DIM 16
#define MAX_VERT 64

int nvert, dimension;
int vertices[MAX_VERT][MAX_DIM];

// Generate vertices of {-1,1}-qube
void gen_cube_vertices(int d){
    dimension = d;
    nvert = (1 << dimension);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++, x >>= 1)
            vertices[i][j] = 2 * (x&1) - 1;
	}	
}

// Generate all vectors with coordinates from {-1,0,1}
void gen_01vectors(int d){
    dimension = d;
    nvert = pow(3,d);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++){
            vertices[i][j] = (x % 3) - 1;
            x /= 3;
        }
	}	
}

// Generate all vectors with coordinates from {-2,-1,1,2}
void gen_12vectors(int d){
    dimension = d;
    nvert = 1 << (2*dimension);
	int i, j, x;
	for (i = 0; i < nvert; i++){
        x = i;
		for (j = 0; j < dimension; j++, x >>= 2)
            vertices[i][j] = ((x&2) - 1) * (1 + (x&1));
	}	
}

// Generate vertices of the simplex with multiplicity mult
void gen_simplex_vertices(int d, int mult){
    dimension = d;
    nvert = (d+1) * mult;
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
void gen_cross_vertices(int d, int mult){
    dimension = d;
    nvert = 2 * d * mult;
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
void gen_cross_edges(int d){
    dimension = d;
    nvert = 2 * d * (d-1);
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


// Программа генерирует несколько примеров диаграмм Гейла
// и сохраняет их в файлы

int main(int argc, char *argv[])
{
	char outfname[64];
    // All {-2,-1,1,2}-vectors for the given dimension
    gen_12vectors(2);
	sprintf (outfname, "%dd%dv.g", dimension, nvert);
    write_gale(outfname);

    // All {-1,0,1}-vectors for the given dimension
    gen_01vectors(3);
	sprintf (outfname, "%dd%dv.g", dimension, nvert);
    write_gale(outfname);
    
    // The vertices of a cube
    gen_cube_vertices(5);
	sprintf (outfname, "%dd%dv-cube.g", dimension, nvert);
    write_gale(outfname);
    
    // The vertices of a simplex with the given multiplicity
    gen_simplex_vertices(10, 2);
	sprintf (outfname, "%dd%dv-simplex.g", dimension, nvert);
    write_gale(outfname);
    
    // The middles of edges of a cross polytope
    gen_cross_edges(5);
	sprintf (outfname, "%dd%dv-cedges.g", dimension, nvert);
    write_gale(outfname);
    
    // The vertices of a cross polytope with the given multiplicity
    gen_cross_vertices(8, 2);
	sprintf (outfname, "%dd%dv-cross.g", dimension, nvert);
    write_gale(outfname);
	return 0;
}
