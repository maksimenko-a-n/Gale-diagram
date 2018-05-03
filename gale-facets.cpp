#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>  // int64_t

#define DIM 8
#define VERT 32
#define LINE_SIZE 1024
#define MAX_NUMBERS 2
#define MAX_VERT 64
#define MAX_FACET 64

int vertices[VERT][DIM]; // The input
int dimension, nvert, nfacets;
int nfacets_for_testing[DIM+1]; // Contains the num of facets with less or equal k vertices
int *curface[VERT]; // The set of tested vertices
int64_t facet_vertex[MAX_FACET+1]; // Incidence matrix
double gauss_data[VERT][DIM+2]; // The memory for matrix
double *matrix[VERT]; // The rows will be swapped many times
double epsilon = 1.0 / (VERT*VERT); // For fabs(x) < epsilon
//double coef[DIM+1]; // Result of the gauss
FILE *outf; // The output file


// Evaluate the number of edges of the incidence matrix 'vertex_facet'
int edges_number(const int vertices, const int facets, const int64_t *vertex_facet)
{
	int v1, v2, j;
    int64_t v1_facets, common_facets;

    int N_edges = 0;
    // Run over all vertices
	for (v1 = 0; v1 < vertices-1; v1++) {
        v1_facets = vertex_facet[v1];
        // Test if [v1,v2] is edge
        for (v2 = v1+1; v2 < vertices; v2++){
            // Find facets for every candidane edge [v1,v2]
            common_facets = v1_facets & vertex_facet[v2];
            // Compare common_facets and vertex_facet[j] for every j != v1 and j != v2
            for (j = 0; j < vertices; j++){
                if ((common_facets & vertex_facet[j]) == common_facets && (j - v1) * (j - v2) != 0) 
                    // if common_facets is contained in vertex_facet[j] and j != v1 and j != v2
                    goto NextEdge;
            }
            N_edges++;
            NextEdge: ;
        }    
    }
    return N_edges;
}

// Use the Gauss method for checking if the rank of the matrix M is equal to ncols
// Return 1 if the rank is equal to ncols, return 0 otherwise
int is_full_rank(int nrows, int ncols, double **M){
    // Gauss method: Forward Elimination
    int step, row, col;
    double diag_entry, y, *pt;
    for (step = 0; step < ncols-1; step++){
        // Find first nonzero element in the column 'step'
        for (row = step; row < nrows; row++){
            diag_entry = M[row][step];
            if (fabs(diag_entry) >= epsilon)
                break;
        }
        // If we cann't find nonzero element
        if (row >= nrows)
            return 0; // Isn't full rank
        if (row != step){// Swap rows
            pt = M[row];
            M[row] = M[step];
            M[step] = pt;
        }    
        // Normalize diag_entry and row 'step'
        if (diag_entry != 1){
            //M[step][step] = 1;
            for (col = step+1; col < ncols; col++)
                M[step][col] /= diag_entry;
        }
        for (row = step + 1; row < nrows; row++){
            y = M[row][step];
            if (fabs(y) >= epsilon){ // subtract step-th and row-th rows
                //M[row][step] = 0;
                for (col = step+1; col < ncols; col++)
                    M[row][col] -= y * M[step][col];
            }
        }
    }
    // For the last step: step = ncols-1
    for (row = step; row < nrows; row++){
        if (fabs(M[row][step]) >= epsilon)
            return 1; // Full rank matrix
    }
    // If we cann't find nonzero element
    return 0; // Isn't full rank
}

// Input is the vertex number from the array vertices
int is_vertex(int vertex){
    int64_t vf = 0 | ((int64_t)1 << vertex);
    int i, nf;
    // Find all the facets which contain the vertex
    int64_t funion = 0;
    for (i = 0, nf = 0; i < nfacets; i++){
        if (vf & facet_vertex[i]){
            funion |= facet_vertex[i];
            nf++;
        }    
    }
    if (nf < nvert - dimension - 1) // Too few common facets
        return 0;
    // Check the affine independence of the vertices in funion
    // Init the matrix
    int nrows = 0;
    for (i = 0, vf = 1; i < nvert; i++, vf <<= 1){
        if (vf & funion){// Add vertex to matrix
            for (int col = 0; col < dimension; col++)
                matrix[nrows][col] = vertices[i][col];
            matrix[nrows][dimension] = 1;
            nrows++;
        }
    }
    if (nrows <= dimension)
        return 0; // Too few rows
    return is_full_rank(nrows, dimension + 1, matrix);
}

// For counting vertices it uses: dimension, nvert, vertices, facet_vertex.
int vertices_number(){
    int num = 0;
    for (int i = 0; i < nvert; i++)
        num += is_vertex(i);
    return num;
}

// Test if face is a facet.
// ncols -- the number of vertices
// Return 0 if it is a facet
// Return 1 if it is not a facet, but can be if we append some vertices
// Return 2 or 3 if it is not a facet and cann't be a facet (after appending any vertices)
int not_facet(int ncols, int **face) {
    int row, col;
    // Init matrix. The columns are vertices (points) from the face
    for (row = 0; row < dimension; row++){
        for (col = 0; col < ncols; col++)
            matrix[row][col] = face[col][row];
        matrix[row][ncols] = 0;
	}
    // The last row in the matrix indicate --- we seak the affine hull
    for (col = 0; col <= ncols; col++)
        matrix[dimension][col] = 1;
        
    // Gauss method
    // Stage 1: Forward Elimination
    int step;
    double *pt;
    double diag_entry, y;
    for (step = 0; step < ncols-1; step++){
        // Find first nonzero element in the column 'step'
        for (row = step; row <= dimension; row++){
            diag_entry = matrix[row][step];
            if (fabs(diag_entry) >= epsilon)
                break;
        }
        // If we cann't find nonzero element
        if (row >= dimension)
            return 2; // Not a facet (singular system)
        // Normalize diag_entry and row 'row'
        if (diag_entry != 1){
            //matrix[row][step] = 1;
            for (col = step+1; col < ncols; col++)
                matrix[row][col] /= diag_entry;
        }
        if (row != step){// Swap rows
            pt = matrix[row];
            matrix[row] = matrix[step];
            matrix[step] = pt;
        }    
        for (row = step + 1; row <= dimension; row++){
            y = matrix[row][step];
            if (fabs(y) >= epsilon){ // subtract step-th and row-th rows
                for (col = step+1; col < ncols; col++)
                    matrix[row][col] -= y * matrix[step][col];
                //matrix[row][step] = 0;
                //matrix[row][ncols] -= 0;
            }
        }
    }
    
    // The last step: step == ncols - 1
    diag_entry = matrix[dimension][step];
    if (fabs(diag_entry) < epsilon)
        return 1; // Has no solution, but may be appended for a good solution
    // Normalize diag_entry and row 'dimension'
    //matrix[dimension][step] = 1;
    matrix[dimension][ncols] /= diag_entry;
    if (dimension != step){ // Swap rows
        pt = matrix[dimension];
        matrix[dimension] = matrix[step];
        matrix[step] = pt;
    }    
    // Test for impossibility of solution
    for (row = step + 1; row <= dimension; row++){
        if (fabs(matrix[row][step]) >= epsilon) // Looks like 0 * x == c, where c != 0.
            return 1; // Has no solution, but can be appended for a good solution
    }
    
    // Stage 2: back substitution
    // All diagonal entries are equal to 1
    double value;
    for (step = ncols - 1; step >= 0; step--){
        value = matrix[step][ncols];
        if (value < epsilon) // if value is nonpositive
            return 3; // Singular or will be singular (if we append some points to it)
        for (row = step - 1; row >= 0; row--)
            matrix[row][ncols] -= matrix[row][step] * value;
    }    
    return 0;
}

// Test of all co-faces with exactly k vertices
// startv -- the current index in 'vertices' (the set of all vertices)
// curnv -- current number of vertices in curface (the tested face)
// curvertexset is the characteristic vector of the set of vertices (negation of curface)
int facets_with_k_vert (int k, int startv, int curnv, int64_t curvertexset){
    // Evaluating the not_facet() for every new vertex is a bad idea
    // The solving of SLAE is expensive
	if (curnv >= k){
        int isnt_facet = not_facet(curnv, curface);
        if (isnt_facet == 0){
            facet_vertex[nfacets] = curvertexset;
           	nfacets++;
            if (nfacets > MAX_FACET) // ATTENTION!!!
                return 1;
        }
        return 0;
	}

    // Add one vertex to the curface and recursively call the faces_with_k_vert()
	int64_t one_bit = 1;
	int endv = nvert - k + curnv;
    // nfacets_for_testing[j] contains the number of facets with less or equal j vertices
    int i, i_max = nfacets_for_testing[curnv < k-1 ? curnv+1 : k-1];
	for (one_bit <<= startv; startv <= endv; one_bit <<= 1, startv++){
        int64_t newset = curvertexset - one_bit;
        // Check if the current set has a subset which is a facet
        // Generally, we have to test not all the facets, but only the ones with smaller number of vertices
        for (i = 0; (i_max - i) * ((newset & facet_vertex[i]) - newset) != 0; i++) ;
        if (i >= i_max){
            curface[curnv] = vertices[startv];
            if (facets_with_k_vert (k, startv+1, curnv+1, newset) == 1)
                return 1;
        }
	}
    return 0;
}

// Find all facets and save them into facet_vertex 
int find_facets (){
	nfacets = 0;
    nfacets_for_testing[0] = 0;
	for (int k = 2; k <= dimension+1; k++){
        nfacets_for_testing[k-1] = nfacets;
		if (facets_with_k_vert (k, 0, 0, ~(int64_t)0) == 1)
            return nfacets;
		//fprintf (outf, "N(%d-cofaces) = %d\n", k, nfaces);
		//printf ("N(%d-cofaces) = %d\n", k, nfaces);
    }
    return nfacets;
}

// Read Gale diagram from the file: dimension, nvert and vertices[][]
int read_gale(char *fname){
	FILE *inf = fopen(fname, "r");
	if (inf == NULL){
		printf ("ERROR: Cann't open the file %s\n", fname);
		return 1;
	}
	//printf ("Open input file %s\n", fname);

	char buffer[LINE_SIZE];
	int length, k;
	char *pch;
	pch = fgets (buffer, LINE_SIZE, inf);
	if ( pch == NULL)
		return 2;

	// Read dimension
	dimension = strtol (buffer, &pch, 10);
	if (dimension < 2 || dimension > DIM){
		printf ("ERROR: dimension = %d, but must be in [2, %d]\n", dimension, DIM);
		return 3;
	}
	printf ("dim = %d; ", dimension);

	// Read dimension
	nvert = strtol (pch, &pch, 10);
	if (nvert < 2 || nvert > VERT){
		printf ("ERROR: number of vertices = %d, but must be in [2, %d]\n", nvert, VERT);
		return 4;
	}
	printf ("vert = %d\n", nvert);
	
	int i, j;
	for (i = 0; i < nvert; i++){
		pch = fgets (buffer, LINE_SIZE, inf);
		if ( pch == NULL){
			printf ("ERROR: unexpected EOF\n");
			return 2;
		}	
		for (j = 0; j < dimension; j++){
			vertices[i][j] = strtol (pch, &pch, 10);
			if (fabs(vertices[i][j]) > MAX_NUMBERS){
				printf ("ERROR with MAX_NUMBERS: in line %d of %s\n", i+2, fname);
				return 4;
			}
		}
	}
	fclose (inf);
		
	return 0;
}	

// Write vertices[][] to file
void write_gale(FILE *f){
	fprintf (f, "Source:\n%d %d\n", dimension, nvert);
	int i, j;
	for (i = 0; i < nvert; i++){
		fprintf (f, "%2d:", i);
		for (j = 0; j < dimension; j++)
			fprintf (f, " %2d", vertices[i][j]);
		fprintf (f, "\n");
	}
}

// Write facets to file
void write_facets(FILE *f){
    fprintf (f, "Facets: %d\n", nfacets);
    int i, j;
    for (i = 0; i < nfacets; i++){
        int64_t mask = facet_vertex[i];
        for (j = 0; j < nvert; j++, mask >>= 1){
            if (mask & 1)
                fprintf (f, " %2d", j);
        }
        fprintf (f, "\n");
    }
}
	
//////////////
//
//    MAIN
//
/////////////


// Input is the name of a file with Gale diagram
// Output is the cofacets

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		return 0;
	}

	if (read_gale(argv[1]))
		return 1;
		
	char outfname[128];
	sprintf (outfname, "%s.out", argv[1]);
	outf = fopen(outfname, "w");
	if (outf == NULL){
		printf ("ERROR: Cann't open file %s\n", outfname);
		return 1;
	}
	printf ("Output is saved in %s\n", outfname);
        
	write_gale(outf);

    // Init matrix pointers for Gauss
    for (int i = 0; i < VERT; i++)
        matrix[i] = gauss_data[i];

	clock_t t;
	t = clock();
    find_facets();
	t = clock() - t; 
    write_facets(outf);
	if (nfacets > MAX_FACET){
        printf ("ERROR: The number of facets is bigger, than %d.\n", MAX_FACET);
        fclose (outf);
        return 0;
    }    
	fprintf (outf, "Facets: %d, %d clicks (%4.3f seconds)\n", nfacets, t, ((float)t)/CLOCKS_PER_SEC);
	printf ("Facets: %d, %d clicks (%4.3f seconds)\n", nfacets, t, ((float)t)/CLOCKS_PER_SEC);
    // Evaluating of the numbers of vertices, edges and ridges is much faster than for the facets
    int num_of_real_vertices = vertices_number();
    printf ("Polytope vertices: %d\n", num_of_real_vertices);
    fprintf (outf, "Polytope vertices: %d\n", num_of_real_vertices);
    if (num_of_real_vertices == nvert){
        // Transpose facet_vertex to vertex_facet
        int64_t vertex_facet[MAX_VERT]; // Incidence matrix
        int64_t v = 1; 
        for (int i = 0; i < nvert; i++, v <<= 1){
            vertex_facet[i] = 0;
            int64_t f = 1; 
            for (int j = 0; j < nfacets; j++, f <<= 1){
                if (facet_vertex[j] & v)
                    vertex_facet[i] |= f;
            }
        }
        int edges = edges_number(nvert, nfacets, vertex_facet);
        int ridges = edges_number(nfacets, nvert, facet_vertex);
        printf ("Edges: %d, Ridges: %d\n", edges, ridges);
        fprintf (outf, "Edges: %d, Ridges: %d\n", edges, ridges);
    }
	fclose (outf);
		
	return 0;
}
