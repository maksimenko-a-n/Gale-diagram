#include "gale-diagram.hpp"

Gale_diagram::Gale_diagram(){
    gauss_epsilon = sqrt(DBL_EPSILON); // For fabs(x) < epsilon
    matrix_data = NULL; 
}

Gale_diagram::~Gale_diagram(){
    free_matrix();
}

void Gale_diagram::init_matrix(int nrows, int ncols){
    free_matrix();
    matrix_data = (double*)malloc(nrows * ncols * sizeof(double));
    if (matrix_data == NULL){
        printf ("ERROR: cann't allocate memory for the matrix_data.\n");
        exit(1);
    }
    for (int i = 0; i < nrows; i++)
        matrix[i] = matrix_data + i * ncols;
}

// Allocate memory for the matrix which is used in Gauss method: not_facet(), is_vertex()
void Gale_diagram::init_Gauss(){
    int nrows = dimension < vertices.size() ? vertices.size() : dimension + 1;
    init_matrix(nrows, dimension+2);
}

void Gale_diagram::free_matrix(){
    if (matrix_data != NULL)
        free(matrix_data);
}

// Read Gale diagram from the file: dimension, nvert and vertices[][]
int Gale_diagram::read(const char *filename){
	FILE *inf = fopen(filename, "r");
	if (inf == NULL){
		printf ("ERROR: Cann't open the file %s\n", filename);
		return 1;
	}
	//printf ("Open input file %s\n", filename);

	char buffer[LINE_SIZE];
	int length, k;
	char *pch;
	pch = fgets (buffer, LINE_SIZE, inf);
	if ( pch == NULL){
		printf ("ERROR: unexpected EOF in %s\n", filename);
		return 2;
    }    

	// Read dimension
	dimension = strtol (buffer, &pch, 10);
	if (dimension < 2 || dimension > MAX_DIM){
		printf ("ERROR in %s: dimension = %d, but must be in [2, %d]\n", filename, dimension, MAX_DIM);
		return 3;
	}
	//printf ("dim = %d; ", dimension);

	// Read dimension
	int nvert = strtol (pch, &pch, 10);
	if (nvert < 2 || nvert > MAX_VERT){
		printf ("ERROR in %s: number of vertices = %d, but must be in [2, %d]\n", filename, nvert, MAX_VERT);
		return 4;
	}
	//printf ("vert = %d\n", nvert);
	
    vertices.clear(); // Remove all old points
	char *new_pch;
	int i, j;
	for (i = 0; i < nvert; i++){
		pch = fgets (buffer, LINE_SIZE, inf);
		if ( pch == NULL){
			printf ("ERROR: unexpected EOF in %s\n", filename);
			return 2;
		}	
        vertices.push_back(vector<int>());
		for (j = 0; j < dimension; j++){
			int x = strtol (pch, &new_pch, 10);
			if (pch == new_pch || fabs(x) > MAX_NUMBERS){
				printf ("ERROR in line %d of %s: coordinates of points must be integers in [%d, %d]\n", i+2, filename, -MAX_NUMBERS, MAX_NUMBERS);
				return 4;
			}
            pch = new_pch;
            vertices[i].push_back(x);
		}
	}
	fclose (inf);
		
	return 0;
}	

// Write vertices[][] to file
void Gale_diagram::write(FILE *outf){
    int nvert = vertices.size();
	fprintf (outf, "Source:\n%d %d\n", dimension, nvert);
	int i, j;
	for (i = 0; i < nvert; i++){
		fprintf (outf, "%2d:", i);
		for (j = 0; j < dimension; j++)
			fprintf (outf, " %2d", vertices[i][j]);
		fprintf (outf, "\n");
	}
}

// Write facets to file
void Gale_diagram::write_facets(FILE *outf){
    int nfacets = facet_vertex.size();
    int nvert = vertices.size();
    fprintf (outf, "Facets %d:\n", nfacets);
    int i, j;
    for (i = 0; i < nfacets; i++){
        int64_t mask = facet_vertex[i];
        for (j = 0; j < nvert; j++, mask >>= 1){
            if (mask & 1)
                fprintf (outf, " %2d", j);
        }
        fprintf (outf, "\n");
    }
}


// Gauss method: Forward Elimination
// Return 1 if the rank of M is smaller than ncols-1
// Otherwise return 0
inline int Gale_diagram::forward_elimination(int nrows, int ncols, double **M){
    int step, row, col;
    double diag_entry, y, *pt;
    for (step = 0; step < ncols-1; step++){
        // Find first nonzero element in the column 'step'
        for (row = step; row < nrows; row++){
            diag_entry = M[row][step];
            if (fabs(diag_entry) >= gauss_epsilon)
                break; // We have found a nonzero element
        }
        // If we cann't find nonzero element
        if (row >= nrows)
            return 1; // Isn't full rank
        // Normalize the nonzero element and the row
        if (diag_entry != 1){
            //M[step][step] = 1;
            for (col = step+1; col < ncols; col++)
                M[row][col] /= diag_entry;
        }
        if (row != step){// Swap rows
            pt = M[row];
            M[row] = M[step];
            M[step] = pt;
        }    
        for (row = step + 1; row < nrows; row++){
            y = M[row][step];
            if (fabs(y) >= gauss_epsilon){ // subtract step-th and row-th rows
                //M[row][step] = 0;
                for (col = step+1; col < ncols; col++)
                    M[row][col] -= y * M[step][col];
            }
        }
    }
    return 0;
}

// Test if curface is a facet.
// ncols -- the number of vertices
// Return 0 if it is a facet
// Return 1 if it is not a facet, but can be if we append some vertices
// Return 2 or 3 if it is not a facet and cann't be a facet (after appending any vertices)
int Gale_diagram::not_facet(int ncols){ //, int64_t x){
    int row, col;
    // Init matrix. The columns are vertices (points) from the curface
    for (row = 0; row < dimension; row++){
        for (col = 0; col < ncols; col++)
            matrix[row][col] = curface[col][row];
        matrix[row][ncols] = 0;
	}
    // The last row in the matrix indicate --- we seak the affine hull
    for (col = 0; col <= ncols; col++)
        matrix[dimension][col] = 1;

    // Gauss method
    // Stage 1: Forward Elimination
    if (forward_elimination(dimension, ncols, matrix) != 0)
        return 2; // Not a facet (singular system)
    int step;
    double y;
    for (step = 0; step < ncols-1; step++){
        y = matrix[dimension][step];
        if (fabs(y) >= gauss_epsilon){ // subtract step-th and dimension-th rows
            for (col = step+1; col < ncols; col++)
                matrix[dimension][col] -= y * matrix[step][col];
        }
    }
    
    // The last step: step == ncols - 1
    double diag_entry = matrix[dimension][step];
    if (fabs(diag_entry) < gauss_epsilon){
        return 1; // Has no solution, but may be appended for a good solution
    }
    // Normalize diag_entry and row 'dimension'
    //matrix[dimension][step] = 1;
    matrix[dimension][ncols] /= diag_entry;
    if (dimension != step){ // Swap rows
        double *pt = matrix[dimension];
        matrix[dimension] = matrix[step];
        matrix[step] = pt;
    }    
    // Test for impossibility of solution
    for (row = step + 1; row <= dimension; row++){
        if (fabs(matrix[row][step]) >= gauss_epsilon) // Looks like 0 * x == c, where c != 0.
            return 1; // Has no solution, but can be appended for a good solution
    }
    
    // Stage 2: back substitution
    // All diagonal entries are equal to 1
    double value;
    for (step = ncols - 1; step >= 0; step--){
        value = matrix[step][ncols];
        if (value < gauss_epsilon){ // if value is nonpositive
            return 3; // Singular or will be singular (if we append some points to it)
        }
        for (row = step - 1; row >= 0; row--)
            matrix[row][ncols] -= matrix[row][step] * value;
    }    
    return 0;
}

// Test of all co-faces with exactly k vertices
// startv -- the current index in 'vertices' (the set of all vertices)
// curnv -- current number of vertices in curface (the tested face)
// curvertexset is the characteristic vector of the set of vertices (negation of curface)
int Gale_diagram::facets_with_k_vert (int k, int startv, int curnv, int64_t curvertexset){
    // Evaluating the not_facet() for every new vertex is a bad idea
    // The solving of SLAE is expensive
	if (curnv >= k){
        int isnt_facet = not_facet(curnv);
        if (isnt_facet == 0){
            facet_vertex.push_back(curvertexset);
            if (facet_vertex.size() >= MAX_FACET) // ATTENTION!!!
                return 1;
        }
        return 0;
	}

    // Add one vertex to the curface and recursively call the faces_with_k_vert()
	int64_t one_bit = 1;
	int endv = vertices.size() - k + curnv;
    // nfacets_for_testing[j] contains the number of facets with less or equal j vertices
    int i, i_max = nfacets_for_testing[curnv < k-1 ? curnv+1 : k-1];
	for (one_bit <<= startv; startv <= endv; one_bit <<= 1, startv++){
        int64_t newset = curvertexset - one_bit;
        // Check if the current set has a subset which is a facet
        // Generally, we have to test not all the facets, but only the ones with smaller number of vertices
        for (i = 0; i < i_max && (newset & facet_vertex[i]) != newset; i++) ;
//        for (i = 0; (i_max - i) * ((newset & facet_vertex[i]) - newset) != 0; i++) ;
        if (i >= i_max){
            curface[curnv] = vertices[startv].data();
            if (facets_with_k_vert (k, startv+1, curnv+1, newset) == 1)
                return 1;
        }
	}
    return 0;
}

// Find all facets and save them into facet_vertex 
int Gale_diagram::find_facets (){
	facet_vertex.clear();
    nfacets_for_testing[0] = 0;
	for (int k = 2; k <= dimension+1; k++){
        nfacets_for_testing[k-1] = facet_vertex.size();
		if (facets_with_k_vert (k, 0, 0, ~(int64_t)0) == 1)
            return facet_vertex.size();
		//fprintf (outf, "N(%d-cofaces) = %d\n", k, nfaces);
		//printf ("N(%d-cofaces) = %d\n", k, nfaces);
    }
    return facet_vertex.size();
}

// Input is the vertex number from the array vertices
int Gale_diagram::is_vertex(int vertex){
    int64_t vf = 1;
    vf <<= vertex;
    int i, nf;
    // Find all the facets which contain the vertex
    int64_t funion = 0;
    int nfacets = facet_vertex.size();
    for (i = 0, nf = 0; i < nfacets; i++){
        if (vf & facet_vertex[i]){
            funion |= facet_vertex[i];
            nf++; // nf is the number of facets which contains the vertex
        }    
    }
    int nvert = vertices.size();
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
    // Gauss method: Forward Elimination
    if (forward_elimination(nrows, dimension + 1, matrix) != 0)
        return 0; // Isn't full rank
    // The last step of Forward Elimination:
    for (int row = dimension; row < nrows; row++){
        if (fabs(matrix[row][dimension]) >= gauss_epsilon)
            return 1; // Full rank matrix
    }
    // If we cann't find nonzero element
    return 0; // Isn't full rank
}

int Gale_diagram::vertices_number(FILE *outf){
    int num = 0;
    vector<int> bad_vertices;
    for (int i = 0; i < vertices.size(); i++){
        if (is_vertex(i))
            num++;
        else
            bad_vertices.push_back(i);
    }
    if (outf != NULL){
        if (bad_vertices.size() > 0){
            fprintf (outf, "Bad vertices:");
            for (int i = 0; i < bad_vertices.size(); i++)
                fprintf (outf, " %2d", bad_vertices[i]);
            fprintf (outf, "\n");
        }
    }
    return num;
}


// Write 0-1 matrix with cols columns (cols <= 64)
void write_matrix (FILE *wfile, const int ncols, const vector<int64_t> &inc_matrix){
	int nrows = inc_matrix.size();
    const int64_t *data = inc_matrix.data();
	fprintf (wfile, "%d %d\n", nrows, ncols);
	for (int i = 0; i < nrows; i++){
		int64_t x = data[i];
		for (int j = 0; j < ncols; j++, x >>= 1)
			fprintf (wfile, " %d", x&1);
		fprintf (wfile, "\n");
	}
	//fprintf (wfile, "end\n");
}


// Evaluate the number of edges of the incidence matrix 'vertex_facet'
int edges_number(const int vertices, const int facets, const vector<int64_t> &vertex_facet, FILE *outf)
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
                    break;
            }
            if (j >= vertices){
                if (outf != NULL)
                    fprintf (outf, " %2d %2d\n", v1, v2);
                N_edges++;
            }    
        }    
    }
    return N_edges;
}

// Evaluate the number of edges of the incidence matrix 'vertex_facet'
// For the BIG number of facets ( > 64)
int edges_number_long(const int vertices, const int facets, const vector< vector<int64_t> > &vertex_facet, FILE *outf)
{
	int v1, v2, j, s;
    int size = (facets - 1) / 64 + 1;
    const int64_t *v1_facets, *v2_facets;
    vector<int64_t> common_facets(size,0);

    int N_edges = 0;
	for (v1 = 0; v1 < vertices-1; v1++) {
        v1_facets = vertex_facet[v1].data();
        for (v2 = v1+1; v2 < vertices; v2++){
            v2_facets = vertex_facet[v2].data();
            for (s = 0; s < size; s++)
                common_facets[s] = v1_facets[s] & v2_facets[s];
            for (j = 0; j < vertices; j++){
                if ((j - v1) * (j - v2) != 0){ 
                    for (s = 0; s < size; s++){
                        if ((common_facets[s] & vertex_facet[j][s]) != common_facets[s])
                            break;
                    }        
                    if (s == size) // common_facets is contained in vertex_facet[j]
                        break; // (v1, v2) is not an edge
                }    
            }
            if (j >= vertices){
                if (outf != NULL)
                    fprintf (outf, " %2d %2d\n", v1, v2);
                N_edges++;
            }    
        }    
    }
    return N_edges;
}

// Transpose facet_vertex to vertex_facet
vector< vector<int64_t> > transpose(const int nfacets, const int nvert, vector<int64_t> &facet_vertex){
    int size = (nfacets - 1) / 64 + 1;
    vector< vector<int64_t> > vertex_facet(nvert, vector<int64_t>(size,0)); // Incidence matrix
    int64_t v = 1; 
    for (int i = 0; i < nvert; i++, v <<= 1){
        //vertex_facet[i] = 0;
        int64_t f = 1; 
        int s = 0;
        for (int j = 0; j < nfacets; j++, f <<= 1){
            if (f == 0){
                s++;
                f = 1;
            }
            if (facet_vertex[j] & v)
                vertex_facet[i][s] |= f;
        }
    }
    return vertex_facet;
}