#include "gale-diagram.hpp"

Gale_diagram::Gale_diagram(){
    gauss_epsilon = sqrt(DBL_EPSILON); // For comparisons like fabs(x) < gauss_epsilon
    matrix_data = NULL; 
}

Gale_diagram::~Gale_diagram(){
    free_matrix(); // Deallocate the memory of matrix
}

// Allocate memory for matrix
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

// Deallocate the memory of matrix
void Gale_diagram::free_matrix(){
    if (matrix_data != NULL)
        free(matrix_data);
}

// Read Gale diagram from the file: dimension, nverts and vertices[][]
int Gale_diagram::read_diagram(FILE *inf){
	char buffer[LINE_SIZE]; // the buffer for reading lines from the file
    // Read the first line
	char *pch = fgets (buffer, LINE_SIZE, inf);
	if ( pch == NULL){
		printf ("Read file ERROR: unexpected EOF\n");
		return 1;
    }    

	// Read dimension
	dimension = strtol (pch, &pch, 10);
	if (dimension < 2 || dimension > MAX_DIM){
		printf ("Read file ERROR: dimension = %d, but must be in [2, %d]\n", dimension, MAX_DIM);
		return 2;
	}

	// Read dimension
	int nverts = strtol (pch, &pch, 10);
	if (nverts < 2 || nverts > MAX_VERT){
		printf ("Read file ERROR: number of vertices = %d, but must be in [2, %d]\n", nverts, MAX_VERT);
		return 3;
	}
	
    vertices.clear(); // Remove all old points
    // Read points from the file
	for (int i = 0; i < nverts; i++){
		pch = fgets (buffer, LINE_SIZE, inf);
		if ( pch == NULL){
			printf ("Read file ERROR: unexpected EOF\n");
			return 1;
		}	
        vertices.push_back(vector<int>());
		for (int j = 0; j < dimension; j++){
            char *new_pch;
			int x = strtol (pch, &new_pch, 10);
			if (pch == new_pch || fabs(x) > MAX_NUMBERS){
				printf ("Read file ERROR: in line %d: coordinates of points must be integers in [%d, %d]\n", i+2, -MAX_NUMBERS, MAX_NUMBERS);
				return 4;
			}
            pch = new_pch;
            vertices[i].push_back(x);
		}
	}
	return 0;
}	

// Read Gale diagram from the file: dimension, nverts and vertices[][]
int Gale_diagram::read(const char *filename){
	FILE *inf = fopen(filename, "r");
	if (inf == NULL){
		printf ("ERROR: Cann't open the file %s\n", filename);
		return 1;
	}
    int err = read_diagram(inf);
	fclose (inf);
	return err;
}	

// Write vertices[][] to file
void Gale_diagram::write(FILE *outf){
    int nverts = vertices.size();
	fprintf (outf, "%d %d\n", dimension, nverts);
	int i, j;
	for (i = 0; i < nverts; i++){
		//fprintf (outf, "%2d:", i);
		for (j = 0; j < dimension; j++)
			fprintf (outf, " %2d", vertices[i][j]);
		fprintf (outf, "\n");
	}
}

// Write facets to file
void Gale_diagram::write_facets(FILE *outf){
    int nfacets = facet_vertex.size();
    int nverts = vertices.size();
    fprintf (outf, "Facets %d:\n", nfacets);
    int i, j;
    for (i = 0; i < nfacets; i++){
        int64_t mask = facet_vertex[i];
        for (j = 0; j < nverts; j++, mask >>= 1){
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

// Test if current_coface is a facet.
// ncols -- the number of vertices
// Return 0 if it is a facet
// Return 1 if it is not a facet, but can be if we append some vertices
// Return 2 or 3 if it is not a facet and cann't be a facet (after appending any vertices)
int Gale_diagram::not_facet(int ncols){ 
    int row, col;
    // Init matrix. The columns are vertices (points) from the current_coface
    for (row = 0; row < dimension; row++){
        for (col = 0; col < ncols; col++)
            matrix[row][col] = current_coface[col][row];
        //matrix[row][col] = -current_coface[col][row];
	}

    // Gauss method
    // Stage 1: Forward Elimination
    if (forward_elimination(dimension, ncols, matrix) != 0)
        return 2; // Not a facet (singular system)
    
    // Test for the impossibility of solution
    int step = ncols - 1;
    for (row = step; row < dimension; row++){
        if (fabs(matrix[row][step]) >= gauss_epsilon) // The equation looks like 0 * x == c, where c != 0.
            return 1; // Has no solution, but can be appended for a good solution
    }
    
    // Stage 2: back substitution
    // All diagonal entries are equal to 1
    double value;
    for (step = ncols - 2; step >= 0; step--){
        value = -matrix[step][ncols-1];
        if (value < gauss_epsilon){ // if value is nonpositive
            return 3; // Singular or will be singular (if we append some points to it)
        }
        for (row = step - 1; row >= 0; row--)
            matrix[row][ncols-1] += matrix[row][step] * value;
    }    
    return 0;
}

// This is a recursive function.
// It enumerates and checks all co-faces with exactly k vertices.
// The one function call adds one vertex from vertices[startv,...,nverts] to the current_coface 
// startv is the current index in 'vertices' (the set of all vertices)
// curnv is the number of vertices in current_coface (the tested coface)
// current_coface is the array of pointers to the current set of vertices
// curvertexset is the characteristic vector of the current face (the complement of current_coface)
int Gale_diagram::facets_with_k_vert (int k, int startv, int curnv, int64_t curvertexset){
    // Evaluating the not_facet() for every new set of vertices is a bad idea
    // The solving of SLAE is expensive
	if (curnv >= k){
        int isnt_facet = not_facet(curnv);
        if (isnt_facet == 0){
            facet_vertex.push_back(curvertexset);
            // If the number of facets is too big, then close the program
            if (facet_vertex.size() >= MAX_FACET)
                return 1;
        }
        return 0;
	}

    // Add one vertex to the current_coface and recursively call the faces_with_k_vert()
	int64_t one_bit = 1;
	int endv = vertices.size() - k + curnv;
    // nfacets_for_testing[j] contains the number of facets with less or equal j vertices
    //int i, i_max = nfacets_for_testing[curnv < k-1 ? curnv+1 : k-1]; // optimization for small num of facets
	for (one_bit <<= startv; startv <= endv; one_bit <<= 1, startv++){
        int64_t newset = curvertexset - one_bit; // the new set of vertices
        // Check whether the newset has a subset which is a facet
        // Generally, we have to test not all the facets, but only the ones with smaller number of vertices
        //for (i = 0; i < i_max && (newset & facet_vertex[i]) != newset; i++) ; // optimization for small num of facets
//        for (i = 0; (i_max - i) * ((newset & facet_vertex[i]) - newset) != 0; i++) ; // Optimization
        //if (i >= i_max){ // optimization for small num of facets
            // Add vertices[startv] to current_coface
            current_coface[curnv] = vertices[startv].data();
            if (facets_with_k_vert (k, startv+1, curnv+1, newset) == 1)
                return 1;
        //} // optimization for small num of facets
	}
    return 0;
}

// Find all facets and save them into facet_vertex 
int Gale_diagram::find_facets (){
	facet_vertex.clear(); // Clear the incidence matrix
	int64_t one_bit = 1;
    //nfacets_for_testing[0] = 0; // optimization for small num of facets
    // Test if there is (0,...,0) among the vertices
    int v, j;
    for (v = 0; v < vertices.size(); v++){
        for (j = 0; j < dimension && vertices[v][j] == 0; j++) ;
        if (j == dimension)
            facet_vertex.push_back((~(int64_t)0) - (one_bit << v));
    }
    // Consistently test every set of k vertices 
	for (int k = 2; k <= dimension+1; k++){
        //nfacets_for_testing[k-1] = facet_vertex.size(); // optimization for small num of facets
		if (facets_with_k_vert (k, 0, 0, ~(int64_t)0) == 1)
            return facet_vertex.size();
		//printf ("Number of %d-cofaces = %d\n", k, facet_vertex.size() - nfacets_for_testing[k-1]);
    }
    return facet_vertex.size();
}

// Input is the vertex number from the array vertices[]
int Gale_diagram::is_vertex(int p){
    int nverts = vertices.size();
    if (nverts < dimension + 2)
        return 0; // Too few vertices
    int64_t vf = 1;
    vf <<= p;
    int i, nf;
    // Find all the facets which contain the vertex
    int64_t intersection = ~(int64_t)0;
    int nfacets = facet_vertex.size();
    for (i = 0, nf = 0; i < nfacets; i++){
        if (vf & facet_vertex[i]){
            intersection &= facet_vertex[i];
            nf++; // nf is the number of facets which contains the vertex
        }    
    }
    if (nf < nverts - dimension - 1)
        return 0; // Too few common facets

    int64_t mask = 1; 
    if (nverts < MAX_VERT)
        mask = (mask << nverts) - 1;
    else
        mask = ~(int64_t)0;
    // Now mask consists of (nverts) ones and (MAX_VERT - nverts) zeroes

    if ((intersection & mask) != vf)
        return 0; // Intersection of facets isn't coincide with the vertex
    // Check the rank of the vertices (without vertex p)
    // Init the matrix
    int nrows = 0;
    for (i = 0; i < nverts; i++){
        if (i != p){// Add vertex to matrix
            for (int col = 0; col < dimension; col++){
                matrix[nrows][col] = vertices[i][col];
            }    
            nrows++;
        }
    }
    // Gauss method: Forward Elimination
    if (forward_elimination(nrows, dimension, matrix) != 0)
        return 0; // Isn't full rank
    // The last step of Forward Elimination:
    for (int row = dimension-1; row < nrows; row++){
        if (fabs(matrix[row][dimension-1]) >= gauss_epsilon)
            return 1; // Full rank matrix
    }
    // If we cann't find nonzero element
    return 0; // Isn't full rank
}

// Is this a diagram of some convex polytope?
// Write the first bad vertex into outf
int Gale_diagram::is_polytope(FILE *outf){
    vector<int> bad_vertices;
    for (int i = 0; i < vertices.size(); i++){
        if (is_vertex(i) == 0){
            if (outf != NULL)
                fprintf (outf, "Bad vertex: %d\n", i);
            return 0;
        }
    }
    return 1;
}

// Transform vertices to output with permutations of coordinates
void Gale_diagram::vector2bin(uint64_t *output, char *permut){
	uint64_t transform[5] = {0,1,3,14,15};
    int i, j, p;
    int nverts = vertices.size();
    for (i = 0; i < nverts; i++)
        output[i] = 0;
    for (j = 0; j < dimension; j++){
        p = permut[j] * 4;
        for (i = 0; i < nverts; i++){
            output[i] |= (transform[vertices[i][j] + 2] << p);
        }
    }
}

// Transform binary input to vertices
void Gale_diagram::bin2vector(uint64_t *input){
	int map[16] = {-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,2};
    int i, j;
    int nverts = vertices.size();
    for (i = 0; i < nverts; i++){
        uint64_t x = input[i];
        for (j = 0; j < dimension; j++, x >>= 4){
            vertices[i][j] = map[x & 0b1111];
        }
    }
}

// Normalize binary input
void Gale_diagram::normalize_bin(uint64_t *input){
    int i, j;
    int nverts = vertices.size();
    for (i = 0; i < nverts; i++){
        uint64_t x = input[i];
        uint64_t mask = 0b1001;
        for (j = 0; j < dimension; j++, x >>= 4, mask <<= 4){
            if ((x & 0b1111) == 0b1100)
                input[i] -= mask;
        }
    }
}


// Find the lexicographically minimal and write it in output
void Gale_diagram::lex_min(uint64_t *output){
    vector< vector<char> > all_permutations;
    all_permutations = gen_permute(dimension);
    int npermut = all_permutations.size();
/*    fprintf (outf, "Permutations:\n");
    for (int i = 0; i < npermut; i++){
        for (int j = 0; j < gale.dimension; j++)
            fprintf (outf, " %d", all_permutations[i][j]);
        fprintf (outf, "\n");
    }
    fprintf (outf, "\n");
*/    
    output[0] = ~(uint64_t)0;
    uint64_t current1[MAX_VERT], current2[MAX_VERT];
    int nverts = vertices.size();
    for (int i = 0; i < npermut; i++){ // Choose permutation of coordinates
        vector2bin(current1, all_permutations[i].data());
        int nreflections = (1 << dimension);
        for (int j = 0; j < nreflections; j++){ // Choose reflection
            uint64_t mask = 0; // mask for reflections
            uint64_t mask2 = 0b1111; // one-bit reflection
            uint32_t mask3 = 1; // one-bit reflection
            int x, k;
            for (x = j, k = 0; k < dimension; k++, x >>= 1, mask2 <<= 4){
                mask |= (x&1) * mask2;
                mask3 |= (mask3 << 4);
            }    
            for (k = 0; k < nverts; k++){ // Realize reflection
                uint32_t cur = current1[k];
                uint32_t y = cur & (cur >> 1);
                y = ((y >> 2) | (~y)) & mask3; // Some magic for zeroes 0b0011
                y |= y << 1;
                y |= y << 2;
                current2[k] = cur ^ (mask & y);
            }    
            int is_compare = 1; // Do compare current2 and output?
            // Bubble sort
            for (k = 0; k < nverts-1; k++){ // Sorting of current2 and comparing with output
                int m_min = k;
                int min = current2[m_min];
                for (int m = k+1; m < nverts; m++){
                    if (min > current2[m]){
                        m_min = m;
                        min = current2[m_min];
                    }    
                }        
                if (is_compare){
                    if (min > output[k])
                        break;
                    if (min < output[k]){
                        output[k] = min;
                        is_compare = 0;
                    }    
                }    
                else
                    output[k] = min;
                // Swap    
                current2[m_min] = current2[k];
                current2[k] = min;
            }
            if (k == nverts-1){
                if ((!is_compare) || current2[k] < output[k])
                    output[k] = current2[k];
            }
        }
    }
    normalize_bin(output);
}

// Count vertices, that are incident at least one facet
int Gale_diagram::vertices_in_facets(){
    int64_t intersection = ~(int64_t)0;
    int nfacets = facet_vertex.size();
    int i;
    for (i = 0; i < nfacets; i++)
        intersection &= facet_vertex[i];
    int nverts = vertices.size();
    int num = 0;
    int64_t funion = ~intersection;
    for (i = 0; i < nverts; i++, funion >>= 1)
        num += funion&1;
    return num;
}

// Count the maximum multiplicity of vertices
// Suppose nverts > 0
int Gale_diagram::max_multiplicity(){
    int nverts = vertices.size(); 
    int max = 1;
    int multiplicity = 1;
    for (int i = 1; i < nverts; i++){
        if (memcmp(vertices[i-1].data(), vertices[i].data(), sizeof(int)*dimension) == 0){
            multiplicity++; // increase multiplicity
            if (multiplicity > max)
                max = multiplicity;
        }
        else{
            multiplicity = 1;
        }
    }
    return max;
}

// Write 0-1 matrix with ncols columns (ncols <= 64) into file wfile
void write_incmatrix (FILE *wfile, const int ncols, const vector<int64_t> &inc_matrix){
	int nrows = inc_matrix.size();
    const int64_t *data = inc_matrix.data();
	fprintf (wfile, "%d %d\n", nrows, ncols);
	for (int i = 0; i < nrows; i++){
		int64_t x = data[i];
		for (int j = 0; j < ncols; j++, x >>= 1)
			fprintf (wfile, " %d", x&1);
		fprintf (wfile, "\n");
	}
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
vector< vector<int64_t> > transpose(const int nfacets, const int nverts, vector<int64_t> &facet_vertex){
    int size = (nfacets - 1) / 64 + 1;
    vector< vector<int64_t> > vertex_facet(nverts, vector<int64_t>(size,0)); // The output matrix
    int64_t v = 1; 
    for (int i = 0; i < nverts; i++, v <<= 1){
        int64_t f = 1; 
        int s = 0;
        for (int j = 0; j < nfacets; j++, f <<= 1){
            if (f == 0){ // New block of facets
                s++;
                f = 1;
            }
            if (facet_vertex[j] & v)
                vertex_facet[i][s] |= f;
        }
    }
    return vertex_facet;
}

void recursive_permute(vector< vector<char> > &array, vector<char> &cur_perm, int d, int step){
    if (step == d-1){
        array.push_back(vector<char>(cur_perm));
        return;
    }
        
    int x = cur_perm[step];
    for(int i = step; i < d; i++){
        cur_perm[step] = cur_perm[i];
        cur_perm[i] = x;
        recursive_permute(array, cur_perm, d, step+1);
        cur_perm[i] = cur_perm[step];
        cur_perm[step] = x;
    }
}

vector< vector<char> > gen_permute(int d){
    int nperm = 1, i, j;
    for (i = 2; i <= d; i++)
        nperm *= i;
    vector< vector<char> > permut_array; // The output
    vector<char> cur_permute(d,0);
    for (i = 0; i < d; i++)
        cur_permute[i] = i;
    recursive_permute(permut_array, cur_permute, d, 0);
    return permut_array;
}
