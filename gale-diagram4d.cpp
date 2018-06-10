#include "gale-diagram4d.hpp"

Gale_diagram::Gale_diagram(){
    gauss_epsilon = sqrt(DBL_EPSILON); // For comparisons like fabs(x) < gauss_epsilon
    // Allocate memory for the matrix which is used in Gauss method: not_facet(), is_vertex()
    for (int i = 0; i < MAX_VERT; i++)
        matrix[i] = matrix_data + i * (DIM+1);
}

// Convert binary format to vertices
void Gale_diagram::convert_bin2vertices(uint8_t *diagram){
    for (int i = 0; i < nverts; i++){
        uint8_t x = diagram[i];
        for (int j = 0; j < DIM; j++, x >>= 2)
            vertices[i][j] = (x & 0b11) - 1;
    }
}	

// Read one diagram from the binary file
int Gale_diagram::read_diagram_bin(FILE *inf){
    uint8_t diagram[MAX_VERT];
    nfacets = 0; // Init the number of facets
    int s = fread (diagram, sizeof(uint8_t), nverts, inf);
    if (s < nverts) return 1;
    convert_bin2vertices(diagram);    
    return 0;
}	

// Write one diagram to the binary file
int Gale_diagram::write_diagram_bin(FILE *outf){
    uint8_t diagram[MAX_VERT];
    for (int i = 0; i < nverts; i++){
        diagram[i] = 0;
        for (int j = 0; j < DIM; j++)
            diagram[i] |= (uint8_t)(vertices[i][j] + 1) << (j*2);
    }
    int s = fwrite (diagram, sizeof(uint8_t), nverts, outf);
    if (s < nverts){
		printf ("Write file ERROR: Cann't write %d bytes\n", nverts);
		return 1;
	}
    return 0;
}	

// Write vertices[][] to file
void Gale_diagram::write(FILE *outf){
	fprintf (outf, "%d %d\n", DIM, nverts);
	int i, j;
	for (i = 0; i < nverts; i++){
		//fprintf (outf, "%2d:", i);
		for (j = 0; j < DIM; j++)
			fprintf (outf, " %2d", (int)(vertices[i][j]));
		fprintf (outf, "\n");
	}
}

// Write vertices[][] to file
void Gale_diagram::write_short(FILE *outf){
	fprintf (outf, "%d %d\n", DIM, nverts);
	int i, j;
	for (i = 0; i < nverts; i++){
		//fprintf (outf, "%2d:", i);
		for (j = -1; j >= 0; j--)
			fprintf (outf, "%d", (int)(vertices[i][j]) + 1);
		fprintf (outf, "\n");
	}
}

// Write facets to file
void Gale_diagram::write_facets(FILE *outf){
    fprintf (outf, "Facets %d:\n", nfacets);
    int i, j;
    for (i = 0; i < nfacets; i++){
        uint32_t mask = facet_vertex[i];
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
    for (row = 0; row < DIM; row++){
        for (col = 0; col < ncols; col++)
            matrix[row][col] = current_coface[col][row];
        //matrix[row][col] = -current_coface[col][row];
	}

    // Gauss method
    // Stage 1: Forward Elimination
    if (forward_elimination(DIM, ncols, matrix) != 0)
        return 2; // Not a facet (singular system)
    
    // Test for the impossibility of solution
    int step = ncols - 1;
    for (row = step; row < DIM; row++){
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
int Gale_diagram::facets_with_k_vert (int k, int startv, int curnv, uint32_t curvertexset){
    // Evaluating the not_facet() for every new set of vertices is a bad idea
    // The solving of SLAE is expensive
	if (curnv >= k){
        int isnt_facet = not_facet(curnv);
        if (isnt_facet == 0){
            // If the number of facets is too big, then close the program
            if (nfacets >= MAX_FACET)
                return 1;
            facet_vertex[nfacets] = curvertexset;
            nfacets++;
        }
        return 0;
	}

    // Add one vertex to the current_coface and recursively call the faces_with_k_vert()
	uint32_t one_bit = 1;
	int endv = nverts - k + curnv;
    // nfacets_for_testing[j] contains the number of facets with less or equal j vertices
    //int i, i_max = nfacets_for_testing[curnv < k-1 ? curnv+1 : k-1]; // optimization for small num of facets
	for (one_bit <<= startv; startv <= endv; one_bit <<= 1, startv++){
        uint32_t newset = curvertexset - one_bit; // the new set of vertices
        // Check whether the newset has a subset which is a facet
        // Generally, we have to test not all the facets, but only the ones with smaller number of vertices
        //for (i = 0; i < i_max && (newset & facet_vertex[i]) != newset; i++) ; // optimization for small num of facets
//        for (i = 0; (i_max - i) * ((newset & facet_vertex[i]) - newset) != 0; i++) ; // Optimization
        //if (i >= i_max){ // optimization for small num of facets
            // Add vertices[startv] to current_coface
            current_coface[curnv] = vertices[startv];
            if (facets_with_k_vert (k, startv+1, curnv+1, newset) == 1)
                return 1;
        //} // optimization for small num of facets
	}
    return 0;
}

// Find all facets with the last vertex and add them to facet_vertex
int Gale_diagram::facets_with_last_vert (){
	uint32_t one_bit = 1;
    int prev_nfacets = nfacets;
    uint32_t curvertexset = ~(uint32_t)0;
    nverts--; // Attention: temporarily the nverts will be less
    curvertexset -= one_bit << nverts; // add one vertex
    current_coface[0] = vertices[nverts];
    // Consistently test every set of k vertices 
	for (int k = 2; k <= DIM+1; k++){
        //nfacets_for_testing[k-1] = facet_vertex.size(); // optimization for small num of facets
		if (facets_with_k_vert (k, 0, 1, curvertexset) == 1){
            nverts++; // Restore nverts before quit
            return nfacets+1;
        }    
        cofacets_num[k] += nfacets - prev_nfacets;
        prev_nfacets = nfacets;
		//printf ("Number of %d-cofaces = %d\n", k, facet_vertex.size() - nfacets_for_testing[k-1]);
    }
    nverts++; // Restore nverts before quit
    return nfacets;
}

// Find all facets and save them into facet_vertex 
int Gale_diagram::find_facets (){
	uint32_t one_bit = 1;
    cofacets_num[0] = 0;
    nfacets = 0;
    int prev_nfacets = 0;
    //nfacets_for_testing[0] = 0; // optimization for small num of facets
    // Test if there is (0,...,0) among the vertices
/*    for (int v = 0; v < nverts; v++){
        int j;
        for (j = 0; j < dimension && vertices[v][j] == 0; j++) ;
        if (j == dimension){
            facet_vertex[nfacets] = (~(uint32_t)0) - (one_bit << v);
            nfacets++;
        }    
    }
    cofacets_num[1] = nfacets - prev_nfacets;
    prev_nfacets = nfacets;*/
    cofacets_num[1] = 0;
    
    // Consistently test every set of k vertices 
    // Attention! Cann't use DIM instead of dimension, when compile with -O2
	for (int k = 2; k <= DIM+1; k++){
        //nfacets_for_testing[k-1] = facet_vertex.size(); // optimization for small num of facets
		if (facets_with_k_vert (k, 0, 0, ~(uint32_t)0) == 1)
            return nfacets+1;
        cofacets_num[k] = nfacets - prev_nfacets;
        prev_nfacets = nfacets;
		//printf ("Number of %d-cofaces = %d\n", k, facet_vertex.size() - nfacets_for_testing[k-1]);
    }
    return nfacets;
}

// Input is the vertex number from the array vertices[]
int Gale_diagram::is_vertex(int p){
    if (nverts < DIM + 2)
        return 0; // Too few vertices
    uint32_t vf = 1;
    vf <<= p;
    int i, nf;
    // Find all the facets which contain the vertex
    uint32_t intersection = ~(uint32_t)0;
    for (i = 0, nf = 0; i < nfacets; i++){
        if (vf & facet_vertex[i]){
            intersection &= facet_vertex[i];
            nf++; // nf is the number of facets which contains the vertex
        }    
    }
    if (nf < nverts - DIM - 1)
        return 0; // Too few common facets

    uint32_t mask = 1; 
    if (nverts < MAX_VERT)
        mask = (mask << nverts) - 1;
    else
        mask = ~(uint32_t)0;
    // Now mask consists of (nverts) ones and (MAX_VERT - nverts) zeroes

    if ((intersection & mask) != vf)
        return 0; // Intersection of facets isn't coincide with the vertex
    // Check the rank of the vertices (without vertex p)
    // Init the matrix
    int nrows = 0;
    for (i = 0; i < nverts; i++){
        if (i != p){// Add vertex to matrix
            for (int col = 0; col < DIM; col++){
                matrix[nrows][col] = vertices[i][col];
            }    
            nrows++;
        }
    }
    // Gauss method: Forward Elimination
    if (forward_elimination(nrows, DIM, matrix) != 0)
        return 0; // Isn't full rank
    // The last step of Forward Elimination:
    for (int row = DIM-1; row < nrows; row++){
        if (fabs(matrix[row][DIM-1]) >= gauss_epsilon)
            return 1; // Full rank matrix
    }
    // If we cann't find nonzero element
    return 0; // Isn't full rank
}

// Is this a diagram of some convex polytope?
// Write the first bad vertex into outf
int Gale_diagram::is_polytope(FILE *outf){
    for (int i = 0; i < nverts; i++){
        if (is_vertex(i) == 0){
            if (outf != NULL)
                fprintf (outf, "Bad vertex: %d\n", i);
            return 0;
        }
    }
    return 1;
}


// Count vertices, that are incident at least one facet
int Gale_diagram::vertices_in_facets(){
    uint32_t intersection = ~(uint32_t)0;
    int i;
    for (i = 0; i < nfacets; i++)
        intersection &= facet_vertex[i];
    int num = 0;
    uint32_t funion = ~intersection;
    for (i = 0; i < nverts; i++, funion >>= 1)
        num += funion&1;
    return num;
}


// Write 0-1 matrix with ncols columns (ncols <= 32) into file wfile
void write_incmatrix (FILE *wfile, const int nrows, const int ncols, const uint32_t *inc_matrix){
	fprintf (wfile, "%d %d\n", nrows, ncols);
	for (int i = 0; i < nrows; i++){
		uint32_t x = inc_matrix[i];
		for (int j = 0; j < ncols; j++, x >>= 1)
			fprintf (wfile, " %d", x&1);
		fprintf (wfile, "\n");
	}
}


// Evaluate the number of edges of the incidence matrix 'vertex_facet'
int edges_number(const int vertices, const int facets, const uint32_t *vertex_facet, FILE *outf)
{
	int v1, v2, j;
    uint32_t v1_facets, common_facets;

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

// Transpose facet_vertex to vertex_facet
void transpose(const int nfacets, const int nverts, const uint32_t *facet_vertex, uint32_t *vertex_facet){
    uint32_t v = 1; 
    for (int i = 0; i < nverts; i++, v <<= 1){
        vertex_facet[i] = 0;
        uint32_t f = 1; 
        for (int j = 0; j < nfacets; j++, f <<= 1){
            if (facet_vertex[j] & v)
                vertex_facet[i] |= f;
        }
    }
}
