#include <stdio.h>
#include <stdlib.h>
//#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <string.h> // memcmp()

#define DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_VERT 20 // The bit length of the type uint32_t
#define MAX_FACET 20 // For the acceleration of the processing
//#define LINE_SIZE 1024 // The maximum len of a line in the input file

using namespace std; 

// Write 0-1 matrix with ncols columns (ncols <= 64)
void write_incmatrix (FILE *wfile, const int nrows, const int ncols, const uint32_t *inc_matrix);
// Find the edges number for the given inc-matrix
// Write the found edges in the file outf (if it isn't NULL)
int edges_number(const int vertices, const int facets, const uint32_t *vertex_facet, FILE *outf = NULL);
// Transpose inc-matrix
void transpose(const int nfacets, const int nverts, const uint32_t *facet_vertex, uint32_t *vertex_facet);

// The Gale diagram (GD)
class Gale_diagram{
public:
//    FILE *outf; // The output file
    int dimension; // The dimension of GD
    int nverts; // The number of points in GD
	double vertices[MAX_VERT][DIM]; // The list of points in the diagram
    uint8_t diagram[MAX_VERT];
    int nfacets; // The number of cofacets in GD
    uint32_t facet_vertex[MAX_FACET]; // Incidence matrix (the result of calculations)
    int cofacets_num[DIM+1]; // cofacets_num[k] is the number of cofacets with k vertices
    
    Gale_diagram();
    ~Gale_diagram(){};
    int read_diagram_bin(FILE *inf); // Read one diagram from the binary file
    int write_diagram_bin(FILE *outf); // Write one diagram to the binary file
    void write(FILE *outf); // Write GD to file
    void write_short(FILE *outf); // Write GD to file
    void write_facets(FILE *outf); // Write the list of the facets (sets of points)
    int not_facet(int ncols); // Check whether curface isn't a facet. Return 0 if it is a facet
    int facets_with_k_vert (int k, int startv, int curnv, uint32_t curvertexset); // Find all facets with k vertices
    int find_facets (); // Find all facets
    int is_vertex(int p); // Can vertices[p] be a vertex of the appropriate polytope?
    int is_polytope(FILE *outf = NULL); // Is it a diagram of a convex polytope?
    int vertices_in_facets(); // Count the total number vertices in cofacets
    int max_multiplicity(); // Count the maximum multiplicity of vertices
private:
    double gauss_epsilon; // For comparisons like fabs(x) < gauss_epsilon
    //int nfacets_for_testing[DIM+1]; // Contains the num of facets with less or equal k vertices
    double *current_coface[MAX_VERT]; // The set of tested vertices
    double *matrix[MAX_VERT]; // For the Gauss method
    double matrix_data[MAX_VERT*(DIM+1)]; // The memory for the matrix
    inline int forward_elimination(int nrows, int ncols, double **M); // First stage of the Gauss method
};
