#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <limits.h>
#include <stdint.h>  // int64_t
#include <float.h>  // DBL_EPSILON

#define MAX_DIM 16 // The maximum dimension of the Gale diagram (GD)
#define MAX_NUMBERS 2 // The maximum abs value of the coordinates of points in GD
#define MAX_VERT 64 // The bit length of the type int64_t
#define MAX_FACET 65536
#define LINE_SIZE 1024 // The maximum len of a line in the input file

using namespace std; 

// Write 0-1 matrix with ncols columns (ncols <= 64)
void write_incmatrix (FILE *wfile, const int ncols, const vector<int64_t> &inc_matrix);
// Find the edges number for the given inc-matrix
// Write the found edges in the file outf (if it isn't NULL)
int edges_number(const int vertices, const int facets, const vector<int64_t> &vertex_facet, FILE *outf = NULL);
int edges_number_long(const int vertices, const int facets, const vector< vector<int64_t> > &vertex_facet, FILE *outf = NULL);
// Transpose inc-matrix
vector< vector<int64_t> > transpose(const int facets, const int vertices, vector<int64_t> &facet_vertex);

// The Gale diagram (GD)
class Gale_diagram{
public:
//    FILE *outf; // The output file
	vector< vector<int> > vertices; // The list of points in the diagram
    int dimension; // The dimension of GD
    vector<int64_t> facet_vertex; // Incidence matrix (the result of calculations)
    
    Gale_diagram();
    ~Gale_diagram();
    int read(const char *filename); // Read GD from file
    void write(FILE *outf); // Write GD to file
    void write_facets(FILE *outf); // Write the list of the facets (sets of points)
    void init_Gauss(); // Allocate memory for the Gauss method. Calls the init_matrix(nrows, dimension+2)
    void init_matrix(int nrows, int ncols); // Allocate memory for matrix
    void free_matrix(); // Free memory of matrix
    int not_facet(int ncols); // Test if curface isn't a facet. Return 0 if it is a facet
    int facets_with_k_vert (int k, int startv, int curnv, int64_t curvertexset); // Find all facets with k vertices
    int find_facets (); // Find all facets
    int is_vertex(int p); // Test the point p from the array vertices
    int vertices_number(FILE *outf = NULL); 
private:
    double gauss_epsilon; // For fabs(x) < epsilon
    int nfacets_for_testing[MAX_DIM+1]; // Contains the num of facets with less or equal k vertices
    int *curface[MAX_VERT]; // The set of tested vertices
    double *matrix[MAX_VERT]; // For the Gauss method
    double *matrix_data; // The memory for the matrix
    inline int forward_elimination(int nrows, int ncols, double **M); // 1 stage of the Gauss method
};
