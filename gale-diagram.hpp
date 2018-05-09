#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <limits.h>
#include <stdint.h>  // int64_t
#include <float.h>  // DBL_EPSILON

#define MAX_DIM 16
#define MAX_VERT 64
#define MAX_NUMBERS 2
#define MAX_FACET 64
#define LINE_SIZE 1024
//#define BIG 4503599627370496*12
//12*15

using namespace std; 

static double epsilon = sqrt(DBL_EPSILON); // For fabs(x) < epsilon
int is_full_rank(int nrows, int ncols, double **M);
int edges_number(const int vertices, const int facets, const vector<int64_t> &vertex_facet);

// The Gale diagram
class Gale_diagram{
public:
//    FILE *outf; // The output file
	vector< vector<int> > vertices; // The list of points in the diagram
    int dimension;
    vector<int64_t> facet_vertex; // Incidence matrix
    
    Gale_diagram();
    ~Gale_diagram();
    int read(const char *filename);
    void write(FILE *outf);
    void write_facets(FILE *outf);
    void init_Gauss();    
    void init_matrix(int nrows, int ncols);    
    void free_matrix();    
    int not_facet(int ncols); // Test if curface isn't facet. Return 0 if it is a facet
    int facets_with_k_vert (int k, int startv, int curnv, int64_t curvertexset); // Find all facets with k vertices
    int find_facets (); // Find all facets
    int is_vertex(int vertex);
    int vertices_number();
private:
    //double epsilon = 1.0 / (VERT*VERT); // For fabs(x) < epsilon
    int nfacets_for_testing[MAX_DIM+1]; // Contains the num of facets with less or equal k vertices
    int *curface[MAX_VERT]; // The set of tested vertices
    double *matrix[MAX_VERT];
    double *matrix_data;
    //int matrix_rows;
};
