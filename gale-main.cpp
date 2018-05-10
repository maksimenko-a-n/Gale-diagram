#include "gale-diagram.hpp"
#include <time.h>

//////////////
//    MAIN
/////////////

// Input is the name of a file with Gale diagram (list of points)
// Output is the cofacets

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		return 0;
	}
    
    Gale_diagram gale; // The Gale diagram
    if(gale.read(argv[1])){
        // Errors in the input file
        return 1;
    }

	char outfname[128]; // The output file name
	sprintf (outfname, "%s.out", argv[1]);
    FILE *outf; // The output file
	outf = fopen (outfname, "w");
	if (outf == NULL){
		printf ("ERROR: Cann't open file %s\n", outfname);
		return 1;
	}
	printf ("Output is saved in %s\n", outfname);
        
	gale.write (outf);

    // Allocate memory for the procedures of finding facets and vertices
    gale.init_Gauss ();  
	clock_t t;
	t = clock();
    int nfacets = gale.find_facets(); // Find facets
	t = clock() - t;
    printf ("Facets = %d (elapsed time: %4.3f sec)\n", nfacets, ((float)t)/CLOCKS_PER_SEC);
	gale.write_facets(outf); // Write facets into the file outf
    int nvert = gale.vertices.size();
    fprintf (outf, "Facets-vertices incidence matrix:\n");
    write_matrix (outf, nvert, gale.facet_vertex);
    // Find the real number of vertices of the appropriate convex polytope
    int num_of_real_vertices = gale.vertices_number(outf);
    //printf ("Vertices = %d\n", num_of_real_vertices);
    //fprintf (outf, "Vertices = %d\n", num_of_real_vertices);

    if (num_of_real_vertices == nvert){
        fprintf (outf, "Ridges:\n");
        int ridges = edges_number(nfacets, nvert, gale.facet_vertex, outf);
        fprintf (outf, "Ridges = %d\n", ridges);
        // Transpose facet_vertex to vertex_facet
        vector< vector<int64_t> > vertex_facet = transpose(nfacets, nvert, gale.facet_vertex);
        fprintf (outf, "Edges:\n");
        int edges = edges_number_long(nvert, nfacets, vertex_facet, outf);
        fprintf (outf, "Edges = %d\n", edges);
        printf ("Edges: %d, Ridges: %d\n", edges, ridges);
    }
    else{
        printf ("This is not convex polytope!\n");
        fprintf (outf, "This is not convex polytope!\n");
    }

	fclose (outf);
		
	return 0;
}
