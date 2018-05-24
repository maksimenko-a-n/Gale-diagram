#include "gale-diagram.hpp"
#include <time.h>

//////////////
//    MAIN
/////////////

/* 
Input is the name of a file with Gale diagram (list of points).
The example of a Gale diagram:
2 4
  1  1
 -1  1
  1 -1
 -1 -1
The number of points (lines in the input file) cann't be greater than MAX_VERT = 64 (the bit size of the type int64_t).
The values of the coordinates must be integers in range [-MAX_NUMBERS, MAX_NUMBERS].
The dimension of the Gale diagram (number of columns in the file) cann't be greater than MAX_DIM = 16.
Output is:
1) the list of sets of points -- facets of the appropriate polytope,
2) the checking whether the Gale diagram corresponds to a convex polytope,
and, if Yes,
3) the graph and the ridge-graph of the polytope.
*/

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		return 0;
	}
    
    Gale_diagram gale; // The Gale diagram
    if(gale.read(argv[1])){
        return 1; // Errors in the input file
    }

	char outfname[256]; // The output file name
	sprintf (outfname, "%s.out", argv[1]);
    FILE *outf = fopen (outfname, "w"); // Open the output file
	if (outf == NULL){
		printf ("ERROR: Cann't open file %s\n", outfname);
		return 1;
	}
	printf ("Output is saved in %s\n", outfname);
	gale.write (outf);
    /*
    uint64_t output[MAX_VERT];
    gale.lex_min(output);
    //int64_t bin_diagram[MAX_VERT];
    //gale.vector2bin(bin_diagram, all_permutations[5].data());
    gale.bin2vector(output);
	gale.write (outf);
	fclose (outf);
    return 0;
    */

    gale.init_Gauss (); // Allocate memory for the procedures of finding facets and vertices
	clock_t t = clock();
    int nfacets = gale.find_facets(); // Find facets
	t = clock() - t;
    printf ("Facets = %d (elapsed time: %4.3f sec)\n", nfacets, ((float)t)/CLOCKS_PER_SEC);
	gale.write_facets(outf); // Write facets into the file outf
    int nverts = gale.vertices.size();
    fprintf (outf, "Facets-vertices incidence matrix:\n");
    write_incmatrix (outf, nverts, gale.facet_vertex);
    // Check whether the Gale diagram corresponds to a convex polytope
    if (gale.is_polytope(outf)){
        t = clock();
        int ridges = edges_number(nfacets, nverts, gale.facet_vertex); // Find ridges
        fprintf (outf, "Ridges = %d\n", ridges);
        // Transpose facet_vertex to vertex_facet
        vector< vector<int64_t> > vertex_facet = transpose(nfacets, nverts, gale.facet_vertex);
        fprintf (outf, "Edges:\n");
        int edges = edges_number_long(nverts, nfacets, vertex_facet, outf); // Find edges
        t = clock() - t;
        fprintf (outf, "Edges = %d\n", edges);
        printf ("Edges: %d, Ridges: %d (elapsed time: %4.3f sec)\n", edges, ridges, ((float)t)/CLOCKS_PER_SEC);
        if (edges * 2 == nverts * (nverts-1))
            printf ("2-Neighborly polytope\n");
        if (ridges * 2 == nfacets * (nfacets-1))
            printf ("Dual 2-neighborly polytope\n");
    }
    else{
        printf ("This is not convex polytope!\n");
        fprintf (outf, "This is not convex polytope!\n");
    }
	fclose (outf);
	return 0;
}
