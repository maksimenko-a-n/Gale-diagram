#include "gale-diagram.hpp"

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
    
    Gale_diagram gale;
    if(gale.read(argv[1])){
        return 1;
    }

	char outfname[128];
	sprintf (outfname, "%s.out", argv[1]);
    FILE *outf; // The output file
	outf = fopen (outfname, "w");
	if (outf == NULL){
		printf ("ERROR: Cann't open file %s\n", outfname);
		return 1;
	}
	printf ("Output is saved in %s\n", outfname);
        
	gale.write (outf);

    gale.init_Gauss ();  
    int nfacets = gale.find_facets ();
    printf ("Facets = %d\n", nfacets);

	gale.write_facets (outf);
    int num_of_real_vertices = gale.vertices_number();
    printf ("Vertices = %d\n", num_of_real_vertices);

    int nvert = gale.vertices.size();
    if (num_of_real_vertices == nvert){
        int ridges = edges_number(nfacets, nvert, gale.facet_vertex);
        // Transpose facet_vertex to vertex_facet
        vector< vector<int64_t> > vertex_facet = transpose(nfacets, nvert, gale.facet_vertex);
        int edges = edges_number_long(nvert, nfacets, vertex_facet);
        printf ("Edges: %d, Ridges: %d\n", edges, ridges);
        //fprintf (outf, "Edges: %d, Ridges: %d\n", edges, ridges);
    }

	fclose (outf);
		
	return 0;
}
