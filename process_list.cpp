#include "gale-diagram.hpp"
#include <string.h> // memcmp()
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
    
	FILE *inf = fopen(argv[1], "r");
	if (inf == NULL){
		printf ("ERROR: Cann't open the file %s\n", argv[1]);
		return 1;
	}
    FILE *outf, *out2f; // Open the output file
    Gale_diagram gale; // The Gale diagram
    int dimension, nverts = 0;
    int newd = 0;
    int min2facets; // Minimum number of facets for a 2-neighborly polytope
    for (int n = 0; ;){
        char buffer[LINE_SIZE]; // the buffer for reading lines from the file
        char *pch = fgets (buffer, LINE_SIZE, inf);
        if ( pch == NULL)
            break;
        const char begin[] = "begin";
        if (memcmp(begin, pch, sizeof(begin)-1) == 0){
            n++;
            if (n % 100 == 0)
                printf (" %d", n);
            int err = gale.read_diagram(inf);
            int d = gale.dimension, nv = gale.vertices.size();
            if (err)
                break;
            if (nverts == 0){
                dimension = d;
                nverts = nv;
                min2facets = 10 * nverts;
                gale.init_Gauss (); // Allocate memory for the procedures of finding facets and vertices
                char outfname[256]; // The output file name
                sprintf (outfname, "%dd%d-.gs", dimension, nverts);
                outf = fopen (outfname, "w"); // Open the output file
                if (outf == NULL){
                    printf ("ERROR: Cann't open file %s\n", outfname);
                    return 1;
                }
                printf ("Output is saved in %s\n", outfname);
                sprintf (outfname, "%dd%d-2ng.gs", dimension, nverts);
                out2f = fopen (outfname, "w"); // Open the output file
                if (out2f == NULL){
                    printf ("ERROR: Cann't open file %s\n", outfname);
                    return 1;
                }
            }
            else{
                if (nv != nverts || d != dimension){
                    printf ("Read file ERROR in the %d-th diagram: ", n);
                    if (nv != nverts)
                        printf ("The number of vertices %d differs from nvert = %d\n", nv, nverts);
                    else
                        printf ("The dim %d differs from the dimension = %d\n", d, dimension);
                    break;
                }    
            }
            // Process
            int nfacets = gale.find_facets(); // Find facets
            int ispoly = gale.is_polytope();
            int ridges = edges_number(nfacets, nverts, gale.facet_vertex); // Find ridges
            vector< vector<int64_t> > vertex_facet = transpose(nfacets, nverts, gale.facet_vertex);
            int edges;
            if (ispoly)    
                edges = edges_number_long(nverts, nfacets, vertex_facet); // Find edges
            int multiplicity = gale.max_multiplicity();
            int freev = nverts - gale.vertices_in_facets();
            //write_incmatrix (outf, nverts, gale.facet_vertex);
            //printf ("vertices in facets = %d\n", gale.vertices_in_facets());
            int max_freev = (nverts <= dimension ? nverts : (2 * dimension - nverts));
            max_freev = (max_freev < 0 ? 0 : max_freev);
            int is_write = 1;
            is_write &= (freev <= max_freev);
            is_write &= (nfacets <= 20);
            //is_write &= (multiplicity < 4);
            //if (nfacets > 0)
            //    is_write &= (ridges * 2 == nfacets * (nfacets - 1));
            if (is_write){
                newd++;
                FILE *wf = outf;
                if (ispoly && (edges*2 == nverts*(nverts-1))){
                    wf = out2f;
                    if (min2facets > nfacets)
                    min2facets = nfacets;
                }    
                fprintf (wf, "begin\n");
                gale.write (wf);
                fprintf (wf, "f = %d", nfacets);
                if (ispoly)
                    fprintf (wf, "; e = %d", edges);
                fprintf (wf, "\nr = %d\n", ridges);
                //fprintf (outf, "v = %d, maxv = %d, multiplicity = %d\n", freev, max_freev, multiplicity);
                fprintf (wf, "end\n");
            }    
        }
    }
    printf ("\nNew number of diagrams = %d\n", newd);
	fclose (inf);
	fclose (outf);
    if (min2facets < 10 * nverts){
        printf ("Min facets for 2-nghb = %d\n", min2facets);
        fprintf (out2f, "Min facets = %d\n", min2facets);
    }    
	fclose (out2f);
	return 0;
}
