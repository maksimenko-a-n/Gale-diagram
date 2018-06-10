#include "gale-diagram.hpp"
#include <string.h> // memcmp()
#include <time.h>

// Get the dimension and the number of vertices from the filename
int get_dv(char *filename, int &dim, int &nv){
    char *pch = filename, *newpch;
    //pch = strchr(filename,'s');
    dim = strtol(pch, &newpch, 10);
    if (pch == newpch)
        return 1;
	if (dim != DIM)
        return 2;
	if (*newpch != 'd')
        return 3;
    pch = newpch + 1;
    nv = strtol(pch, &newpch, 10);
    if (pch == newpch)
        return 4;
	if (nv < 1 || nv > MAX_VERT)
        return 5;
    return 0;
}

void process_one_diagram(FILE *outf, FILE *out2f, Gale_diagram &gale, int &newd, int &min2facets){
    int dimension = gale.dimension;
    int old_nverts = gale.nverts;
    uint8_t *diagram = gale.diagram;
    int M = pow(3, dimension);
    uint8_t zero = 0;
    for (int k = 0; k < dimension; k++)
        zero |= ((uint8_t)1 << (k*2));
    // is_old_vertices[i] == 0 <=>  i is new vertex and is not zero
    char is_old_vertices[1 << 8];
    memset (is_old_vertices, 0, 1 << 8);
    is_old_vertices[zero] = 1;
    for (int i = 0; i < old_nverts; i++)
        is_old_vertices[diagram[i]] = 1;
    // Add one vertex
    gale.nverts = old_nverts + 1;
    for (int i = 0; i < M; i++){
        diagram[old_nverts] = 0;
        int x = i;
        for (int j = 0; j < dimension; j++, x /= 3){
            uint8_t y = x % 3;
            diagram[old_nverts] |= (y << (j*2));
            gale.vertices[old_nverts][j] = y - 1;
        }    
        if (is_old_vertices[diagram[old_nverts]])
            continue;
        // Process    
        int nfacets = gale.find_facets(); // Find facets
        int ridges = edges_number(nfacets, gale.nverts, gale.facet_vertex); // Find ridges
        int ispoly = gale.is_polytope();
        int edges = 0;
        if (ispoly){    
            uint32_t vertex_facet[MAX_VERT];
            transpose(nfacets, gale.nverts, gale.facet_vertex, vertex_facet);
            edges = edges_number(gale.nverts, nfacets, vertex_facet); // Find edges
        }    
        int freev = gale.nverts - gale.vertices_in_facets();
        int max_freev = (gale.nverts <= dimension ? gale.nverts : (2 * dimension - gale.nverts));
        max_freev = (max_freev < 0 ? 0 : max_freev);
        int is_write = 1;
        is_write &= (freev <= max_freev);
        is_write &= (nfacets <= 20);
        int big = gale.cofacets_num[dimension+1] + gale.cofacets_num[dimension];
        is_write &= (big <= 0);
        is_write &= (gale.cofacets_num[2] == 0);
        if (is_write){
            newd++;
            if (ispoly && (edges*2 == gale.nverts*(gale.nverts-1))){
                if (min2facets > nfacets)
                    min2facets = nfacets;
                fprintf (out2f, "b\n");
                gale.write (out2f);
                fprintf (out2f, "f=%d", nfacets);
                if (ispoly)
                    fprintf (out2f, "; e=%d", edges);
                fprintf (out2f, "\nr=%d\n", ridges);
                //fprintf (out2f, "v = %d, maxv = %d, multiplicity = %d\n", freev, max_freev, multiplicity);
                fprintf (out2f, "e\n");
                //write_incmatrix (out2f, nfacets, gale.nverts, gale.facet_vertex);
            }
            else
                gale.write_diagram_bin (outf);
        }
    }
    gale.nverts = old_nverts;
}

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
The number of points (lines in the input file) cann't be greater than MAX_VERT = 20 (the bit size of the type uint32_t).
The values of the coordinates must be integers in range [-1, 1].
The dimension of the Gale diagram (number of columns in the file) cann't be greater than DIM = 4.
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
		return 1;
	}
    
    int dimension, nverts;
    int get_value = get_dv(argv[1], dimension, nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
        return 2;
    }
    printf ("dimension = %d, nverts = %d\n", dimension, nverts);
    
	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%dp.gb", dimension, nverts+1);
	FILE *outf = fopen(outfname, "wb");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 4;
	}
    sprintf (outfname, "%dd%d-2ng.gs", dimension, nverts+1);
	FILE *out2f = fopen(outfname, "w");
	if (out2f == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 5;
	}

    Gale_diagram gale; // The Gale diagram
    gale.dimension = dimension;
    gale.nverts = nverts;
    int newd = 0;
    int min2facets = 10 * nverts; // Minimum number of facets for a 2-neighborly polytope
	clock_t t = clock();
    for (int n = 0; ;n++){
        if (n % 1000 == 0)
            printf (" %d", n);
        if (gale.read_diagram_bin(inf))
            break;
        // Process
        int nfacets = gale.find_facets(); // Find facets
        int ispoly = gale.is_polytope();
        process_one_diagram(outf, out2f, gale, newd, min2facets);
    }
	t = clock() - t;
    printf ("\nElapsed time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
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
