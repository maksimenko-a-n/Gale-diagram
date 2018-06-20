#include <stdio.h>
#include <string.h> // memcmp()
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <time.h>

#define DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_VERT 20 // The bit length of the type uint32_t

using namespace std; 

int nverts;

// Get the dimension and the number of vertices from the filename
int get_dv(char *filename, int &nv){
    char *pch = filename, *newpch;
    //pch = strchr(filename,'s');
    int dim = strtol(pch, &newpch, 10);
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

int compare_vertices (const void * a, const void * b)
{
  return memcmp(*(uint8_t**)a, *(uint8_t**)b, sizeof(uint8_t) * nverts);
}


//////////////
//    MAIN
/////////////

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		return 1;
	}
    
    int get_value = get_dv(argv[1], nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
        return 2;
    }

	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    unsigned long long inputlen = ftell(inf); // The length of the input file
    const unsigned long long MemorySize = 1E10;
    int nthreads = ((inputlen-1) / MemorySize) + 1;
    printf ("%d threads\n", nthreads);
    inputlen /=  nverts * sizeof(uint8_t); // The total number of diagrams
    unsigned long long one_file_size = inputlen / nthreads; // The num of diagrams in one new file
    //fseek(inf, (nthreads - rank - 1) * one_file_size * (nverts * sizeof(uint8_t)), SEEK_SET); // Go to begin of block
    fseek(inf, 0, SEEK_SET); // Go to begin of block

    vector<uint8_t *> all_diagrams;
    for (int rank = 0; rank < nthreads; rank++){
        char outfname[256]; // The output file name
        sprintf (outfname, "%da%d", nverts, rank);
        FILE *outf = fopen(outfname, "wb");
        if (outf == NULL){
            printf ("Write file ERROR: Cann't open the file %s\n", outfname);
            return 4;
        }
        printf ("Write in %s\n", outfname);

        for (unsigned long long n = 0; n < one_file_size; n++){
            uint8_t *data = (uint8_t *)malloc(sizeof(uint8_t)*(nverts));
            if (data == NULL){
                printf ("ERROR: Out of memory\n");
                break;
            }
            int read_size = fread (data, sizeof(uint8_t), nverts, inf);
            if(read_size < nverts){ // The end of input file
                free(data);
                break;
            }
            all_diagrams.push_back(data);
        }

        if (all_diagrams.size() == 0){
            fclose (outf);
            return 0;
        }
    
        qsort (all_diagrams.data(), all_diagrams.size(), sizeof(uint8_t *), compare_vertices);

        uint8_t **pt = all_diagrams.data();
        fwrite (*pt, sizeof(uint8_t), nverts, outf);
        int writen = 1;
        for (int i = 1; i < all_diagrams.size(); i++, pt++){
            if (compare_vertices(pt, pt+1)){
                fwrite (*(pt+1), sizeof(uint8_t), nverts, outf);
                writen++;
            }    
            free(*pt);
            *pt = NULL;
        }
        free(*pt);
        *pt = NULL;
        all_diagrams.clear();
        fclose (outf);
    }
    fclose (inf);
	return 0;
}
