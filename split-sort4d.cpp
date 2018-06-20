#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memcmp()
#include <stdint.h>  // uint32_t
//#include <time.h>

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
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    const unsigned long long MemorySize = 16E9; // RAM size
    //const unsigned long long MemorySize = 4294967296;
    //const unsigned long long MemorySize = 4E7; // RAM size
    int nfiles = ((inputlen-1) / MemorySize) + 1; // The number of parts
    printf ("%d files will be written\n", nfiles);
    inputlen /=  nverts * sizeof(uint8_t); // The total number of diagrams
    unsigned long long one_file_size = inputlen / nfiles; // The num of diagrams in one new file
    //fseek(inf, (nfiles - rank - 1) * one_file_size * (nverts * sizeof(uint8_t)), SEEK_SET); // Go to begin of block

    printf ("Allocate memory\n");
    uint8_t **all_diagrams = (uint8_t **)malloc(sizeof(uint8_t *) * one_file_size);
    if (all_diagrams == NULL){
        printf ("ERROR: Cann't allocate memory for all_diagrams\n");
        return 4;
    }
    for (unsigned long i = 0; i < one_file_size; i++){
        uint8_t *data = (uint8_t *)malloc(sizeof(uint8_t)*(nverts));
        if (data == NULL){
            printf ("ERROR: Out of memory on %d step\n", i);
            return 5;
        }
        all_diagrams[i] = data;
    }        
    
    for (int part = 0; part < nfiles; part++){
        char outfname[256]; // The output file name
        sprintf (outfname, "%da%d", nverts, part);
        FILE *outf = fopen(outfname, "wb");
        if (outf == NULL){
            printf ("Write file ERROR: Cann't open the file %s\n", outfname);
            break;
        }

        printf ("%d-th block: Read", part);
        unsigned long num;
        for (num = 0; num < one_file_size; num++){
            if(fread (all_diagrams[num], sizeof(uint8_t), nverts, inf) < nverts)
                break; // The end of input file
        }
        if (num == 0)  break;
    
        printf (", Sort");
        qsort (all_diagrams, num, sizeof(uint8_t *), compare_vertices);

        printf (", Write in %s", outfname);
        uint8_t **pt = all_diagrams;
        fwrite (*pt, sizeof(uint8_t), nverts, outf);
        //unsigned long writen = 1;
        for (unsigned long i = 1; i < num; i++, pt++){
            if (compare_vertices(pt, pt+1)){
                fwrite (*(pt+1), sizeof(uint8_t), nverts, outf);
                //writen++;
            }    
        }
        fclose (outf);
        printf ("\n");
    }
    fclose (inf);
    printf ("Deallocate memory");
    for (unsigned long i = 0; i < one_file_size; i++)
        free(all_diagrams[i]);
    free(all_diagrams);
	return 0;
}
