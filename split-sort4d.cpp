#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memcmp()
#include <stdint.h>  // uint32_t
#include <time.h>

#define DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_VERT 20 // The bit length of the type uint32_t
#define BLOCK_SIZE 65536 // The size of memory block

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

 	FILE *logf = fopen("log.txt", "a");
	if (logf == NULL){
		printf ("Write file ERROR: Cann't open the file 'log'\n");
		return 3;
	}
    // Print the start time in log-file
    time_t rawtime;
    time ( &rawtime );
    fprintf (logf, "\n%sSplit and sort %s\n", ctime (&rawtime), argv[1]);
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    unsigned long long inputlen = ftell(inf); // The length of the input file
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    const unsigned long long MemorySize = 12E9; // RAM size
    //const unsigned long long MemorySize = 4294967296;
    //const unsigned long long MemorySize = 4E7; // RAM size
    int nfiles = ((inputlen-1) / MemorySize) + 1; // The number of parts
    fprintf (logf, "%d files will be written\n", nfiles);
    inputlen /=  nverts * sizeof(uint8_t); // The total number of diagrams
    unsigned long long one_file_size = (inputlen - 1) / nfiles + 1; // The num of diagrams in one new file
    unsigned long long nblocks = (one_file_size - 1) / BLOCK_SIZE + 1;
    one_file_size = nblocks * BLOCK_SIZE; // file size is multiple of BLOCK_SIZE

    fprintf (logf, "Allocate memory");
    fflush (logf);
	clock_t begt = clock();
	clock_t prevt = begt;
    uint8_t **all_diagrams = (uint8_t **)malloc(sizeof(uint8_t *) * one_file_size);
    if (all_diagrams == NULL){
        fprintf (logf, "ERROR: Cann't allocate memory for all_diagrams\n");
        fclose (logf);
        return 4;
    }
    uint8_t **data = (uint8_t **)malloc(sizeof(uint8_t *) * nblocks);
    for (unsigned long i = 0; i < nblocks; i++){
        data[i] = (uint8_t *)malloc(sizeof(uint8_t) * nverts * BLOCK_SIZE);
        if (data[i] == NULL){
            fprintf (logf, "ERROR: Out of memory on %d step\n", i);
            fclose (logf);
            return 5;
        }
    }

    
    clock_t curt;
    for (int part = 0; part < nfiles; part++){
        curt = clock();
        fprintf (logf, " (%2d sec)\n%d-th file: Read", (curt - prevt)/CLOCKS_PER_SEC, part);
        prevt = curt;
        fflush (logf);

        
        unsigned long num = 0;
        int read_size = BLOCK_SIZE;
        for (unsigned long row = 0; row < nblocks && read_size == BLOCK_SIZE; row++){
            uint8_t *currec = data[row];
            read_size = fread (currec, sizeof(uint8_t), BLOCK_SIZE * nverts, inf) / nverts;
            for (unsigned long col = 0; col < read_size; col++){
                all_diagrams[num] = currec;
                num++;
                currec += nverts;
            }
        }
        if (num == 0)  break;
    
        curt = clock();
        fprintf (logf, " %d (%2d sec), Sort", num, (curt - prevt)/CLOCKS_PER_SEC);
        prevt = curt;
        fflush (logf);

        qsort (all_diagrams, num, sizeof(uint8_t *), compare_vertices);
        curt = clock();
        fprintf (logf, " (%2d sec),", (curt - prevt)/CLOCKS_PER_SEC);
        prevt = curt;

        char outfname[256]; // The output file name
        sprintf (outfname, "%da%d", nverts, part);
        FILE *outf = fopen(outfname, "wb");
        if (outf == NULL){
            fprintf (logf, "Write file ERROR: Cann't open the file %s\n", outfname);
            break;
        }

        fprintf (logf, " Write in %s", outfname);
        fflush (logf);

        uint8_t **pt = all_diagrams;
        fwrite (*pt, sizeof(uint8_t), nverts, outf);
        unsigned long written = 1;
        for (unsigned long i = 1; i < num; i++, pt++){
            if (compare_vertices(pt, pt+1)){
                fwrite (*(pt+1), sizeof(uint8_t), nverts, outf);
                written++;
            }    
        }
        fprintf (logf, " %d", written);
        fclose (outf);
    }
    fclose (inf);
    curt = clock();
    fprintf (logf, " (%2d sec)\nDeallocate memory", (curt - prevt)/CLOCKS_PER_SEC);
    prevt = curt;
    fflush (logf);
    // Free memory
    for (unsigned long i = 0; i < nblocks; i++)
        free(data[i]);
    free(data);
    free(all_diagrams);
    curt = clock();
    fprintf (logf, " (%2d sec)\nElapsed time: %3d sec\n", (curt - prevt)/CLOCKS_PER_SEC, (curt - begt)/CLOCKS_PER_SEC);
    fclose (logf);
	return 0;
}
