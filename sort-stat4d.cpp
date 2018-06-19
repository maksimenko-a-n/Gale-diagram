#include <stdio.h>
#include <string.h> // memcmp()
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <time.h>

#define DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_VERT 32 // The bit length of the type uint32_t

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

	FILE *logf = fopen("log", "a");
	if (logf == NULL){
		printf ("Write file ERROR: Cann't open the file 'log'\n");
		return 1;
	}
    // Print the start time in log-file
    time_t rawtime;
    time ( &rawtime );
    fprintf (logf, "\n%sSort-stat %s\n", ctime (&rawtime), argv[1]);
    
    int get_value = get_dv(argv[1], nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
        return 2;
    }
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    unsigned long inputlen = ftell(inf)/nverts; // The total number of diagrams
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    fprintf (logf, "Input: %d diagrams\n", inputlen);
    fflush(logf);

    const uint32_t size = 256*256*(DIM+1);
    uint32_t all_diagrams[size];
    memset(all_diagrams, 0, size*sizeof(uint32_t));
    uint8_t diagram[MAX_VERT];
	clock_t begt = clock();
    uint64_t counter = 0;
    while (fread (diagram, sizeof(uint8_t), nverts, inf) == nverts){
        int first = diagram[0];
        uint32_t n;
        for (n = 0; n < DIM; n++, first >>= 2)
            if (first == 0)  break;
        all_diagrams[(n*256 + diagram[1])*256 + diagram[2]] += 1;
        counter++;
/*        if (counter == 10000 || counter % 100000 == 0){
            fprintf (logf, " Read %d", counter);
            if (counter != 0){
                clock_t curt = clock();
                int sec = (curt - begt)*(inputlen - counter) / (counter*CLOCKS_PER_SEC);
                int min = sec / 60;
                int hour = min / 60;
                fprintf (logf, " (left %d:%02d:%02d)", hour, min%60, sec%60);
            }
            fflush (logf);
        }    */
    }
	fclose (inf);
	clock_t t = clock() - begt;
    fprintf (logf, "Read time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
    fflush (logf);

    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%d.stat", DIM, nverts);
    fprintf (logf, "Write in %s\n", outfname);
	FILE *outf = fopen(outfname, "w");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 4;
	}
    fprintf (logf, "Totally = %d\n", inputlen);
    fflush (logf);
    
    uint32_t max = 0;
    uint32_t writen = 0;
    for (uint32_t i = 0; i < size; i++){
        if (all_diagrams[i] > 0){
            writen++;
            uint32_t x = i;
            fprintf (outf, "%d:%03d:%03d  %d\n", x >> 16, (x >> 8) & 255, x & 255, all_diagrams[i]);
            if (all_diagrams[i] > max)
                max = all_diagrams[i];
        }    
    }
    fprintf (outf, "Writen %d\n", writen);
    fprintf (outf, "MAX = %d\n", max);
	fclose (outf);
	fclose (logf);
	return 0;
}
