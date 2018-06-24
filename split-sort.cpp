#include "main_params.hpp"

const unsigned long BLOCK_SIZE = 65536; // The number of diagrams read from a file at a time
const unsigned long long MAX_MEMORY_SIZE = 2E9; // The maximum number of entries in the array all_diagrams
//const unsigned long long MEMORY_SIZE = 12E9; // RAM size
//const unsigned long long MEMORY_SIZE = 5E9; // RAM size
//const unsigned long long MemorySize = 4294967296;
//const unsigned long long MemorySize = 4E7; // RAM size

int nverts; // The length of one diagram (is used in compare_diagrams)

// Extract the dimension and the number of vertices from the filename
int get_dv(char *filename, int &nv){
    char *pch = filename, *newpch;
    //pch = strchr(filename,'s');
    int dim = strtol(pch, &newpch, 10);
    if (pch == newpch) return 1;
	if (dim != DIM) return 2;
	if (*newpch != 'd') return 3;
    pch = newpch + 1;
    nv = strtol(pch, &newpch, 10);
    if (pch == newpch) return 4;
	if (nv < 1 || nv > MAX_VERT) return 5;
    return 0;
}

int compare_diagrams (const void * a, const void * b){
    return memcmp(*(Tpoint**)a, *(Tpoint**)b, sizeof(Tpoint) * nverts);
}

//////////////
//    MAIN
/////////////

// Split huge file, sort parts and save them in several files
// Only different diagrams will be saved
int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: splitsort [file name]\n");
		return 1;
	}
    
    int get_value = get_dv(argv[1], nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
        return 2;
    }

 	FILE *logf = fopen("log", "a");
	if (logf == NULL){
		printf ("Write file ERROR: Cann't open the file 'log'\n");
		return 3;
	}
    // Print the start time in log-file
    time_t rawtime;
    time ( &rawtime );
    fprintf (logf, "\n%sSplit and sort %s\n", ctime (&rawtime), argv[1]);
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

    // Open the input file
	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    long long inputlen = ftell(inf); // The length of the input file
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    inputlen /=  nverts * sizeof(Tpoint); // The total number of diagrams
    fprintf (logf, "Input: %lld diagrams\n", inputlen);
    int nfiles = ((inputlen - 1) / MAX_MEMORY_SIZE) + 1; // The number of parts for splitting
    fprintf (logf, "%d files will be written\n", nfiles);
    //unsigned long long one_file_size = (inputlen - 1) / nfiles + 1; // The num of diagrams in one new file
    //unsigned long long nblocks = (one_file_size - 1) / BLOCK_SIZE + 1; // The number of blocks in one file
    unsigned long long nblocks = (inputlen - 1) / BLOCK_SIZE + 1; // Total number of blocks
    nblocks = (nblocks - 1) / nfiles + 1; // The number of blocks in one file
    unsigned long long one_file_size = nblocks * BLOCK_SIZE; // The num of diagrams in one new file is multiple of BLOCK_SIZE

    fprintf (logf, "Allocate memory");
    fflush (logf);
	clock_t begt = clock();
	clock_t prevt = begt;
    // Allocate memory for the array of pointers to diagrams
    Tpoint **all_diagrams = (Tpoint **)malloc(sizeof(Tpoint *) * one_file_size);
    if (all_diagrams == NULL){
        fprintf (logf, "ERROR: Cann't allocate memory for all_diagrams\n");
        fclose (logf);
        return 4;
    }
    // The input file will be read into data by blocks
    Tpoint **data = (Tpoint **)malloc(sizeof(Tpoint *) * nblocks);
    for (unsigned long i = 0; i < nblocks; i++){
        data[i] = (Tpoint *)malloc(sizeof(Tpoint) * nverts * BLOCK_SIZE);
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
        
        // Read
        unsigned long num = 0;
        int read_size = BLOCK_SIZE;
        for (unsigned long row = 0; row < nblocks && read_size == BLOCK_SIZE; row++){
            Tpoint *currec = data[row];
            read_size = fread (currec, sizeof(Tpoint), BLOCK_SIZE * nverts, inf) / nverts;
            for (unsigned long col = 0; col < read_size; col++){
                all_diagrams[num] = currec;
                num++;
                currec += nverts;
            }
        }
        if (num == 0)  break;
    
        curt = clock();
        fprintf (logf, " %lu (%2d sec), Sort", num, (curt - prevt)/CLOCKS_PER_SEC);
        prevt = curt;
        fflush (logf);

        // Sort
        qsort (all_diagrams, num, sizeof(Tpoint *), compare_diagrams);
        curt = clock();
        fprintf (logf, " (%2d sec),", (curt - prevt)/CLOCKS_PER_SEC);
        prevt = curt;

        // Open the new file for output
        char outfname[256];
        sprintf (outfname, "%da%d", nverts, part);
        FILE *outf = fopen(outfname, "wb");
        if (outf == NULL){
            fprintf (logf, "Write file ERROR: Cann't open the file %s\n", outfname);
            break;
        }

        fprintf (logf, " Write in %s", outfname);
        fflush (logf);

        // Write only different diagrams
        Tpoint **pt = all_diagrams;
        fwrite (*pt, sizeof(Tpoint), nverts, outf);
        unsigned long written = 1;
        for (unsigned long i = 1; i < num; i++, pt++){
            if (compare_diagrams(pt, pt+1)){ // If pt+1 differs from pt
                fwrite (*(pt+1), sizeof(Tpoint), nverts, outf);
                written++;
            }    
        }
        fprintf (logf, " %lu", written);
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
