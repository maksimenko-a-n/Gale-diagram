#include <stdio.h>
#include <string.h> // memcmp()
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <time.h>
#include <mpi.h>

#define I_AM_FREE_TAG 2
#define TO_WORK_TAG 4
#define FOR_PROCESSING_TAG 8
#define RESULT_TAG 16
#define RESULT_READY_TAG 32
#define BUNCH 50

#define DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_NUMBERS 1 // The maximum abs value of the coordinates of points in GD
#define MAX_VERT 32 // The bit length of the type uint32_t
#define MAX_FACET 32 // To avoid the memory overfull
#define LINE_SIZE 1024 // The maximum len of a line in the input file

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

// Read one diagram from the binary file
int read_diagram_bin(FILE *inf, int nv, uint8_t *diagram){
    int s = fread (diagram, sizeof(uint8_t), nv, inf);
    if (s < nv){
		//printf ("Read file ERROR: Unexpected End Of File\n");
		return 1;
	}
    return 0;
}	

// Write one diagram to the binary file
int write_diagram_bin(FILE *outf, int nv, uint8_t *diagram){
    int s = fwrite (diagram, sizeof(uint8_t), nv, outf);
    if (s < nv){
		printf ("Write file ERROR: Cann't write %d bytes\n", nv);
		return 2;
	}
    return 0;
}	


void recursive_permute(vector< vector<char> > &array, vector<char> &cur_perm, int step){
    if (step == DIM-1){
        array.push_back(vector<char>(cur_perm));
        return;
    }
        
    int x = cur_perm[step];
    for(int i = step; i < DIM; i++){
        cur_perm[step] = cur_perm[i];
        cur_perm[i] = x;
        recursive_permute(array, cur_perm, step+1);
        cur_perm[i] = cur_perm[step];
        cur_perm[step] = x;
    }
}

vector< vector<char> > gen_permute(){
    int nperm = 1, i, j;
    for (i = 2; i <= DIM; i++)
        nperm *= i;
    vector< vector<char> > permut_array; // The output
    vector<char> cur_permute(DIM,0);
    for (i = 0; i < DIM; i++)
        cur_permute[i] = i;
    recursive_permute(permut_array, cur_permute, 0);
    return permut_array;
}

// Transform vertices to output with permutations of coordinates
void permute_coord(uint8_t *output, uint8_t *input, int nv, char *permut){
    for (int i = 0; i < nv; i++){
        output[i] = 0;
        uint8_t x = input[i];
        for (int j = 0; j < DIM; j++, x >>= 2){
            int p = permut[j] * 2;
            output[i] |= ((x & 0b11) << p);
        }    
    }    
}

// Find the lexicographically minimal and write it in output
void lex_min(uint8_t *output, uint8_t *input, int nv, vector< vector<char> > *all_permutations){
    int npermut = all_permutations->size();
    output[0] = ~(uint8_t)0;
    uint8_t current1[MAX_VERT], current2[MAX_VERT];
    for (int i = 0; i < npermut; i++){ // Choose permutation of coordinates
        permute_coord(current1, input, nv, (*all_permutations)[i].data());
        int nreflections = (1 << DIM);
/*        printf ("Perm:");
        for (int j = 0; j < dim; j++)
            printf ("%2d", all_permutations[i][j]);
        printf ("\n");
*/        for (int j = 0; j < nreflections; j++){ // Choose reflection
            uint8_t mask = 0; // mask for reflections
            uint8_t mask2 = 0b10; // one-bit reflection
            uint8_t mask3 = 1; // mask for zeroes
            int x, k;
//            printf ("Refl:");
            for (x = j, k = 0; k < DIM; k++, x >>= 1, mask2 <<= 2){
                mask |= (x&1) * mask2;
//                printf ("%2d", x&1);
                mask3 |= (mask3 << 2);
            }    
//            printf (": ");
            for (k = 0; k < nv; k++){ // Realize reflection
                uint8_t cur = current1[k];
                uint8_t y;
                y = ((cur >> 1) | (~cur)) & mask3; // Some magic for zeroes 0b01
                y |= y << 1;
                current2[k] = cur ^ (mask & y);
            }    
/*            int map[16] = {-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,2};
            for (k = 0; k < nv; k++){
                uint32_t y = current2[k];
                for (int m = 0; m < dim; m++, y >>= 4)
                    printf (" %2d", map[y & 0b1111]);
                printf (" %03X;", current2[k]);
            }
*/            int is_compare = 1; // Do compare current2 and output?
            // Bubble sort
            for (k = 0; k < nv-1; k++){ // Sorting of current2 and comparing with output
                int m_min = k;
                int min = current2[m_min];
                for (int m = k+1; m < nv; m++){
                    if (min > current2[m]){
                        m_min = m;
                        min = current2[m_min];
                    }    
                }        
                if (is_compare){
                    if (min > output[k])
                        break;
                    if (min < output[k]){
                        output[k] = min;
                        is_compare = 0;
//                        printf (" best");
                    }    
                }    
                else
                    output[k] = min;
                // Swap    
                current2[m_min] = current2[k];
                current2[k] = min;
            }
            if (k == nv-1){
//                printf (" cur = %03X; out = %03X;", current2[k], output[k]);
                if ((!is_compare) || current2[k] < output[k]){
                    output[k] = current2[k];
//                    if (is_compare)
//                        printf (" best");
                }    
            }
//            printf ("\n");
        }
    }
    //normalize_bin(output, dim, nv);
};

int compare_diagrams (const void * a, const void * b)
{
  return memcmp(*(uint8_t**)a, *(uint8_t**)b, sizeof(uint8_t) * nverts);
}


//////////////
//    MAIN
/////////////

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv); // Init MPI
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		MPI_Finalize();
		return 1;
	}
    
    int get_value = get_dv(argv[1], nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
		MPI_Finalize();
        return 2;
    }

    uint8_t diagram[BUNCH*MAX_VERT], sorted[BUNCH*MAX_VERT];
    
	// Get the rank of the curent thread and the number of all threads
	int rank, nthreads;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
	MPI_Status status;
    MPI_Request request;

	if (rank != 0) {
        vector< vector<char> > all_permutations;
        all_permutations = gen_permute();
		int size = 0;
		MPI_Isend(&size, 1, MPI_INT, 0, RESULT_READY_TAG, MPI_COMM_WORLD, &request);
		while(1) {
            int size;
			MPI_Recv(&size, 1, MPI_INT, 0, TO_WORK_TAG, MPI_COMM_WORLD, &status);
			if(size < 1) break;
			MPI_Recv(diagram, size*nverts*sizeof(uint8_t), MPI_BYTE, 0, FOR_PROCESSING_TAG, MPI_COMM_WORLD, &status);
            for (int i = 0; i < size; i++)
                lex_min(sorted + i*nverts, diagram + i*nverts, nverts, &all_permutations);  // The main process
			MPI_Isend(&size, 1, MPI_INT, 0, RESULT_READY_TAG, MPI_COMM_WORLD, &request);
            MPI_Send(sorted, size*nverts*sizeof(uint8_t), MPI_BYTE, 0, RESULT_TAG, MPI_COMM_WORLD);
		}
		printf ("%d's thread fin  ", rank);
		MPI_Finalize();
		return 0;
	}
    
	FILE *logf = fopen("log", "a");
	if (logf == NULL){
		printf ("Write file ERROR: Cann't open the file 'log'\n");
		MPI_Finalize();
		return 1;
	}
    // Print the start time in log-file
    time_t rawtime;
    time ( &rawtime );
    fprintf (logf, "\n%sLex sort %s\nOn %d threads\n", ctime (&rawtime), argv[1], nthreads);
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

    // Open the input file
	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		MPI_Finalize();
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    unsigned long inputlen = ftell(inf)/nverts; // The total number of diagrams
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    fprintf (logf, "Input: %d diagrams\n", inputlen);
    // Open the output file
    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%dlex.gb", DIM, nverts);
    fprintf (logf, "Write in %s\n", outfname);
	FILE *outf = fopen(outfname, "wb");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		MPI_Finalize();
		return 4;
	}
	clock_t begt = clock();
	long long sent = 0; // amount of input
	long long received = 0; // amount of output
//	int rcv_tag = MPI_ANY_TAG;
	long long to_slaves = 0, from_slaves = 1 - nthreads, end_of_file = 0;
	while(!end_of_file || to_slaves != from_slaves) {
        int size = 0;
		MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, RESULT_READY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG != RESULT_READY_TAG){
			printf ("ERROR: Unexpected status.MPI_TAG = %d\n", status.MPI_TAG);
			break;
		}
        from_slaves++;
		if (size > 0){ // Receive results from worker
            MPI_Recv(sorted, size*nverts*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
            fwrite (sorted, sizeof(uint8_t), size*nverts, outf);
            received += size;
            fflush(outf);
            if (received%1000000 == 0){
                fprintf (logf, " Received %d (%d calls)", received, from_slaves);
            }
        }
        size = 0;
        if (!end_of_file){
            // Read diagram from file and send to a thread
            int read_size = fread (diagram, sizeof(uint8_t), BUNCH*nverts, inf);
            size = read_size / nverts;
            if(size < BUNCH) // The end of input file
                end_of_file = 1;
        /*  for (size = 0; size < BUNCH; size++){
                int read_size = fread (diagram + size*nverts, sizeof(uint8_t), nverts, inf);
                if(read_size < nverts){ // The end of input file
                    end_of_file = 1;
                    break;
                }    
            }  */
        }
            
        // Send diagram to worker
        MPI_Isend(&size, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD, &request);
        if (size <= 0)
            continue;
		
		MPI_Send(diagram, size*nverts*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, FOR_PROCESSING_TAG, MPI_COMM_WORLD);
		to_slaves++;
		sent += size;
        if (sent == 50000 || sent % 1000000 == 0){
            fprintf (logf, " Sent %d", sent);
            if (sent != 0){
                clock_t curt = clock();
                int sec = (curt - begt)*(inputlen - sent) / (sent*CLOCKS_PER_SEC);
                int min = sec / 60;
                int hour = min / 60;
                fprintf (logf, " (left %d:%02d:%02d)", hour, min%60, sec%60);
            }
            fflush (logf);
            fflush (outf);
        }    
    }
	clock_t t = clock() - begt;
    fprintf (logf, "\nElapsed time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
	fclose (inf);
	fclose (outf);
	fclose (logf);
    printf ("Root finalized\n");
	MPI_Finalize();
	return 0;
}
