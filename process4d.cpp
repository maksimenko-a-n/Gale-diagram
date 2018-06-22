#include "gale-diagram4d.hpp"
#include <string.h> // memcmp()
#include <time.h>
//#include <unistd.h> // sleep()
#include <vector> 
#include <mpi.h>

#define I_AM_FREE_TAG 2
#define TO_WORK_TAG 4
#define FOR_PROCESSING_TAG 8
#define RESULT_TAG 16
#define RESULT_READY_TAG 32
//#define BUNCH 64

// Extract the dimension and the number of vertices from the filename
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
void lex_min(uint8_t *output, uint8_t *input, int size, vector< vector<char> > *all_permutations){
    int npermut = all_permutations->size();
    output[0] = ~(uint8_t)0;
    uint8_t current1[MAX_VERT], current2[MAX_VERT];
    for (int i = 0; i < npermut; i++){ // Choose permutation of coordinates
        permute_coord(current1, input, size, (*all_permutations)[i].data());
        int nreflections = (1 << DIM);
        for (int j = 0; j < nreflections; j++){ // Choose reflection
            uint8_t mask = 0; // mask for reflections
            uint8_t mask2 = 0b10; // one-bit reflection
            uint8_t mask3 = 1; // mask for zeroes
            int x, k;
            for (x = j, k = 0; k < DIM; k++, x >>= 1, mask2 <<= 2){
                mask |= (x&1) * mask2;
                mask3 |= (mask3 << 2);
            }    
            for (k = 0; k < size; k++){ // Realize reflection
                uint8_t cur = current1[k];
                uint8_t y;
                y = ((cur >> 1) | (~cur)) & mask3; // Some magic for zeroes 0b01
                y |= y << 1;
                current2[k] = cur ^ (mask & y);
            }    
            int is_compare = 1; // Do compare current2 and output?
            // Bubble sort
            for (k = 0; k < size-1; k++){ // Sorting of current2 and comparing with output
                int m_min = k;
                int min = current2[m_min];
                for (int m = k+1; m < size; m++){
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
                    }    
                }    
                else
                    output[k] = min;
                // Swap    
                current2[m_min] = current2[k];
                current2[k] = min;
            }
            if (k == size-1){
                if ((!is_compare) || current2[k] < output[k]){
                    output[k] = current2[k];
                }    
            }
        }
    }
};

void process_one_diagram(int old_nverts, uint8_t *diagram, uint8_t *buffer, int &num, Gale_diagram &gale, vector< vector<char> > *all_permutations){
    num = 0;
    gale.nverts = old_nverts;
    gale.convert_bin2vertices(diagram);
    int old_nfacets = gale.find_facets(); // Find facets
    int ispoly = gale.is_polytope();
    //if (old_nfacets >= MAX_FACET || (ispoly && (old_nfacets > MAX_FACET-2)))
    //    return;
    // Copy the data from the gale
    int old_cofacets_num[DIM+2];
    memcpy(old_cofacets_num, gale.cofacets_num, sizeof(int)*(DIM+2));
    //if (old_cofacets_num[2] != 0 || old_cofacets_num[4] != 0 || old_cofacets_num[5] != 0)
    //    return;
    
    //printf ("Process: nv = %d, nf = %d\n", old_nverts, old_nfacets);

    uint8_t zero = 0; // The zero-point
    for (int k = 0; k < DIM; k++)
        zero |= ((uint8_t)1 << (k*2));
    // is_old_vertices[i] == 0  <=>  i is a new vertex and is not zero
    char is_old_vertices[1 << 8];
    memset (is_old_vertices, 0, 1 << 8);
    is_old_vertices[zero] = 1;
    for (int i = 0; i < old_nverts; i++)
        is_old_vertices[diagram[i]] = 1;

    // Add one vertex
    gale.nverts = old_nverts + 1;
    uint8_t *cur_diagram = buffer;
    for (int i = 0; i < Nall; i++){ // Look over all possible points
        diagram[old_nverts] = 0; 
        int x = i;
        for (int j = 0; j < DIM; j++, x /= 3){
            uint8_t y = x % 3;
            diagram[old_nverts] |= (y << (j*2));
            gale.vertices[old_nverts][j] = y - 1;
        }    
        if (is_old_vertices[diagram[old_nverts]])
            continue; // If we have choosen an old point or zero
        // Process
        gale.nfacets = old_nfacets; // Init the number of facets
        memcpy(gale.cofacets_num, old_cofacets_num, sizeof(int)*(DIM+2));
        int nfacets = gale.facets_with_last_vert (); // Find cofacets only with the new point
        int ridges = edges_number(nfacets, gale.nverts, gale.facet_vertex); // Find ridges
        int ispoly = gale.is_polytope();
        int edges = 0;
        if (ispoly){ // Find the number of edges
            uint32_t vertex_facet[MAX_VERT];
            transpose(nfacets, gale.nverts, gale.facet_vertex, vertex_facet);
            edges = edges_number(gale.nverts, nfacets, vertex_facet); // Find edges
        }
        int is_2neighborly = (edges*2 == old_nverts*(old_nverts+1));
        int is_dual2neighborly = (ispoly && (ridges*2 == nfacets*(nfacets-1)));
        int freev = gale.nverts - gale.vertices_in_facets(); // The number of vertices not in facets
        int max_freev = (gale.nverts <= DIM ? gale.nverts : (2 * DIM - gale.nverts));
        //max_freev = (max_freev < 0 ? 0 : max_freev);
        max_freev = (max_freev < 1 ? 1 : max_freev);
        int is_write = 1; // Is this diagram good for our purposes?
        is_write = (nfacets <= MAX_FACET);
        is_write &= (freev <= max_freev);
        if (!is_2neighborly)
            is_write &= (ispoly && nfacets < MAX_FACET-1) || (nfacets < MAX_FACET-2);
        is_write &= ispoly; // Attention!
        int big = gale.cofacets_num[DIM] + gale.cofacets_num[DIM+1];
        //is_write &= (big <= small);
        //is_write &= (big <= 0);
        //is_write &= (gale.cofacets_num[2] == 0);
        if (is_write){
            lex_min(cur_diagram+1, diagram, old_nverts+1, all_permutations);  // Normalize points in the diagram
            //memcpy(cur_diagram+1, diagram, (old_nverts+1)*sizeof(uint8_t));
            cur_diagram[old_nverts+2] = nfacets;
            cur_diagram[old_nverts+3] = edges;
            cur_diagram[old_nverts+4] = ridges;
            cur_diagram[0] = 0;
            if (is_2neighborly)
                cur_diagram[0] = 1;
            if (is_dual2neighborly)
                cur_diagram[0] = 2;
            num++;
            cur_diagram += old_nverts+5;
        }
    }
}

int write_results(FILE *outf, FILE *out2f, FILE *out2ftxt, int nverts, uint8_t *buffer, int num, long long &savednum, int &min2facets, int &minf, int &maxf, long long &polynum, long long *distribution){
    uint8_t *diagram = buffer;
    for (int n = 0; n < num; n++, diagram += nverts + 4){
        int nfacets = diagram[nverts+1];
        int edges = diagram[nverts+2];
        if (diagram[0] != 1){ // Will be used for the next generation
            minf = minf > nfacets ? nfacets : minf;
            maxf = maxf < nfacets ? nfacets : maxf;
            if (edges > 0) polynum++; // The total number of polytopes
            int s = fwrite (diagram+1, sizeof(uint8_t), nverts, outf);
            if (s < nverts){
                printf ("Write file ERROR: Cann't write %d bytes\n", nverts);
                return 1;
            }
            savednum++; // The number of saved diagrams
            distribution[nfacets] += 1;
        }    
        if (diagram[0] != 0){ // Write special polytopes 
            // Write in the binary format
            int s = fwrite (diagram+1, sizeof(uint8_t), nverts, out2f);
            if (s < nverts){
                printf ("Write file ERROR: Cann't write %d bytes\n", nverts);
                return 2;
            }
            fflush (out2f);
             // Write in the readable format
            int ridges = diagram[nverts+3];
            if (edges*2 == nverts*(nverts-1)){
                // For a 2-neighborly diagram
                if (min2facets > nfacets)
                    min2facets = nfacets;
                fprintf (out2ftxt, "b 2\n");
            }    
            if (edges > 0 && (ridges*2 == nfacets*(nfacets-1))){
                // For a dual 2-neighborly diagram
                fprintf (out2ftxt, "b dual\n");
            }    
            for (int i = 1; i <= nverts; i++){
                uint8_t x = diagram[i];
                for (int j = 0; j < DIM; j++, x >>= 2)
                    fprintf (out2ftxt, " %2d", (x & 0b11) - 1);
                fprintf (out2ftxt, "\n");
            }
            fprintf (out2ftxt, "f=%d; e=%d\nr=%d\ne\n", nfacets, edges, ridges);
            //fprintf (out2ftxt, "v = %d, maxv = %d, multiplicity = %d\n", freev, max_freev, multiplicity);
            //write_incmatrix (out2ftxt, nfacets, gale.nverts, gale.facet_vertex);
            fflush (out2ftxt);
        }
    }
    return 0;
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
	MPI_Init(&argc, &argv); // Init MPI
	if (argc < 2){
		printf("Usage: proc [file name]\n");
		MPI_Finalize();
		return 1;
	}
   
    int nverts;
    // Extract the dimension and the number of vertices from the filename
    int get_value = get_dv(argv[1], nverts);
    if (get_value){
        printf ("Wrong file name: get_dv() error %d\n", get_value);
		MPI_Finalize();
        return 2;
    }

    uint8_t diagram[MAX_VERT];
	uint8_t buffer[Nall * (nverts+5)]; // Every block has (nverts+5) bytes
    // The first byte -- type: ordinary diagram (0) or special one (1)
    // [1, nverts] -- points; nverts+1 -- nfacets; 
    // nverts+2 -- edges (or 0 if it is not a polytope); nverts+3 -- ridges
    
	// Get the rank of the curent thread and the number of all threads
	int rank, nthreads;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    int rcv_get;
	MPI_Status status;
    MPI_Request request;

	if (rank != 0) {
        Gale_diagram gale; // The Gale diagram
        vector< vector<char> > all_permutations;
        all_permutations = gen_permute();
		int size = 0;
		MPI_Isend(&size, 1, MPI_INT, 0, RESULT_READY_TAG, MPI_COMM_WORLD, &request);
		while(1) {
			MPI_Recv(&rcv_get, 1, MPI_INT, 0, TO_WORK_TAG, MPI_COMM_WORLD, &status);
			if(rcv_get < 1) break;
			MPI_Recv(diagram, rcv_get*nverts*sizeof(uint8_t), MPI_BYTE, 0, FOR_PROCESSING_TAG, MPI_COMM_WORLD, &status);
            process_one_diagram(nverts, diagram, buffer, size, gale, &all_permutations); // The main process
			MPI_Isend(&size, 1, MPI_INT, 0, RESULT_READY_TAG, MPI_COMM_WORLD, &request);
            if (size > 0)
                MPI_Send(buffer, size*(nverts+5)*sizeof(uint8_t), MPI_BYTE, 0, RESULT_TAG, MPI_COMM_WORLD);
		}
		//printf ("%d's thread fin  ", rank);
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
    fprintf (logf, "\n%sAdd and Proc %s\nOn %d threads\n", ctime (&rawtime), argv[1], nthreads);
    //struct tm * timeinfo = localtime ( &rawtime );
    //fprintf (logf, "\n%s\nAdd and Proc %s\n", asctime (timeinfo), argv[1]);
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

    // Open the input file
	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		MPI_Finalize();
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    long long inputlen = ftell(inf)/nverts; // The total number of diagrams
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    fprintf (logf, "Input: %d diagrams\n", inputlen);
    // Open the output file
    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%dp.gb", DIM, nverts+1);
    fprintf (logf, "Write in %s\n", outfname);
	FILE *outf = fopen(outfname, "wb");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		MPI_Finalize();
		return 4;
	}
    // Open the binary output file for 2-neighborly polytopes
    sprintf (outfname, "%dd%d-2ng.gb", DIM, nverts+1);
	FILE *out2f = fopen(outfname, "wb");
	if (out2f == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		MPI_Finalize();
		return 5;
	}
    // Open the txt-file for 2-neighborly polytopes
    sprintf (outfname, "%dd%d-2ng.gs", DIM, nverts+1);
	FILE *out2ftxt = fopen(outfname, "w");
	if (out2ftxt == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		MPI_Finalize();
		return 5;
	}

    long long savednum = 0;
    int min2facets = MAX_FACET+1; // Minimum number of facets for a 2-neighborly polytope
    int minf = MAX_FACET+1, maxf = 0; // Minimum and maximum number of facets
    long long polynum = 0; // The number of saved polytopes
    long long distribution[MAX_FACET+1]; // distribution[i] -- the num of diagrams with i facets
    memset(distribution, 0, (MAX_FACET+1)*sizeof(long long));
	clock_t begt = clock();
	// READ AND EVALUATE
	long long sent = 0; // amount of input
	long long received = 0; // amount of output
//	int rcv_tag = MPI_ANY_TAG;
	long long to_slaves = 0, from_slaves = 1 - nthreads, end_of_file = 0;
	while(!end_of_file || to_slaves != from_slaves) {
		MPI_Recv(&rcv_get, 1, MPI_INT, MPI_ANY_SOURCE, RESULT_READY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG != RESULT_READY_TAG){
			printf ("ERROR: Unexpected status.MPI_TAG = %d\n", status.MPI_TAG);
			break;
		}
        from_slaves++;
		if (rcv_get > 0){ // Receive results from worker
            MPI_Recv(buffer, rcv_get*(nverts+5)*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
            write_results(outf, out2f, out2ftxt, nverts+1, buffer, rcv_get, savednum, min2facets, minf, maxf, polynum, distribution);
            received += rcv_get;
            fflush(outf);
        }
        if (from_slaves%10000000 == 0)
            fprintf (logf, " Received %d (%d calls)", received, from_slaves);
        if (end_of_file){
            int send = 0;
            MPI_Isend(&send, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD, &request);
            continue;
        }    
		// Read diagram from file and send to a thread
		int size;
		size = fread (diagram, sizeof(uint8_t), nverts, inf);
		if(size < nverts){ // The end of input file
            size = 0;
            MPI_Isend(&size, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD, &request);
			end_of_file = 1;
            continue;
		}
        // Send diagram to worker
        size = 1;
        MPI_Isend(&size, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD, &request);
		MPI_Send(diagram, size*nverts*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, FOR_PROCESSING_TAG, MPI_COMM_WORLD);
		to_slaves++;
		sent += size;
        if (sent == 100000 || sent % 10000000 == 0){
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
	fclose (inf);
	fclose (outf);
    fprintf (logf, "\nElapsed time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
    fprintf (logf, "New number of diagrams = %d\n", savednum);
    fprintf (logf, "min facets = %d, max facets = %d, polytopes num = %d\n", minf, maxf, polynum);
    fprintf (logf, "Facets distr:");
    for (int i = 0; i <= MAX_FACET; i++)
        fprintf (logf, " d[%d] = %d,", i, distribution[i]);
    fprintf (logf, "\n");
    //printf ("New number of diagrams = %d\n", savednum);
    if (min2facets <= MAX_FACET){
        fprintf (logf, "Min facets for 2-nghb = %d\n", min2facets);
        fprintf (out2ftxt, "Min facets = %d\n", min2facets);
    }    
	fclose (out2ftxt);
	fclose (logf);
    //printf ("Root finalized\n");
	MPI_Finalize();
	return 0;
}
