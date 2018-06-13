#include "gale-diagram4d.hpp"
#include <string.h> // memcmp()
#include <time.h>
#include <unistd.h> // sleep()
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

void process_one_diagram(int old_nverts, uint8_t *diagram, uint8_t *buffer, int &num, Gale_diagram &gale){
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
        int big = gale.cofacets_num[DIM] + gale.cofacets_num[DIM+1];
        //is_write &= (big <= small);
        //is_write &= (big <= 0);
        //is_write &= (gale.cofacets_num[2] == 0);
        if (is_write){
            memcpy(cur_diagram+1, diagram, (old_nverts+1)*sizeof(uint8_t));
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

int write_results(FILE *outf, FILE *out2f, int nverts, uint8_t *buffer, int num, int &newd, int &min2facets, int &minf, int &maxf, unsigned long int &polynum){
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
            newd++; // The number of new diagrams
        }    
        if (diagram[0] != 0){ // Write in the readable format
            int ridges = diagram[nverts+3];
            if (edges*2 == nverts*(nverts-1)){
                // For a 2-neighborly diagram
                if (min2facets > nfacets)
                    min2facets = nfacets;
                fprintf (out2f, "b 2\n");
            }    
            if (edges > 0 && (ridges*2 == nfacets*(nfacets-1))){
                // For a dual 2-neighborly diagram
                fprintf (out2f, "b dual\n");
            }    
            for (int i = 1; i <= nverts; i++){
                uint8_t x = diagram[i];
                for (int j = 0; j < DIM; j++, x >>= 2)
                    fprintf (out2f, " %2d", (x & 0b11) - 1);
                fprintf (out2f, "\n");
            }
            fprintf (out2f, "f=%d; e=%d\nr=%d\ne\n", nfacets, edges, ridges);
            //fprintf (out2f, "v = %d, maxv = %d, multiplicity = %d\n", freev, max_freev, multiplicity);
            //write_incmatrix (out2f, nfacets, gale.nverts, gale.facet_vertex);
            fflush (out2f);
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

	if (rank != 0) {
        Gale_diagram gale; // The Gale diagram
		while(1) {
			MPI_Request request;
			MPI_Isend(&rank, 1, MPI_INT, 0, I_AM_FREE_TAG, MPI_COMM_WORLD, &request);
			MPI_Recv(&rcv_get, 1, MPI_INT, 0, TO_WORK_TAG, MPI_COMM_WORLD, &status);
			if(rcv_get < 1){
				//if (rcv_get < 0)
					break;
				//else{	
				//	sleep(1); // Delay
				//	continue;	
				//}	
			}	

			MPI_Recv(diagram, rcv_get*nverts*sizeof(uint8_t), MPI_BYTE, 0, FOR_PROCESSING_TAG, MPI_COMM_WORLD, &status);
			int size;
            process_one_diagram(nverts, diagram, buffer, size, gale); // The main process
			MPI_Send(&size, 1, MPI_INT, 0, RESULT_READY_TAG, MPI_COMM_WORLD);
            if (size > 0)
                MPI_Send(buffer, size*(nverts+5)*sizeof(uint8_t), MPI_BYTE, 0, RESULT_TAG, MPI_COMM_WORLD);
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
    fprintf (logf, "\n%sAdd and Proc %s\nOn %d threads\n", ctime (&rawtime), argv[1], nthreads);
    //struct tm * timeinfo = localtime ( &rawtime );
    //fprintf (logf, "\n%s\nAdd and Proc %s\n", asctime (timeinfo), argv[1]);
    fprintf (logf, "dim = %d, nverts = %d\n", DIM, nverts);

    // Open the input file
	FILE *inf = fopen(argv[1], "rb");
	if (inf == NULL){
		printf ("Read file ERROR: Cann't open the file %s\n", argv[1]);
		return 3;
	}
    fseek(inf, 0, SEEK_END); // Go to end of file
    unsigned long inputlen = ftell(inf)/nverts; // The total number of diagrams
    fseek(inf, 0, SEEK_SET); // Go to begin of file
    fprintf (logf, "Input: %d diagrams\n", inputlen);
    // Open the output file
    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%dp.gb", DIM, nverts+1);
    fprintf (logf, "Write in %s\n", outfname);
	FILE *outf = fopen(outfname, "wb");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 4;
	}
    // Open the output file for 2-neighborly polytopes
    sprintf (outfname, "%dd%d-2ng.gs", DIM, nverts+1);
	FILE *out2f = fopen(outfname, "w");
	if (out2f == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 5;
	}

    int savednum = 0;
    int min2facets = MAX_FACET+1; // Minimum number of facets for a 2-neighborly polytope
    int minf = MAX_FACET+1, maxf = 0; // Minimum and maximum number of facets
    unsigned long int  polynum = 0; // The number of saved polytopes
	clock_t begt = clock();
	// READ AND EVALUATE
	long long sent = 0; // amount of input
	long long recieved = 0; // amount of output
	int rcv_tag = MPI_ANY_TAG;
	long to_slaves = 0, from_slaves = 0, end_of_file = 0;
	while(!end_of_file || to_slaves != from_slaves) {
		MPI_Recv(&rcv_get, 1, MPI_INT, MPI_ANY_SOURCE, rcv_tag, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == I_AM_FREE_TAG){
			if (status.MPI_SOURCE != rcv_get)
				printf ("WARNING: status.MPI_SOURCE != process_num\n");
			
			// Read string from file and send to a thread
			int size;
			size = fread (diagram, sizeof(uint8_t), nverts, inf);
			if(size < nverts)
				size = 0;
			else
				size = 1;

			if(size == 0){
				MPI_Send(&size, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD);
				rcv_tag = RESULT_READY_TAG;
				end_of_file = 1;
			}
			else{
				MPI_Send(&size, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD);
				MPI_Send(diagram, size*nverts*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, FOR_PROCESSING_TAG, MPI_COMM_WORLD);
				to_slaves++;
				sent += size;
                if (sent == 10000 || sent % 100000 == 0){
                    fprintf (logf, " Sent %d", sent);
                    if (sent != 0){
                        clock_t curt = clock();
                        int sec = (curt - begt)*(inputlen - sent) / (sent*CLOCKS_PER_SEC);
                        int min = sec / 60;
                        int hour = min / 60;
                        fprintf (logf, " (left %d:%02d:%02d)", hour, min%60, sec%60);
                    }
                    fflush (logf);
                }    
			}
        }    
		else {
			if (status.MPI_TAG != RESULT_READY_TAG){
				printf ("ERROR: Unexpected status.MPI_TAG = %d\n", status.MPI_TAG);
				break;
			}
            from_slaves++;
			if (rcv_get > 0){	
                MPI_Recv(buffer, rcv_get*(nverts+5)*sizeof(uint8_t), MPI_BYTE, status.MPI_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
                write_results(outf, out2f, nverts+1, buffer, rcv_get, savednum, min2facets, minf, maxf, polynum);
                recieved += rcv_get;
                fflush(outf);
                if (recieved%1000000 == 0){
                    fprintf (logf, " Recieve %d (%d calls)", recieved, from_slaves);
                }
            }
            if (end_of_file){
                MPI_Request request;
                int send = 0;
                MPI_Isend(&send, 1, MPI_INT, status.MPI_SOURCE, TO_WORK_TAG, MPI_COMM_WORLD, &request);
            }    
		}
    }
	clock_t t = clock() - begt;
    fprintf (logf, "\nElapsed time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
    fprintf (logf, "New number of diagrams = %d\n", savednum);
    fprintf (logf, "min facets = %d, max facets = %d, polytopes num = %d\n", minf, maxf, polynum);
    //printf ("New number of diagrams = %d\n", savednum);
	fclose (inf);
	fclose (outf);
    if (min2facets <= MAX_FACET){
        fprintf (logf, "Min facets for 2-nghb = %d\n", min2facets);
        fprintf (out2f, "Min facets = %d\n", min2facets);
    }    
	fclose (out2f);
	fclose (logf);
    printf ("Root finalized\n");
	MPI_Finalize();
	return 0;
}
