#include <stdio.h>
#include <string.h> // memcmp()
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <time.h>

#define MAX_DIM 4 // The maximum dimension of the Gale diagram (GD)
#define MAX_NUMBERS 1 // The maximum abs value of the coordinates of points in GD
#define MAX_VERT 32 // The bit length of the type uint32_t
#define MAX_FACET 32 // To avoid the memory overfull
#define LINE_SIZE 1024 // The maximum len of a line in the input file

using namespace std; 

int dimension, nverts;

// Get the dimension and the number of vertices from the filename
int get_dv(char *filename, int &dim, int &nv){
    char *pch = filename, *newpch;
    //pch = strchr(filename,'s');
    dim = strtol(pch, &newpch, 10);
    if (pch == newpch)
        return 1;
	if (dim < 2 || dim > MAX_DIM)
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

// Read Gale diagram from the file: dimension, nverts and vertices[][]
int read_diagram(FILE *inf, int &dim, int &nv, uint8_t *diagram){
	char buffer[LINE_SIZE]; // the buffer for reading lines from the file
    // Read the first line
	char *pch = fgets (buffer, LINE_SIZE, inf);
	if ( pch == NULL){
		printf ("Read file ERROR: unexpected EOF\n");
		return 1;
    }    

	// Read dimension
	dim = strtol (pch, &pch, 10);
	if (dim < 2 || dim > MAX_DIM){
		printf ("Read file ERROR: dimension = %d, but must be in [2, %d]\n", dim, MAX_DIM);
		return 2;
	}

	// Read dimension
	nv = strtol (pch, &pch, 10);
	if (nv < 1 || nv > MAX_VERT){
		printf ("Read file ERROR: number of vertices = %d, but must be in [1, %d]\n", nv, MAX_VERT);
		return 3;
	}
	
    // Read points from the file
	for (int i = 0; i < nv; i++){
		pch = fgets (buffer, LINE_SIZE, inf);
		if ( pch == NULL){
			printf ("Read file ERROR: unexpected EOF\n");
			return 1;
		}	
        diagram[i] = 0;
		for (int j = 0; j < dim; j++){
            char *new_pch;
			int x = strtol (pch, &new_pch, 10);
			if (pch == new_pch || fabs(x) > MAX_NUMBERS){
				printf ("Read file ERROR: in line %d: coordinates of points must be integers in [%d, %d]\n", i+2, -MAX_NUMBERS, MAX_NUMBERS);
				return 4;
			}
            pch = new_pch;
            diagram[i] |= ((x + 1) << (j*2));
		}
	}
	return 0;
}	

// Read Gale diagram from the file: dimension, nverts and vertices[][]
int read_diagram_short(FILE *inf, int &dim, int &nv, uint8_t *diagram){
	char buffer[LINE_SIZE]; // the buffer for reading lines from the file
    // Read the first line
	char *pch = fgets (buffer, LINE_SIZE, inf);
	if ( pch == NULL){
		printf ("Read file ERROR: unexpected EOF\n");
		return 1;
    }    

	// Read dimension
	dim = strtol (pch, &pch, 10);
	if (dim < 2 || dim > MAX_DIM){
		printf ("Read file ERROR: dimension = %d, but must be in [2, %d]\n", dim, MAX_DIM);
		return 2;
	}

	// Read dimension
	nv = strtol (pch, &pch, 10);
	if (nv < 1 || nv > MAX_VERT){
		printf ("Read file ERROR: number of vertices = %d, but must be in [1, %d]\n", nv, MAX_VERT);
		return 3;
	}
	
    // Read points from the file
	for (int i = 0; i < nv; i++){
		pch = fgets (buffer, LINE_SIZE, inf);
		if ( pch == NULL){
			printf ("Read file ERROR: unexpected EOF\n");
			return 1;
		}	
        char *new_pch;
        int x = strtol (pch, &new_pch, 3);
		if (pch == new_pch){
			printf ("Read file ERROR: in line %d: coordinates of points must be integers in [%d, %d]\n", i+2, -MAX_NUMBERS, MAX_NUMBERS);
			return 4;
		}
        pch = new_pch;
        diagram[i] = 0;
		for (int j = 0; j < dim; j++){
            diagram[i] |= (((x % 3) + 1) << (j*2));
            x /= 3;
		}
	}
	return 0;
}	

// Write vertices[][] to file
void write_diagram(FILE *outf, int dim, int nv, uint8_t *diagram){
	fprintf (outf, "b\n%d %d\n", dim, nv);
	for (int i = 0; i < nv; i++){
        uint8_t x = diagram[i];
		for (int j = 0; j < dim; j++, x >>= 2)
			fprintf (outf, " %2d", (int)(x & 0b11) - 1);
		fprintf (outf, "\n");
	}
	fprintf (outf, "e\n");
}


void recursive_permute(vector< vector<char> > &array, vector<char> &cur_perm, int d, int step){
    if (step == d-1){
        array.push_back(vector<char>(cur_perm));
        return;
    }
        
    int x = cur_perm[step];
    for(int i = step; i < d; i++){
        cur_perm[step] = cur_perm[i];
        cur_perm[i] = x;
        recursive_permute(array, cur_perm, d, step+1);
        cur_perm[i] = cur_perm[step];
        cur_perm[step] = x;
    }
}

vector< vector<char> > gen_permute(int d){
    int nperm = 1, i, j;
    for (i = 2; i <= d; i++)
        nperm *= i;
    vector< vector<char> > permut_array; // The output
    vector<char> cur_permute(d,0);
    for (i = 0; i < d; i++)
        cur_permute[i] = i;
    recursive_permute(permut_array, cur_permute, d, 0);
    return permut_array;
}

// Transform vertices to output with permutations of coordinates
void permute_coord(uint8_t *output, uint8_t *input, int dim, int nv, char *permut){
    for (int i = 0; i < nv; i++){
        output[i] = 0;
        uint8_t x = input[i];
        for (int j = 0; j < dim; j++, x >>= 2){
            int p = permut[j] * 2;
            output[i] |= ((x & 0b11) << p);
        }    
    }    
}

// Find the lexicographically minimal and write it in output
void lex_min(uint8_t *output, uint8_t *input, int dim, int nv){
    vector< vector<char> > all_permutations;
    all_permutations = gen_permute(dim);
    int npermut = all_permutations.size();
    output[0] = ~(uint8_t)0;
    uint8_t current1[MAX_VERT], current2[MAX_VERT];
    for (int i = 0; i < npermut; i++){ // Choose permutation of coordinates
        permute_coord(current1, input, dim, nv, all_permutations[i].data());
        int nreflections = (1 << dim);
/*        printf ("Perm:");
        for (int j = 0; j < dim; j++)
            printf ("%2d", all_permutations[i][j]);
        printf ("\n");
*/        for (int j = 0; j < nreflections; j++){ // Choose reflection
            uint8_t mask = 0; // mask for reflections
            uint8_t mask2 = 0b10; // one-bit reflection
            uint8_t mask3 = 1; // one-bit reflection
            int x, k;
//            printf ("Refl:");
            for (x = j, k = 0; k < dim; k++, x >>= 1, mask2 <<= 2){
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
    sprintf (outfname, "%dd%ds.gb", dimension, nverts);
	FILE *outf = fopen(outfname, "wb");
	if (outf == NULL){
		printf ("Write file ERROR: Cann't open the file %s\n", outfname);
		return 4;
	}
    vector<uint8_t *> all_diagrams;
    uint8_t diagram[MAX_VERT];
	clock_t t = clock();
    for (int n = 0; ;n++){
        if (n % 1000 == 0)
            printf (" %d", n);
        if (read_diagram_bin(inf, nverts, diagram))
            break;
        
        uint8_t *data = (uint8_t *)malloc(sizeof(uint8_t)*(nverts));
        if (data == NULL){
            printf ("ERROR: Out of memory\n");
            break;
        }
        lex_min(data, diagram, dimension, nverts);
        all_diagrams.push_back(data);
    }
	t = clock() - t;
    printf ("\nElapsed time: %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);
	fclose (inf);

//    nverts++;
    printf ("Totally = %d\n", all_diagrams.size());
    if (all_diagrams.size() == 0){
        fclose (outf);
        return 0;
    }
    
    printf ("Sort time:");
	t = clock();
    qsort (all_diagrams.data(), all_diagrams.size(), sizeof(uint8_t *), compare_vertices);
	t = clock() - t;
    printf (" %4.3f sec\n", ((float)t)/CLOCKS_PER_SEC);

    uint8_t **pt = all_diagrams.data();
    write_diagram_bin(outf, nverts, *pt);
    //write_diagram(outf, dimension, nverts, *pt);
    int writen = 1;
    for (int i = 1; i < all_diagrams.size(); i++, pt++){
        if (compare_vertices(pt, pt+1)){
            write_diagram_bin(outf, nverts, *(pt+1));
            //write_diagram(outf, dimension, nverts, *(pt+1));
            writen++;
        }    
        free(*pt);
        *pt = NULL;
    }
    free(*pt);
    *pt = NULL;
    printf ("Writen %d\n", writen);
	fclose (outf);
	return 0;
}
