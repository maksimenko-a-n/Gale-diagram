#include <stdio.h>
#include <string.h> // memcmp()
#include <stdlib.h>
#include <vector>  
#include <math.h>
#include <stdint.h>  // uint32_t
#include <float.h>  // DBL_EPSILON
#include <time.h>

#define MAX_DIM 8 // The maximum dimension of the Gale diagram (GD)
#define MAX_NUMBERS 2 // The maximum abs value of the coordinates of points in GD
#define MAX_VERT 32 // The bit length of the type uint32_t
#define MAX_FACET 32 // To avoid the memory overfull
#define LINE_SIZE 1024 // The maximum len of a line in the input file

using namespace std; 

int dimension, nverts;
const uint32_t transform[5] = {0,0b1,0b11,0b1110,0b1111};
const int map[16] = {-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,2};

// Read Gale diagram from the file: dimension, nverts and vertices[][]
int read_diagram(FILE *inf, int &dim, int &nv, uint32_t *diagram){
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
            diagram[i] |= (transform[x + 2] << (j*4));
		}
	}
	return 0;
}	

// Write vertices[][] to file
void write_diagram(FILE *outf, int dim, int nv, uint32_t *diagram){
	fprintf (outf, "begin\n%d %d\n", dim, nv);
	for (int i = 0; i < nv; i++){
        uint32_t x = diagram[i];
		for (int j = 0; j < dim; j++, x >>= 4)
			fprintf (outf, " %2d", map[x & 0b1111]);
		fprintf (outf, "\n");
	}
	fprintf (outf, "end\n");
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
void permute_coord(uint32_t *output, uint32_t *input, int dim, int nv, char *permut){
    for (int i = 0; i < nv; i++){
        output[i] = 0;
        uint32_t x = input[i];
        for (int j = 0; j < dim; j++, x >>= 4){
            int p = permut[j] * 4;
            output[i] |= ((x & 0b1111) << p);
        }    
    }    
}

// Normalize binary input
void normalize_bin(uint32_t *input, int dim, int nv){
    int i, j;
    for (i = 0; i < nv; i++){
        uint32_t x = input[i];
        uint32_t mask = 0b1001;
        for (j = 0; j < dim; j++, x >>= 4, mask <<= 4){
            if ((x & 0b1111) == 0b1100)
                input[i] -= mask;
        }
    }
}

// Find the lexicographically minimal and write it in output
void lex_min(uint32_t *output, uint32_t *input, int dim, int nv){
    vector< vector<char> > all_permutations;
    all_permutations = gen_permute(dim);
    int npermut = all_permutations.size();
    output[0] = ~(uint32_t)0;
    uint32_t current1[MAX_VERT], current2[MAX_VERT];
    for (int i = 0; i < npermut; i++){ // Choose permutation of coordinates
        permute_coord(current1, input, dim, nv, all_permutations[i].data());
        int nreflections = (1 << dim);
/*        printf ("Perm:");
        for (int j = 0; j < dim; j++)
            printf ("%2d", all_permutations[i][j]);
        printf ("\n");
*/        for (int j = 0; j < nreflections; j++){ // Choose reflection
            uint32_t mask = 0; // mask for reflections
            uint32_t mask2 = 0b1111; // one-bit reflection
            uint32_t mask3 = 1; // one-bit reflection
            int x, k;
//            printf ("Refl:");
            for (x = j, k = 0; k < dim; k++, x >>= 1, mask2 <<= 4){
                mask |= (x&1) * mask2;
//                printf ("%2d", x&1);
                mask3 |= (mask3 << 4);
            }    
//            printf (": ");
            for (k = 0; k < nv; k++){ // Realize reflection
                uint32_t cur = current1[k];
                uint32_t y = cur & (cur >> 1);
                y = ((y >> 2) | (~y)) & mask3; // Some magic for zeroes 0b0011
                y |= y << 1;
                y |= y << 2;
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
  return memcmp(*(uint32_t**)a, *(uint32_t**)b, sizeof(uint32_t) * nverts);
}


//////////////
//    MAIN
/////////////

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gale [file name]\n");
		return 0;
	}
    
	FILE *inf = fopen(argv[1], "r");
	if (inf == NULL){
		printf ("ERROR: Cann't open the file %s\n", argv[1]);
		return 1;
	}
	char outfname[256]; // The output file name
	sprintf (outfname, "%s.out", argv[1]);
    FILE *outf = fopen (outfname, "w"); // Open the output file
	if (outf == NULL){
		printf ("ERROR: Cann't open file %s\n", outfname);
		return 1;
	}
	printf ("Output is saved in %s\n", outfname);


	char buffer[LINE_SIZE]; // the buffer for reading lines from the file
    const char begin[] = "begin";
	char *pch;
    uint32_t diagram[MAX_VERT];
    vector<uint32_t *> all_diagrams;
    uint32_t *data;
    nverts = 0;
    for (int n = 0; ;){
        pch = fgets (buffer, LINE_SIZE, inf);
        if ( pch == NULL)
            break;
        if (memcmp(begin, pch, sizeof(begin)-1) == 0){
            n++;
            printf (" %2d", n);
            int d, nv;
            int err = read_diagram(inf, d, nv, diagram);
            if (err)
                break;
            if (nverts == 0){
                dimension = d;
                nverts = nv;
            }
            else{
                if (nv != nverts || d != dimension){
                    printf ("Read file ERROR in the %d-th diagram: ", n);
                    if (nv != nverts)
                        printf ("The number of vertices %d differs from nvert = %d\n", nv, nverts);
                    else
                        printf ("The dim %d differs from the dimension = %d\n", d, dimension);
                    break;
                }    
            }
            // Copy data from diagram to all_diagrams[i]
            int M = pow(3, dimension);
            for (int i = 0; i < M; i++){
                diagram[nverts] = 0;
                int x = i;
                int s = 0;
                for (int j = 0; j < dimension; j++){
                    s |= (x % 3) - 1;
                    diagram[nverts] |= (transform[(x % 3) + 1] << (j*4));
                    x /= 3;
                }
                if (s != 0){
                    data = (uint32_t *)malloc(sizeof(uint32_t)*(nverts+1));
                    if (data == NULL){
                        printf ("ERROR: Out of memory\n");
                        break;
                    }
                    lex_min(data, diagram, dimension, nverts+1);
                    //memcpy(data, diagram, sizeof(uint32_t)*nverts);
                    all_diagrams.push_back(data);
                }
            }
        }
    }
	fclose (inf);

    printf ("Before sort\n");
    nverts += 1;
    qsort (all_diagrams.data(), all_diagrams.size(), sizeof(uint32_t *), compare_vertices);

    uint32_t **pt = all_diagrams.data();
    write_diagram(outf, dimension, nverts, *pt);
    int writen = 1;
    for (int i = 1; i < all_diagrams.size(); i++, pt++){
        if (compare_vertices(pt, pt+1)){
            write_diagram(outf, dimension, nverts, *(pt+1));
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
