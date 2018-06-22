#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memcmp()
#include <stdint.h>  // uint32_t

#define MAX_VERT 20
#define DIM 4
const uint8_t maxvalue = 255;

// Get the dimension and the number of vertices from the filename
int get_dv(char *filename, int &nv, char &letter){
    char *pch = filename, *newpch;
    nv = strtol(pch, &newpch, 10);
    if (pch == newpch)  return 1;
	if (nv < 2 || nv > MAX_VERT)  return 2;
	letter = *newpch + 1;
    return 0;
}

// Read one diagram from f to diagram
// If cann't, then put maxvalue into diagram[0]
void read_diagram(FILE *f, uint8_t *diagram, int nvert){
    if (fread(diagram, sizeof(uint8_t), nvert, f) < nvert)  diagram[0] = maxvalue;
}

int merge(const char *infile1, const char *infile2, const char *outfile){
   char *pch;
	int nvert = strtol(infile1, &pch, 10);
    if (pch == infile1 || nvert < 2 || nvert > MAX_VERT){
        printf ("Wrong input file name\n");
        return 2;
    }
 		
	FILE *ina, *inb, *out;
	ina = fopen(infile1, "rb");
	inb = fopen(infile2, "rb");
	out = fopen(outfile, "wb");
	if (ina == NULL || inb == NULL || out == NULL){
		printf ("ERROR: Cann't open one of files\n");
		return 1;
	}
	printf ("Merge %s and %s. Result is saved in %s.\n", infile1, infile2, outfile);

	uint8_t a[MAX_VERT], b[MAX_VERT];
    // The first byte is never equal to maxvalue
    // a[0] == maxvalue is equivalent to the end of file ina
    read_diagram(ina, a, nvert);
    read_diagram(inb, b, nvert);
	
	while(a[0] < maxvalue || b[0] < maxvalue) {
		int cmp_result = memcmp(a, b, nvert*sizeof(uint8_t));
		if(cmp_result > 0) {
			fwrite(b, sizeof(uint8_t), nvert, out);
            read_diagram(inb, b, nvert);
		}
		else { // if(cmp_result <= 0)
			fwrite(a, sizeof(uint8_t), nvert, out);
            read_diagram(ina, a, nvert);
			if(cmp_result == 0)
                read_diagram(inb, b, nvert);
		}
	}
	
	fclose(ina);
	fclose(inb);
	fclose(out);
    return 0;
}

int merge_all(int nverts, char &letter, int &nfiles){
    char in1[32], in2[32], out[32];
    int half = nfiles/2;
    for (int i = 0; i < half; i++){
        sprintf (in1, "%d%c%d", nverts, letter, i);
        sprintf (in2, "%d%c%d", nverts, letter, half + i);
        sprintf (out, "%d%c%d", nverts, letter+1, i);
        //printf ("From %s and %s to %s\n", in1, in2, out);
        merge(in1, in2, out);    
    }
    if (nfiles % 2 == 1){
        sprintf (in1, "%d%c%d", nverts, letter, half*2);
        sprintf (out, "%d%c%d", nverts, letter+1, half);
        printf ("Rename %s to %s\n", in1, out);
        rename (in1, out);
        nfiles = half+1;
    }
    else
        nfiles = half;
    letter++;
    return 0;
}

int main(int argc, char *argv[])
{
	if(argc < 2){
		printf ("Usage: merge [nvertices] [nfiles]\n");
		return 0;
	}
    
    char *pch = argv[1], *newpch;
    int nverts = strtol(pch, &newpch, 10);
    if (pch == newpch)  return 1;
	if (nverts < 2 || nverts > MAX_VERT)  return 2;
    pch = argv[2];
    int nfiles = strtol(pch, &newpch, 10);
    if (pch == newpch)  return 1;
	if (nfiles < 2 || nfiles > 100)  return 2;
	char letter = 'a';
    
    while (nfiles > 1)
        merge_all(nverts, letter, nfiles);
    return 0;
}
