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

int main(int argc, char *argv[])
{
	if(argc < 3){
		printf ("Usage: merge [infile1] [infile2] [union]\n");
		return 0;
	}	

    char *pch;
	int nvert = strtol(argv[1], &pch, 10);
    if (pch == argv[1] || nvert < 2 || nvert > MAX_VERT){
        printf ("Wrong input file name\n");
        return 2;
    }
 		
	FILE *ina, *inb, *out;
	ina = fopen(argv[1], "rb");
	inb = fopen(argv[2], "rb");
	out = fopen(argv[3], "wb");
	if (ina == NULL || inb == NULL || out == NULL){
		printf ("ERROR: Cann't open one of files\n");
		return 1;
	}
	printf ("Merge %s and %s\n", argv[1], argv[2]);

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
