#include "main_params.hpp"


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
// If cann't, then put MAX_Tpoint into diagram[0]
inline void read_diagram(FILE *f, Tpoint *diagram, int size){
    if (fread(diagram, sizeof(Tpoint), size, f) < size)  diagram[0] = MAX_Tpoint;
}

// Write diagram to file f
inline void write_diagram(FILE *f, Tpoint *diagram, int size){
    fwrite(diagram, sizeof(Tpoint), size, out);
}

int merge(const char *infile1, const char *infile1, const char *outfile){
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
	printf ("Merge %s and %s. Result is saved in %s.\n", infile1, infile1, outfile);

	Tpoint a[MAX_VERT], b[MAX_VERT];
    // The first byte is never equal to MAX_Tpoint
    // a[0] == MAX_Tpoint is equivalent to the end of file ina
    read_diagram(ina, a, nvert);
    read_diagram(inb, b, nvert);
	
	while(a[0] < MAX_Tpoint || b[0] < MAX_Tpoint) {
		int cmp_result = memcmp(a, b, nvert * sizeof(Tpoint));
		if(cmp_result > 0) {
            write_diagram(out, b, nvert);
            read_diagram(inb, b, nvert);
		}
		else { // if(cmp_result <= 0)
            write_diagram(out, a, nvert);
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

int main(int argc, char *argv[])
{
	if(argc < 3){
		printf ("Usage: merge [infile1] [infile2] [union]\n");
		return 0;
	}	

    return merge(argv[1], argv[2], argv[3]);
}
