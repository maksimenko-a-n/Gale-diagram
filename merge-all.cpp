#include "main_params.hpp"

// Read one diagram from f to diagram
// If cann't, then put MAX_Tpoint into diagram[0]
inline void read_diagram(FILE *f, Tpoint *diagram, int size){
    if (fread(diagram, sizeof(Tpoint), size, f) < size)  diagram[0] = MAX_Tpoint;
}

// Write diagram to file f
inline void write_diagram(FILE *f, Tpoint *diagram, int size){
    fwrite(diagram, sizeof(Tpoint), size, out);
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

	Tpoint a[MAX_VERT], b[MAX_VERT];
    // The first byte is never equal to MAX_Tpoint
    // a[0] == MAX_Tpoint is equivalent to the end of file ina
    read_diagram(ina, a, nvert);
    read_diagram(inb, b, nvert);
	
	while(a[0] < MAX_Tpoint || b[0] < MAX_Tpoint) {
		int cmp_result = memcmp(a, b, nvert*sizeof(Tpoint));
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

// Merge all sorted files [nv]a0, ..., [nv]a[nfiles-1] in one file
// Because files are huge, it will be done in log(nfiles) stages
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
	if (nfiles < 2 || nfiles > 128)  return 2;
	char letter = 'a'; // The indicator of a stage
    
    while (nfiles > 1)
        merge_all(nverts, letter, nfiles);
    return 0;
}
