#include "main_params.hpp"

//////////////
//    MAIN
/////////////

int main(int argc, char *argv[])
{
	if (argc < 2){
		printf("Usage: gen [dimension]\n");
		return 0;
	}
    
    char *pch = argv[1];
    int dimension = strtol (pch, &pch, 10);
	if (dimension != DIM){
		printf ("Read file ERROR: dimension = %d, but must be %d\n", dimension, DIM);
		return 1;
	}
    
    char outfname[256]; // The output file name
    sprintf (outfname, "%dd%dv.gb", dimension, 1);
    FILE *outf = fopen (outfname, "wb"); // Open the output file
    if (outf == NULL){
        printf ("ERROR: Cann't open file %s\n", outfname);
        return 2;
    }
    printf ("Output is saved in %s\n", outfname);
    
    Tpoint w = 0, x = 1;
    for (int i = 0; i < dimension; i++, x <<= 2){
        fwrite (&w, sizeof(Tpoint), 1, outf);
        w += x;
    }
	fclose (outf);
	return 0;
}
