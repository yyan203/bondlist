// This is an example c++ (really c-style c++) code for bond list creation for lammps readdata, for Dreiding potential, the list is required by Moltemplate to get 3body and 4body interaction
//
// The basic components are essentially the same as the perl script lmps.analysis.pl. However, note the differences in implementing and speed of these two languages.

//Tasks
//(1) Setup the neighbors of each cell (Expand function SimulationBox::SetupCell())
//(2) Throw atoms in Each Cell (Write function SimulationBox::DistributeAtoms())
//(3) Calculate RDF using LinkCell

//Memory check in CCNI
// valgrind --leak-check=full ./md.cell dump.melt haha 1 3

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include <vector>

#include "simulationbox.h"
#include "readdata.h"


#define VERBOSE 1
#define MAXDUMP 5

////////////////////
// Main functions //
////////////////////

int main(int argc, char **argv)
{
    char inputfilename[100], outputfilename[100];

    SimulationBox **mysystem;
    Atom **myatom; //convenient pointer
    float Lx, Ly, Lz;
    float cellsize;

    /************* Command line parameter processing ************/

    if(argc != 7) {
        // Provide the correct syntax of this analysis code
        printf("Correct syntax: md.analysis.cxx inputfile outputfile cellsize BoxLengthx Boxlengthy Boxlengthz\n");
        exit(0);
    }
    else {
        //Note that argv[0] is "md.analysis.cxx"
        sscanf(argv[1], "%s", inputfilename);
        sscanf(argv[2], "%s", outputfilename);
        sscanf(argv[3], "%f", &cellsize);
        sscanf(argv[4], "%f", &Lx);
        sscanf(argv[5], "%f", &Ly);
        sscanf(argv[6], "%f", &Lz);
    }

    //Initialize Simulation Boxes
    mysystem = new SimulationBox *[MAXDUMP];
    for(i=0;i<MAXDUMP;i++) {mysystem[i] = new SimulationBox();}
    printf("Initialization of Simulation Boxes has been finished! %d\n",i);

    //Read xyz
    readxyz(inputfilename, mysystem, cellsize, Lx, Ly, Lz);

    //Bond analysis

    //Output results
    ofp = fopen(outputfilename, "w");
}
