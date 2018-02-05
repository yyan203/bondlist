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
#include <iostream>
#include <fstream>
#include <vector>
#include "bond.h"
//#include "simulationbox.h"
#include "readdata.h"



#define VERBOSE 1
#define MAXDUMP 1

////////////////////
// Main functions //
////////////////////

int main(int argc, char **argv)
{

    // This code is used to calculate bonds from xyz file,
    // The output is used as input for Moltemplate to further generate Angle, Dihedral, Improper angles.

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
    for(int i=0;i<MAXDUMP;i++) {mysystem[i] = new SimulationBox();}
    printf("Initialization of Simulation Boxes has been finished!\n");

    //Read xyz
       readxyz(inputfilename, mysystem, cellsize, Lx, Ly, Lz);

    //Bond analysis

    //Output results
    mysystem[0]->PrintInfo();
    mysystem[0]->myatom[0]->PrintInfo();
    mysystem[0]->myatom[271]->PrintInfo();
    Bonds allbonds(mysystem, 0);

    int bondnum = allbonds.getbondnum();
    printf("bondnum=%d\n",bondnum);
    std::ofstream ofp;
    ofp.open(outputfilename);
    std::vector<int> bond;
    for (auto i = 0;i < bondnum; i++){
        bond = allbonds.getbond(i);
        //ofp << "$bond:" << bond[2] << "  $atom:" << bond[0] + 1 << "  $atom:" << bond[1] + 1 << "\n"; // bond type + atomi + atomj
        ofp <<  bond[2] << "  " << bond[0] + 1 << "  " << bond[1] + 1 << "\n"; // bond type + atomi + atomj
    }
    ofp.close();

}
