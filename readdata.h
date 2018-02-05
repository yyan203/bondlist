// yongjian Feb 3 2018
#include <stdio.h>
#include <stdlib.h>
#include <string>
//#include "simulationbox.h"
#define ONELINEMAX 2000
#define MAXDUMP 5

// read xyz data
   void readxyz( char inputfilename[], SimulationBox ** mysystem, const float cellsize,
                            const float Lx, const float Ly, const float Lz){

    /************** Import the data file ********************/
    FILE *ifp;
    ifp = std::fopen(inputfilename, "r");
    if(ifp == NULL) {
        printf("File %s does not exist!\n", inputfilename);
        exit(0);
    }
    Atom **myatom;
    int totalframes;
    int i=0,j=0,k=0;
    int il=0,jl=0,kl=0;
    int iframe = 0; //input number of frames
    float xlo, xhi, ylo, yhi, zlo, zhi;
    int id;
    char type[2];
    int current_timestep, nmols;
    char templine[ONELINEMAX], str1[100], str2[100];
    float ix, iy, iz, ivx, ivy, ivz, ifx, ify, ifz;

    while(!feof(ifp)) {
            if(fgets(templine, ONELINEMAX, ifp) == NULL) {
                break;
            } // ITEM: NUMBER OF ATOMS  && also check if end of line is reached

            sscanf(templine, "%ld", &nmols);  //number of atoms

            fgets(templine, ONELINEMAX, ifp); // Comment line

            mysystem[iframe]->origin[0] = 0.0;
            mysystem[iframe]->origin[1] = 0.0;
            mysystem[iframe]->origin[2] = 0.0;

            mysystem[iframe]->L[0] = Lx;
            mysystem[iframe]->L[1] = Ly;
            mysystem[iframe]->L[2] = Lz;

            //initialize myatom
            if(nmols>0) {
                mysystem[iframe]->nmols = nmols;
                mysystem[iframe]->myatom = new Atom *[nmols];
                myatom = mysystem[iframe]->myatom;
                for(auto i=0; i < nmols; i++) {
                    myatom[i] = new Atom();
                }
                printf("Initialization of myatom has been successfully finished!\n");
            }
            else {
                printf("Error in nmols (=%d) \n", nmols);
                exit(0);
            }

            //Be aware that id may start from 1
            //id is the TRUE identification of atoms
            //your i index is NOT the indentification of atoms

            printf("Importing %d atoms...\n", nmols);
            for(auto i=0; i < nmols; i++) {
                fgets(templine, ONELINEMAX, ifp);
                sscanf(templine, "%s %f %f %f", &type, &ix, &iy, &iz);
                myatom[i]->id = i + 1; // if (id==4000){printf("%d\n",id);}
                myatom[i]->type_ += type[0];
                if (type[1] == 'n') myatom[i]->type_ += type[1];
                myatom[i]->x[0] = ix;
                myatom[i]->x[1] = iy;
                myatom[i]->x[2] = iz;
                if (ix < mysystem[iframe]->origin[0]) mysystem[iframe]->origin[0] = ix;
                if (iy < mysystem[iframe]->origin[1]) mysystem[iframe]->origin[1] = iy;
                if (iz < mysystem[iframe]->origin[2]) mysystem[iframe]->origin[2] = iz;
            }
            //Setup the LinkCell
            mysystem[iframe]->SetupCell(cellsize);
            printf("%d\n",iframe);

            //Throw in atoms, Coding here
            float x1=0,x2=0,x3=0;
            float originx,originy,originz;
            float cellsizex,cellsizey,cellsizez;
            int p1,p2,p3,memberid;
            originx= mysystem[iframe]->origin[0]; printf("originx=%f\n",originx);
            originy= mysystem[iframe]->origin[1];
            originz= mysystem[iframe]->origin[2];
            cellsizex=(float)(mysystem[iframe]->L[0]/mysystem[iframe]->ncell[0]); printf("cellsizex=%f\n",cellsizex);
            cellsizey=(float)(mysystem[iframe]->L[1]/mysystem[iframe]->ncell[1]);
            cellsizez=(float)(mysystem[iframe]->L[2]/mysystem[iframe]->ncell[2]);
            //loop of every atoms of one frame in dump file
            for (i=0;i<nmols;i++) {
                x1=myatom[i]->x[0]; // if (i==400){printf("%d\n",i);}
                x2=myatom[i]->x[1];
                x3=myatom[i]->x[2];
                //according to the atom's position, allocate them to the right cell. Since cellsizex is float number, to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell.
                p1=(int)((x1-originx)/cellsizex);
                if (p1==mysystem[iframe]->ncell[0]){p1=p1-1;}
                p2=(int)((x2-originy)/cellsizey);
                if (p2==mysystem[iframe]->ncell[1]){p2=p2-1;}
                p3=(int)((x3-originz)/cellsizez);
                if (p3==mysystem[iframe]->ncell[2]){p3=p3-1;}
                memberid=mysystem[iframe]->mycell[p1][p2][p3]->nmember;
                //allocate the atom's index to the right cell
                mysystem[iframe]->mycell[p1][p2][p3]->member[memberid] = i;
                //increase the nmember of cell when a new atom is added into it.
                mysystem[iframe]->mycell[p1][p2][p3]->nmember = mysystem[iframe]->mycell[p1][p2][p3]->nmember + 1;
            }
            //Increment
            iframe++;
            if(iframe>=MAXDUMP) {
                printf("Too many dumps (%d), increase MAXDUMP!\n", iframe);
                exit(0);
            }
    }
    totalframes = iframe;
    printf("Total number of frames is %d\n", totalframes);
    fclose(ifp);
}

/////////////////
/************** Import the data file ********************/
//Read lata4olivia data
/*
void readlata4olivia(char inputfilename[], SimulationBox ** mysystem, const float cellsize){

    FILE *ifp;
    ifp = std::fopen(inputfilename, "r");

    if(ifp==NULL) {
        printf("File %s does not exist!\n", inputfilename);
        exit(0);
    }
    //initialize Simulation Boxes

    Atom **myatom;
    int totalframes;
    int iframe=0; //input number of frames
    float xlo, xhi, ylo, yhi, zlo, zhi;
    int id;
    char type[2];
    long current_timestep, nmols;
    char templine[ONELINEMAX], str1[100], str2[100];
    float ix, iy, iz, ivx, ivy, ivz, ifx, ify, ifz;

    while(!feof(ifp)) {
        fgets(templine, ONELINEMAX, ifp); // TIMESTEP line
        sscanf(templine, "%s %s", str1, str2);

        //    printf("For the first line: (%s) (%s) (%s)\n", templine, str1, str2);
        if(strcmp(str1, "ITEM:")==0) {
            fgets(templine, ONELINEMAX, ifp); //actual timesteps
            sscanf(templine, "%ld", &current_timestep);
            fgets(templine, ONELINEMAX, ifp); // ITEM: NUMBER OF ATOMS
            fgets(templine, ONELINEMAX, ifp); // nmols
            sscanf(templine, "%ld", &nmols);  //number of atoms

            fgets(templine, ONELINEMAX, ifp); // ITEM: BOX BOUND
            fgets(templine, ONELINEMAX, ifp); // xlo xhi
            sscanf(templine, "%f %f", &xlo, &xhi);
            fgets(templine, ONELINEMAX, ifp); // ylo yhi
            sscanf(templine, "%f %f", &ylo, &yhi);
            fgets(templine, ONELINEMAX, ifp); // zlo zhi
            sscanf(templine, "%f %f", &zlo, &zhi);
            fgets(templine, ONELINEMAX, ifp); //empty line

            mysystem[iframe]->origin[0] = xlo;
            mysystem[iframe]->origin[1] = ylo;
            mysystem[iframe]->origin[2] = zlo;

            mysystem[iframe]->L[0] = xhi - xlo;
            mysystem[iframe]->L[1] = yhi - ylo;
            mysystem[iframe]->L[2] = zhi - zlo;

            //initialize myatom
            if(nmols>0) {
                mysystem[iframe]->nmols = nmols;
                mysystem[iframe]->myatom = new Atom *[nmols];
                myatom = mysystem[iframe]->myatom;
                for(int i=0;i<nmols;i++) {
                    myatom[i] = new Atom();
                }
                printf("Initialization of myatom has been successfully finished! \n");
            }
            else {
                printf("Error in nmols (=%d) \n", nmols);
                exit(0);
            }

            //Be aware that id may start from 1
            //id is the TRUE identification of atoms
            //your i index is NOT the indentification of atoms

            printf("Importing %d atoms...\n", nmols);
            for(int i=0;i<nmols;i++) {
                fgets(templine, ONELINEMAX, ifp);
                sscanf(templine, "%d %s %f %f %f %f %f %f %f %f %f",
                       &id, &type, &ix, &iy, &iz, &ivx, &ivy, &ivz, &ifx, &ify, &ifz);
                myatom[i]->id = id; // if (id==4000){printf("%d\n",id);}
                myatom[i]->type_ = type[0] + type[1];
                myatom[i]->x[0] = ix;
                myatom[i]->x[1] = iy;
                myatom[i]->x[2] = iz;
            }
            //Setup the LinkCell
            mysystem[iframe]->SetupCell(cellsize);
            printf("%d\n",iframe);

//////////////////////////////////////////////////////
//////Yongjian Coded the following part///////////////
//////////////////////////////////////////////////////

            //Throw in atoms, Coding here
            float x1=0,x2=0,x3=0;
            float originx,originy,originz;
            float cellsizex,cellsizey,cellsizez;
            int p1,p2,p3,memberid;
            originx= mysystem[iframe]->origin[0]; printf("originx=%f\n",originx);
            originy= mysystem[iframe]->origin[1];
            originz= mysystem[iframe]->origin[2];
            cellsizex=(float)(mysystem[iframe]->L[0]/mysystem[iframe]->ncell[0]); printf("cellsizex=%f\n",cellsizex);
            cellsizey=(float)(mysystem[iframe]->L[1]/mysystem[iframe]->ncell[1]);
            cellsizez=(float)(mysystem[iframe]->L[2]/mysystem[iframe]->ncell[2]);
            //loop of every atoms of one frame in dump file
            for (int i=0;i<nmols;i++) {
                x1=myatom[i]->x[0]; // if (i==400){printf("%d\n",i);}
                x2=myatom[i]->x[1];
                x3=myatom[i]->x[2];
                //according to the atom's position, allocate them to the right cell. Since cellsizex is float number, to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell.
                p1=(int)((x1-originx)/cellsizex);
                if (p1==mysystem[iframe]->ncell[0]){p1=p1-1;}
                p2=(int)((x2-originy)/cellsizey);
                if (p2==mysystem[iframe]->ncell[1]){p2=p2-1;}
                p3=(int)((x3-originz)/cellsizez);
                if (p3==mysystem[iframe]->ncell[2]){p3=p3-1;}
                memberid=mysystem[iframe]->mycell[p1][p2][p3]->nmember;
//printf ("p1=%d p2=%d p3=%d\n",p1,p2,p3);
                //allocate the atom's index to the right cell
                mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=i;
                //mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=myatom[i]->id;
//printf("mycell%d%d%d->member[%d]:%d\n",p1,p2,p3,memberid,mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]);
//if (p1==0 && p2==0 && p3==0){printf("mycell-000 has atom(id): %d\n",mysystem[iframe]->mycell[p1][p2][p3]->member[mysystem[iframe]->mycell[p1][p2][p3]->nmember]);}
                //increase the nmember of cell when a new atom is added into it.
                mysystem[iframe]->mycell[p1][p2][p3]->nmember=mysystem[iframe]->mycell[p1][p2][p3]->nmember+1;
//if(i==400){printf("nmember of mycell[%d][%d][%d] is %d\n",p1,p2,p3,mysystem[iframe]->mycell[p1][p2][p3]->nmember);}
            }

//for (i=0;i<7;i++) {printf("mycell[0][%d][0] has nmember:%d\n",i,mysystem[iframe]->mycell[0][i][0]->nmember);}


            //Increment
            iframe++;
            if(iframe>=MAXDUMP) {
                printf("Too many dumps (%d), increase MAXDUMP!\n", iframe);
                exit(0);
            }
        }
        else {
            if(!feof(ifp))
                printf("Not the right syntax, end of file.\n");
            else
                printf("Finish reading in input file\n");
        }
    }
    totalframes = iframe;
    printf("Total number of frames is %d\n", totalframes);
    //  printf("Atom %f\n", mysystem[0]->myatom[1599]->x[0]);
    //  printf("atom: %f %f %f\n", mysystem[0]->myatom[1599]->x[0], mysystem[0]->myatom[1599]->x[1], mysystem[0]->myatom[1599]->x[2]);

    //mysystem[0]->PrintInfo();

    fclose(ifp);

}
*/
//deallocate memory
/*
for(j=0;j<totalframes;j++) {
    myatom = mysystem[j]->myatom;
    for(i=0;i<mysystem[j]->nmols;i++) {
        delete myatom[i];
    }
    delete [] myatom;
}
for(j=0;j<MAXDUMP;j++) {
    delete mysystem[j];
}
delete [] mysystem;
 */
