#include <vector>

/////////////////////////////////////////////
//Do RDF calculations, Coding here
/////////////////////////////////////////////



/*
class RDF {
public:
    std::vector<int> neighbours;
private:
    std::vector<int> getRDF() { return }

    float lx = mysystem[iframe]->L[0], ly = mysystem[iframe]->L[1], lz = mysystem[iframe]->L[2];
    float dr = 0.1, dx, dy, dz, dist, area, angle;
    int nbins = int(180 / dr) + 1, index;
    int bin[nbins];
    long anglenumber = 0;
    float g[nbins];
//          printf("Done 1\n");
//initialize the array!
    for(
    i = 0;
    i<nbins;
    i++) {
        bin[i] = 0;
        g[i] = 0.0;
    }

//for (iframe=0, iframe<MAXDUMP,iframe++) //loop of different atoms in different cell;
    int cellx, celly, cellz, nmember, imember; //cell index
    int ncellx, ncelly, ncellz; //neighbour cell index
    int ncell[3], iii, ii, c, inmember;
    int comparei, comparej;
    int itype, jtype;

    for(iii = 0;iii<3;iii++) { ncell[iii] = mysystem[iframe]->ncell[iii]; } //initialize array ncell[3]

//loop among different cell
    for (cellx = 0;cellx<ncell[0];cellx++){
        for (celly = 0; celly < ncell[1]; celly++) {
            for (cellz = 0; cellz < ncell[2]; cellz++) {
                nmember = mysystem[iframe]->mycell[cellx][celly][cellz]->nmember;

//loop among all the atoms contained in one cell

                for (imember = 0; imember < nmember; imember++) {
                    neighbours.clear();
//loop among all neighbour cells
                    for (ii = 0; ii < 27; ii++) {
                        c = mysystem[iframe]->mycell[cellx][celly][cellz]->nbr[ii];

                        ncellx = (int) (c / (ncell[1] * ncell[2]));
                        ncelly = (int) (fmod((int) (c / ncell[2]), ncell[1]));
                        ncellz = (int) (fmod(c, ncell[2]));
//printf("mysystem[%d]->mycell[%d][%d][%d]: ncellx:%d ncelly:%d ncellz:%d\n",iframe,cellx,celly,cellz,ncellx,ncelly,ncellz);
                        inmember = mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->nmember;
//loop among all the neighbour's atoms
                        for (j = 0; j < inmember; j++) {
                            comparei = mysystem[iframe]->mycell[cellx][celly][cellz]->member[imember]; // get atom index
                            comparej = mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->member[j]; // get atom index
                            itype = myatom[comparei]->type; // get type
                            jtype = myatom[comparej]->type; // get type
//printf("comparei:%d  comparej:%d\n",comparei,comparej);
//if (comparei==1066 && comparej==1067){printf("atom pair 1066-1067 has been included in the loops\n1066 in cell %d %d %d\n1067 in cell %d %d %d\n",cellx,celly,cellz,ncellx,ncelly,ncellz);}
//Avoid double counting of pair comparei-comparej
                            if (itype == speciesA && jtype == speciesB && comparei != comparej) {
//if (comparei<comparej && itype==speciesA && jtype==speciesB ) {
//printf("comparei:%d  comparej:%d\n",comparei,comparej);
                                float xi, yi, zi, xj, yj, zj;
                                xi = myatom[comparei]->x[0];
                                xj = myatom[comparej]->x[0];
                                yi = myatom[comparei]->x[1];
                                yj = myatom[comparej]->x[1];
                                zi = myatom[comparei]->x[2];
                                zj = myatom[comparej]->x[2];
                                dx = xi - xj; //printf("dx=%f\n",dx);
                                dy = yi - yj;
                                dz = zi - zj;
                                while (dx > lx / 2.0) { dx = dx - lx; }
                                while (dx < -lx / 2.0) { dx = dx + lx; }
                                while (dy > ly / 2.0) { dy = dy - ly; }
                                while (dy < -ly / 2.0) { dy = dy + ly; }
                                while (dz > lz / 2.0) { dz = dz - lz; }
                                while (dz < -lz / 2.0) { dz = dz + lz; }
                                dist = sqrt(dx * dx + dy * dy + dz *
                                                                dz);// printf("%d-%d:type%d-%d:dist %f ",myatom[comparei]->id,myatom[comparej]->id,itype,jtype,dist);
//if(myatom[comparei]->id==10 && myatom[comparej]->id==11 && myatom[comparej]->id==12) printf("10-11-12 dist=%f",dist);
//if(myatom[comparei]->id==10 && myatom[comparej]->id==12 && myatom[comparej]->id==11) printf("10-11-12 dist=%f",dist);
                                if (dist <= bondlength)
                                    neighbours.push_back(comparej); //printf("atomID %d\n",comparej);
                            }
                        }
                    }

// loop all neighbour bonds atoms and calculate bond angles
//if(neighbours.size()>0)printf("neighbours.size()=%d\n",neighbours.size());
                    for (il = 0; il < neighbours.size(); il++)
                        for (jl = il + 1; jl < neighbours.size(); jl++) {
                            angle = bondangle(myatom[neighbours[il]], myatom[comparei], myatom[neighbours[jl]], lx, ly,
                                              lz);
                            index = int(angle / dr); //printf("angle=%f,index=%d ",angle,index);
                            if (index <= nbins) {
                                anglenumber++;
                                bin[index]++;//printf("bin%d=%d\n",index,bin[index]);
                            }

                        }


                }
            }
        }
    }


    printf( "This is frame:(%d) \n", iframe);
    float ndensity=float(nmols_speciesB/(lx*ly*lz));
    float norm=(1.0/nmols_speciesA);
//printf("nmols:%d lx:%f ly:%f ly:%f ndensity:%f\n",nmols,lx,ly,lz,ndensity);
    printf ("anglenumber=%d\n",anglenumber);
    for (index=0;index<nbins;index++){
        angle=index*dr;
        if(index==0){g[index]=0;}
        else {
//area=4*PI*dist*dist*dr;g[index]=bin[index]/(area*ndensity)*norm;
            g[index]=float(bin[index])/float(anglenumber);
//printf("%d %d %f \n",anglenumber,bin[index], g[index]);
        }
    }

//The simulation box information is in mysystem[iframe]. For instance mysystem[iframe]->L[0] will be X-dimension of the system
//The atom information is in myatom. For instance myatom[i]->x[0] is the x-coordinate of myatom

};

*/

/*
fprintf(ofp, "\n#This is frame: (%d) \n#index  angle  angle-distribution\n", iframe-1);
for (index=0;index<nbins;index++)
{
dist=index*dr;
if(g[index]>0.0) fprintf (ofp, "%d %f %f \n",index,dist,g[index]);
}

//output
for(j=0;j<totalframes;j++) {
for(i=0;i<mysystem[j]->nmols;i++) {
if(mysystem[j]->myatom[i]->id==inputparameter)
fprintf(ofp, "%d %f %f %f\n", j,
mysystem[j]->myatom[i]->x[0],
mysystem[j]->myatom[i]->x[1],
mysystem[j]->myatom[i]->x[2]
);
}
}
*/
