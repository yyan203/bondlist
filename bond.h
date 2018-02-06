#include <map>
#include <unordered_map>
#include "math.h"
#include "simulationbox.h"
#include <string>
#include <vector>

#define PI 3.1415926

// Fun1
   // define bond length within which two atoms form a bond
std::unordered_map<std::string, float > bondlen{
  {"Zn-N" , 1.970000 },
  {"N-Zn" , 1.970000 },
  {"H-C"  , 1.020000 },
  {"C-H"  , 1.020000 },
  {"C-N"  , 1.340000 },
  {"N-C"  , 1.340000 },
  {"C-C"  , 1.390000 },
};
// new
std::unordered_map<std::string, int > bondtype{
        {"Zn-N" , 1},
        {"N-Zn" , 1},
        {"H-C"  , 2},
        {"C-H"  , 2},
        {"C-N"  , 3},
        {"N-C"  , 3},
        {"C-C"  , 4}
};

class Bonds {
private:
    int bondnum_;
    std::vector<std::vector<int>> bonds_;

public:
    Bonds (SimulationBox **mysystem, int frame = 0 ){
        bondnum_ = 0;
        findbonds(mysystem, frame);
    };

    float distance(const Atom * A, const Atom * B, const SimulationBox * system) const {
        float xi, yi, zi, xj, yj, zj;
        float dx, dy, dz, lx, ly, lz;
        lx = system->L[0];
        ly = system->L[1];
        lz = system->L[2];
        xi = A->x[0];
        xj = B->x[0];
        yi = A->x[1];
        yj = B->x[1];
        zi = A->x[2];
        zj = B->x[2];
        dx = xi - xj; dy = yi - yj; dz = zi - zj;
        while (dx > lx / 2.0) { dx = dx - lx; }
        while (dx < -lx / 2.0) { dx = dx + lx; }
        while (dy > ly / 2.0) { dy = dy - ly; }
        while (dy < -ly / 2.0) { dy = dy + ly; }
        while (dz > lz / 2.0) { dz = dz - lz; }
        while (dz < -lz / 2.0) { dz = dz + lz; }
        return sqrt(dx * dx + dy * dy + dz * dz);
    }
    // check if two atoms form a bond, if no return 0, else return bond type 1, 2, 3, ...
    int ifbonds (const Atom * A, const Atom * B, const SimulationBox * system) const {
        std::string bondname;
        bondname += A->type_;
        bondname += "-";
        bondname += B->type_;
        auto got = bondtype.find(bondname);
        if (got == bondtype.end()){return 0;}
        if (this->distance(A, B, system) < 1.20 * bondlen[bondname]) {  // 1.15 is good enough for zif62
            // if (this->distance(A, B, system) < 1.18 * bondlen[bondname]) {
            return bondtype[bondname];
        }
        return 0;
    }

    void findbonds (SimulationBox **mysystem, int frame){
        Atom **myatom = mysystem[frame]->myatom;
        int nmols = mysystem[frame]->nmols;
        float x1=0,x2=0,x3=0;
        float originx,originy,originz;
        int cellx, celly, cellz, nmember, imember; //cell index
        float cellsizex,cellsizey,cellsizez;
        int ncellx, ncelly, ncellz; //neighbour cell index
        int ncell[3], iii, ii, c, inmember;
        int comparei, comparej;
        Atom *iatom, *jatom;
        int p1,p2,p3,memberid;

        // dimension of mycell
        for(int d=0;d<3;d++) {
            ncell[d] = mysystem[frame]->ncell[d];
        }
        originx= mysystem[frame]->origin[0];
        originy= mysystem[frame]->origin[1];
        originz= mysystem[frame]->origin[2];
        cellsizex=(float)(mysystem[frame]->L[0]/mysystem[frame]->ncell[0]);
        cellsizey=(float)(mysystem[frame]->L[1]/mysystem[frame]->ncell[1]);
        cellsizez=(float)(mysystem[frame]->L[2]/mysystem[frame]->ncell[2]);

        // loop all atoms
        for (int i = 0; i < nmols; i++){
            x1=myatom[i]->x[0]; // if (i==400){printf("%d\n",i);}
            x2=myatom[i]->x[1];
            x3=myatom[i]->x[2];
            //according to the atom's position, allocate them to the right cell.
            //Since cellsizex is float number, to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell.
            p1=(int)((x1-originx)/cellsizex);
            if (p1==mysystem[frame]->ncell[0]){p1=p1-1;}
            p2=(int)((x2-originy)/cellsizey);
            if (p2==mysystem[frame]->ncell[1]){p2=p2-1;}
            p3=(int)((x3-originz)/cellsizez);
            if (p3==mysystem[frame]->ncell[2]){p3=p3-1;}
            memberid = mysystem[frame]->mycell[p1][p2][p3]->nmember;

            // loop all neighbor cells
            for (ii = 0; ii < 27; ii++) {
                c = mysystem[frame]->mycell[p1][p2][p3]->nbr[ii];
                ncellx = (int) (c / (ncell[1] * ncell[2]));
                ncelly = (int) (fmod((int) (c / ncell[2]), ncell[1]));
                ncellz = (int) (fmod(c, ncell[2]));
                inmember = mysystem[frame]->mycell[ncellx][ncelly][ncellz]->nmember; // get num of atom member
            //loop among all the neighbour's atoms
                for (int j = 0; j < inmember; j++) {
                    //comparei = mysystem[frame]->mycell[cellx][celly][cellz]->member[imember]; // get atom index
                    comparej = mysystem[frame]->mycell[ncellx][ncelly][ncellz]->member[j]; // get atom index
                    iatom = myatom[i];        // get atom object
                    jatom = myatom[comparej]; // get atom object
            //Avoid double counting of pair comparei-comparej
                    if (i < comparej) {
                       int type = this->ifbonds(iatom, jatom, mysystem[frame]);
                       if ( type > 0) {
                            std::vector<int> bond;
                            bond.push_back(i);
                            bond.push_back(comparej);
                            bond.push_back(type);
                            this->bonds_.push_back(bond);
                            this->bondnum_ = this->bondnum_ + 1;
                        }
                    }
                }
            }

        }  // loop all atoms

    }  // findbonds

    int get_bondtype(const std::string& s)
    {
      return  bondtype[s];
    }

    std::vector<int> getbond(int index) const
    {
        if (index >= bondnum_) { printf("Not right\n");exit(0);}
        return  bonds_[index];
    }
    int getbondnum() const
    {
        return bondnum_;
    }

    float get_bondlen(const std::string& s)
    {
        return  bondlen[s];
    }


    // Fun2
       // calculate angle between A-B-C
       //float bondangle( Atom& A,  Atom& B,  Atom& C)

    float bondangle(const Atom * A, const Atom * B, const Atom * C, const float &lx, const float &ly, const float &lz)
    {
      float angle=0;
      float v1x,v1y,v1z;
      float v2x,v2y,v2z;
      v1x=A->x[0]-B->x[0];
      v1y=A->x[1]-B->x[1];
      v1z=A->x[2]-B->x[2];
      v2x=C->x[0]-B->x[0];
      v2y=C->x[1]-B->x[1];
      v2z=C->x[2]-B->x[2];
      while (v1x>lx/2.0) {v1x=v1x-lx;}
      while (v1x<-lx/2.0){v1x=v1x+lx;}
      while (v1y>ly/2.0) {v1y=v1y-ly;}
      while (v1y<-ly/2.0){v1y=v1y+ly;}
      while (v1z>lz/2.0) {v1z=v1z-lz;}
      while (v1z<-lz/2.0){v1z=v1z+lz;}
      while (v2x>lx/2.0) {v2x=v2x-lx;}
      while (v2x<-lx/2.0){v2x=v2x+lx;}
      while (v2y>ly/2.0) {v2y=v2y-ly;}
      while (v2y<-ly/2.0){v2y=v2y+ly;}
      while (v2z>lz/2.0) {v2z=v2z-lz;}
      while (v2z<-lz/2.0){v2z=v2z+lz;}
      angle=acos((v1x*v2x+v1y*v2y+v1z*v2z)/(sqrt(v1x*v1x+v1y*v1y+v1z*v1z)*sqrt(v2x*v2x+v2y*v2y+v2z*v2z)))*180.0/PI;
  if(angle<0) angle+=180;
  return angle;
  }

}; //class Bonds
