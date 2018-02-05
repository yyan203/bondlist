#include "atom.h"
#define MaxCellAtoms 500
#define MaxCellNbr 100

// define cell and SimulationBox
class Cell {   // this is about the cell formed by dividing the box 
public:
  int member[MaxCellAtoms];
  int nmember;
  int nbr[MaxCellNbr];
  int nnbr;
  float origin[3];
  float size[3];
  Cell() {nmember=0;nnbr=0;};
};


class SimulationBox { 
public:
  Atom **myatom;
  Cell ****mycell;
  int ncell[3]; //dimension of the cells
  float L[3];
  float origin[3]; //origin of the simulation box, default will be (0,0,0)
  int nmols;
  void SetupCell(float cellsize) {    
#if VERBOSE ==1 
    printf("Initiating setupcell\n");
#endif
    int i,j,k,p,q,r,pi,qi,ri,l;
    if(nmols==0) {
      printf("Simulation Box is empty!\n");
      exit(0);
    }
    for(int d=0;d<3;d++) {
      ncell[d] = (int)(L[d]/cellsize)-1; //make it a little larger than necessary
      //if(ncell[d]<4) {ncell[d]=1;} //if the cells are too few, there are no advantages of doing link-cell.
    }
#if VERBOSE ==1 
    printf("cell[%d,%d,%d]\n", ncell[0],ncell[1],ncell[2]);
#endif
    //if(ncell[0]==1||ncell[1]==1||ncell[2]==1){printf ("Warning!!! Cells are too few to use linkcell!!! Please decrese cellsize.\n"); exit(0);}
    mycell = new Cell ***[ncell[0]];
    for(i=0;i<ncell[0];i++) {
      mycell[i] = new Cell **[ncell[1]];
      for(j=0;j<ncell[1];j++) {
	mycell[i][j] = new Cell *[ncell[2]];
	for(k=0;k<ncell[2];k++) {
	  mycell[i][j][k] = new Cell();
	}
      }
    }

    //Setup NBR cells, Coding here
///////Yongjian code following/////////////
   //loop of cells
  for(i=0;i<ncell[0];i++) {
     for(j=0;j<ncell[1];j++) {
       for(k=0;k<ncell[2];k++) {
//printf("k has a value: %d\n",k);
          //loop of neighbour cells
          mycell[i][j][k]->nnbr=27;l=0; 
        for (pi=i-1;pi<i+2;pi++){ 
           for (qi=j-1;qi<j+2;qi++){
              for (ri=k-1;ri<k+2;ri++){
//printf("ri has a value: %d\n",ri);
         //convert neighbour cells outside the range (periodic boundary condition) 
         p=pi;q=qi;r=ri;l=(pi-i+1)*9+(qi-j+1)*3+(ri-k+1);     
      	   if (l<0 || l>26) {printf( "l:%d is wrong",l); exit(0);}
           if (p==-1){p=ncell[0]-1;}
           if (p==ncell[0]){p=0;}  //if (p==-1||p==ncell[0]){p=abs(abs(p)-ncell[0]);}
           if (q==-1){q=ncell[1]-1;}
           if (q==ncell[1]){q=0;}  //if (q==-1||q==ncell[1]){q=abs(abs(q)-ncell[1]);} 
           if (r==-1){r=ncell[2]-1;} 
           if (r==ncell[2]){r=0;}  //if (r==-1||r==ncell[2]){r=abs(abs(r)-ncell[2]);}
          //convert every neighbour's vector index into a scalor index 
            mycell[i][j][k]->nbr[l]=p*ncell[1]*ncell[2]+q*ncell[2]+r;
//printf("l has a value: %d\nnbr[l] has a value:%d\ni=%d j=%d k=%d\n",l,mycell[i][j][k]->nbr[l],i,j,k);

} } }
//printf("mycell[i][j][k] has %d neighbour\ni j k are:%d %d %d\n", mycell[i][j][k]->nnbr,i,j,k);     
   }
  }
 }
printf("i j k are  %d %d %d\n", i,j,k);
};     
////////Yongjian code ends////////////////

  SimulationBox() {ncell[0]=ncell[1]=ncell[2]=0; nmols=0;mycell=NULL;};
  
  ~SimulationBox() {
    int i,j,k;
    if(mycell!=NULL) {
#if VERBOSE == 1
      printf("Now destroy mycell.\n");
#endif
      for(i=0;i<ncell[0];i++) {
	for(j=0;j<ncell[1];j++) {
	  for(k=0;k<ncell[2];k++) {
	    delete mycell[i][j][k];
	  }
	  delete [] mycell[i][j];
	}
	delete [] mycell[i];
      }      
      delete [] mycell;
    }
  };
  void PrintInfo() {
    printf("SimulationBox: ncell(%d,%d,%d), L(%f,%f,%f), origin(%f,%f,%f),nmols(%d)\n",
	   ncell[0],ncell[1],ncell[2],
	   L[0],L[1],L[2],
	   origin[0],origin[1],origin[2],
	   nmols);
  };
};
