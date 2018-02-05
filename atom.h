#include <string>
#include <cstring>
#include <iostream>

# define atom class
class Atom {  // this is about the information of every single atom
public:
  int id;
  std::string type_; // e.g. Zn, C, N, H
  int index; // e.g. the 1 in Zn1,  the 10 in N10
  float x[3];
  float v[3];
  float f[3];
  void PrintInfo() {
    std::cout << "Atom: id["<< id <<"] type[" << type_ << "] x["<< x[0] << ","<< x[1]<< ","<<x[2]<< "]\n";
    // printf("Atom: id[%d] type[%s] x[%f,%f,%f]\n", id, type_.c_str(),x[0],x[1],x[2]);
  }
};

