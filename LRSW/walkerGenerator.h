#include <iostream>
#include <sstream>
#include <fstream>
//#include "ewald2d.h" 
#include "ewald.h"

class walkerGenerator{
private:
  simpleLattice* Lattice;
  std::vector<double> latticeConstant, latticeLength;
  std::vector<int> cellSize;
  int N, numVert, numCell, dim;
  std::string latticename, interaction;
  double sigma;
  std::string absolutePATH;

public:
  walkerGenerator(simpleLattice* Lattice_, std::string interaction_, double sigma_, std::string absolutePATH_){
    Lattice = Lattice_;
    latticeConstant = Lattice->latticeConstant;
    latticeLength = Lattice->latticeLength;
    cellSize = Lattice->cellSize;
    N = Lattice->num_sites();
    numVert = Lattice->numVert;
    numCell = Lattice->numCell;
    latticename = Lattice->latticename;
    dim = Lattice->dim;
    interaction = interaction_;
    sigma = sigma_;
    absolutePATH = absolutePATH_;
  }

  void binaryfile(){
    std::vector<double> bondprob;
    bondprob.resize(N*numVert,0.0);
    double Jtot=0;

    calcBondprob(bondprob,Jtot);

    looper::random_choice walkerChoice(bondprob);

    std::string walkerFN, JtotFN;
    walkerFN = walkerFilename();
    JtotFN = JtotFilename();
    std::ofstream Jtotfile;
    Jtotfile.open(JtotFN.c_str());
    Jtotfile.precision(16);
    Jtotfile << Jtot << std::endl;
    alps::OXDRFileDump odump(walkerFN);
    walkerChoice.save(odump);
    return;
  }

  looper::random_choice instance(double& Jtot){
    std::vector<double> bondprob;
    bondprob.resize(N*numVert,0.0);
    Jtot=0;

    calcBondprob(bondprob,Jtot);

    looper::random_choice walkerChoice(bondprob);
    return walkerChoice;
  }

  std::string walkerFilename(){
    std::ostringstream filename;
    filename << absolutePATH << "walkerTable/";
    filename << dim <<"D"<< latticename <<"_";
    for(int i=0 ; i<dim ; i++){
      filename << "L" << i+1 << "-" << cellSize[i]<<"_";
    }
    for(int i=0 ; i<dim ; i++){
      filename << "a" << i+1 << "-" << latticeConstant[i] <<"_";
    }
    if(interaction == "LRI"){
      filename << "sigma-" << sigma;
    }
    else{
      filename << interaction;
    }
    filename << ".xdr" ;
    std::string tmp;
    tmp = filename.str();
    return tmp;
  }

  std::string JtotFilename(){
    std::ostringstream filename;
    filename << absolutePATH << "Jtot/";
    filename << dim <<"D"<< latticename <<"_";
    for(int i=0 ; i<dim ; i++){
      filename << "L" << i+1 << "-" << cellSize[i]<<"_";
    }
    for(int i=0 ; i<dim ; i++){
      filename << "a" << i+1 << "-" << latticeConstant[i] <<"_";
    }
    if(interaction == "LRI"){
      filename << "sigma-" << sigma;
    }
    else{
      filename << interaction;
    }
    filename << ".txt" ;
    std::string tmp;
    tmp = filename.str();
    return tmp;
  }


private:
  void calcBondprob(std::vector<double>& bondprob, double& Jtot){
    ewald ewaldsum(latticeLength, sigma,4,4);//for the system with anisotropy, autohnumax must be implemented.
//    ewald ewaldsum(latticeLength, sigma);
    std::vector<double> origin;
    origin.resize(dim,0.0);
    for(int j=0 ; j<numVert ; j++){
      origin = Lattice->coordinate(j);
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int i=0 ; i<N ; i++){
        if(i!=j){ // to avoid self bonding
          std::vector<double> ci = Lattice->coordinate(i);
          double J;
          if(interaction == "LRI"){
            J = ewaldsum.ewaldsum(origin,ci);
          }
          else if(interaction == "MF"){
            J = 1.0/N;
          }
          else if(interaction == "SR"){
            if(Lattice->discriminateNearest(j,i)){ J = 1; }
            else{ J =0; }
          }
          else if(interaction == "free"){
            J = 0;
          }
          else{
            std::cout <<  "Interaction name \"" << interaction << "\" does not exist." << std::endl;
            J = 0;
          }
          Jtot += J;
          bondprob[i+N*j] = J;
          if(i==N-1){
            std::cout << i <<"/"<< N <<"  "<< origin[0] <<"  "<< origin[1] <<"  "<< ci[0] <<"  "<< ci[1] <<"  "<< J << std::endl;
          }
        }
      }
    }
    Jtot *= numCell;
    Jtot /= 2.0; //  /2.0 exists to counter double count.
  }
};



