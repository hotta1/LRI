#include <alps/parapack/parapack.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <boost/assign/list_of.hpp>
#include"../../LRSW/RFGenerator.h"
#include <boost/random.hpp>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[]){
  std::string absolutepath = argv[1];
  int num_clones = atoi(argv[2]);
  std::string latticename = argv[3];
  double fieldintensity = atof(argv[4]);
  double rfstddev = atof(argv[5]);
  double rho = atof(argv[6]);
  int dim = atoi(argv[7]);
  int L = atoi(argv[8]);
  double T = atof(argv[9]);

  simpleLattice *Lattice;
  if(latticename=="square lattice" || latticename=="square_lattice"){ Lattice = new squareLattice(L, dim); }
  else if(latticename == "triangular lattice" || latticename=="triangular_lattice"){ Lattice = new triangularLattice(L); }
  else{ std::cout << "Aveilable latticename is square_lattice or triangular_lattice" << std::endl; std::exit(1); }

  RFGenerator RFGen;
  RFGen.generate(Lattice, fieldintensity, rfstddev, rho, T, absolutepath, num_clones);

  for(int i = 0 ; i < num_clones ; ++i){
    std::vector<double> RF = RFGen.load(Lattice, fieldintensity, rfstddev, rho, T, i);
    for(int j=0 ; j < 4 ; ++j){
      std::cout << RF[j];
    }
    std::cout << std::endl;
  }
  return 0;
}

