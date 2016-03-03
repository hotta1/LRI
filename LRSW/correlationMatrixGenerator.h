//#include "ewald.h"

class correlationMatrixGenerator{
private:
  std::vector<std::vector<double> > corMatrix;
  simpleLattice* Lattice;
  int N;
//  double rho;

public:
  correlationMatrixGenerator(simpleLattice* Lattice_){
    Lattice = Lattice_;
    N = Lattice->num_sites();
    corMatrix.resize(N);
    for(int i=0 ; i<N ; ++i){ corMatrix[i].resize(N); }

  }

  std::vector<std::vector<double> > algebra(double rho){
    //ewald ewaldsum(Lattice->latticeLength, rho);
    ewald ewaldsum(Lattice->latticeLength, rho,4,4);
    std::vector<double> origin;
    origin.resize(Lattice->dim, 0.0);
    for(int i=0 ; i<Lattice->numVert ; i++){
      origin = Lattice->coordinate(i);
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int j=0 ; j<N ; j++){
        if(i==j){
          corMatrix[i][j] = 1.0;
          //std::cout <<corMatrix[i][j]<<"  ";
        }
        else{
          std::vector<double> ci = Lattice->coordinate(j);
          //corMatrix[i][j] = ewaldsum.ewaldsum(origin,ci);
          corMatrix[i][j] = ewaldsum.nearest(origin,ci)/3.0;
          //std::cout <<corMatrix[i][j]<<"  ";
        }
      }
      //std::cout << std::endl;
    }
    for(int i = Lattice->numVert ; i<N ; ++i){
      for(int j=0 ; j<N ; j++){
        int effectivei = i;
        int effectivej = j;
        Lattice->shiftRel(effectivei,effectivej);
        corMatrix[i][j] = corMatrix[effectivei][effectivej];
      }
    }
    return corMatrix;
  }
};

