#include <alps/parapack/parapack.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <boost/assign/list_of.hpp>
#include <boost/random.hpp>
#include <unistd.h>

class RFGenerator{
public:
  RFGenerator(){}

  void generate(simpleLattice *Lattice, double fieldintensity, double rfstddev, double rho, double T, int num_clones, std::string absolutepath, std::string binarymatrixfile, double cormatnormalization, boost::mt19937 &mt){
    boost::normal_distribution<double> dist(fieldintensity, rfstddev);
//    correlatedDistributionGenerator<boost::normal_distribution<double> > corDistGen(Lattice, dist, cormatnormalization);
//    corDistGen.eigen(absolutepath, rho, binarymatrixfile);
    correlatedDistributionGenerator<boost::normal_distribution<double> > corDistGen(Lattice, dist);
    corDistGen.eigen(rho, cormatnormalization, absolutepath, binarymatrixfile);
    #ifdef _OPENMP
    #pragma omp parallel
    if(omp_get_thread_num()==0){ std::cout << "This program is done with OPENMP. omp_num_threads="<< omp_get_num_threads() <<", omp_max_threads="<< omp_get_max_threads() <<std::endl; }
    #pragma omp for
    #endif
    for(int i=1 ; i < num_clones+1 ; ++i){
      std::vector<double> extField;
      extField = corDistGen.generate(mt);
      std::string extFieldFN = extFieldFilename(Lattice, fieldintensity, rfstddev, cormatnormalization, rho, T, i);
      alps::OXDRFileDump odump(extFieldFN);
      odump << extField;
    }
    return;
  }

  std::vector<double> load(simpleLattice *Lattice, double fieldintensity, double rfstddev, double cormatnormalization, double rho, double T, int clone_id){
    std::vector<double> extField;
    std::ifstream extFieldFile;
    std::string extFieldFN = extFieldFilename(Lattice, fieldintensity, rfstddev, cormatnormalization, rho, T, clone_id);
    extFieldFile.open(extFieldFN.c_str());
    if(!extFieldFile){ std::cout << "error: " << extFieldFN << "does not exist." << std::endl; std::exit(1); }
    alps::IXDRFileDump idump(extFieldFN);
    idump >> extField;
    return extField;
  }

private:
  std::string extFieldFilename(simpleLattice *Lattice, double fieldintensity, double rfstddev, double cormatnormalization, double rho, double T, int clone_id){
    std::ostringstream filename;
    filename << get_current_directory();
    filename <<"/"<< Lattice->dim <<"D"<< Lattice->latticename <<"_";
    for(int i=0 ; i<Lattice->dim ; i++){ filename << "L" << i+1 << "-" << Lattice->cellSize[i]<<"_"; }
    //for(int i=0 ; i<Lattice->dim ; i++){ filename << "a" << i+1 << "-" << Lattice->latticeConstant[i] <<"_"; }
    filename << "rfave-" << fieldintensity;
    filename << "_rfstddev-" << rfstddev;
    filename << "_normalization-" << cormatnormalization;
    filename << "_rho-" << rho;
    filename << "_T-" << T;
    filename << "_ID-" << clone_id;
    filename << ".xdr" ;
    std::string tmp;
    tmp = filename.str();
    return tmp;
  }

  inline std::string get_current_directory(){
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    return cwd;
  }
};

