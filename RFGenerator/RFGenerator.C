#include <alps/parapack/parapack.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/random.hpp>
#include <unistd.h>
#include <boost/algorithm/string.hpp>
#include<boost/timer/timer.hpp>
#include "../LRSW/ewald.h"
#include "../LRSW/simpleLattice.h"
#include "../LRSW/correlatedDistributionGenerator.h"
#include "../LRSW/RFGenerator.h"

std::string replaceString(const std::string target, const std::string from, const std::string to);

int main(){
  std::string absolutepath;
  int num_clones;
  std::string latticename;
  std::string binarymatrixfile="on";
  double fieldintensity;
  double rfstddev;
  double rho;
  int dim;
  double transmatnormalization = -1.0;
  std::vector<int> L;
  std::vector<std::vector<double> > T;

  std::ifstream ipfile("LRSW_params");
  std::string line;
  while(ipfile){
    getline(ipfile,line);
    std::vector<std::string> strvec;
    boost::algorithm::split( strvec, line, boost::algorithm::is_space());
    if(strvec.size()>2){
      if(strvec[0]=="BINARYRFFILE"&&strvec[2]!="\"on\""){ std::cout <<"BINARYRFFILE==off. Run without binaryrffile."<<std::endl; return 0; }
      if(strvec[0]=="RANDOMFIELD"&&strvec[2]!="\"on\""){ std::cout <<"RANDOMFIELD==off. Run without binaryrffile."<<std::endl; return 0; }
      if(strvec[0]=="RFCORRELATION"&&strvec[2]!="\"on\""){ std::cout <<"RFCORRELATION==off. Run without binaryrffile."<<std::endl; return 0; }
      if(strvec[0]=="TRANSMATNORMALIZATION"){ transmatnormalization = std::atof(strvec[2].c_str()); }
      if(strvec[0]=="BINARYMATRIXFILE"){ binarymatrixfile = strvec[2]; }
      if(strvec[0]=="ABSOLUTEPATH"){ absolutepath = strvec[2]; }
      if(strvec[0]=="NUM_CLONES"){ num_clones = std::atoi(strvec[2].c_str()); }
      if(strvec[0]=="LATTICE"){ latticename = strvec[2]; }
      if(strvec[0]=="FIELDINTENSITY"){ fieldintensity = std::atof(strvec[2].c_str()); }
      if(strvec[0]=="RFSTDDEV"){ rfstddev = std::atof(strvec[2].c_str()); }
      if(strvec[0]=="RHO"){ rho = std::atof(strvec[2].c_str()); }
      if(strvec[0]=="DIMENSION"){ dim = std::atoi(strvec[2].c_str()); }
      if(strvec[0]=="L"){ 
        L.push_back(std::atoi(strvec[2].c_str()));
        std::vector<double> tmp;
        T.push_back(tmp);
      }
      if(strvec[0]=="{" && strvec[1]=="T"){ 
        T[L.size()-1].push_back(std::atof(strvec[3].c_str()));
      }
    }
  }

  absolutepath = replaceString(absolutepath, "\"","");
  latticename = replaceString(latticename, "\"","");
  binarymatrixfile = replaceString(binarymatrixfile, "\"","");
  rho = static_cast<double>(dim) - rho;
  if(transmatnormalization < 0.0){ transmatnormalization = static_cast<double>(dim) + 2.0; }

  boost::mt19937 mt(std::time(0));

  for(int i=0 ; i<T.size() ; ++i){
    for(int j=0 ; j<T[i].size() ; ++j){
      boost::timer::cpu_timer timer;
      simpleLattice *Lattice;
      if(latticename=="square lattice" || latticename=="square_lattice"){ Lattice = new squareLattice(L[i], dim); }
      else if(latticename == "triangular lattice" || latticename=="triangular_lattice"){ Lattice = new triangularLattice(L[i]); }
      else{ std::cout << "Aveilable latticename is square_lattice or triangular_lattice" << std::endl; std::exit(1); }
    
      RFGenerator RFGen;
      RFGen.generate(Lattice, fieldintensity, rfstddev, rho, T[i][j], num_clones, absolutepath, binarymatrixfile, transmatnormalization, mt);
      std::cout << "Finished creating random field. It takes " << timer.format() << std::endl;
      delete Lattice;
    }
  }

  return 0;
}

std::string replaceString(const std::string target, const std::string from, const std::string to) {
  std::string result = target;
  std::string::size_type pos = 0;
  while(pos = result.find(from, pos), pos != std::string::npos) {
    result.replace(pos, from.length(), to);
    pos += to.length();
  }
  return result;
}
