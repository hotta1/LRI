#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <boost/math/special_functions/gamma.hpp>

class ewald{
private:
  double kappa;
  int numax, hmax;
  std::vector<double> L;
  double sigma;
  int dim;

public:
  ewald(std::vector<double> L_, double sigma_){
    L = L_;
    dim = L.size();
    sigma = sigma_;
    double Lsum=0;
    for(int i=0 ; i<dim ; i++){Lsum += L[i];}
    kappa = 2.0*dim/Lsum; // verification is needed...
    autohnumax();
//    numax = 3*dim;// verification is needed...
//    hmax = 3*dim;// verification is needed...
  }
  ewald(std::vector<double> L_, double sigma_, int numax_, int hmax_){
    L = L_;
    dim = L.size();
    sigma = sigma_;
    double Lsum=0;
    for(int i=0 ; i<dim ; i++){Lsum += L[i];}
    kappa = 2.0*dim/Lsum; // verification is needed...
    numax = numax_;
    hmax = hmax_;
  }
  ewald(std::vector<double> L_, double sigma_, double kappa_){
    L = L_;
    dim = L.size();
    sigma = sigma_;
    kappa = kappa_;
    autohnumax();
//    numax = 3*dim;// verification is needed...
//    hmax = 3*dim;// verification is needed...
  }
  ewald(std::vector<double> L_, double sigma_, int numax_, int hmax_, double kappa_){
    L = L_;
    dim = L.size();
    sigma = sigma_;
    kappa = kappa_;
    numax = numax_;
    hmax = hmax_;
  }

private:
  void autohnumax(){
    std::vector<double> r_j, r_k;
    r_j.resize(dim, 0.0);
    r_k.resize(dim);
    for(int i=0 ; i<dim ; ++i){ r_k[i] = L[i]/3.0; }
    hmax = 0;
    numax = 0;
    double phi1checkOld=1.0, phi1checkNew=0.0, phi2checkOld=1.0, phi2checkNew=0.0;
    std::vector<double> realVecCheck, reciVecCheck;
    realVecCheck.resize(dim,0.0);
    reciVecCheck.resize(dim,0.0);
    while(phi1checkOld!=phi1checkNew){
    //while(fabs(phi1checkOld - phi1checkNew) < 1.0e-30){
      phi1checkOld = phi1checkNew;
      numax += 1;
      int minele = minElement(L);// minimum element of L has the most influence to Ewald summation
      for(int i=0 ; i<dim ; i++){
        if(i==minele){ realVecCheck[i] = L[i]*numax + r_k[i] - r_j[i]; }
        else{          realVecCheck[i] = r_k[i] - r_j[i]; }
      }
      double absr;
      absr = absolute(realVecCheck);
      phi1checkNew += 1.0/(pow(absr,sigma)*tgamma(sigma/2.0))*boost::math::tgamma(sigma/2.0,kappa*kappa*absr*absr);
      std::cout << phi1checkNew<<" "<<phi1checkOld<<" "<<phi1checkNew-phi1checkOld<<std::endl;
    }
    while(phi2checkOld!=phi2checkNew){
    //while(fabs(phi2checkOld - phi2checkNew) < 10e-30){
      if(hmax%2==0){phi2checkOld = phi2checkNew;}//Elements of the Ewald summation in wavespace converges with fluctuation. These elements are capable of almost zero value and/or negative-sign. To avoid excessive small cutoff, sum up two elements and average out the fluctuation.
      hmax += 1;
      int maxele = maxElement(L);// maximum element of L(minimum element of reciprocal vector) has the most linfluence to Ewald summation
      reciVecCheck[maxele] = static_cast<double>(hmax)/L[maxele];
      double absk, integral;
      absk = absolute(reciVecCheck);
      integral = pow(M_PI*absk,sigma-static_cast<double>(dim))/2.0*negativetgamma(-(sigma-static_cast<double>(dim))/2.0, M_PI*M_PI*absk*absk/kappa/kappa);
      phi2checkNew += 2.0*pow(M_PI,static_cast<double>(dim)/2.0)/tgamma(sigma/2.0)*std::cos(2.0*M_PI*(innerProduct(reciVecCheck,subtraction(r_j,r_k))))*integral;
      std::cout <<phi2checkNew<<" "<<phi2checkOld<<" "<<phi2checkNew-phi2checkOld<< std::endl;
    }
    std::cout << numax <<" "<< hmax << std::endl;
/*    if(autohnumax){//by a convergence of elements
      hmax = 0;
      numax = 0;
      double phi1check=1.0, phi2check=1.0;
      std::vector<double> realVecCheck, reciVecCheck;
      realVecCheck.resize(dim,0.0);
      reciVecCheck.resize(dim,0.0);
      while(phi1check>1.0e-300){
      //while(phi1check!=0.0){
        numax += 1;
        int minele = minElement(L);// minimum element of L has the most influence to Ewald summation
        for(int i=0 ; i<dim ; i++){
          if(i==minele){ realVecCheck[i] = L[i]*numax + r_k[i] - r_j[i]; }
          else{          realVecCheck[i] = r_k[i] - r_j[i]; }
        }
        double absr;
        absr = absolute(realVecCheck);
        phi1check = 1.0/(pow(absr,sigma)*tgamma(sigma/2.0))*boost::math::tgamma(sigma/2.0,kappa*kappa*absr*absr);
//        std::cout <<phi1check<<std::endl;
      }
      while(fabs(phi2check)>1.0e-300){
      //while(phi2check!=0.0){
        hmax += 1;
        int maxele = maxElement(L);// maximum element of L(minimum element of reciprocal vector) has the most linfluence to Ewald summation
        reciVecCheck[maxele] = static_cast<double>(hmax)/L[maxele];
        double absk, integral;
        absk = absolute(reciVecCheck);
        integral = pow(M_PI*absk,sigma-static_cast<double>(dim))/2.0*negativetgamma(-(sigma-static_cast<double>(dim))/2.0, M_PI*M_PI*absk*absk/kappa/kappa);
        phi2check = 2.0*pow(M_PI,static_cast<double>(dim)/2.0)/tgamma(sigma/2.0)*std::cos(2.0*M_PI*(innerProduct(reciVecCheck,subtraction(r_j,r_k))))*integral;
//        std::cout <<phi2check<<std::endl;
      }
//      std::cout << numax <<" "<< hmax << std::endl;
    }*/
    return;
  }
  template<typename type>
  inline int minElement(std::vector<type> x){
    type min=x[0];
    int minEleNum=0;
    for(int i=1 ; i<x.size() ; i++){
      if(min > x[i]){
        min = x[i];
        minEleNum = i;
      }
    }
    return minEleNum;
  }

public:
  double nearest(std::vector<double> r_k, std::vector<double> r_j){
    if((r_k.size()!=dim)||(r_j.size()!=dim)){std::cout << "error: dimension of r_k and r_j must be same as that of L" << std::endl; return 0;}
    double PHI=0;
    std::vector<double> realVector;
    realVector.resize(dim);
    recursiveNearest(PHI,realVector,r_k,r_j,0);
    return PHI;
  }
private:
  void recursiveNearest(double& PHI, std::vector<double> realVector, std::vector<double> r_k, std::vector<double> r_j, int currentdim){
    for(int i=-1 ; i<=1 ; i++){
      realVector[currentdim] = L[currentdim]*i + r_k[currentdim] - r_j[currentdim];
      if(currentdim==dim-1){
        if(PHI < 1.0/pow(absolute(realVector),sigma)){ PHI = 1.0/pow(absolute(realVector),sigma);}
      }
      else{
        recursiveNearest(PHI,realVector,r_k,r_j,currentdim+1);
      }
    }
  }

public:
  double simplesum(std::vector<double> r_k, std::vector<double> r_j, int max){
    if((r_k.size()!=dim)||(r_j.size()!=dim)){std::cout << "error: dimension of r_k and r_j must be same as that of L" << std::endl; return 0;}
    double PHI=0;
    std::vector<double> realVector;
    realVector.resize(dim);
    recursiveSimplesum(PHI,realVector,r_k,r_j,0,max);
    return PHI;
  }
private:
  void recursiveSimplesum(double& PHI, std::vector<double> realVector, std::vector<double> r_k, std::vector<double> r_j, int currentdim, int max){
    for(int i=-max ; i<=max ; i++){
      realVector[currentdim] = L[currentdim]*i + r_k[currentdim] - r_j[currentdim];
      if(currentdim==dim-1){
        PHI += 1.0/pow(absolute(realVector),sigma);
      }
      else{
        recursiveSimplesum(PHI,realVector,r_k,r_j,currentdim+1,max);
      }
    }
  }

public:
  double ewaldsum(std::vector<double> r_k, std::vector<double> r_j){
    if((r_k.size()!=dim)||(r_j.size()!=dim)){std::cout << "error: dimension of r_k and r_j must be same as that of L" << std::endl; return 0;}
    double phi1=0, phi2=0;
    std::vector<double> realVector, reciprocalVector;
    realVector.resize(dim);
    reciprocalVector.resize(dim);
    recursiveRealspace(phi1,realVector,r_k,r_j,0);
    recursiveWavespace(phi2,reciprocalVector,r_k,r_j,0);
    phi2 += 2.0*pow(M_PI,static_cast<double>(dim)/2.0)*pow(kappa,sigma-static_cast<double>(dim))/(tgamma(sigma/2.0)*(sigma-static_cast<double>(dim)));
    for(int i=0 ; i<dim ; i++){
      phi2 /= L[i];
    }
    return phi1+phi2;
  }

private:
  void recursiveRealspace(double& phi1, std::vector<double> realVector, std::vector<double> r_k, std::vector<double> r_j, int currentdim){
    for(int nu=-numax ; nu<=numax ; nu++){
      realVector[currentdim] = L[currentdim]*nu + r_k[currentdim] - r_j[currentdim];
      if(currentdim==dim-1){
        double absr;
        absr = absolute(realVector);
        phi1 += 1.0/(pow(absr,sigma)*tgamma(sigma/2.0))*boost::math::tgamma(sigma/2.0,kappa*kappa*absr*absr);
      }
      else{
        recursiveRealspace(phi1, realVector, r_k, r_j, currentdim+1);
      }
    }
  }

  void recursiveWavespace(double& phi2, std::vector<double> reciprocalVector, std::vector<double> r_k, std::vector<double> r_j, int currentdim){
    for(int h=-hmax ; h<=hmax ; h++){
      reciprocalVector[currentdim] = static_cast<double>(h)/L[currentdim];
      if(currentdim==dim-1){
        double absk, integral;
        absk = absolute(reciprocalVector);
//        if(absk>=0.000000000001){//It has some suspicion on reliability?(Is double absk exact zero?)
        if(absk!=0){//It has some suspicion on reliability?(Is double absk exact zero?)
          integral = pow(M_PI*absk,sigma-static_cast<double>(dim))/2.0*negativetgamma(-(sigma-static_cast<double>(dim))/2.0, M_PI*M_PI*absk*absk/kappa/kappa);
          phi2 += 2.0*pow(M_PI,static_cast<double>(dim)/2.0)/tgamma(sigma/2.0)*std::cos(2.0*M_PI*(innerProduct(reciprocalVector,subtraction(r_j,r_k))))*integral;
        }
      }
      else{
        recursiveWavespace(phi2, reciprocalVector, r_k, r_j, currentdim+1);
      }
    }
  }

  inline double absolute(double x){ return std::sqrt(x*x); }
  inline double absolute(double x, double y){ return std::sqrt(x*x +y*y); }
  inline double absolute(double x, double y, double z){ return std::sqrt(x*x+y*y+z*z); }
  template<typename type>
  inline type absolute(std::vector<type> x){
    type sum=0;
    for(int i=0 ; i<x.size() ; i++){ sum += x[i]*x[i]; }
    return std::sqrt(sum);
  }
  template<typename type>
  inline type innerProduct(std::vector<type> x, std::vector<type> y){
    type sum=0;
    if(x.size()!=y.size()){std::cout << "error: inner product between x and y cannot be defined since these vector have differ dimension" << std::endl; return sum;}
    for(int i=0 ; i<x.size() ; i++){ sum += x[i]*y[i]; }
    return sum;
  }
  template<typename type>
  inline std::vector<type> summation(std::vector<type> x, std::vector<type> y){
    std::vector<type> ans;
    ans.resize(x.size());
    if(x.size()!=y.size()){std::cout << "error: summation between x and y cannot be defined since these vector have differ dimension" << std::endl; return ans;}
    for(int i=0 ; i<x.size() ; i++){ ans[i] = x[i] + y[i]; }
    return ans;
  }
  template<typename type>
  inline std::vector<type> summation(std::vector<type> x, std::vector<type> y, std::vector<type> z){
    std::vector<type> ans;
    ans.resize(x.size());
    if(x.size()!=y.size()){std::cout << "error: summation between x, y and z cannot be defined since these vector have differ dimension" << std::endl; return ans;}
    for(int i=0 ; i<x.size() ; i++){ ans[i] = x[i] + y[i] + z[i]; }
    return ans;
  }
  template<typename type>
  inline std::vector<type> subtraction(std::vector<type> x, std::vector<type> y){
    std::vector<type> ans;
    ans.resize(x.size());
    if(x.size()!=y.size()){std::cout << "error: subtraction between x and y cannot be defined since these vector have differ dimension" << std::endl; return ans;}
    for(int i=0 ; i<x.size() ; i++){ ans[i] = x[i] - y[i]; }
    return ans;
  }
  template<typename type>
  inline int maxElement(std::vector<type> x){
    type max=x[0];
    int maxEleNum=0;
    for(int i=1 ; i<x.size() ; i++){
      if(max < x[i]){
        max = x[i];
        maxEleNum = i;
      }
    }
    return maxEleNum;
  }
  inline double negativetgamma(double a, double x){
    if(a>0.0){ return boost::math::tgamma(a,x); }
    else if(a==0.0){ return -1.0*boost::math::expint(-x); }
    else{ return -1.*pow(x,a)*exp(-x)/a + negativetgamma(a+1.,x)/a; }
  }

};



