#include<stdlib.h>
#include<boost/random.hpp>
#include<vector>
#include<iostream>
#include<fstream>
#include<sstream>
#include<alps/parapack/filelock.h>
#include <alps/parapack/worker.h>

extern "C"{ 
//  void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *, int, int); 
  void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *); 
  void dpotrf_(char *, int *, double *, int *, int *);
  void dgetrf_(int *, int *, double *, int *, int *, int *);
};


template<class randomDistribution> class correlatedDistributionGenerator{
public:
  correlatedDistributionGenerator(simpleLattice* Lattice_,  randomDistribution distribution, std::string absolutePATH_):mt(std::time(0)){
    Lattice = Lattice_;
    N = Lattice->num_sites();
    corMatrix.resize(N);
    for(int i=0 ; i<N ; ++i){ corMatrix[i].resize(N); }
    absolutePATH = absolutePATH_;
    dist = distribution;
  }

  std::vector<std::vector<double> > corAlgebra(double rho_){
    rho = rho_;
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
        }
        else{
          std::vector<double> ci = Lattice->coordinate(j);
          //corMatrix[i][j] = ewaldsum.ewaldsum(origin,ci);
          corMatrix[i][j] = ewaldsum.nearest(origin,ci)/3.0;
        }
        //std::cout <<corMatrix[i][j]<<"  ";
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
    corType = "algebra";
    return corMatrix;
  }

  void eigen(std::string binarymatrixfile){
    if(binarymatrixfile=="on"){
      std::string transMatrixFN;
      transMatrixFN = transMatrixFilename("EIGEN");
      std::string lockFN;
      lockFN = transMatrixFN;
      lockFN += ".lck";
      boost::filesystem::path lockfile(lockFN);
      alps::filelock lock(lockfile);
      lock.lock();
      std::ifstream transMatrixFile;
      transMatrixFile.open(transMatrixFN.c_str());
      if( !transMatrixFile ){
        calcEigenTrans();
        alps::OXDRFileDump odump(transMatrixFN);
        odump << transMatrix;
        std::cout << "saving " << transMatrixFN <<std::endl;
      }
      lock.release();
      transMatrixFile.close();  //reopen files for the case file couldn't be detected.
      transMatrixFile.open(transMatrixFN.c_str());
      alps::IXDRFileDump idump(transMatrixFN);
      idump >> transMatrix;
      std::cout << "loading " << transMatrixFN <<std::endl;
    }
    else{
      calcEigenTrans();
    }
    return;
  }

  void cholesky(std::string binarymatrixfile){
    if(binarymatrixfile=="on"){
      std::string transMatrixFN;
      transMatrixFN = transMatrixFilename("CHOLESKY");
      std::string lockFN;
      lockFN = transMatrixFN;
      lockFN += ".lck";
      boost::filesystem::path lockfile(lockFN);
      alps::filelock lock(lockfile);
      lock.lock();
      std::ifstream transMatrixFile;
      transMatrixFile.open(transMatrixFN.c_str());
      if( !transMatrixFile ){
        calcCholeskyTrans();
        alps::OXDRFileDump odump(transMatrixFN);
        odump << transMatrix;
        std::cout << "saving " << transMatrixFN <<std::endl;
      }
      lock.release();
      transMatrixFile.close();  //reopen files for the case file couldn't be detected.
      transMatrixFile.open(transMatrixFN.c_str());
      alps::IXDRFileDump idump(transMatrixFN);
      idump >> transMatrix;
      std::cout << "loading " << transMatrixFN <<std::endl;
    }
    else{
      calcCholeskyTrans();
    }
    return;
  }

  std::vector<double> generate(alps::rng_helper::generator_type &gen){
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(gen);}
    sample = product(sample,transMatrix);
    return sample;
  }

  std::vector<double> genind(alps::rng_helper::generator_type &gen){
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(gen);}
    return sample;
  }

  std::vector<double> generate(){
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(mt);}
    sample = product(sample,transMatrix);
    return sample;
  }

  std::vector<double> genind(){
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(mt);}
    return sample;
  }


protected:
  double rho;
  std::string corType;
  std::string absolutePATH;
  std::vector<std::vector<double> > corMatrix;
  simpleLattice* Lattice;
  int N;
  boost::mt19937 mt;
  randomDistribution dist;
  std::vector<std::vector<double> > transMatrix;

  std::vector<std::vector<double> > calcEigenTrans(){
    std::vector<std::vector<double> > EVector;
    EVector.resize(N);
    for(int i=0; i<N ; ++i){ EVector[i].resize(N); }
    std::vector<double> EValue;
    EValue.resize(N);
    
    dsyev(corMatrix, EVector, EValue);
    std::vector<double> sqrtEValue(N);
    for(int i=0 ; i<N ; ++i){
      if(EValue[i] < 0.0){ std::cout <<"error: eigen value of correlation matrix must be positive"<<std::endl; std::exit(1); }
      sqrtEValue[i] = std::sqrt(EValue[i]);
    }
    transMatrix.resize(N);
    for(int i=0 ; i<N ; ++i){transMatrix[i].resize(N);}
    for(int i=0 ; i<N ; ++i){
      for(int j=0 ; j<N ; ++j){
        transMatrix[i][j] = sqrtEValue[i]*EVector[i][j];
      }
    }
    std::cout << "creating eigen transMatrix" <<std::endl;
    return transMatrix;
  }

  std::vector<std::vector<double> > calcCholeskyTrans(){
    dpotrf(corMatrix, transMatrix);
    std::cout << "creating cholesky transMatrix" <<std::endl;
    return transMatrix;
  }

  int dsyev(std::vector<std::vector<double> > &inMatrix, std::vector<std::vector<double> > &EVector, std::vector<double> &EValue){
    int N = inMatrix.size();
    int info;
    double work[4*N];
    double *U, *E;
    U = (double *)malloc(sizeof(double)*N*N);
    E = (double *)malloc(sizeof(double)*N);
    
    for( int i=0; i<N; i++ ){
      for( int j=0; j<N; j++ ){
        U[i*N+j] = inMatrix[j][i];
      }
    }
    char JOBZ ='V', UPLO='U';
    int LWORK=N*4;
    
    dsyev_( &JOBZ, &UPLO, &N, U, &N, E, work, &LWORK, &info);
    
    for( int i=0; i<N; i++ ){
      for( int j=0; j<N; j++ ){
        EVector[i][j] = U[i*N+j];//E[i] is ith eigen vector
      }
    }
    
    for( int i=0 ; i<N ; i++){EValue[i] = E[i];}
    
    return info;
  }

  int dpotrf(std::vector<std::vector<double> > &inMatrix, std::vector<std::vector<double> > &outMatrix){
    int N = inMatrix.size();
    int info;
    double *U;
    U = (double *)malloc(sizeof(double)*N*N);
    
    for( int i=0; i<N; i++ ){
      for( int j=0; j<N; j++ ){
        U[i*N+j] = inMatrix[j][i];
      }
    }
    char UPLO='U';
    
    dpotrf_(&UPLO, &N, U, &N, &info);

    outMatrix.resize(N);
    for(int i=0 ; i<N ; i++){ outMatrix[i].resize(N,0.0); }
    
    if(UPLO=='U'){
      for( int i=0; i<N; i++ ){
        for( int j=i; j<N; j++ ){
          outMatrix[i][j] = U[i*N+j];
        }
      }
    }
    else if(UPLO=='L'){
      for( int i=0; i<N; i++ ){
        for( int j=0; j<=i; j++ ){
          outMatrix[i][j] = U[i*N+j];
        }
      }
    }
    return info;
  }

  int dgetrf(std::vector<std::vector<double> > &inMatrix, std::vector<std::vector<double> > &outMatrixU, std::vector<std::vector<double> > &outMatrixL){
    int N = inMatrix.size();
    int info;
    double *U;
    U = (double *)malloc(sizeof(double)*N*N);
    
    for( int i=0; i<N; i++ ){
      for( int j=0; j<N; j++ ){
        U[i*N+j] = inMatrix[j][i];
      }
    }
    int *IPIV;
    IPIV = (int *)malloc(sizeof(int)*N);

    dgetrf_(&N, &N, U, &N, IPIV, &info);

    outMatrixU.resize(N);
    for(int i=0 ; i<N ; i++){ outMatrixU[i].resize(N,0.0); }
    outMatrixL.resize(N);
    for(int i=0 ; i<N ; i++){ outMatrixL[i].resize(N,0.0); }

    for( int i=0; i<N; i++ ){
      for( int j=i; j<N; j++ ){
        outMatrixU[i][j] = U[i*N+j];
      }
    }
    for( int i=0; i<N; i++ ){
      for( int j=0; j<=i; j++ ){
        if(i==j){outMatrixL[i][j] = 1;}
        else{outMatrixL[i][j] = U[i*N+j];}
      }
    }
    return info;
  }

  template<class Type>
  std::vector<std::vector<Type> > product(std::vector<std::vector<Type> > A, std::vector<std::vector<Type> > B){
    for(int i=0 ; i<A.size() ; ++i){
      if( A[i].size() != B.size() ){ std::cout << "Number of row A and column B must be same." <<std::endl; exit(0);}
    }
    std::vector<std::vector<Type> > prd;
    prd.resize(A.size());
    for(int i=0 ; i<A.size() ; ++i){ prd[i].resize(B[i].size(),0.0); }
    for(int i=0 ; i<A.size() ; ++i){
      for(int j=0 ; j<B[i].size() ; ++j){
        for(int k=0 ; k<A[i].size() ; ++k){ prd[i][j] += A[i][k]*B[k][j]; }
      }
    }
    return prd;
  }
  template<class Type>
  std::vector<Type> product(std::vector<std::vector<Type> > A, std::vector<Type> B){
    for(int i=0 ; i<A.size() ; ++i){
      if( A[i].size() != B.size() ){ std::cout << "Number of row A and column B must be same." <<std::endl; exit(0);}
    }
    std::vector<Type > prd;
    prd.resize(A.size(),0.0);
    for(int i=0 ; i<A.size() ; ++i){
      for(int k=0 ; k<B.size() ; ++k){ prd[i] = prd[i] + A[i][k]*B[k]; }
    }
    return prd;
  }
  template<class Type>
  std::vector<Type> product(std::vector<Type> A, std::vector<std::vector<Type> > B){
    for(int i=0 ; i<A.size() ; ++i){
      if( A.size() != B.size() ){ std::cout << "Number of column A and row B must be same." <<std::endl; exit(0);}
    }
    std::vector<Type > prd;
    prd.resize(A.size(),0.0);
    for(int i=0 ; i<B[0].size() ; ++i){
      for(int k=0 ; k<B.size() ; ++k){ prd[i] = prd[i] + A[k]*B[k][i]; }
    }
    return prd;
  }

  std::string transMatrixFilename(std::string decompositionType){
    std::ostringstream filename;
    filename << absolutePATH << "transMatrix/";
    filename << Lattice->dim <<"D"<< Lattice->latticename <<"_";
    for(int i=0 ; i<Lattice->dim ; i++){
      filename << "L" << i+1 << "-" << Lattice->cellSize[i]<<"_";
    }
    for(int i=0 ; i<Lattice->dim ; i++){
      filename << "a" << i+1 << "-" << Lattice->latticeConstant[i] <<"_";
    }
    if(corType == "algebra"){
      filename << "rho-" << rho;
    }
    else{
    }
    filename << decompositionType;
    filename << ".xdr" ;
    std::string tmp;
    tmp = filename.str();
    return tmp;
  }

};

