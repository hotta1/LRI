#include<stdlib.h>
#include<boost/random.hpp>
#include<vector>
#include<iostream>
#include<fstream>
#include<sstream>

extern "C"{ 
//  void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *, int, int); 
  void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *); 
  void dpotrf_(char *, int *, double *, int *, int *);
  void dgetrf_(int *, int *, double *, int *, int *, int *);
};

template<class randomDistribution> class correlatedDistribution{
public:
  //correlatedDistribution():mt(1L){} //L means long type
  //correlatedDistribution():mt(std::time(0)){}
  correlatedDistribution(){}

  std::vector<double> generate(alps::rng_helper::generator_type &gen){
    //std::random_device rnd;
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(gen);}
    sample = product(sample,transMatrix);
    return sample;
  }

  std::vector<double> genind(alps::rng_helper::generator_type &gen){
    //std::random_device rnd;
    std::vector<double> sample(N);
    for(int i=0 ; i<N ; ++i){ sample[i] = dist(gen);}
    return sample;
  }

protected:
//  boost::mt19937 mt;
  randomDistribution dist;
  std::vector<std::vector<double> > transMatrix;
  int N;
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
    
//    dsyev_( &JOBZ, &UPLO, &N, U, &N, E, work, &LWORK, &info,1,1);
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

/*  for(int i=0 ; i<N ; ++i){
    for(int j=0 ; j<N ; ++j){
      std::cout << outMatrix[i][j] <<" " ;
    }
    std::cout << std::endl;
  }*/

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

/*    for(int i=0 ; i<N ; ++i){
      for(int j=0 ; j<N ; ++j){
        std::cout << outMatrixU[i][j] <<" " ;
      }
      std::cout << std::endl;
    }
    for(int i=0 ; i<N ; ++i){
      for(int j=0 ; j<N ; ++j){
        std::cout << outMatrixL[i][j] <<" " ;
      }
      std::cout << std::endl;
    }*/

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
};

//the case using inheritance with template class, "this->" operator is needed to use member of a base class.
template<class randomDistribution> class corDistCholesky : public correlatedDistribution<randomDistribution>{
public:
  corDistCholesky(std::vector<std::vector<double> > corMatrix, randomDistribution distribution){
    this->dist = distribution;
    this->N = corMatrix.size();
    for(int i=0 ; i<this->N ; i++){ 
      if(corMatrix[i].size() != this->N){
        std::cout<<"error: corMatrix must be square matrix."<<std::endl;
        exit(0);
      }
    }
    this->dpotrf(corMatrix, this->transMatrix);

    std::vector<std::vector<double> > tmp;
    tmp.resize(this->N);
    for(int i=0 ; i<this->N ; i++){tmp[i].resize(this->N);}
    for(int i=0; i<this->N; i++ ){
      for(int j=0; j<this->N; j++ ){
        tmp[i][j] = this->transMatrix[j][i];
      }
    }

/*    for(int i=0 ; i<this->N ; ++i){
      for(int j=0 ; j<this->N ; ++j){
        std::cout << tmp[i][j] <<" " ;
      }
      std::cout << std::endl;
    }
    tmp = this->product(tmp,this->transMatrix);
    for(int i=0 ; i<this->N ; ++i){
      for(int j=0 ; j<this->N ; ++j){
        std::cout << tmp[i][j] <<" " ;
      }
      std::cout << std::endl;
    }*/

  }
};

template<class randomDistribution> class corDistEigen : public correlatedDistribution<randomDistribution>{
public:
  corDistEigen(std::vector<std::vector<double> > corMatrix, randomDistribution distribution){
    this->dist = distribution;
    this->N = corMatrix.size();
    for(int i=0 ; i<this->N ; i++){ 
      if(corMatrix[i].size() != this->N){
        std::cout<<"error: corMatrix must be square matrix."<<std::endl;
        exit(0);
      }
    }
    std::vector<std::vector<double> > EVector;
    EVector.resize(this->N);
    for(int i=0; i<this->N ; ++i){ EVector[i].resize(this->N); }
    std::vector<double> EValue;
    EValue.resize(this->N);
    
    this->dsyev(corMatrix, EVector, EValue);
    std::vector<double> sqrtEValue(this->N);
    for(int i=0 ; i<this->N ; ++i){
      if(EValue[i] < 0.0){ std::cout <<"error: eigen value of correlation matrix must be positive"<<std::endl; std::exit(1); }
      sqrtEValue[i] = std::sqrt(EValue[i]);
    }
  //  for(int i=0 ; i<this->N ; ++i){sqrtEValue[i] = std::sqrt(fabs(EValue[i]));}
    this->transMatrix.resize(this->N);
    for(int i=0 ; i<this->N ; ++i){this->transMatrix[i].resize(this->N);}
    for(int i=0 ; i<this->N ; ++i){
      for(int j=0 ; j<this->N ; ++j){
        this->transMatrix[i][j] = sqrtEValue[i]*EVector[i][j];
      }
    }
//    std::vector<std::vector<double> > sqrtDiagEValue(this->N);
//    for(int i=0 ; i<this->N ; ++i){sqrtDiagEValue[i].resize(this->N,0);}
//    for(int i=0 ; i<this->N ; ++i){sqrtDiagEValue[i][i] = std::sqrt(EValue[i]);}
//   this->transMatrix=this->product(sqrtDiagEValue,EVector);

/*    std::cout <<"tst"<<std::endl;
    for(int i=0 ; i<this->N ; ++i){
      for(int j=0 ; j<this->N ; ++j){
        std::cout << EVector[i][j] <<"  " ;
      }
      std::cout << std::endl;
    }
    std::cout <<"tst"<<std::endl;
    for(int i=0 ; i<this->N ; ++i){ std::cout << EValue[i]<<"  "; }
    std::cout <<std::endl;
    std::cout <<"tst"<<std::endl;
    int x=1;
    std::vector<double> R,L;
    R = this->product(corMatrix,EVector[x]);
    L.resize(this->N);
    for(int i=0 ; i<this->N ; ++i){ L[i] = EValue[x]*EVector[x][i]; }
    for(int i=0 ; i<this->N ; ++i){ std::cout << R[i] <<" "; }
    std::cout <<std::endl;
    for(int i=0 ; i<this->N ; ++i){ std::cout << L[i] <<" "; }
    std::cout <<std::endl;
    for(int i=0 ; i<this->N ; ++i){ std::cout << R[i]-L[i] <<" "; }
    std::cout <<std::endl<<std::endl;*/
  }
};

