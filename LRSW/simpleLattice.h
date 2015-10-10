class simpleLattice{
public:
  simpleLattice(){}
  virtual ~simpleLattice(){}

  std::string latticename;
  std::vector<int> cellSize; //The last vector element is equal to a number of vertices(numVert).
  int N,dim;
  std::vector<double> latticeConstant, latticeLength;
  int numCell, numVert;
  std::vector<int> numNearest;

  int num_sites(){ return N; }
  int num_nearest(int i){ 
    i%=numVert;
    return numNearest[i]; 
  }

  virtual std::vector<double> coordinate(int n)=0;

  void shiftAbs(int& n1, int& n2){ //translate n2 to n2+n1
    std::vector<int> n1_coord, n2_coord;
    n1_coord.resize(dim+1); //+1(for 3d system, n_coord[3]) are number of site in the unit cell.
    n2_coord.resize(dim+1); //for 3d system, n = n_coord[0]*cell_axes[1]*cell_axes[2]*numVert + n_coord[1]*cell_axes[2]*numVert + n_coord[2]*numVert + n_coord[3] 

    n1_coord = vectorRepresentation(n1);
    n2_coord = vectorRepresentation(n2);

    for(int i=0 ; i<dim ; i++){
      n2_coord[i] += n1_coord[i];
      n2_coord[i] %= cellSize[i];
    }

    // n1 == siteRepresentation(n1_coord);
    n2 = siteRepresentation(n2_coord);
    return;
  }

  int shiftCoordinate(int n, std::vector<int> x){ // x[x],x[y]... are relative coordinate to shift. On the other hand, x[vert] is absolute coordinate.
    std::vector<int> n_coord;
    n_coord = vectorRepresentation(n);
    for(int i=0 ; i<dim+1 ; i++){
      if(i!=dim){n_coord[i] += x[i];}
      else{      n_coord[i] = x[i]; }
      while(n_coord[i]<0){ n_coord[i] += cellSize[i]; }
      n_coord[i] %= cellSize[i];
    }
    int nout;
    nout = siteRepresentation(n_coord);
    return nout;
  }

  virtual int nearestSite(int n, int x)=0;//x is xth nearest site. it has double count.

  virtual bool discriminateNearest(int n1, int n2)=0;//it has double cout.

protected:
  int siteRepresentation(std::vector<int> n_coord){
    int n=0;
    for(int i=0 ; i<dim+1 ; i++){
      int A=1;
      for(int j=i+1 ; j<dim+1 ; j++){
        A *= cellSize[j];
      }
      n += n_coord[i]*A;
    }
    return n;
  }

  std::vector<int> vectorRepresentation(int n){
    n %= N;
    std::vector<int> n_coord;
    n_coord.resize(dim+1);
    for(int i=0 ; i<dim+1 ; i++){
      int A=1, B=1;
      for(int j=i ; j<dim+1 ; j++){
        A *= cellSize[j];
      }
      for(int j=i+1 ; j<dim+1 ; j++){
        B *= cellSize[j];
      }
      n_coord[i] = n%A/B;
    }
    return n_coord;
  }
};


class squareLattice : public simpleLattice{ //for arbitrary dimension and lattice constant
public:
  squareLattice(){}

  squareLattice(int cellSize_, int dimension){
    dim = dimension;
    latticename = "square";
    cellSize.resize(dim,cellSize_);
    numVert = 1;
    numNearest.resize(numVert,2*dim);
    cellSize.push_back(numVert);//last element of cellSize is numVert
    latticeConstant.resize(dim,1.0);
    latticeLength.resize(dim);
    numCell=1; //numCell = N/numVert
    for(int i=0 ; i<dim ; i++){
      latticeLength[i] = latticeConstant[i]*cellSize[i];
      numCell *= cellSize[i];
    }
    N = numCell*numVert;
  }
  squareLattice(std::vector<int> cellSize_, std::vector<double> latticeConstant_ ){
    if(cellSize_.size()!=latticeConstant_.size()){std::cout<<"error: size of cellSize_ and latticeConstant_ must be same"<<std::endl;}
    dim = cellSize_.size();
    latticename = "square";
    cellSize = cellSize_;
    numVert = 1;
    cellSize.push_back(numVert);//last element of cellSize is numVert
    latticeConstant = latticeConstant_;
    latticeLength.resize(dim);
    numCell=1; //numCell = N/numVert
    for(int i=0 ; i<dim ; i++){
      latticeLength[i] = latticeConstant[i]*cellSize[i];
      numCell *= cellSize[i];
    }
    N = numCell*numVert;
  }

  virtual ~squareLattice(){}

  virtual std::vector<double> coordinate(int n){
    std::vector<int> n_coord;
    n_coord = vectorRepresentation(n);
    std::vector<double> nout;
    nout.resize(dim);
    for(int i=0 ; i<dim ; i++){
      nout[i] = latticeConstant[i]*n_coord[i];
    }
    return nout;
  }

  virtual int nearestSite(int n, int x){//x is xth nearest site. it has double count.
    std::vector<int> nearest;
    nearest.resize(dim+1,0);
    if( (x<0) || (x>=2*dim) ){ std::cout << "second argument of nearestSite(n,x) must be 0<=x<2*dim"<< std::endl; exit(0); }
    if(x<dim){ nearest[x] = 1; }
    else{      nearest[x%dim] = -1; }

    return shiftCoordinate(n, nearest);
  }

  virtual bool discriminateNearest(int n1, int n2){//it has double cout.
    std::vector< std::vector<int> > nearest;
    nearest.resize(2*dim);
    for(int i=0 ; i<2*dim ; i++){
      nearest[i].resize(dim+1,0);
    }
    for(int i=0 ; i<dim ; i++){
      nearest[2*i][i] = 1;
      nearest[2*i+1][i] = -1;
    }
    for(int i=0 ; i<2*dim ; i++){
      if( n2 == shiftCoordinate(n1, nearest[i]) ){ return true; }
    }
    return false;
  }
  private:
};


class triangularLattice : public simpleLattice{
  public:
  triangularLattice(){}
  
  triangularLattice(int cellSize_){
    dim = 2;
    latticename = "triangular";
    cellSize.resize(dim);
    cellSize[0]=cellSize_;
    cellSize[1]=cellSize_/2;
    numVert = 2;
    numNearest.resize(numVert,6);
    cellSize.push_back(numVert);
    latticeConstant.resize(dim);
    latticeConstant[0]=1.0;
    latticeConstant[1]=sqrt(3.0);
    latticeLength.resize(dim);
    numCell=1; //numCell = N/numVert
    for(int i=0 ; i<dim ; i++){
      latticeLength[i] = latticeConstant[i]*cellSize[i];
      numCell *= cellSize[i];
    }
    N = numCell*numVert;
  }

  triangularLattice(std::vector<int> cellSize_, std::vector<double> latticeConstant_){
    if(!(cellSize_.size()==2 && latticeConstant_.size()==2)){std::cout<<"error: size of cellSize_ and latticeConstant_ must be 2"<<std::endl;}
    dim = 2;
    latticename = "triangular";
    cellSize = cellSize_;
    numVert = 2;
    cellSize.push_back(numVert);
    latticeConstant = latticeConstant_;
    latticeLength.resize(dim);
    numCell=1; //numCell = N/numVert
    for(int i=0 ; i<dim ; i++){
      latticeLength[i] = latticeConstant[i]*cellSize[i];
      numCell *= cellSize[i];
    }
    N = numCell*numVert;
  }

  virtual ~triangularLattice(){}

  virtual std::vector<double> coordinate(int n){
    std::vector<int> n_coord;
    n_coord = vectorRepresentation(n);
    std::vector<double> nout;
    nout.resize(dim);
    for(int i=0; i<dim ; i++){
      nout[i] = latticeConstant[i]*n_coord[i];
      if(n_coord[dim]==1){
        nout[i] += latticeConstant[i]/2.0;
      }
    }
    return nout;
  }

  virtual int nearestSite(int n, int x){//x is xth nearest site(from 0 to num_nearest-1). it has double count.
    std::vector<int> n_coord;
    n_coord = vectorRepresentation(n);
    std::vector<int> nearest;
    nearest.resize(dim+1,0);
    if(n_coord[dim+1-1]==0){
      if(x==0){
        nearest = boost::assign::list_of(0)(0)(1);
      }
      else if(x==1){
        nearest = boost::assign::list_of(1)(0)(0);
      }
      else if(x==2){
        nearest = boost::assign::list_of(0)(-1)(1);
      }
      else if(x==3){
        nearest = boost::assign::list_of(-1)(-1)(1);
      }
      else if(x==4){
        nearest = boost::assign::list_of(-1)(0)(0);
      }
      else if(x==5){
        nearest = boost::assign::list_of(-1)(0)(1);
      }
      else{
        std::cout << "error in simpleLattice.h:nearestSite"<<std::endl;
      }
    }
    else if(n_coord[dim+1-1]==1){
      if(x==0){
        nearest = boost::assign::list_of(1)(1)(0);
      }
      else if(x==1){
        nearest = boost::assign::list_of(1)(0)(1);
      }
      else if(x==2){
        nearest = boost::assign::list_of(1)(0)(0);
      }
      else if(x==3){
        nearest = boost::assign::list_of(0)(0)(0);
      }
      else if(x==4){
        nearest = boost::assign::list_of(-1)(0)(1);
      }
      else if(x==5){
        nearest = boost::assign::list_of(0)(1)(0);
      }
      else{
        std::cout << "error in simpleLattice.h:nearestSite"<<std::endl;
      }
    }
    return shiftCoordinate(n, nearest);
  }

  virtual bool discriminateNearest(int n1, int n2){
    std::vector< std::vector< std::vector<int> > > nearest;
    nearest.resize(numVert);
    for(int i=0 ; i<numVert ; i++){
      nearest[i].resize(6);
      for(int j=0 ; j<6 ; j++){
        nearest[i][j].resize(dim+1,0);
      }
    }
    nearest[0][0]=boost::assign::list_of( 0)( 0)( 1);
    nearest[0][1]=boost::assign::list_of( 1)( 0)( 0);
    nearest[0][2]=boost::assign::list_of( 0)(-1)( 1);
    nearest[0][3]=boost::assign::list_of(-1)(-1)( 1);
    nearest[0][4]=boost::assign::list_of(-1)( 0)( 0);
    nearest[0][5]=boost::assign::list_of(-1)( 0)( 1);
    nearest[1][0]=boost::assign::list_of( 1)( 1)( 0);
    nearest[1][1]=boost::assign::list_of( 1)( 0)( 1);
    nearest[1][2]=boost::assign::list_of( 1)( 0)( 0);
    nearest[1][3]=boost::assign::list_of( 0)( 0)( 0);
    nearest[1][4]=boost::assign::list_of(-1)( 0)( 1);
    nearest[1][5]=boost::assign::list_of( 0)( 1)( 0);
    if(n1%numVert==0){
      for(int i=0 ; i<6 ; i++){
        if( n2 == shiftCoordinate(n1, nearest[0][i]) ){ return true; }
      }
    }
    else if(n1%numVert==1){
      for(int i=0 ; i<6 ; i++){
        if( n2 == shiftCoordinate(n1, nearest[1][i]) ){ return true; }
      }
    }
    return false;
  }

  private:
};


