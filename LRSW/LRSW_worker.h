#include <alps/parapack/worker.h>
#include <vector>
#include <boost/random.hpp>
//#include <looper/random_choice.h>
#include "random_choice.h"
#include <boost/math/special_functions/erf.hpp>   //for erfc()
#include <alps/parapack/filelock.h>
#include <boost/assign/list_of.hpp> //for boost::assign::list_of(x)(y)
#include "simpleLattice.h"//libraries must be included above header files.
#include "walkerGenerator.h"
//#include <omp.h>

//class LRSW_worker : public alps::parapack::lattice_mc_worker<> {
class LRSW_worker : public alps::parapack::mc_worker {
private:
//  typedef alps::parapack::lattice_mc_worker<> super_type;
  typedef alps::parapack::mc_worker super_type;

public:
  LRSW_worker(alps::Parameters const& params) : super_type(params) {
    T = params.value_or_default("T", 2.2);
    mcs = 0;
    MCSTEP = params.value_or_default("SWEEPS", 1 << 15);
    MCTHRM = params.value_or_default("THERMALIZATION", MCSTEP >> 3);

    epsilon = 0.0;
    int numclones;
    numclones = params.value_or_default("NUM_CLONES", 1); 
    std::string absolutePATH;
    absolutePATH = params.value_or_default("ABSOLUTEPATH", "/home/hotta/LRI/LRSW/");
/*    std::string genwalkermode;
    genwalkermode = params.value_or_default("GENERATEWALKERMODE", "off");*/
    std::string binarywalkerfile;
    binarywalkerfile = params.value_or_default("BINARYWALKERFILE", "on");
    swapconfiguration = params.value_or_default("SWAPCONFIGURATION", "on");
    interaction = params.value_or_default("INTERACTION", "LRI");
    normalization = params.value_or_default("NORMALIZATION", "off");
    latticename = params.value_or_default("LATTICE", "square lattice");
    a = params.value_or_default("a", 1.0);
    L = params.value_or_default("L", 64);
    dim = params.value_or_default("dimension", 2);
    sigma = params.value_or_default("sigma", 1.0);
    sigma += static_cast<double>(dim);
    if(latticename=="square lattice" || latticename=="square_lattice"){ Lattice = new squareLattice(L, dim); }
    else if(latticename == "triangular lattice" || latticename=="triangular_lattice"){ Lattice = new triangularLattice(L); }
    numVert = Lattice->numVert;
    numCell = Lattice->numCell;
    N = Lattice->num_sites();
    spin.resize(N, 1); // spin configuration
    sz = N;
//    histogramPartition = 16384;
    if( N <= 128 ){histogramPartition = 2.0/static_cast<double>(N);}
    else{               histogramPartition = 2.0/static_cast<double>(128);}
/*    if( N*16 <= MCSTEP*numclones ){histogramPartition = 2.0/static_cast<double>(N);}
    else{                          histogramPartition = 16.0/static_cast<double>(MCSTEP*numclones);}*/

    cluster.resize(N);  // create clusters by union find
    numDescendantNode.resize(N,1); // includes itself
    walkerGenerator genWalker(Lattice, interaction, sigma, absolutePATH);
/*    generateWalker2d genWalker(Lattice2d, interaction, sigma, absolutePATH);*/
/*    if(genwalkermode=="on"){
      makeAllSquareWalkerTable(absolutePATH);
    }*/

    if(interaction=="nearest" && dim<5){
      std::cout << "This program is worked by SW" << std::endl;
    }
    else{
      std::cout << "This program is worked by FT" << std::endl;
      if(binarywalkerfile=="on"){
      //generate file
        std::string walkerFN, JtotFN;
        walkerFN = genWalker.walkerFilename();
        JtotFN = genWalker.JtotFilename();
        std::string lockFN;
        lockFN = walkerFN;
        lockFN += ".lck";
        boost::filesystem::path lockfile(lockFN);
        alps::filelock lock(lockfile);
        lock.lock();
        std::ifstream Jtotfile, walkerfile;
        walkerfile.open(walkerFN.c_str());
        Jtotfile.open(JtotFN.c_str());
        if( !walkerfile || !Jtotfile ){
          genWalker.binaryfile();
        }
        lock.release();
        walkerfile.close();  //reopen files for the case file couldn't be detected.
        Jtotfile.close();
        walkerfile.open(walkerFN.c_str());
        Jtotfile.open(JtotFN.c_str());
        alps::IXDRFileDump idmup(walkerFN);
        walkerChoice.load(idmup);
        Jtotfile >> Jtot;
      }
      else{
        walkerChoice = genWalker.instance(Jtot);
      }
      if(normalization == "off"){
        lambdatot = 2.0*Jtot/T;
      }
      else{
        lambdatot = N/T;
      }
      pois = boost::poisson_distribution<>(lambdatot);
    }
  }

  virtual ~LRSW_worker() {delete Lattice;}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::RealObservable("Magnetization")
        << alps::RealObservable("Magnetization^2")
        << alps::RealObservable("Magnetization^3")
        << alps::RealObservable("Magnetization^4")
        << alps::RealObservable("Magnetization^5")
        << alps::RealObservable("Magnetization^6")
        << alps::RealObservable("Magnetization^7")
        << alps::RealObservable("Magnetization^8")
        << alps::RealObservable("Magnetization^9")
        << alps::RealObservable("Magnetization^10")
        << alps::RealObservable("|Magnetization|")
        << alps::RealObservable("Magnetization^2 by graph")
        << alps::RealObservable("Magnetization^4 by graph")
        << alps::RealObservable("Magnetization^6 by graph")
        << alps::RealObservable("Magnetization^8 by graph")
        << alps::RealObservable("CorrelationFunction1_L/2")
        << alps::RealObservable("CorrelationFunction1_L/4")
        << alps::RealObservable("CorrelationFunction2_L/2")
        << alps::RealObservable("CorrelationFunction2_L/4")
        << alps::RealObservable("CorrelationFunction3_L/2")
        << alps::RealObservable("CorrelationFunction3_L/4")
        << alps::RealObservable("Temperature")
        << alps::RealObservable("Number of Sites")
        << alps::RealObservable("CorrelationFunction_1")
        << alps::RealObservable("CorrelationFunction_16")
        << alps::RealObservable("CorrelationFunction_64")
        << alps::RealObservable("CorrelationFunction_256")
        << alps::HistogramObservable<double>("histogram", -1.0, 1.0, histogramPartition)
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Energy_normalized")
        << alps::RealObservable("Energy*Magnetization^2")
        << alps::RealObservable("Energy*Magnetization^4")
        << alps::RealObservable("Energy*Magnetization^6");
//        << alps::RealObservable("Second Moment");
  }

  bool is_thermalized() const { return mcs >= MCTHRM; }
  double progress()
  const { return 1.0 * mcs / (MCTHRM + MCSTEP); }

  void run(alps::ObservableSet& obs) {
    ++mcs;
    //initialize bond flags
    for(int i=0 ; i<N ; i++){
      cluster[i]=i;
      numDescendantNode[i]=1;
    }

    //make clusters
    int ktot=0;
    if( interaction=="nearest" && dim<5 ){
      double prob = 1.0 - std::exp(-2.0/T);
      if(normalization=="on"){
        int sumNumNearest = 0;
        for(int i=0 ; i<numVert ; i++){ sumNumNearest += Lattice->num_nearest(i); }
        prob = 1.0 - std::exp(-2.0*numVert/sumNumNearest/2.0/T);// /2.0 exists to make double count
      }
      for(int i=0 ; i<numCell ; i++){
        for(int j=0 ; j<numVert ; j++){
        int n1 = i*numVert + j;
          for(int k=0 ; k < Lattice->num_nearest(j) ; k++){
            int n2 = Lattice->nearestSite(n1,k);
            if(random_01()<prob && spin[n1]==spin[n2]){
              ktot++;
              unionfind(n1,n2);
            }
          }
        }
      }
    }
    else{
      int Ktot = pois(generator_01());
      for(int i=0 ; i < Ktot ; i++){
        int X, abs_n1, abs_n2, n1, n2;
        X = walkerChoice(engine());// X = N*n2_coord[vert] + n1
        abs_n1 = X%N;
        abs_n2 = (int)(numCell*random_01())*numVert;
        Lattice->shiftAbs(abs_n1, abs_n2 ,n1 ,n2);
        n2 += X/N;
        if( spin[n1]==spin[n2] ){
          ktot++;
          unionfind(n1,n2);
        }
      }
    }

    //flip spins and get clustersize_X
    double dclustersize_2=0, dclustersize_4=0, dclustersize_6=0, dclustersize_8=0;
    for(int i=0 ;i<N ; i++){
      if(cluster[i] == i){
        if( 2*random_01() < 1 ){spin[i] *= -1;}
        int size=0;
        size = numDescendantNode[i];
        dclustersize_2 += powint(size/static_cast<double>(N),2);
        dclustersize_4 += powint(size/static_cast<double>(N),4);
        dclustersize_6 += powint(size/static_cast<double>(N),6);
        dclustersize_8 += powint(size/static_cast<double>(N),8);
      }
    }
    for(int i=0 ; i<N ; i++){
      if(cluster[i] != i){
        while(cluster[cluster[i]] != cluster[i]){
          numDescendantNode[cluster[i]] -= numDescendantNode[i];
          cluster[i] = cluster[cluster[i]];
        }
        spin[i] = spin[cluster[i]];
      }
    }

    int correlation1_2, correlation1_4, correlation2_2, correlation2_4, correlation3_2, correlation3_4, correlation_1, correlation_16, correlation_64, correlation_256, second_moment=0;
    sz = 0;
    correlation1_2 = 0;
    correlation1_4 = 0;
    correlation2_2 = 0;
    correlation2_4 = 0;
    correlation3_2 = 0;
    correlation3_4 = 0;
    correlation_1 = 0;
    correlation_16 = 0;
    correlation_64 = 0;
    correlation_256 = 0;
    for(int i=0 ; i<N ; i++){
      sz += spin[i] ;
      int n1_2, n1_4, n2_2, n2_4, n3_2, n3_4, n_1, n_16, n_64, n_256;
      std::vector<int> shift;
      shift = boost::assign::list_of(L/4)(0)(0);
      n1_2 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(L/8)(0)(0);
      n1_4 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(0)(L/4)(0);
      n2_2 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(0)(L/8)(0);
      n2_4 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(L/4)(L/4)(0);
      n3_2 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(L/8)(L/8)(0);
      n3_4 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(1)(0)(0);
      n_1  = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(16)(0)(0);
      n_16 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(64)(0)(0);
      n_64 = Lattice->shiftCoordinate(i, shift);
      shift = boost::assign::list_of(256)(0)(0);
      n_256 = Lattice->shiftCoordinate(i, shift);
      if(cluster[i] == cluster[n1_2]){ correlation1_2++; }
      if(cluster[i] == cluster[n1_4]){ correlation1_4++; }
      if(cluster[i] == cluster[n2_2]){ correlation2_2++; }
      if(cluster[i] == cluster[n2_4]){ correlation2_4++; }
      if(cluster[i] == cluster[n3_2]){ correlation3_2++; }
      if(cluster[i] == cluster[n3_4]){ correlation3_4++; }
      if(cluster[i] == cluster[n_1]){  correlation_1++; }
      if(cluster[i] == cluster[n_16]){ correlation_16++; }
      if(cluster[i] == cluster[n_64]){ correlation_64++; }
      if(cluster[i] == cluster[n_256]){ correlation_256++; }
//      second_moment += L/(2.0*M_PI)*std::sqrt(spin[i]/spin[n_1]-1.0);
    }

    double dsz = sz / static_cast<double>(N);
    double dcorrelation1_2 = correlation1_2/static_cast<double>(N);
    double dcorrelation1_4 = correlation1_4/static_cast<double>(N);
    double dcorrelation2_2 = correlation2_2/static_cast<double>(N);
    double dcorrelation2_4 = correlation2_4/static_cast<double>(N);
    double dcorrelation3_2 = correlation3_2/static_cast<double>(N);
    double dcorrelation3_4 = correlation3_4/static_cast<double>(N);
    double dcorrelation_1  = correlation_1/static_cast<double>(N);
    double dcorrelation_16 = correlation_16/static_cast<double>(N);
    double dcorrelation_64 = correlation_64/static_cast<double>(N);
    double dcorrelation_256 = correlation_256/static_cast<double>(N);
//    double dsecond_moment = second_moment/static_cast<double>(N);
    obs["Magnetization"] << dsz;
    obs["Magnetization^2"] << dsz * dsz;
    obs["Magnetization^3"] << dsz * dsz * dsz;
    obs["Magnetization^4"] << dsz * dsz * dsz * dsz;
    obs["Magnetization^5"] << dsz * dsz * dsz * dsz * dsz;
    obs["Magnetization^6"] << dsz * dsz * dsz * dsz * dsz * dsz;
    obs["Magnetization^7"] << dsz * dsz * dsz * dsz * dsz * dsz * dsz;
    obs["Magnetization^8"] << dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz;
    obs["Magnetization^9"] << dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz;
    obs["Magnetization^10"] << dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz * dsz;
    obs["|Magnetization|"] << std::abs(dsz);
    double m2g, m4g, m6g, m8g;
    m2g = dclustersize_2;
    m4g = 3*dclustersize_2*dclustersize_2 - 2*dclustersize_4;
    m6g = 15*dclustersize_2*dclustersize_2*dclustersize_2 - 30*dclustersize_4*dclustersize_2 + 16*dclustersize_6;
    m8g = 105*dclustersize_2*dclustersize_2*dclustersize_2*dclustersize_2 +448*dclustersize_2*dclustersize_6 -420*dclustersize_2*dclustersize_2*dclustersize_4 +140*dclustersize_4*dclustersize_4 -272*dclustersize_8;
    obs["Magnetization^2 by graph"] << m2g;
    obs["Magnetization^4 by graph"] << m4g;
    obs["Magnetization^6 by graph"] << m6g;
    obs["Magnetization^8 by graph"] << m8g;
    obs["CorrelationFunction1_L/2"] << dcorrelation1_2;
    obs["CorrelationFunction1_L/4"] << dcorrelation1_4;
    obs["CorrelationFunction2_L/2"] << dcorrelation2_2;
    obs["CorrelationFunction2_L/4"] << dcorrelation2_4;
    obs["CorrelationFunction3_L/2"] << dcorrelation3_2;
    obs["CorrelationFunction3_L/4"] << dcorrelation3_4;
    obs["CorrelationFunction_1"] << dcorrelation_1;
    obs["CorrelationFunction_16"] << dcorrelation_16;
    obs["CorrelationFunction_64"] << dcorrelation_64;
    obs["CorrelationFunction_256"] << dcorrelation_16;
    obs["Temperature"] << T;
    obs["Number of Sites"] << static_cast<double>(N);//observable must be real number
    obs["histogram"] << dsz;
    double energy = (1.0+epsilon)*Jtot-T*ktot;
    obs["Energy"] << energy;
    obs["Energy^2"] << energy*energy;
    obs["Energy_normalized"] << energy/Jtot;//(1.0+epsilon) - T*ktot/Jtot;
    obs["Energy*Magnetization^2"] << energy*m2g;
    obs["Energy*Magnetization^4"] << energy*m4g;
    obs["Energy*Magnetization^6"] << energy*m6g;
//    obs["Second Moment"] << dsecond_moment;
    return;
  }


  void save(alps::ODump& dp) const { 
    if(swapconfiguration=="off"){}
    else{dp << mcs << spin;}
  }
  void load(alps::IDump& dp) {
    if(swapconfiguration=="off"){}
    else{dp >> mcs >> spin;}
  }

private:
  double T;
  int mcs;
  int MCSTEP;
  int MCTHRM;
  int N, dim; // number of lattice sites
  std::vector<int> spin; // spin configuration
  int sz;

  double Jtot, epsilon;
  int numVert, numCell;
  std::string latticename;
  std::string interaction, normalization, swapconfiguration;
  double histogramPartition;
  double sigma, a;
  int L;
  std::vector<double> latticeLength;
  simpleLattice* Lattice;
  std::vector<int> cluster; //for memorize cluster
  std::vector<int> numDescendantNode;
  double lambdatot;
  boost::poisson_distribution<> pois;
  looper::random_choice walkerChoice; //make a function to choose at assigned probability by walker's method
  void unionfind(int n1, int n2){
    while(cluster[n1]!=cluster[cluster[n1]]){
      numDescendantNode[cluster[n1]] -= numDescendantNode[n1];
      cluster[n1] = cluster[cluster[n1]];
    }
    while(cluster[n2]!=cluster[cluster[n2]]){
      numDescendantNode[cluster[n2]] -= numDescendantNode[n2];
      cluster[n2] = cluster[cluster[n2]];
    }
    if(cluster[n1]!=cluster[n2]){
      if(numDescendantNode[cluster[n1]] < numDescendantNode[cluster[n2]]){ 
        numDescendantNode[cluster[n2]] += numDescendantNode[cluster[n1]];
        cluster[cluster[n1]] = cluster[n2];
      }
      else{
        numDescendantNode[cluster[n1]] += numDescendantNode[cluster[n2]];
        cluster[cluster[n2]] = cluster[n1];
      }
    }
    return;
  }

  template<typename type>
  type powint(type x, int n){
    type ans=1;
    if(n>=0){
      for(int i=0 ; i<n ; i++){ ans*=x; }
    }
    else{
      for(int i=0 ; i>n ; i--){ ans/=x; }
    }
    return ans;
  }


  void makeAllSquareWalkerTable(std::string absolutePATH){
    double a_;
    std::string interaction_;
    a_ = 1.0;
    interaction_ = "LRI";
    for(int i=3 ; i<15 ; i++){
      int L_;
      L_ = 1 << i;
      for(int j=8 ; j<20 ; j++){
        simpleLattice* Lattice_;
        Lattice_ = new squareLattice(L_, a_);
        double sigma_ = j*0.25;
//        if( ( sigma_ < 2.0+delta && sigma_ > 2.0-delta ) || ( sigma_ < 4.0+delta && sigma_ > 4.0-delta ) || ( sigma_ < 6.0+delta && sigma_ > 6.0-delta ) )
        if(j%8==0){
          double sigma_even;
          sigma_even = sigma_ + 0.0001;
          walkerGenerator genWalker(Lattice_, interaction_, sigma_even, absolutePATH);
          genWalker.binaryfile();
          std::cout << sigma_even << std::endl;
        }
        else{
          walkerGenerator genWalker(Lattice_, interaction_, sigma_, absolutePATH);
          genWalker.binaryfile();
          std::cout << sigma_ << std::endl;
        }
        delete Lattice_;
      }
    }
    return;
  }
};


class LRSW_evaluator : public alps::parapack::simple_evaluator {
public:
  LRSW_evaluator(alps::Parameters const&) {}
  void evaluate(alps::ObservableSet& obs) const {
    alps::RealObsevaluator m2 = obs["Magnetization^2 by graph"];
    alps::RealObsevaluator m4 = obs["Magnetization^4 by graph"];
    alps::RealObsevaluator m6 = obs["Magnetization^6 by graph"];
    alps::RealObsevaluator m8 = obs["Magnetization^8 by graph"];
    double dev[8] = {0.214002, 0.18885, 0.0938762, 0.262266, 0.257176, 0.276397, 0.162642, 0.221176};
//    alps::RealObsevaluator m2 = obs["Magnetization^2"];
//    alps::RealObsevaluator m4 = obs["Magnetization^4"];
//    alps::RealObsevaluator m6 = obs["Magnetization^6"];
//    alps::RealObsevaluator m8 = obs["Magnetization^8"];

    alps::RealObsevaluator binder("Binder Ratio of Magnetization");
    binder = m2*m2 / m4;
    obs.addObservable(binder);
    alps::RealObsevaluator binder2("Binder Ratio of Magnetization 2");
    binder2 = m2*m2*m2 / m6;
    obs.addObservable(binder2);
    alps::RealObsevaluator binder3("Binder Ratio of Magnetization 3");
    binder3 = m2*m2*m2*m2 / m8;
    obs.addObservable(binder3);
    alps::RealObsevaluator binder4("Binder Ratio of Magnetization 4");
    binder4 = m4*m4 / m8;
    obs.addObservable(binder4);
    alps::RealObsevaluator binder5("Binder Ratio of Magnetization 5");
    binder5 = m4*m2 / m6;
    obs.addObservable(binder5);
    alps::RealObsevaluator binder6("Binder Ratio of Magnetization 6");
    binder6 = m6*m2 / m8;
    obs.addObservable(binder6);
    alps::RealObsevaluator binder7("Binder Ratio of Magnetization 7");
    binder7 = m2*m2*m4 / m8;
    obs.addObservable(binder7);
    alps::RealObsevaluator binder8("Binder Ratio of Magnetization 8");
    binder8 = m4*m4 / (m6*m2);
    obs.addObservable(binder8);

    alps::RealObsevaluator binderinv("Binder Ratio of Magnetization Inverse");
    binderinv = m4 / (m2*m2);
    obs.addObservable(binderinv);
    alps::RealObsevaluator combinderinv("Combined Binder Inverse");
    combinderinv = binder/0.456947 + binderinv*0.456947;
    obs.addObservable(combinderinv);

    alps::RealObsevaluator binderinv2("Binder Ratio of Magnetization Inverse 2");
    binderinv2 = m6 / (m4*m2);
    obs.addObservable(binderinv2);
    alps::RealObsevaluator combinderinv2("Combined Binder Inverse 2");
    combinderinv2 = binder5*3.0 + binderinv2/3.0;
    obs.addObservable(combinderinv2);

    alps::RealObsevaluator combinderinv3dNN("Combined Binder Inverse 3dNN");
    combinderinv3dNN = binder/0.635 + binderinv*0.635;
    obs.addObservable(combinderinv3dNN);

    alps::RealObsevaluator combinderinv2dNN("Combined Binder Inverse 2dNN");
    combinderinv2dNN = binder/0.856216 + binderinv*0.856216;
    obs.addObservable(combinderinv2dNN);

    alps::RealObsevaluator combinderopt123("Combined Binder Optimized 123");
    combinderopt123 = 27.5241*binder - 42.2937*binder2 + 22.3375*binder3;
    obs.addObservable(combinderopt123);
    alps::RealObsevaluator combinderopt456("Combined Binder Optimized 456");
    combinderopt456 = 8.61608*binder4 + 52.7925*binder5 - 57.2968*binder6;
    obs.addObservable(combinderopt456);
    alps::RealObsevaluator combinderopt148("Combined Binder Optimized 148");
    combinderopt148 = 25.3101*binder - 13.8139*binder4 - 8.10894*binder8;
    obs.addObservable(combinderopt148);
    alps::RealObsevaluator combinderopt125("Combined Binder Optimized 125");
    combinderopt125 = 42.1668*binder + 2.77314*binder2 - 37.1244*binder5;
    obs.addObservable(combinderopt125);
    alps::RealObsevaluator combinder12("Combined Binder 12");
    combinder12 = binder - dev[0]/dev[1]*binder2;
    obs.addObservable(combinder12);
    alps::RealObsevaluator combinder13("Combined Binder 13");
    combinder13 = binder - dev[0]/dev[2]*binder3;
    obs.addObservable(combinder13);
    alps::RealObsevaluator combinder14("Combined Binder 14");
    combinder14 = binder - dev[0]/dev[3]*binder4;
    obs.addObservable(combinder14);
    alps::RealObsevaluator combinder15("Combined Binder 15");
    combinder15 = binder - dev[0]/dev[4]*binder5;
    obs.addObservable(combinder15);
    alps::RealObsevaluator combinder16("Combined Binder 16");
    combinder16 = binder - dev[0]/dev[5]*binder6;
    obs.addObservable(combinder16);
    alps::RealObsevaluator combinder17("Combined Binder 17");
    combinder17 = binder - dev[0]/dev[6]*binder7;
    obs.addObservable(combinder17);
    alps::RealObsevaluator combinder18("Combined Binder 18");
    combinder18 = binder - dev[0]/dev[7]*binder8;
    obs.addObservable(combinder18);
    alps::RealObsevaluator combinder23("Combined Binder 23");
    combinder23 = binder2 - dev[1]/dev[2]*binder3;
    obs.addObservable(combinder23);
    alps::RealObsevaluator combinder24("Combined Binder 24");
    combinder24 = binder2 - dev[1]/dev[3]*binder4;
    obs.addObservable(combinder24);
    alps::RealObsevaluator combinder25("Combined Binder 25");
    combinder25 = binder2 - dev[1]/dev[4]*binder5;
    obs.addObservable(combinder25);
    alps::RealObsevaluator combinder26("Combined Binder 26");
    combinder26 = binder2 - dev[1]/dev[5]*binder6;
    obs.addObservable(combinder26);
    alps::RealObsevaluator combinder27("Combined Binder 27");
    combinder27 = binder2 - dev[1]/dev[6]*binder7;
    obs.addObservable(combinder27);
    alps::RealObsevaluator combinder28("Combined Binder 28");
    combinder28 = binder2 - dev[1]/dev[7]*binder8;
    obs.addObservable(combinder28);
    alps::RealObsevaluator combinder34("Combined Binder 34");
    combinder34 = binder3 - dev[2]/dev[3]*binder4;
    obs.addObservable(combinder34);
    alps::RealObsevaluator combinder35("Combined Binder 35");
    combinder35 = binder3 - dev[2]/dev[4]*binder5;
    obs.addObservable(combinder35);
    alps::RealObsevaluator combinder36("Combined Binder 36");
    combinder36 = binder3 - dev[2]/dev[5]*binder6;
    obs.addObservable(combinder36);
    alps::RealObsevaluator combinder37("Combined Binder 37");
    combinder37 = binder3 - dev[2]/dev[6]*binder7;
    obs.addObservable(combinder37);
    alps::RealObsevaluator combinder38("Combined Binder 38");
    combinder38 = binder3 - dev[2]/dev[7]*binder8;
    obs.addObservable(combinder38);
    alps::RealObsevaluator combinder45("Combined Binder 45");
    combinder45 = binder4 - dev[3]/dev[4]*binder5;
    obs.addObservable(combinder45);
    alps::RealObsevaluator combinder46("Combined Binder 46");
    combinder46 = binder4 - dev[3]/dev[5]*binder6;
    obs.addObservable(combinder46);
    alps::RealObsevaluator combinder47("Combined Binder 47");
    combinder47 = binder4 - dev[3]/dev[6]*binder7;
    obs.addObservable(combinder47);
    alps::RealObsevaluator combinder48("Combined Binder 48");
    combinder48 = binder4 - dev[3]/dev[7]*binder8;
    obs.addObservable(combinder48);
    alps::RealObsevaluator combinder56("Combined Binder 56");
    combinder56 = binder5 - dev[4]/dev[5]*binder6;
    obs.addObservable(combinder56);
    alps::RealObsevaluator combinder57("Combined Binder 57");
    combinder57 = binder5 - dev[4]/dev[6]*binder7;
    obs.addObservable(combinder57);
    alps::RealObsevaluator combinder58("Combined Binder 58");
    combinder58 = binder5 - dev[4]/dev[7]*binder8;
    obs.addObservable(combinder58);
    alps::RealObsevaluator combinder67("Combined Binder 67");
    combinder67 = binder6 - dev[5]/dev[6]*binder7;
    obs.addObservable(combinder67);
    alps::RealObsevaluator combinder68("Combined Binder 68");
    combinder68 = binder6 - dev[5]/dev[7]*binder8;
    obs.addObservable(combinder68);
    alps::RealObsevaluator combinder78("Combined Binder 78");
    combinder78 = binder7 - dev[6]/dev[7]*binder8;
    obs.addObservable(combinder78);

    alps::RealObsevaluator combinder12inv("Combinde Binder 12 Inverse");
    combinder12inv = 1.0/combinder12;
    obs.addObservable(combinder12inv);
    alps::RealObsevaluator comcombinder12inv("Combined Combined Binder 12 Inverse");
    comcombinder12inv = combinder12/0.2843448 + combinder12inv*0.2843448;
    obs.addObservable(comcombinder12inv);
/*    alps::RealObsevaluator combinder("Combined Binder");
    combinder = binder - 0.214002/0.18885*binder2;
    obs.addObservable(combinder);
    alps::RealObsevaluator combinder2("Combined Binder 2");
    combinder2 = binder - 0.214002/0.0938762*binder3;
    obs.addObservable(combinder2);
    alps::RealObsevaluator combinder3("Combined Binder 3");
    combinder3 = binder2 - 0.18885/0.0938762*binder3;
    obs.addObservable(combinder3);
    alps::RealObsevaluator combinder4("Combined Binder 4");
    combinder4 = binder2 - 0.18885/0.262266*binder4;
    obs.addObservable(combinder4);
    alps::RealObsevaluator combinder5("Combined Binder 5");
    combinder5 = binder2 - 0.18885/0.257176*binder5;
    obs.addObservable(combinder5);
    alps::RealObsevaluator combinder6("Combined Binder 6");
    combinder6 = binder2 - 0.18885/0.276397*binder6;
    obs.addObservable(combinder6);
    alps::RealObsevaluator combinder7("Combined Binder 7");
    combinder7 = binder2 - 0.18885/0.162642*binder7;
    obs.addObservable(combinder7);
    alps::RealObsevaluator combinder8("Combined Binder 8");
    combinder8 = binder2 - 0.18885/0.221176*binder8;
    obs.addObservable(combinder8);*/


    alps::RealObsevaluator T = obs["Temperature"];
    alps::RealObsevaluator energy = obs["Energy"];
    alps::RealObsevaluator eneM2 = obs["Energy*Magnetization^2"];
    alps::RealObsevaluator eneM4 = obs["Energy*Magnetization^4"];
    alps::RealObsevaluator eneM6 = obs["Energy*Magnetization^6"];
    alps::RealObsevaluator derbinder1K("Derivative of Binder1 by K");
    derbinder1K = 2*m2*(energy*m2 - eneM2)/m4 - m2*m2*(energy*m4-eneM4)/(m4*m4);
    obs.addObservable(derbinder1K);
    alps::RealObsevaluator derbinder2K("Derivative of Binder2 by K");
    derbinder2K = ((energy*m2 - eneM2)*3*m2*m2*m6 - m2*m2*m2*(energy*m6 - eneM6))/(m6*m6);
    obs.addObservable(derbinder2K);
    alps::RealObsevaluator derbinder1T("Derivative of Binder1 by T");
    derbinder1T = -1*derbinder1K/(T*T);
    obs.addObservable(derbinder1T);
    alps::RealObsevaluator derbinder2T("Derivative of Binder2 by T");
    derbinder2T = -1*derbinder2K/(T*T);
    obs.addObservable(derbinder2T);
    alps::RealObsevaluator derbinder1invK("Derivative of Binder1 inverse by K");
    derbinder1invK = (energy*m4-eneM4)/(m2*m2) - 2*m4*(energy*m2 - eneM2)/(m2*m2*m2);
    obs.addObservable(derbinder1invK);
    alps::RealObsevaluator dercombinderK("Derivative of Combined Binder by K");
    dercombinderK = derbinder1K - 0.214002/0.18885*derbinder2K;
    obs.addObservable(dercombinderK);
    alps::RealObsevaluator dercombinderT("Derivative of Combined Binder by T");
    dercombinderT = (-1)*dercombinderK/(T*T);
    obs.addObservable(dercombinderT);
    alps::RealObsevaluator dercombinderinvK("Derivative of Combined Binder inverse by K");
    dercombinderinvK = derbinder1K/0.456947 + derbinder1invK*0.456947;
    obs.addObservable(dercombinderinvK);
    alps::RealObsevaluator dercombinderinvT("Derivative of Combined Binder inverse by T");
    dercombinderinvT = (-1)*dercombinderinvK/(T*T);
    obs.addObservable(dercombinderinvT);

    alps::RealObsevaluator mabs = obs["|Magnetization|"];
    alps::RealObsevaluator N = obs["Number of Sites"];
    alps::RealObsevaluator energy2 = obs["Energy^2"];
    alps::RealObsevaluator magsus("Magnetic Susceptibility");
    magsus = (m2-mabs*mabs)/T*N; // (<M^2>-<M>^2)/T/V = (<m^2>-<m>^2)*V^2/T/V = (<m^2>-<m>^2)/T*V 
    obs.addObservable(magsus);
    alps::RealObsevaluator specheat("Specific Heat");
    specheat = (energy2-energy*energy)/(T*T*N);// (<E^2>-<E>^2)/(T*T)/V
    obs.addObservable(specheat);

    alps::RealObsevaluator C1_2 = obs["CorrelationFunction1_L/2"];
    alps::RealObsevaluator C1_4 = obs["CorrelationFunction1_L/4"];
    alps::RealObsevaluator C2_2 = obs["CorrelationFunction2_L/2"];
    alps::RealObsevaluator C2_4 = obs["CorrelationFunction2_L/4"];
    alps::RealObsevaluator C3_2 = obs["CorrelationFunction3_L/2"];
    alps::RealObsevaluator C3_4 = obs["CorrelationFunction3_L/4"];
    alps::RealObsevaluator CorRatio1("Correlation Ratio 1");
    alps::RealObsevaluator CorRatio2("Correlation Ratio 2");
    alps::RealObsevaluator CorRatio3("Correlation Ratio 3");
    CorRatio1 = C1_2 / C1_4;
    CorRatio2 = C2_2 / C2_4;
    CorRatio3 = C3_2 / C3_4;
    obs.addObservable(CorRatio1);
    obs.addObservable(CorRatio2);
    obs.addObservable(CorRatio3);
  }
};

