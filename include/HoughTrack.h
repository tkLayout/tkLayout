#ifndef HOUGH_TRACK_H
#define HOUGH_TRACK_H
/*
// nice but not good for sparse data!!
//
template<class T, const int D>
class Histo {
  Histo<T, D-1>* bins_[];
  const int nbins_;
  const std::string label_;

public:
  Histo(int nbins, std::string label = "") : nbins_(nbins), label_(label), bins_(NULL) {}

  Histo<T, D>& operator()(int nbins, std::string label = "") { 
    if (!bins_) { 
      bins_ = new *Histo<T, D-1>[nbins_]; 
      std::fill(bins_, bins_+nbins_, Histo<T, D-1>(nbins, label)); 
    } else {
      for (int i = 0; i < nbins_; i++) 
        (*bins_[i])(nbins, label);
    }
    return *this;
  }

  bool inited() const { return bins_ != NULL ? bins_[0]->inited() : false; }

  Histo<T, D-1>& operator[](int index) { return *bins_[index]; }


  ~Histo() { delete[] bins_; }
};


template<class T>
class Histo<1> {
  T bins_[];
  const int nbins_;
  const std::label label_;
public:
  Histo(int nbins, std::string label = "") : nbins_(nbins), label_(label), bins_(new T[nbins]) {}

  T& operator[](int index) { return bins_[index]; }

  bool inited() const { return bins_ != NULL; }

  int nbins() const { return nbins_; }

  ~Histo() { delete[] bins_; }
  //Histo& operator()(int, std::string) { return *this; } // no op to break recursion
};
/////


template<class T, const int D>
class Histo {
public:
  class Index {
    int i_[D];
  public:
    Index(int i[D]) { std::copy(i, i+D, i_); }
    bool operator<(const Index& other) const {
      for(int i=0; i<D; i++) {
        if ((*this)[i] < other[i]) return true;
      }
      return false;
    }
    int operator[](int i) { return i_[i]; }
  };
private:
  std::map<Index, T> bins_;
public:
  void fill(const Index& index, const T& weight) { bins_[index] += w; }
  const T& get(const Index& index) const { return bins_.at(index); }
};
*/


#include <math.h>
#include <iterator>
#include <map>
#include <fstream>
#include <stdint.h>


#include <TStyle.h>
#include <TH3.h>
#include <TRandom.h>

#include <Histo.h>
#include <TrackShooter.h>
#include <ptError.h>
#include <global_funcs.h>



// Simple sparsified histogram
template<class T, int D>
class SparseMatrix {
  std::map<int, SparseMatrix<T, D-1> > elem_;
public:
  SparseMatrix<T, D-1>& operator[](int index) { return elem_[index]; }
};

template<class T>
class SparseMatrix<T, 1> {
  std::map<int, T> elem_;
public:
  T& operator[](int index) { return elem_[index]; }
};

// -----------------



struct TracksP { // holder struct for TTree export
  std::vector<unsigned>* eventn;
  std::vector<unsigned>* trackn;
  std::vector<double>* eta;
  std::vector<double>* phi0;
  std::vector<double>* z0;
  std::vector<double>* pt;
  std::vector<char>* nhits;
  TracksP() : eventn(0), trackn(0), eta(0), phi0(0), z0(0), pt(0), nhits(0) {}
};

struct HitsP { // holder struct for TTree export
  std::vector<double> *glox, *gloy, *gloz;
  std::vector<double> *locx, *locy;
  std::vector<float> *pterr, *hitprob;
  std::vector<float> *deltas;
  std::vector<char> *cnt, *z, *rho, *phi;
  HitsP() : glox(0), gloy(0), gloz(0), locx(0), locy(0), pterr(0), hitprob(0), deltas(0), cnt(0), z(0), rho(0), phi(0) {} 
};



class SmartBin {
  uint8_t count_;
  uint8_t eventid_;
  uint16_t hitmask_;
public:
  SmartBin(int count = 0, uint8_t eventid = 0, uint16_t hitmask = 0) : count_(count), eventid_(eventid), hitmask_(hitmask) {}
  SmartBin& operator+=(const SmartBin& other) {
    if (eventid_ != other.eventid_) {
      eventid_ = other.eventid_;
      hitmask_ = 0; 
    }
    if (!(hitmask_ & other.hitmask_)) {
      count_ += other.count_;
      hitmask_ |= other.hitmask_;
    }
    return *this;
  }
  SmartBin& operator=(const SmartBin& other) {
    count_ = other.count_;
    eventid_ = other.eventid_;
    hitmask_ = other.hitmask_;
    return *this;
  }
  uint8_t count() const { return count_; }
  uint8_t eventid() const { return eventid_; }
  uint16_t hitmask() const { return hitmask_; }

  operator int() const { return count_; }
};



class HoughTrack {

  SparseMatrix<ModuleData, 4> mods_;
  //TH2I* histo_;
  Histo<4, SmartBin> histo_;

  TRandom3 die_;

  double rectangularSmear(double mean, double sigma, int nsteps, int step);  

  double calcPhi0(double x, double y, double pt);
  double calcTheta(double z, double r, double z0, double pt);
  void processHit(int evid, int hitid, double x, double y, double z, double rho, double pt, double ptError, double yres);
  void loadGeometryData(TFile* infile);

  enum { H_K = 0, H_PHI0 = 1, H_Z0 = 2, H_THETA = 3 };
public:                     //     invPt,phi0,  z0, theta
  HoughTrack() : histo_((int[])   {1000,  1000, 100, 1000}, 
                        (double[]){-0.5, -3.14, -70.5,   0.},
                        (double[]){ 0.5,  3.14,  69.5, 3.14}),
                 die_(0xcafebabe) {}
  void processTree(std::string filename, long int startev, long int howmany);
  ~HoughTrack();
};



#endif
