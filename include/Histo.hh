#ifndef HISTO_H
#define HISTO_H

#include <limits>
#include <iterator>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

/*
template<class T>
class Seqv {
protected:
  std::vector<T> elems_;
public:
  Seqv(const T& elem) { elems_.push_back(elem); }
  Seqv<T>& operator()(const T& elem) { elems_.push_back(elem); return *this; }
  T& operator[](int i) { return elems_[i]; }
  const T& operator[](int i) const { return elems_[i]; }
  size_t size() const { return elems_.size(); }
};
*/

template<int N, class T>
class Seq {
protected:
  T elems_[N];
  int nexti_;
  typedef T (&arrayref)[N];
public:
  Seq() : nexti_(0) {}
  Seq(const T& elem) : nexti_(0) { elems_[nexti_++] = elem; }
  Seq<N, T>& operator()(const T& elem) { elems_[nexti_++] = elem; return *this; }
  T& operator[](int i) { return elems_[i]; }
  const T& operator[](int i) const { return elems_[i]; }
  size_t size() const { return N; }
  operator arrayref() { return elems_; } 
};

template<int N, class T> Seq<N, T> seq(const T& elem) { return Seq<N, T>(elem); }



template<int N, class T>
struct ArrayCompare {
  static bool less(const T l[], const T r[]) {
    return ArrayCompare<1, T>::less(&l[0], &r[0]) || (!ArrayCompare<1, T>::less(&r[0], &l[0]) && ArrayCompare<N-1, T>::less(&l[1], &r[1]));
  }
  static bool equal(const T l[], const T r[]) {
    return ArrayCompare<1, T>::equal(&l[0], &r[0]) && ArrayCompare<N-1, T>::equal(&l[1], &r[1]);
  }
};

template<class T> 
struct ArrayCompare<1, T> { 
  static bool less(const T l[], const T r[]) {
    return l[0] < r[0];
  }
  static bool equal(const T l[], const T r[]) {
    return l[0] == r[0];
  }
};


template<int N, class T, T (*MinGetter)() = std::numeric_limits<T>::min, T (*MaxGetter)() = std::numeric_limits<T>::max>
class BinKey {
  T elems_[N];
public:
  typedef T ElemType;
  //T& operator[](int i) { return elems_[i]; }
  //const T& operator[](int i) const { return elems_[i]; }
  void set(int i, T k) { elems_[i] = k; }
  T at(int i) const { return elems_[i]; }
  size_t size() const { return N; }
  bool operator<(const BinKey<N, T>& other) const { return ArrayCompare<N, T>::less(BinKey<N, T>::elems_, other.elems_); }
  bool operator!=(const BinKey<N, T>& other) const { return !ArrayCompare<N, T>::equal(BinKey<N, T>::elems_, other.elems_); }
  bool underflow() const { return std::count(elems_, elems_+N, MinGetter()); } 
  bool overflow() const { return std::count(elems_, elems_+N, MaxGetter()); }
  static T min() { return MinGetter(); }
  static T max() { return MaxGetter(); }
};


template<int N, class T, class B = BinKey<N, unsigned> >
class Histo {
public:
  enum { Dimensions = N };
  typedef T BinType;
  typedef B InternalBinKey;
  typedef BinKey<N, double> ExportableBinKey;

  typedef std::pair<ExportableBinKey, T> exportable_iterator_element;
  typedef std::pair<InternalBinKey, T> internal_iterator_element;
  typedef typename std::map<InternalBinKey, T>::const_iterator internal_iterator;

  template<class U, class ExportableIteratorElement = typename U::exportable_iterator_element, class InternalIterator = typename U::internal_iterator>
  class ConstIterator : public std::iterator<std::input_iterator_tag, ExportableIteratorElement> { // this is by definition a const iterator
    typedef ConstIterator<U, ExportableIteratorElement, InternalIterator> exportable_iterator;
    typedef exportable_iterator_element* pointer;
    typedef exportable_iterator_element& reference;
    const U& histo_;
    InternalIterator binit_;
    exportable_iterator_element current_;
  public:
    ConstIterator(const InternalIterator& binit, const U& histo) : histo_(histo), binit_(binit) {}
    ConstIterator(const exportable_iterator& other) : histo_(other.histo_), binit_(other.binit_) {}
    exportable_iterator& operator++() { ++binit_; return *this; }
    exportable_iterator operator++(int) { exportable_iterator tmp(*this); operator++(); return tmp; }
    bool operator==(const exportable_iterator& other) { return binit_==other.binit_; }
    bool operator!=(const exportable_iterator& other) { return binit_!=other.binit_; } 
    reference operator*() {
      return (current_ = exportable_iterator_element(histo_.makeExportableBinKey(binit_->first), binit_->second)); 
    }
    pointer operator->() {
      return &(current_ = exportable_iterator_element(histo_.makeExportableBinKey(binit_->first), binit_->second));
    }
  };

  typedef ConstIterator<Histo<N, T, B> > const_iterator;

  template<int M, class H>
    class Indexer {
      template<int O, class I> friend class Indexer;
      Indexer<M-1, H> sub;
      H& h_;
      typename H::InternalBinKey& k_;
      typename H::InternalBinKey& key() { return sub.key(); }
    public:
      Indexer(H& h, typename H::InternalBinKey& k) : sub(h, k) , h_(h), k_(sub.key()) {}
      //Indexer<M-1, H> operator[](double x) {
      //  k_[H::Dimensions-M] = h_.coordToKey(x, H::Dimensions-M);
      //  return Indexer<M-1, H>(h_, k_);
     // }
      Indexer<M-1, H>& operator[](double x) {
        k_.set(H::Dimensions-M, h_.coordToKey(x, H::Dimensions-M));
        return sub;
      }
    };

  template<class H>
    class Indexer<0, H> {
      H& h_;
      template<int O, class I> friend class Indexer;
      typename H::InternalBinKey k_;
      typename H::InternalBinKey& key() { return k_; }
      //const typename H::InternalBinKey& key() const { return k_; }
    public:
      Indexer(H& h, typename H::InternalBinKey& k) : h_(h), k_(k) {}
      operator const typename H::BinType&() { return h_.bins_[k_]; }
      Indexer<0, H>& operator+=(const typename H::BinType& w) { 
        typename H::BinType& bin = h_.bins_[k_];
        bin += w; h_.minMax(bin);
        return *this;
      }
      Indexer<0, H>& operator=(const typename H::BinType& w) {
        typename H::BinType& bin = h_.bins_[k_];
        bin = w; h_.minMax(bin);
        return *this;
      }
      
      bool overflow() const { return k_.overflow(); }
      bool underflow() const { return k_.underflow(); }

      //operator const typename H::BinType&() const { return h_.bins_.at(k_); }
    };

protected:
  friend class ConstIterator<Histo<N, T, B> >;

  std::map<InternalBinKey, T> bins_;
  InternalBinKey currentBin_;

  int nbins_[N];
  double lo_[N], hi_[N];
  T min_, max_;
  
  ExportableBinKey makeExportableBinKey(const InternalBinKey& bk) const {
    ExportableBinKey ebk;
    for (int i=0; i < N; i++) {
      double binw = (hi_[i] - lo_[i])/nbins_[i];
      ebk.set(i, (bk.at(i) != bk.min() ? (bk.at(i) != bk.max() ? (bk.at(i)-1)*binw + lo_[i] + binw/2. : ebk.max()) : ebk.min())); 
      // new bin key returns the center of the bin (to avoid falling in the prev/next bin due to rounding errors)
      // underflow or overflow status is transferred in the bin key
    }
    return ebk;
  }

  InternalBinKey makeInternalBinKey(const ExportableBinKey& ebk) const {
    InternalBinKey ibk;
    for (int i=0; i < N; i++) {
      ibk.set(i, coordToKey(ebk[i], i));
    }
    return ibk;
  }

  int coordToKey(double value, int index) {
    int key = (value - lo_[index])/((hi_[index] - lo_[index])/nbins_[index]) + 1; // bin keys start at 1 (0 might be reserved for underflow)
    //return key;
    return key >= 1 ? (key <= nbins_[index] ? key : InternalBinKey::max()) : InternalBinKey::min();  // overflow bin
  }


  T& findBin(double coords[N]) {
    InternalBinKey key;
    for (int i=0; i < N; i++) {
      key.set(i, coordToKey(coords[i], i));
    }
    return bins_[key];
  }

  void setBinning(int k, int nbins, double lo, double hi) {
    nbins_[k] = nbins;
    lo_[k] = lo;
    hi_[k] = hi;
  }

  void minMax(const T& bin) {
    max_ = bin > max_ ? bin : max_;
    min_ = bin < min_ ? bin : min_;
  }

  Histo() {}

public:
  Histo(const int nbins[N], const double lo[N], const double hi[N]) {
    for (int i=0; i<N; i++) {
      nbins_[i] = nbins[i];
      lo_[i] = lo[i];
      hi_[i] = hi[i];
    }
    min_ = std::numeric_limits<T>::max();
    max_ = T(0);
  }
//  Histo(const Seq<N, int>& nbins, const Seq<N, double>& lo, const Seq<N, double>& hi) 
//    : nbins_(nbins), lo_(lo), hi_(hi), 
//      min_(std::numeric_limits<T>::max()), max_(T(0)) {}

  Histo(std::istream& is) {
    deserialize(is);
  }

  void fill(double coords[N], const T& weight = T(1)) {
    T& bin = findBin(coords);
    bin += weight;
    minMax(bin);
  }

  T get(double coords[N]) { return findBin(coords); }

  InternalBinKey getKey(double coords[N]) {
    InternalBinKey key;
    for (int i=0; i < N; i++) {
      key.set(i, coordToKey(coords[i], i));
    }
    return key;
  }

  template<int M> Histo<M, T, B>* fold(int indices[M]) const {
    Histo<M, T, B>* folded = new Histo<M, T>();
    for (int i=0; i < M; i++) {
      int index = indices[i];
      folded.setBinning(i, nbins_[index], lo_[index], hi_[index]);
    }

    for (typename std::map<InternalBinKey, T>::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
      typename Histo<M, T>::InternalBinKey key;
      for (int i=0; i < M; i++) key.set(i, it->first.at(indices[i]));
      folded.bins_[key] += it->second;
    }

    return folded;
  }

  int getNbins(int k) const { return nbins_[k]; }
  double getLo(int k) const { return lo_[k]; }
  double getHi(int k) const { return hi_[k]; }
  double getWbins(int k) const { return (hi_[k]-lo_[k])/nbins_[k]; }

  T minimumValue() const { return min_; } // unset bins don't count 
  T maximumValue() const { return max_; }

  size_t size() const { return bins_.size(); }

  const_iterator begin() const { return const_iterator(bins_.begin(), *this); }
  const_iterator end() const { return const_iterator(bins_.end(), *this); }

  void clear() {
    bins_.clear();
    min_ = std::numeric_limits<T>::max;
    max_ = T(0);
  }

  void suppressZeros(const T& threshold = T(0)) {
    for (typename std::map<InternalBinKey, T>::iterator it = bins_.begin(); it != bins_.end(); it++) {
      if ((-threshold <= it->second && it->second <= threshold) || it->first.overflow() || it->first.underflow()) bins_.erase(it);
    }
  }

  void serialize(std::ostream& out) {
    for (int i=0; i<N; i++) out << nbins_[i] << (i<N-1 ? " " : "\r\n");
    for (int i=0; i<N; i++) out << lo_[i] << (i<N-1 ? " " : "\r\n");
    for (int i=0; i<N; i++) out << hi_[i] << (i<N-1 ? " " : "\r\n");
    out << min_ << " " << max_ << std::endl;
    out << bins_.size() << std::endl;
    for (typename std::map<InternalBinKey, T>::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
      for (int i=0; i<N; i++) out << it->first.at(i) << " ";
      out << it->second << std::endl;
    }
  }

  bool deserialize(std::istream& in) {
    for (int i=0; i<N; i++) in >> nbins_[i];
    for (int i=0; i<N; i++) in >> lo_[i];
    for (int i=0; i<N; i++) in >> hi_[i];
    in >> min_ >> max_;
    size_t size;
    in >> size;
    size_t i = 0;
    for (; i < size && !in.eof(); i++) {
      InternalBinKey key;
      for (int j=0; j<N; j++) { typename InternalBinKey::ElemType tmp; in >> tmp; key.set(j, tmp); }
      in >> bins_[key];
    }
    if (i < size) return false;

    return true;
  }


  Indexer<N-1, Histo<N,T,B> > operator[](double x) {
    InternalBinKey k; k.set(0, coordToKey(x,0));
    return Indexer<N-1, Histo<N,T,B> >(*this, k);
  }

  /*Indexer<0, Histo<N,T> > operator[](const InternalBinKey& k) {
    currentBin_ = k;
    return Indexer<0, Histo<N,T> >(*this, currentBin_);
  }*/

  Indexer<0, Histo<N,T,B> > operator[](const ExportableBinKey& k) {
    currentBin_ = makeInternalBinKey(k);
    return Indexer<0, Histo<N,T,B> >(*this, currentBin_);
  }
};



template<class H> void toTH1(H& histo, TH1& thisto, int k = 0) {
  thisto.SetBins(histo.getNbins(k), histo.getLo(k), histo.getHi(k));
  for (typename H::const_iterator it = histo.begin(); it != histo.end(); ++it) {
    if (!it->first.overflow() && !it->first.underflow()) thisto.Fill(it->first.at(k), it->second);
  }
}

template<class H> void toTH2(H& histo, TH2& thisto, int k = 0, int l = 1) {
  thisto.SetBins(histo.getNbins(k), histo.getLo(k), histo.getHi(k),
                 histo.getNbins(l), histo.getLo(l), histo.getHi(l));
  for (typename H::const_iterator it = histo.begin(); it != histo.end(); ++it) {
    if (!it->first.overflow() && !it->first.underflow()) thisto.Fill((*it).first.at(k), it->first.at(l), it->second);
  }
}

template<class H> void toTH3(H& histo, TH3& thisto, int k = 0, int l = 1, int m = 2) {
  thisto.SetBins(histo.getNbins(k), histo.getLo(k), histo.getHi(k),
                 histo.getNbins(l), histo.getLo(l), histo.getHi(l),
                 histo.getNbins(m), histo.getLo(m), histo.getHi(m));
  for (typename H::const_iterator it = histo.begin(); it != histo.end(); ++it) {
    if (!it->first.overflow() && !it->first.underflow()) thisto.Fill(it->first.at(k), it->first.at(l), it->first.at(m), it->second);
  }
}


#endif
