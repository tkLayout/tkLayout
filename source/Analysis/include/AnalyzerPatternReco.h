/*
 * AnalyzerPatternReco.h
 *
 *  Created on: 28. 11. 2016
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERPATTERNRECO_H_
#define INCLUDE_ANALYZERPATTERNRECO_H_

#include <AnalyzerUnit.h>
#include "global_constants.h"
#include <map>
#include <vector>

// Fwd declaration
class Detector;
class IrradiationMap;
class RootWSite;
class TProfile;

/*
 * @class AnalyzerPatternReco
 * Analyze tracker pattern recognition capabilities for various combination of tracker layers/discs, tagged with given name: pixel (pixel only),
 * strip (strip only), fwd (forward only), etc. Two approaches are studied: inside-out and outside-in approach. The quantities used to qualify the
 * pattern reco capabilities are meant for primary tracks only rather than for secondary tracks, as one always starts in the innermost layer (+ IP contraint)
 * or outermost layer in case of out-in approach. The extrapolation part of pattern recognition uses formulae calculated in parabolic approximation and
 * particle fluences (occupancies) scaled to given pile-up scenario. The fluences are read in from an external file and are assumed to be simulated by Fluka
 * for given p-p collision energies & rough detector (material) setup.
 *
 * Unique module name defined as "AnalyzerPatternReco".
 */
class AnalyzerPatternReco : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerPatternReco(const Detector& detector);

  //! Destructor
  virtual ~AnalyzerPatternReco();

  //! Initialize - mostly histograms & other containers
  //! @return True if OK
  bool init(int nTracks);

  //! Analyzer pattern recognition capabilities and fill histogram (if init OK)
  //! @return True if OK
  bool analyze();

  //! Visualize pattern recognition results (if init & analysis OK) -> add html page with collected created histograms
  //! @return True if OK
  bool visualize(RootWSite& webSite);

  //! Get number of used tracks
  int getNSimTracks() const { return m_nTracks;}

 private:

  //! Check that a file can be opened
  bool checkFile(const std::string& fileName, const std::string& filePath);

  int    m_nTracks; //!< Number of simulated tracks to be used in the analysis
  double m_etaMin;  //!< Minimum eta value;
  double m_etaMax;  //!< Maximum eta value

  IrradiationMap* m_chargedMap;   //!< Charged particles fluxes

  std::map<std::string, std::vector<TProfile*>> m_hisPtHitD0ProjInOut;    //!< InOut approach: D0 projection @ ith+3 measurement plane at given eta for set of pt
  std::map<std::string, std::vector<TProfile*>> m_hisPHitD0ProjInOut;     //!< InOut approach: D0 projection @ ith+3 measurement plane at given eta for set of p
  std::map<std::string, std::vector<TProfile*>> m_hisPtHitZ0ProjInOut;    //!< InOut approach: Z0 projection @ ith+3 measurement plane at given eta for set of pt
  std::map<std::string, std::vector<TProfile*>> m_hisPHitZ0ProjInOut;     //!< InOut approach: Z0 projection @ ith+3 measurement plane at given eta for set of p
  std::map<std::string, std::vector<TProfile*>> m_hisPtHitProbContamInOut;//!< InOut approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of pt
  std::map<std::string, std::vector<TProfile*>> m_hisPHitProbContamInOut; //!< InOut approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of p

  std::vector<TProfile*> m_hisPtFakeProbInOut; //! InOut approach: Bkg contamination probability accumulated across eta for set of pT
  std::vector<TProfile*> m_hisPFakeProbInOut;  //! InOut approach: Bkg contamination probability accumulated across eta for set of p
  std::vector<TProfile*> m_hisPtFakeProbOutIn; //! OutIn approach: Bkg contamination probability accumulated across eta for set of pT
  std::vector<TProfile*> m_hisPFakeProbOutIn;  //! OutIn approach: Bkg contamination probability accumulated across eta for set of p

  const int    c_nBins;

}; // Class

#endif /* INCLUDE_ANALYZERPATTERNRECO_H_ */
