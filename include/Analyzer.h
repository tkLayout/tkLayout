// 
// File:   Analyzer.h
// Author: ndemaio
//
// Created on May 19, 2009, 1:58 PM
//

/**
 * @file Analyzer.h
 * @brief This class takes care of analysing the material budget
 */

#ifndef _ANALYZER_H
#define	_ANALYZER_H

#define MY_RANDOM_SEED 0xcaffe

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <hit.hh>
#include <module.hh>
#include <ModuleCap.h>
#include <InactiveElement.h>
#include <InactiveSurfaces.h>
#include <MaterialBudget.h>
#include <TProfile.h>
#include <TGraph.h>
#include "TRandom3.h"

namespace insur {
    /**
     * A warning that may occur during processing
     */
    static const std::string msg_module_warning = "Warning: tracker module with undefined subdetector type found.";
    
    /**
     * Two comparison functions for <i>std::pair<int, int></i> entries.
     */
    bool compareIntPairFirst(std::pair<int, int> p, std::pair<int, int> q);
    bool compareIntPairSecond(std::pair<int, int> p, std::pair<int, int> q);

    /**
     * @class SummaryTable
     * @brief A generic object to build summary tables
     */
    class SummaryTable {
    public:
      SummaryTable() {};
      void setCell(const int row, const int column, std::string content) { summaryTable[std::make_pair(row,column)]=content;};
      std::string getCell(const int row, const int column) { return summaryTable[make_pair(row,column)];};
      std::map<std::pair<int, int>, std::string>& getContent() { return summaryTable; };
    private:
      std::map<std::pair<int, int>, std::string> summaryTable;
    };

    /**
     * @class profileBag
     * @brief A bag of profiles sorted by variable, scope and track's pt
     */
    class profileBag {
    public:
      static const int RhoProfile;
      static const int PhiProfile;
      static const int DProfile;
      static const int CtgthetaProfile;
      static const int Z0Profile;
      static const int PProfile;
      static const int IdealProfile;
      static const int RealProfile;
      static const int TriggerProfile;
      static const int StandardProfile;
      std::map<double, TGraph>& getProfiles(const int& attribute);
      int clearTriggerProfiles();
      int clearStandardProfiles();
      static int buildAttribute(bool ideal, bool isTrigger);
    private:
      std::map<int, std::map<double, TGraph> > graphMap_;
      int clearProfiles(const int& attributeMask);
    };

    /**
     * @class Analyzer
     * @brief This class analyses the properties of a given <i>MaterialBudget</i> instance with respect to eta.
     *
     * It simulates a series of tracks that start at the origin (z = 0), maintain a fixed value of PI / 2 for phi and cover
     * an eta range from 0 to the maximal eta found in the provided geometry. Each volume hit by a track contributes
     * its radiation and interaction lengths to a grand total for that track. Those grand totals, recorded by eta, are
     * stored in a series of histograms that give a complete profile of the expected interaction of the tracker itself with
     * the particles that pass through it.
     */
    class Analyzer {
    public:
        Analyzer();
        virtual ~Analyzer() {}
        TH1D& getHistoModulesBarrelsR() { return ractivebarrel; }
        TH1D& getHistoModulesBarrelsI() { return iactivebarrel; }
        TH1D& getHistoModulesEndcapsR() {return ractiveendcap; }
        TH1D& getHistoModulesEndcapsI() { return iactiveendcap; }
        TH1D& getHistoServicesBarrelsR() { return rserfbarrel; }
        TH1D& getHistoServicesBarrelsI() { return iserfbarrel; }
        TH1D& getHistoServicesEndcapsR() { return rserfendcap; }
        TH1D& getHistoServicesEndcapsI() { return iserfendcap; }
        TH1D& getHistoSupportsBarrelsR() { return rlazybarrel; }
        TH1D& getHistoSupportsBarrelsI() { return ilazybarrel; }
        TH1D& getHistoSupportsEndcapsR() { return rlazyendcap; }
        TH1D& getHistoSupportsEndcapsI() { return ilazyendcap; }
        TH1D& getHistoSupportsBarrelTubesR() { return rlazybtube; }
        TH1D& getHistoSupportsBarrelTubesI() { return ilazybtube; }
        TH1D& getHistoSupportsTubesR() { return rlazytube; }
        TH1D& getHistoSupportsTubesI() { return ilazytube; }
        TH1D& getHistoSupportsUserDefinedR() { return rlazyuserdef; }
        TH1D& getHistoSupportsUserDefinedI() { return ilazyuserdef; }
        TH1D& getHistoBarrelsAllR() { return rbarrelall; }
        TH1D& getHistoBarrelsAllI() { return ibarrelall; }
        TH1D& getHistoEndcapsAllR() { return rendcapall; }
        TH1D& getHistoEndcapsAllI() { return iendcapall; }
        TH1D& getHistoModulesAllR() { return ractiveall; }
        TH1D& getHistoModulesAllI() { return iactiveall; }
        TH1D& getHistoServicesAllR() { return rserfall; }
        TH1D& getHistoServicesAllI() { return iserfall; }
        TH1D& getHistoSupportsAllR() { return rlazyall; }
        TH1D& getHistoSupportsAllI() { return ilazyall; }
        TH1D& getHistoExtraServicesR() { return rextraservices; }
        TH1D& getHistoExtraServicesI() { return iextraservices; }
        TH1D& getHistoExtraSupportsR() { return rextrasupports; }
        TH1D& getHistoExtraSupportsI() { return iextrasupports; }
        TH1D& getHistoGlobalR() { return rglobal; }
        TH1D& getHistoGlobalI() { return iglobal; }
        TH2D& getHistoIsoR() { return isor; }
        TH2D& getHistoIsoI() { return isoi; }
        TH2D& getHistoMapRadiation();
        TH2D& getHistoMapInteraction();
        //std::vector<Track>& getTracks() { return tv; } // useless ?! remove !
        std::map<double, TGraph>& getRhoProfiles(bool ideal, bool isTrigger);
        std::map<double, TGraph>& getPhiProfiles(bool ideal, bool isTrigger);
        std::map<double, TGraph>& getDProfiles(bool ideal, bool isTrigger);
        std::map<double, TGraph>& getCtgThetaProfiles(bool ideal, bool isTrigger);
        std::map<double, TGraph>& getZ0Profiles(bool ideal, bool isTrigger);
        std::map<double, TGraph>& getPProfiles(bool ideal, bool isTrigger);
	profileBag& getProfileBag() { return myProfileBag; }
        virtual void analyzeMaterialBudget(MaterialBudget& mb, std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL);
        virtual void analyzeMaterialBudgetTrigger(MaterialBudget& mb, std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL);
	virtual void analyzeTrigger(MaterialBudget& mb, std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL);
	void analyzeGeometry(Tracker& tracker, int nTracks = 1000); // TODO: why virtual?
	void computeBandwidth(Tracker& tracker);
	void createGeometryLite(Tracker& tracker);
	TH2D& getMapPhiEta() { return mapPhiEta; }
        TCanvas& getEtaProfileCanvas() {return etaProfileCanvas;};
        TH1D& getHitDistribution() {return hitDistribution;};
        TProfile& getTotalEtaProfile() {return totalEtaProfile;};
	TGraph& getPowerDensity() {return powerDensity;};
        std::vector<TProfile>& getTypeEtaProfiles() {return typeEtaProfile;};
        std::vector<TObject> getSavingVector();
	TCanvas* getGeomLite() {if (geomLiteCreated) return geomLite; else return NULL; };
	TCanvas* getGeomLiteXY() {if (geomLiteXYCreated) return geomLiteXY; else return NULL; };
	TCanvas* getGeomLiteYZ() {if (geomLiteYZCreated) return geomLiteYZ; else return NULL; };
	TCanvas* getGeomLiteEC() {if (geomLiteECCreated) return geomLiteEC; else return NULL; };
	TH1D& getChanHitDistribution() { return chanHitDistribution; };
	TH1D& getBandwidthDistribution() { return bandwidthDistribution; };
	TH1D& getBandwidthDistributionSparsified() {return bandwidthDistributionSparsified; } ;
	int getGeometryTracksUsed() {return geometryTracksUsed; };
	int getMaterialTracksUsed() {return materialTracksUsed; };
	// Hadrons
	TGraph& getHadronTotalHitsProfile() {return hadronTotalHitsProfile;};
	TGraph& getHadronAverageHitsProfile() {return hadronAverageHitsProfile;};
	std::vector<double>& getHadronNeededHitsFraction() {return hadronNeededHitsFraction;};
	std::vector<TGraph>& getHadronGoodTracksFraction() { return hadronGoodTracksFraction; };


	static std::vector<double> average(TGraph& myGraph, std::vector<double> cuts);

	static const double ZeroHitsRequired;
	static const double OneHitRequired;

	void computeWeightSummary(MaterialBudget& mb);
	std::map<std::string, SummaryTable>& getBarrelWeightSummary() { return barrelWeights;};
	std::map<std::string, SummaryTable>& getEndcapWeightSummary() { return endcapWeights;};
	std::map<std::string, SummaryTable>& getBarrelWeightComponentSummary() { return barrelComponentWeights;};
	std::map<std::string, SummaryTable>& getEndcapWeightComponentSummary() { return endcapComponentWeights;};
	std::map<std::string, double>& getTypeWeigth() { return typeWeight; };
    protected:
        /**
         * @struct Cell
         * @brief This struct contains the cumulative radiation and interaction lengths over the area described by r and z.
         * @param rlength The cumulative radiation length
         * @param ilength The cumulative interaction length
         * @param rmin The minimal radius value of the cell
         * @param rmax The maximal radius value of the cell
         * @param etamin The minimal eta value of the cell
         * @param etamax The maximal eta value of the cell
         */
        struct Cell { double rlength; double ilength; double rmin; double rmax; double etamin; double etamax; };
        std::vector<std::vector<Cell> > cells;
        TH1D ractivebarrel, ractiveendcap, rserfbarrel, rserfendcap, rlazybarrel, rlazyendcap, rlazybtube, rlazytube, rlazyuserdef;
        TH1D iactivebarrel, iactiveendcap, iserfbarrel, iserfendcap, ilazybarrel, ilazyendcap, ilazybtube, ilazytube, ilazyuserdef;
        TH1D rbarrelall, rendcapall, ractiveall, rserfall, rlazyall;
        TH1D ibarrelall, iendcapall, iactiveall, iserfall, ilazyall;
        TH1D rextraservices, rextrasupports;
        TH1D iextraservices, iextrasupports;
        TH1D rglobal, iglobal;
        TH2D isor, isoi;
	TH2D mapRadiation, mapInteraction;
	TH2I mapRadiationCount, mapInteractionCount;
	TH2D mapRadiationCalib, mapInteractionCalib;
	TH2D mapPhiEta;
	TCanvas etaProfileCanvas;
	TCanvas* geomLite; bool geomLiteCreated;
	TCanvas* geomLiteXY; bool geomLiteXYCreated;
	TCanvas* geomLiteYZ; bool geomLiteYZCreated;
	TCanvas* geomLiteEC; bool geomLiteECCreated;
	TH1D chanHitDistribution;
	TH1D bandwidthDistribution;
	TH1D bandwidthDistributionSparsified;

	std::map<std::string, SummaryTable> barrelWeights;
	std::map<std::string, SummaryTable> endcapWeights;
	std::map<std::string, SummaryTable> barrelComponentWeights;
	std::map<std::string, SummaryTable> endcapComponentWeights;
	std::map<std::string, double> typeWeight;


	TH1D hitDistribution;
        //std::vector<Track> tv;  // remove ?
        //std::vector<Track> tvIdeal; // remove ?
        //std::map<double, TGraph> rhoprofiles, phiprofiles, dprofiles, ctgThetaProfiles, z0Profiles, pProfiles;
        //std::map<double, TGraph> rhoprofilesIdeal, phiprofilesIdeal, dprofilesIdeal, ctgThetaProfilesIdeal, z0ProfilesIdeal, pProfilesIdeal;
        //std::map<double, TGraph> rhoprofilesTrigger, phiprofilesTrigger, dprofilesTrigger, ctgThetaProfilesTrigger, z0ProfilesTrigger, pProfilesTrigger;
        //std::map<double, TGraph> rhoprofilesTriggerIdeal, phiprofilesTriggerIdeal, dprofilesTriggerIdeal, ctgThetaProfilesTriggerIdeal, z0ProfilesTriggerIdeal, pProfilesTriggerIdeal;
	profileBag myProfileBag;
        //std::vector<Track> triggerTv; // remove ?
        //std::vector<Track> triggerTvIdeal; //remove ?
	
	// Hadrons
	TGraph hadronTotalHitsProfile;
	TGraph hadronAverageHitsProfile;
	std::vector<double> hadronNeededHitsFraction;
	std::vector<TGraph> hadronGoodTracksFraction;

	TGraph powerDensity;
        TProfile totalEtaProfile;
        std::vector<TProfile> typeEtaProfile;
  
        std::vector<TObject> savingGeometryV; // Vector of ROOT objects to be saved
        std::vector<TObject> savingMaterialV; // Vector of ROOT objects to be saved

	Material findAllHits(MaterialBudget& mb, MaterialBudget* pm, 
			     double& eta, double& theta, double& phi, Track& track);

	void computeDetailedWeights(std::vector<std::vector<ModuleCap> >& tracker, std::map<std::string, SummaryTable>& weightTables, bool byMaterial);
        virtual Material analyzeModules(std::vector<std::vector<ModuleCap> >& tr,
                                                                                          double eta, double theta, double phi, Track& t, bool isPixel = false);
        virtual Material findHitsModules(std::vector<std::vector<ModuleCap> >& tr,
					 double eta, double theta, double phi, Track& t, bool isPixel = false);
        virtual Material findHitsModuleLayer(std::vector<ModuleCap>& layer, double eta, double theta, double phi, Track& t, bool isPixel = false);

        virtual Material findModuleLayerRI(std::vector<ModuleCap>& layer, double eta, double theta, double phi, Track& t, bool isPixel = false);
        virtual Material analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta, double theta,
								  Track& t, MaterialProperties::Category cat = MaterialProperties::no_cat, bool isPixel = false);
        virtual Material findHitsInactiveSurfaces(std::vector<InactiveElement>& elements, double eta, double theta,
						  Track& t, bool isPixel = false);

	void calculateProfiles(std::vector<double>& p,
			       std::vector<Track>& trackVector,
			       profileBag& aProfileBag,
			       int profileAttributes);
	//std::map<double, TGraph>& thisRhoProfiles,
	//std::map<double, TGraph>& thisPhiProfiles,
	//std::map<double, TGraph>& thisDProfiles,
	//std::map<double, TGraph>& thisCtgThetaProfiles,
	//std::map<double, TGraph>& thisZ0Profiles,
	//std::map<double, TGraph>& thisPProfiles);

        void clearMaterialBudgetHistograms();
        void clearTriggerPerformanceHistograms();
        void clearGeometryHistograms();
        void clearCells();
        void setHistogramBinsBoundaries(int bins, double min, double max);
        void setTriggerHistogramBinsBoundaries(int bins, double min, double max);
        void setCellBoundaries(int bins, double minr, double maxr, double minz, double maxz);
        void fillCell(double r, double eta, double theta, Material mat);
        void fillMapRT(const double& r, const double& theta, const Material& mat);
        void fillMapRZ(const double& r, const double& z, const Material& mat);
        void transformEtaToZ();
    private:
        // A random number generator
	TRandom3 myDice; 
        int findCellIndexR(double r);
        int findCellIndexEta(double eta);
	int createResetCounters(Tracker& tracker, std::map <std::string, int> &modTypes);
	std::pair <XYZVector, double > shootDirection(double minEta, double maxEta);
	ModuleVector trackHit(const XYZVector& origin, const XYZVector& direction, ModuleVector* properModules);
	void resetTypeCounter(std::map<std::string, int> &modTypes);
	double diffclock(clock_t clock1, clock_t clock2);
	Color_t colorPicker(std::string);
	std::map<std::string, Color_t> colorPickMap;
	Color_t lastPickedColor;
	int geometryTracksUsed;
	int materialTracksUsed;
    };
}
#endif	/* _ANALYZER_H */

