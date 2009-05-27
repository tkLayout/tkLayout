/**
 * @file Analyzer.cc
 * @brief
 */

#include <Analyzer.h>
namespace insur {
    // public
    void Analyzer::analyzeMaterialBudget(MaterialBudget& mb, int etaSteps, int phiSteps) {
        int nTracks, nScans;
        double etaStep, phiStep, etaMax, phi, eta, theta;
        if (analysed) {
            clearHistograms();
            analysed = !analysed;
        }
        // find etaMax
        etaMax = findEtaMax(mb);
        // prepare etaStep, phiStep, nTracks, nScans
        etaStep = (2 * etaMax) / (double)etaSteps;
        phiStep = 360.0 / (double)phiSteps;
        nTracks = (2 * etaMax) / etaStep + 1;
        nScans = 360 / phiStep;
        // reset the number of bins and the histogram boundaries (-etaMax to etaMax) for all histograms
        setHistogramBinsBoundaries(nTracks, -1 * etaMax, etaMax);
        // loop over nScans (phi range [0, 2Pi) )
        for (int i_phi = 0; i_phi < nScans; i_phi++) {
            phi = i_phi * phiStep;
            //      loop over nTracks (eta range [-etaMax, etaMax])
            for (int i_eta = 0; i_eta < nTracks; i_eta++) {
                std::pair<double, double> tmp;
                eta = i_eta * etaStep - etaMax;
                theta = 2 * atan(pow(E, -1 * eta));
                std::cout << "Phi = " << phi << ", Theta = " << theta << ", Eta = " << eta << std::endl;
                //      active volumes, barrel
                tmp = analyzeModules(mb.getBarrelModuleCaps(), theta, phi);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                ractivebarrel.Fill(eta, tmp.first);
                iactivebarrel.Fill(eta, tmp.second);
                rbarrelall.Fill(eta, tmp.first);
                ibarrelall.Fill(eta, tmp.second);
                ractiveall.Fill(eta, tmp.first);
                iactiveall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      active volumes, endcap
                tmp = analyzeModules(mb.getEndcapModuleCaps(), theta, phi);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                ractiveendcap.Fill(eta, tmp.first);
                iactiveendcap.Fill(eta, tmp.second);
                rendcapall.Fill(eta, tmp.first);
                iendcapall.Fill(eta, tmp.second);
                ractiveall.Fill(eta, tmp.first);
                iactiveall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      services, barrel
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), eta, theta);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rserfbarrel.Fill(eta, tmp.first);
                iserfbarrel.Fill(eta, tmp.second);
                rbarrelall.Fill(eta, tmp.first);
                ibarrelall.Fill(eta, tmp.second);
                rserfall.Fill(eta, tmp.first);
                iserfall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      services, endcap
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), eta, theta);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rserfendcap.Fill(eta, tmp.first);
                iserfendcap.Fill(eta, tmp.second);
                rendcapall.Fill(eta, tmp.first);
                iendcapall.Fill(eta, tmp.second);
                rserfall.Fill(eta, tmp.first);
                iserfall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      supports, barrel
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::b_sup);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rlazybarrel.Fill(eta, tmp.first);
                ilazybarrel.Fill(eta, tmp.second);
                rbarrelall.Fill(eta, tmp.first);
                ibarrelall.Fill(eta, tmp.second);
                rlazyall.Fill(eta, tmp.first);
                ilazyall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      supports, endcap
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::e_sup);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rlazyendcap.Fill(eta, tmp.first);
                ilazyendcap.Fill(eta, tmp.second);
                rendcapall.Fill(eta, tmp.first);
                iendcapall.Fill(eta, tmp.second);
                rlazyall.Fill(eta, tmp.first);
                ilazyall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      supports, tubes
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::t_sup);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rlazytube.Fill(eta, tmp.first);
                ilazytube.Fill(eta, tmp.second);
                rlazyall.Fill(eta, tmp.first);
                ilazyall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
                //      supports, user defined
                tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::u_sup);
                tmp.first = tmp.first / (double)nScans;
                tmp.second = tmp.second / (double)nScans;
                rlazyuserdef.Fill(eta, tmp.first);
                ilazyuserdef.Fill(eta, tmp.second);
                rlazyall.Fill(eta, tmp.first);
                ilazyall.Fill(eta, tmp.second);
                rglobal.Fill(eta, tmp.first);
                iglobal.Fill(eta, tmp.second);
            }
        }
        analysed = true;
    }
    
    // protected
    std::pair<double, double> Analyzer::analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double theta, double phi) {
        std::cout << "analyzeModules()" << std::endl;
        std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
        std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
        std::pair<double, double> res, tmp;
        res.first = 0.0;
        res.second = 0.0;
        while (iter != guard) {
            tmp = findModuleLayerRI(*iter, theta, phi);
            res.first = res.first + tmp.first;
            res.second = res.second + tmp.second;
            iter++;
        }
        std::cout << "analyzeModules(): done." << std::endl;
        return res;
    }
    
    std::pair<double, double> Analyzer::findModuleLayerRI(std::vector<ModuleCap>& layer, double theta, double phi) {
        std::cout << "findModuleLayerRI()" << std::endl;
        std::vector<ModuleCap>::iterator iter = layer.begin();
        std::vector<ModuleCap>::iterator guard = layer.end();
        std::pair<double, double> res, tmp;
        XYZVector origin, direction;
        Polar3DVector dir;
        double distance;
        res.first = 0.0;
        res.second = 0.0;
        dir.SetCoordinates(1, theta, phi);
        direction = dir;
        while (iter != guard) {
            if ((iter->getModule().getSubdetectorType() == Module::Barrel) ||
                    (iter->getModule().getSubdetectorType() == Module::Endcap)) {
                distance = iter->getModule().trackCross(origin, direction);
                if (distance > 0) {
                    tmp.first = iter->getRadiationLength();
                    tmp.second = iter->getInteractionLength();
                    if (iter->getModule().getSubdetectorType() == Module::Barrel) {
                        tmp.first = tmp.first / sin(theta);
                        tmp.second = tmp.second / sin(theta);
                    }
                    else {
                        tmp.first = tmp.first / cos(theta);
                        tmp.second = tmp.second / cos(theta);
                    }
                    // TODO: incorporate phi in calculation
                }
            }
            else std::cout << msg_module_warning << std::endl;
            iter++;
        }
        std::cout << "findModuleLayerRI(): done." << std::endl;
        return res;
    }
    
    std::pair<double, double> Analyzer::analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
            double theta, MaterialProperties::Category cat) {
        std::vector<InactiveElement>::iterator iter = elements.begin();
        std::vector<InactiveElement>::iterator guard = elements.end();
        std::pair<double, double> res, tmp;
        double s = 0.0;
        res.first = 0.0;
        res.second = 0.0;
        while (iter != guard) {
            if ((cat == MaterialProperties::no_cat) || (cat == iter->getCategory())) {
                tmp = iter->getEtaMinMax();
                if ((tmp.first < eta) && (tmp.second > eta)) {
                    if (iter->isVertical()) {
                        if (cat == MaterialProperties::u_sup) {
                            s = iter->getZLength() / cos(theta);
                            if (s > (iter->getRWidth() / sin(theta))) s = iter->getRWidth() / sin(theta);
                            res.first = res.first + iter->getRadiationLength() * s / iter->getZLength();
                            res.second = res.second + iter->getInteractionLength() * s / iter->getZLength();
                        }
                        else {
                            res.first = res.first + iter->getRadiationLength() / cos(theta);
                            res.second = res.second + iter->getInteractionLength() / cos(theta);
                        }
                    }
                    else {
                        if (cat == MaterialProperties::u_sup) {
                            s = iter->getZLength() / sin(theta);
                            if (s > (iter->getRWidth() / cos(theta))) s = iter->getRWidth() / cos(theta);
                            res.first = res.first + iter->getRadiationLength() * s / iter->getZLength();
                            res.second = res.second + iter->getInteractionLength() * s / iter->getZLength();
                        }
                        else {
                            res.first = res.first + iter->getRadiationLength() / sin(theta);
                            res.second = res.second + iter->getInteractionLength() / sin(theta);
                        }
                    }
                }
            }
            iter++;
        }
        return res;
    }
    
    void Analyzer::clearHistograms() {
        ractivebarrel.Reset();
        ractivebarrel.SetNameTitle("ractivebarrels", "Barrel Modules Radiation Length");
        ractiveendcap.Reset();
        ractiveendcap.SetNameTitle("ractiveendcap", "Endcap Modules Radiation Length");
        rserfbarrel.Reset();
        rserfbarrel.SetNameTitle("rserfbarrel", "Barrel Services Radiation Length");
        rserfendcap.Reset();
        rserfendcap.SetNameTitle("rserfendcap", "Endcap Services Radiation Length");
        rlazybarrel.Reset();
        rlazybarrel.SetNameTitle("rlazybarrel", "Barrel Supports Radiation Length");
        rlazyendcap.Reset();
        rlazyendcap.SetNameTitle("rlazyendcap", "Endcap Supports Radiation Length");
        rlazytube.Reset();
        rlazytube.SetNameTitle("rlazytube", "Support Tubes Radiation Length");
        rlazyuserdef.Reset();
        rlazyuserdef.SetNameTitle("rlazyuserdef", "Userdefined Supports Radiation Length");
        iactivebarrel.Reset();
        iactivebarrel.SetNameTitle("iactivebarrel", "Barrel Modules Interaction Length");
        iactiveendcap.Reset();
        iactiveendcap.SetNameTitle("iactiveendcap", "Endcap Modules Interaction Length");
        iserfbarrel.Reset();
        iserfbarrel.SetNameTitle("iserfbarrel", "Barrel Services Interaction Length");
        iserfendcap.Reset();
        iserfendcap.SetNameTitle("iserfendcap", "Endcap Services Interaction Length");
        ilazybarrel.Reset();
        ilazybarrel.SetNameTitle("ilazybarrel", "Barrel Supports Interaction Length");
        ilazyendcap.Reset();
        ilazyendcap.SetNameTitle("ilazyendcap", "Endcap Supports Interaction Length");
        ilazytube.Reset();
        ilazytube.SetNameTitle("ilazytube", "Support Tubes Interaction Length");
        ilazyuserdef.Reset();
        ilazyuserdef.SetNameTitle("ilazyuserdef", "Userdefined Supports Interaction Length");
        // composite
        rbarrelall.Reset();
        rbarrelall.SetNameTitle("rbarrelall", "Barrel Radiation Length");
        rendcapall.Reset();
        rendcapall.SetNameTitle("rendcapall", "Endcap Radiation Length");
        ractiveall.Reset();
        ractiveall.SetNameTitle("ractiveall", "Modules Radiation Length");
        rserfall.Reset();
        rserfall.SetNameTitle("rserfall", "Services Radiation Length");
        rlazyall.Reset();
        rlazyall.SetNameTitle("rlazyall", "Supports Radiation Length");
        ibarrelall.Reset();
        ibarrelall.SetNameTitle("ibarrelall", "Barrel Interaction Length");
        iendcapall.Reset();
        iendcapall.SetNameTitle("iendcapall", "Endcap Interaction Length");
        iactiveall.Reset();
        iactiveall.SetNameTitle("iactiveall", "Modules Interaction Length");
        iserfall.Reset();
        iserfall.SetNameTitle("iserfall", "Services Interaction Length");
        ilazyall.Reset();
        ilazyall.SetNameTitle("ilazyall", "Supports Interaction Length");
        // global
        rglobal.Reset();
        rglobal.SetNameTitle("rglobal", "Overall Radiation Length");
        iglobal.Reset();
        iglobal.SetNameTitle("iglobal", "Overall Interaction Length");
    }
    
    void Analyzer::setHistogramBinsBoundaries(int bins, double min, double max) {
        ractivebarrel.SetBins(bins, min, max);
        ractiveendcap.SetBins(bins, min, max);
        rserfbarrel.SetBins(bins, min, max);
        rserfendcap.SetBins(bins, min, max);
        rlazybarrel.SetBins(bins, min, max);
        rlazyendcap.SetBins(bins, min, max);
        rlazytube.SetBins(bins, min, max);
        rlazyuserdef.SetBins(bins, min, max);
        iactivebarrel.SetBins(bins, min, max);
        iactiveendcap.SetBins(bins, min, max);
        iserfbarrel.SetBins(bins, min, max);
        iserfendcap.SetBins(bins, min, max);
        ilazybarrel.SetBins(bins, min, max);
        ilazyendcap.SetBins(bins, min, max);
        ilazytube.SetBins(bins, min, max);
        ilazyuserdef.SetBins(bins, min, max);
        rbarrelall.SetBins(bins, min, max);
        rendcapall.SetBins(bins, min, max);
        ractiveall.SetBins(bins, min, max);
        rserfall.SetBins(bins, min, max);
        rlazyall.SetBins(bins, min, max);
        ibarrelall.SetBins(bins, min, max);
        iendcapall.SetBins(bins, min, max);
        iactiveall.SetBins(bins, min, max);
        iserfall.SetBins(bins, min, max);
        ilazyall.SetBins(bins, min, max);
        rglobal.SetBins(bins, min, max);
        iglobal.SetBins(bins, min, max);
    }
    
    double Analyzer::findEtaMax(MaterialBudget& mb) {
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<InactiveElement>::iterator iter, guard;
        std::pair<double, double> etaMinMax;
        double eta;
        etaMinMax = mb.getTracker().getEtaMinMax();
        eta = (fabs(etaMinMax.first) < fabs(etaMinMax.second)) ? fabs(etaMinMax.second) : fabs(etaMinMax.first);
        guard = is.getBarrelServices().end();
        for (iter = is.getBarrelServices().begin(); iter != guard; iter++) {
            etaMinMax = iter->getEtaMinMax();
            if (fabs(etaMinMax.first) > eta) eta = fabs(etaMinMax.first);
            if (fabs(etaMinMax.second) > eta) eta = fabs(etaMinMax.second);
        }
        guard = is.getEndcapServices().end();
        for (iter = is.getEndcapServices().begin(); iter != guard; iter++) {
            etaMinMax = iter->getEtaMinMax();
            if (fabs(etaMinMax.first) > eta) eta = fabs(etaMinMax.first);
            if (fabs(etaMinMax.second) > eta) eta = fabs(etaMinMax.second);
        }
        guard = is.getSupports().end();
        for (iter = is.getSupports().begin(); iter != guard; iter++) {
            etaMinMax = iter->getEtaMinMax();
            if (fabs(etaMinMax.first) > eta) eta = fabs(etaMinMax.first);
            if (fabs(etaMinMax.second) > eta) eta = fabs(etaMinMax.second);
        }
        return eta;
    }
    
}
