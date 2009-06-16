/**
 * @file Analyzer.cc
 * @brief This is the implementation of the class that analyses a material budget
 */

#include <Analyzer.h>
namespace insur {
    // public
    /**
     * The main analysis function provides a frame for the scan in eta, defers summing up the radiation
     * and interaction lengths for each volume category to subfunctions, and sorts those results into the
     * correct histograms.
     * @param mb A reference to the instance of <i>MaterialBudget</i> that is to be analysed
     * @param etaSteps The number of wedges in the fan of tracks covered by the eta scan
     */
    void Analyzer::analyzeMaterialBudget(MaterialBudget& mb, int etaSteps) {
        int nTracks;
        double etaStep, eta, theta, phi;
        clearHistograms();
        // prepare etaStep, phiStep, nTracks, nScans
        etaStep = etaMax / (double)(etaSteps - 1);
        nTracks = etaSteps;
        // reset the number of bins and the histogram boundaries (-etaMax to etaMax) for all histograms
        setHistogramBinsBoundaries(nTracks, 0, etaMax);
        // used fixed phi
        phi = PI / 2.0;
        //      loop over nTracks (eta range [0, etaMax])
        for (int i_eta = 0; i_eta < nTracks; i_eta++) {
            std::pair<double, double> tmp;
            eta = i_eta * etaStep;
            theta = 2 * atan(pow(E, -1 * eta));
            //      active volumes, barrel
            tmp = analyzeModules(mb.getBarrelModuleCaps(), theta, phi);
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
            rlazytube.Fill(eta, tmp.first);
            ilazytube.Fill(eta, tmp.second);
            rlazyall.Fill(eta, tmp.first);
            ilazyall.Fill(eta, tmp.second);
            rglobal.Fill(eta, tmp.first);
            iglobal.Fill(eta, tmp.second);
            //      supports, user defined
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::u_sup);
            rlazyuserdef.Fill(eta, tmp.first);
            ilazyuserdef.Fill(eta, tmp.second);
            rlazyall.Fill(eta, tmp.first);
            ilazyall.Fill(eta, tmp.second);
            rglobal.Fill(eta, tmp.first);
            iglobal.Fill(eta, tmp.second);
        }
        for (int r = 2; r <= isor.GetNbinsY(); r++) {
            for (int z = 1; z <= isor.GetNbinsX(); z++) {
                isor.SetBinContent(z, r, isor.GetBinContent(z, r) + isor.GetBinContent(z, r - 1));
                isoi.SetBinContent(z, r, isor.GetBinContent(z, r) + isor.GetBinContent(z, r - 1));
            }
        }
    }
    
    // protected
    /**
     * The layer-level analysis function for modules forms the frame for sending a single track through the active modules.
     * It loops through all layers and adds up the results returned from the analysis of each of them.
     * @param tr A reference to the <i>ModuleCap</i> vector of vectors that sits on top of the tracker modules
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @return The summed up radiation and interaction lengths for the given track, bundled into a <i>std::pair</i>
     */
    std::pair<double, double> Analyzer::analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double theta, double phi) {
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
        return res;
    }
    
    /**
     * The module-level analysis function loops through all modules of a given layer, checking for collisions with the given track.
     * If one is found, the radiation and interaction lengths are scaled with respect to theta, then summed up into a grand total,
     * which is returned. As phi is fixed at the moment and the tracks hit the modules orthogonally with respect to it, it is so far
     * not used to scale the results further.
     * @param layer A reference to the <i>ModuleCap</i> vector linking the collection of material properties to the current layer
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @return The scaled and summed up radiation and interaction lengths for the given layer and track, bundled into a <i>std::pair</i>
     */
    std::pair<double, double> Analyzer::findModuleLayerRI(std::vector<ModuleCap>& layer, double theta, double phi) {
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
            if ((iter->getModule().getMaxZ() > 0) &&
                    (theta > iter->getModule().getMinTheta()) && (theta < iter->getModule().getMaxTheta())) {
                if ((iter->getModule().getSubdetectorType() == Module::Barrel) ||
                        (iter->getModule().getSubdetectorType() == Module::Endcap)) {
                    distance = iter->getModule().trackCross(origin, direction);
                    if (distance > 0) {
                        double r, z;
                        tmp.first = iter->getRadiationLength();
                        tmp.second = iter->getInteractionLength();
                        if (iter->getModule().getSubdetectorType() == Module::Barrel) {
                            tmp.first = tmp.first / sin(theta);
                            tmp.second = tmp.second / sin(theta);
                            r = iter->getModule().getMeanPoint().Rho();
                            if (theta != PI / 2.0) z = r / tan(theta);
                            else z = 0.0;
                        }
                        else {
                            tmp.first = tmp.first / cos(theta);
                            tmp.second = tmp.second / cos(theta);
                            z = iter->getModule().getMeanPoint().Z();
                            r = z * tan(theta);
                        }
                        res.first = res.first + tmp.first;
                        res.second = res.second + tmp.second;
                        isor.Fill(z, r, res.first);
                        isoi.Fill(z, r, res.second);
                    }
                }
                else std::cout << msg_module_warning << std::endl;
            }
            iter++;
        }
        return res;
    }
    
    /**
     * The analysis function for inactive volumes loops through the given vector of elements, checking for collisions with
     * the given track. If one is found, the radiation and interaction lengths are scaled with respect to theta, then summed
     * up into a grand total, which is returned. As all inactive volumes are symmetric with respect to rotation around the
     * z-axis, the track angle phi is not necessary.
     * @param elements A reference to the collection of inactive surfaces that is to be checked for collisions with the track
     * @param eta The pseudorapidity value of the current track
     * @param theta The track angle in the yz-plane
     * @param cat The category of inactive surfaces that need to be considered within the collection; none if the function is to look at all of them
     * @return The scaled and summed up radiation and interaction lengths for the given collection of elements and track, bundled into a <i>std::pair</i>
     */
    std::pair<double, double> Analyzer::analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
            double theta, MaterialProperties::Category cat) {
        std::vector<InactiveElement>::iterator iter = elements.begin();
        std::vector<InactiveElement>::iterator guard = elements.end();
        std::pair<double, double> res, tmp;
        double s = 0.0;
        res.first = 0.0;
        res.second = 0.0;
        while (iter != guard) {
            if (((iter->getZOffset() + iter->getZLength()) > 0)
                    && ((cat == MaterialProperties::no_cat) || (cat == iter->getCategory()))) {
                tmp = iter->getEtaMinMax();
                if ((tmp.first < eta) && (tmp.second > eta)) {
                    double r, z;
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
                        z = iter->getZOffset() + iter->getZLength() / 2.0;
                        r = z * tan(theta);
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
                        r = iter->getInnerRadius() + iter->getRWidth() / 2.0;
                        if (theta != PI / 2.0) z = r / tan(theta);
                        else z = 0.0;
                    }
                    isor.Fill(z, r, res.first);
                    isoi.Fill(z, r, res.second);
                }
            }
            iter++;
        }
        return res;
    }
    
    /**
     * This convenience function resets and empties all histograms so they are ready for a new round of analysis.
     */
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
        // isolines
        isor.Reset();
        isor.SetNameTitle("isor", "Radiation Length Contours");
        isoi.Reset();
        isoi.SetNameTitle("isoi", "Interaction Length Contours");
    }
    
    /**
     * This convenience function sets the number of bins and the lower and upper range for their contents for
     * each of the available histograms.
     */
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
        isor.SetBins(z_bins, 0.0, max_length / 2.0, r_bins, inner_radius - volume_width, outer_radius);
        isoi.SetBins(z_bins, 0.0, max_length / 2.0, r_bins, inner_radius - volume_width, outer_radius);
    } 
}
