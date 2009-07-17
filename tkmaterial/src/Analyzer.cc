/**
 * @file Analyzer.cc
 * @brief This is the implementation of the class that analyses a material budget
 */

#include <Analyzer.h>
namespace insur {
    // public
    /**
     * A comparison function for the first elements in two pairs of integers.
     * @param p The first pair
     * @param q The second pair
     * @return True if <i>p.first</i> is smaller than <i>q.first</i>, false otherwise
     */
    bool compareIntPairFirst(std::pair<int, int> p, std::pair<int, int> q) {
        return (p.first < q.first);
    }
    
    /**
     * A comparison function for the second elements in two pairs of integers.
     * @param p The first pair
     * @param q The second pair
     * @return True if <i>p.second</i> is smaller than <i>q.second</i>, false otherwise
     */
    bool compareIntPairSecond(std::pair<int, int> p, std::pair<int, int> q) {
        return (p.second < q.second);
    }
    
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
        clearCells();
        // prepare etaStep, phiStep, nTracks, nScans
        if (etaSteps > 1) etaStep = etaMax / (double)(etaSteps - 1);
        else etaStep = etaMax;
        nTracks = etaSteps;
        // reset the number of bins and the histogram boundaries (0.0 to etaMax) for all histograms, recalculate the cell boundaries
        setHistogramBinsBoundaries(nTracks, 0.0, etaMax);
        setCellBoundaries(nTracks, 0.0, outer_radius + volume_width, 0.0, etaMax);
        // used fixed phi
        phi = PI / 2.0;
        //      loop over nTracks (eta range [0, etaMax])
        for (int i_eta = 0; i_eta < nTracks; i_eta++) {
            std::pair<double, double> tmp;
            eta = i_eta * etaStep;
            theta = 2 * atan(pow(E, -1 * eta));
            //      active volumes, barrel
            tmp = analyzeModules(mb.getBarrelModuleCaps(), eta, theta, phi);
            ractivebarrel.Fill(eta, tmp.first);
            iactivebarrel.Fill(eta, tmp.second);
            rbarrelall.Fill(eta, tmp.first);
            ibarrelall.Fill(eta, tmp.second);
            ractiveall.Fill(eta, tmp.first);
            iactiveall.Fill(eta, tmp.second);
            rglobal.Fill(eta, tmp.first);
            iglobal.Fill(eta, tmp.second);
            //      active volumes, endcap
            tmp = analyzeModules(mb.getEndcapModuleCaps(), eta, theta, phi);
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
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, MaterialProperties::o_sup);
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
        // integration over eta
        for (unsigned int i = 0; i < cells.size(); i++) {
            for (unsigned int j = 1; j < cells.at(i).size(); j++) {
                cells.at(i).at(j).rlength = cells.at(i).at(j).rlength + cells.at(i).at(j - 1).rlength;
                cells.at(i).at(j).ilength = cells.at(i).at(j).ilength + cells.at(i).at(j - 1).ilength;
            }
        }
        // transformation from (eta, r) to (z, r) coordinates
        transformEtaToZ();
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
    std::pair<double, double> Analyzer::analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double eta, double theta, double phi) {
        std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
        std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
        std::pair<double, double> res, tmp;
        res.first = 0.0;
        res.second = 0.0;
        while (iter != guard) {
            tmp = findModuleLayerRI(*iter, eta, theta, phi);
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
    std::pair<double, double> Analyzer::findModuleLayerRI(std::vector<ModuleCap>& layer, double eta, double theta, double phi) {
        std::vector<ModuleCap>::iterator iter = layer.begin();
        std::vector<ModuleCap>::iterator guard = layer.end();
        std::pair<double, double> res, tmp, offset;
        XYZVector origin, direction;
        Polar3DVector dir;
        double distance, r;
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
                        r = distance * sin(theta);
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
                        fillCell(r, eta, tmp.first, tmp.second);
                        res.first = res.first + tmp.first;
                        res.second = res.second + tmp.second;
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
        std::pair<double, double> res, tmp, offset;
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
                        z = iter->getZOffset() + iter->getZLength() / 2.0;
                        r = z * tan(theta);
                        if (cat == MaterialProperties::u_sup) {
                            s = iter->getZLength() / cos(theta);
                            if (s > (iter->getRWidth() / sin(theta))) s = iter->getRWidth() / sin(theta);
                            if (iter->track()) {
                                res.first = res.first + iter->getRadiationLength() * s / iter->getZLength();
                                res.second = res.second + iter->getInteractionLength() * s / iter->getZLength();
                                fillCell(r, eta, iter->getRadiationLength() * s / iter->getZLength(), iter->getInteractionLength() * s / iter->getZLength());
                            }
                            else {
                                rextrasupports.Fill(eta, iter->getRadiationLength() * s / iter->getZLength());
                                iextrasupports.Fill(eta, iter->getInteractionLength() * s / iter->getZLength());
                            }
                        }
                        else {
                            if (iter->track()) {
                                res.first = res.first + iter->getRadiationLength() / cos(theta);
                                res.second = res.second + iter->getInteractionLength() / cos(theta);
                                fillCell(r, eta, iter->getRadiationLength() / cos(theta), iter->getInteractionLength() / cos(theta));
                            }
                            else {
                                if ((iter->getCategory() == MaterialProperties::b_ser)
                                      || (iter->getCategory() == MaterialProperties::e_ser)) {
                                    rextraservices.Fill(eta, iter->getRadiationLength() / cos(theta));
                                    iextraservices.Fill(eta, iter->getInteractionLength() / cos(theta));
                                }
                                else if ((iter->getCategory() == MaterialProperties::b_sup)
                                             || (iter->getCategory() == MaterialProperties::e_sup)
                                             || (iter->getCategory() == MaterialProperties::o_sup)
                                             || (iter->getCategory() == MaterialProperties::t_sup)) {
                                    rextrasupports.Fill(eta, iter->getRadiationLength() / cos(theta));
                                    iextrasupports.Fill(eta, iter->getInteractionLength() / cos(theta));
                                }
                            }
                        }
                    }
                    else {
                        r = iter->getInnerRadius() + iter->getRWidth() / 2.0;
                        if (cat == MaterialProperties::u_sup) {
                            s = iter->getZLength() / sin(theta);
                            if (s > (iter->getRWidth() / cos(theta))) s = iter->getRWidth() / cos(theta);
                            if (iter->track()) {
                                res.first = res.first + iter->getRadiationLength() * s / iter->getZLength();
                                res.second = res.second + iter->getInteractionLength() * s / iter->getZLength();
                                fillCell(r, eta, iter->getRadiationLength() * s / iter->getZLength(), iter->getInteractionLength() * s / iter->getZLength());
                            }
                            else {
                                rextrasupports.Fill(eta, iter->getRadiationLength() * s / iter->getZLength());
                                iextrasupports.Fill(eta, iter->getInteractionLength() * s / iter->getZLength());
                            }
                        }
                        else {
                            if (iter->track()) {
                                res.first = res.first + iter->getRadiationLength() / sin(theta);
                                res.second = res.second + iter->getInteractionLength() / sin(theta);
                                fillCell(r, eta, iter->getRadiationLength() / sin(theta), iter->getInteractionLength() / sin(theta));
                            }
                            else {
                                if ((iter->getCategory() == MaterialProperties::b_ser)
                                      || (iter->getCategory() == MaterialProperties::e_ser)) {
                                    rextraservices.Fill(eta, iter->getRadiationLength() / sin(theta));
                                    iextraservices.Fill(eta, iter->getInteractionLength() / sin(theta));
                                }
                                else if ((iter->getCategory() == MaterialProperties::b_sup)
                                             || (iter->getCategory() == MaterialProperties::e_sup)
                                             || (iter->getCategory() == MaterialProperties::o_sup)
                                             || (iter->getCategory() == MaterialProperties::t_sup)) {
                                    rextrasupports.Fill(eta, iter->getRadiationLength() / sin(theta));
                                    iextrasupports.Fill(eta, iter->getInteractionLength() / sin(theta));
                                }
                            }
                        }
                    }
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
        // outside tracking volume
        rextraservices.Reset();
        rextraservices.SetNameTitle("rextraservices", "Services Outside Tracking Volume: Radiation Length");
        rextrasupports.Reset();
        rextrasupports.SetNameTitle("rextrasupports", "Supports Outside Tracking Volume: Radiation Length");
        iextraservices.Reset();
        iextraservices.SetNameTitle("iextraservices", "Services Outside Tracking Volume: Interaction Length");
        iextrasupports.Reset();
        iextrasupports.SetNameTitle("iextrasupports", "Supports Outside Tracking Volume: Interaction Length");
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
     * This convenience function sets all values in the internal array <i>cells</i> to zero.
     */
    void Analyzer::clearCells() {
        for (unsigned int i = 0; i < cells.size(); i++) {
            cells.at(i).clear();
        }
        cells.clear();
    }
    
    /**
     * This convenience function sets the number of bins and the lower and upper range for their contents for
     * each of the available histograms.
     * @param bins The number of bins in each 1D histogram
     * @param min The minimal eta value that should be plotted
     * @param max the maximal eta value that should be plotted
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
        rextraservices.SetBins(bins, min, max);
        rextrasupports.SetBins(bins, min, max);
        iextraservices.SetBins(bins, min, max);
        iextrasupports.SetBins(bins, min, max);
        rglobal.SetBins(bins, min, max);
        iglobal.SetBins(bins, min, max);
        isor.SetBins(bins, 0.0, max_length, bins / 2, 0.0, outer_radius + volume_width);
        isoi.SetBins(bins, 0.0, max_length, bins / 2, 0.0, outer_radius + volume_width);
    }
    /**
     * This convenience function sets the number of bins and the lower and upper range for their contents for
     * each of the cells that make up the basis for the 2D histograms.
     * @param bins The number of bins in eta; the number of bins in r will be half that
     * @param minr The minimum radius that will be considered for tracking
     * @param maxr The maximum radius that will be considered for tracking
     * @param mineta The minimum value of eta that will be considered for tracking
     * @param maxeta The maximum value of eta that will be considered for tracking
     */
    void Analyzer::setCellBoundaries(int bins, double minr, double maxr, double mineta, double maxeta) {
        double rstep, etastep;
        rstep = 2 * (maxr - minr) / bins;
        etastep = (maxeta - mineta) / bins;
        Cell c;
        c.rlength = 0.0;
        c.ilength = 0.0;
        cells.resize(bins);
        for (unsigned int i = 0; i < cells.size(); i++) if (bins != 1) cells.at(i).resize(bins / 2, c);
        for (unsigned int i = 0; i < cells.size(); i++) {
            for (unsigned int j = 0; j < cells.at(i).size(); j++) {
                cells.at(i).at(j).rmin = minr + j * rstep;
                cells.at(i).at(j).rmax = minr + (j+1) * rstep;
                cells.at(i).at(j).etamin = mineta + i * etastep;;
                cells.at(i).at(j).etamax = mineta + (i+1) * etastep;
            }
        }
    }
    
    /**
     * This function assigns the local radiation and interaction lengths of a detected hit to their position in the
     * (eta, r) space.
     * @param r The radius at which the hit was detected
     * @param eta The eta value of the current track
     * @param rl The local radiation length
     * @param il The local interaction length
     */
    void Analyzer::fillCell(double r, double eta, double rl, double il) {
        int rindex, etaindex;
        if (cells.size() > 0) {
            for (rindex = 0; (unsigned int) rindex < cells.at(0).size(); rindex++) {
                if ((cells.at(0).at(rindex).rmin <= r) && (cells.at(0).at(rindex).rmax > r)) break;
            }
            if ((unsigned int) rindex < cells.at(0).size()) {
                for (etaindex = 0; (unsigned int) etaindex < cells.size(); etaindex++) {
                    if ((cells.at(etaindex).at(rindex).etamin <= eta) && (cells.at(etaindex).at(rindex).etamax > eta)) break;
                }
                if ((unsigned int) etaindex < cells.size()) {
                    cells.at(etaindex).at(rindex).rlength = cells.at(etaindex).at(rindex).rlength + rl;
                    cells.at(etaindex).at(rindex).ilength = cells.at(etaindex).at(rindex).ilength + il;
                }
            }
        }
    }
    
    /**
     * The integrated radiation and interaction lengths in (eta, r) are converted to (z, r) coordinates and stored in
     * <i>isor</i> and <i>isoi</i> in this function.
     */
    void Analyzer::transformEtaToZ() {
        int size_z, size_r, rindex, etaindex;
        double z, r, eta, z_max, z_min, r_max, r_min, z_c, r_c;
        size_z = isor.GetNbinsX();
        size_r = isor.GetNbinsY();
        z_max = isor.GetXaxis()->GetXmax();
        z_min = isor.GetXaxis()->GetXmin();
        r_max = isor.GetYaxis()->GetXmax();
        r_min = isor.GetYaxis()->GetXmin();
        z_c = (z_max - z_min) / (2 * size_z);
        r_c = (r_max - r_min) / (2 * size_r);
        for (int i = 1; i <= size_z; i++) {
            z = z_min + 2 * (i - 1) * z_c + z_c;
            for (int j = 1; j <= size_r; j++) {
                r = r_min + 2 * (j - 1) * r_c + r_c;
                eta = -log(tan(atan(r / z) / 2.0));
                etaindex = findCellIndexEta(eta);
                rindex = findCellIndexR(r);
                if ((etaindex >= 0) && (rindex >= 0)) isor.Fill(z, r, cells.at(etaindex).at(rindex).rlength);
            }
        }
        size_z = isoi.GetNbinsX();
        size_r = isoi.GetNbinsY();
        z_max = isoi.GetXaxis()->GetXmax();
        z_min = isoi.GetXaxis()->GetXmin();
        r_max = isoi.GetYaxis()->GetXmax();
        r_min = isoi.GetYaxis()->GetXmin();
        z_c = (z_max - z_min) / (2 * size_z);
        r_c = (r_max - r_min) / (2 * size_r);
        for (int i = 0; i < size_z; i++) {
            z = z_min + 2 * i * z_c + z_c;
            for (int j = 0; j < size_r; j++) {
                r = r_min + 2 * j * r_c + r_c;
                eta = -log(tan(atan(r / z) / 2.0));
                etaindex = findCellIndexEta(eta);
                rindex = findCellIndexR(r);
                if ((etaindex >= 0) && (rindex >= 0)) isoi.Fill(z, r, cells.at(etaindex).at(rindex).ilength);
            }
        }
    }
    
    // private
    /**
     * This is a convenience function that converts a radius to an index into the second dimension of <i>cells</i>.
     * @param r The given radius value
     * @return The corresponding index into the vector
     */
    int Analyzer::findCellIndexR(double r) {
        int index = -1;
        if (r >= 0) {
            if (cells.size() > 0) {
                index = 0;
                while (index < (int)cells.at(0).size()) {
                    if (cells.at(0).at(index).rmax < r) index++;
                    else break;
                }
                if (index == (int)cells.at(0).size()) index = -1;
                return index;
            }
        }
        return index;
    }
    
    /**
     * This is a convenience function that converts an eta value to an index into the first dimension of <i>cells</i>.
     * @param eta The given eta value
     * @return the corresponding index into the vector
     */
    int Analyzer::findCellIndexEta(double eta) {
        int index = -1;
        if (eta >= 0) {
            index = 0;
            while (index < (int)cells.size()) {
                if (cells.at(index).size() > 0) {
                    if (cells.at(index).at(0).etamax < eta) index++;
                    else break;
                }
                else {
                    index = -1;
                    break;
                }
            }
            if (index == (int)cells.size()) index = -1;
            return index;
        }
        return index;
    }
}
