/**
 * @file Usher.cc
 * @brief This is the base class implementation for the algorithm that places the inactive elements around an existing <i>Tracker</i> object of active modules.
 */

#include <Usher.h>
namespace insur {
    // public
    /**
     * This is the function that provides a frame for the steps that are necessary to build up the inactive surfaces around
     * a given collection of active modules. The tracker object is analysed first, then the information is used to create the
     * inactive surfaces on the z+ side. In a last step, those surfaces are mirrored around the origin of z to create the
     * complete set of volumes.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @param printstatus A flag that turns the command line summary at the end of processing on or off
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrange(Tracker& tracker, InactiveSurfaces& is, std::string geomfile, bool printstatus) {
        TrackerIntRep tintrep;
        is.setUp(tintrep.analyze(tracker));
        if (is.isUp()) is = arrangeUp(tintrep, is, geomfile);
        else is = arrangeDown(tintrep, is, geomfile);
        is = mirror(tintrep, is);
        if (printstatus) print(tintrep, is, false);
        return is;
    }
    
    // protected
    /**
     * This function provides a frame for the steps that are necessary to to build the inactive surfaces on the z+
     * side for the UP configuration.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrangeUp(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        std::cout << "Arranging UP configuration...";
        is = servicesUp(tracker, is);
        is = supportsAll(tracker, is, geomfile);
        std::cout << "done." << std::endl;
        return is;
    }
    
    /**
     * This function provides a frame for the steps that are necessary to to build the inactive surfaces on the z+
     * side for the DOWN configuration.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrangeDown(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        std::cout << "Arranging DOWN configuration...";
        is = servicesDown(tracker, is);
        is = supportsAll(tracker, is, geomfile);
        std::cout << "done." << std::endl;
        return is;
    }
    
    /**
     * This is the function that performs the final mirror operation once the z+ side has been completed.
     * It takes care of tracking the correct feeder and neighbour volumes in the mirror images as well.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that the mirror operations will be applied to
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::mirror(TrackerIntRep& tracker, InactiveSurfaces& is) {
        std::cout << "Mirroring barrel services...";
        // number of barrel service volumes that need to be reflected
        unsigned int half = is.getBarrelServices().size();
        // barrel service loop
        for (unsigned int i = 0; i < half; i++) {
            InactiveElement& blueprint = is.getBarrelServicePart(i);
            // vertical services => complex neighbourhood
            if (blueprint.isVertical()) {
                InactiveRing ring = mirrorRing(blueprint);
                // feeder is module layer
                if (ring.getFeederType() == InactiveElement::tracker) {
                    // find out if feeder is short layer
                    bool short_layer = false;
                    std::list<std::pair<int, double> >::const_iterator iter = tracker.shortBarrelsList().begin();
                    std::list<std::pair<int, double> >::const_iterator guard = tracker.shortBarrelsList().end();
                    while (iter != guard) {
                        if (blueprint.getFeederIndex() == tracker.realIndexLayer(iter->first)) {
                            short_layer = true;
                            break;
                        }
                        else iter++;
                    }
                    // feeder index of short layer must be adjusted by 1
                    if (short_layer) ring.setFeederIndex(ring.getFeederIndex() - 1);
                }
                // feeder is another service
                else if (ring.getFeederType() == InactiveElement::barrel) {
                    ring.setFeederIndex(is.getBarrelServices().size() - 1);
                }
                // neighbour is barrel service
                if (blueprint.getNeighbourType() == InactiveElement::barrel) {
                    ring.setNeighbourIndex(is.getBarrelServices().size() - i + blueprint.getNeighbourIndex());
                }
                // neighbour is endcap service
                else if (blueprint.getNeighbourType() == InactiveElement::endcap) {
                    ring.setNeighbourIndex(2 * tracker.totalDiscs() - blueprint.getNeighbourIndex() - 1);
                }
                // append the new volume to the list of barrel services
                is.addBarrelServicePart(ring);
            }
            // horizontal services => simple neighbourhood
            else {
                InactiveTube tube = mirrorTube(blueprint);
                // find out if feeder is short layer
                bool short_layer = false;
                std::list<std::pair<int, double> >::const_iterator iter = tracker.shortBarrelsList().begin();
                while (iter != tracker.shortBarrelsList().end()) {
                    if (blueprint.getFeederIndex() == iter->first) {
                        short_layer = true;
                        break;
                    }
                    else iter++;
                }
                // feeder index of short layer must be adjusted by 1
                if (short_layer) tube.setFeederIndex(tube.getFeederIndex() - 1);
                // append the new volume to the list of barrel services
                is.addBarrelServicePart(tube);
            }
        }
        std::cout << "done." << std::endl << "Mirroring endcap services...";
        // number of endcap service volumes that need to be reflected
        half = is.getEndcapServices().size();
        // endcap service loop
        for (unsigned int i = 0; i < half; i++) {
            InactiveElement& blueprint = is.getEndcapServicePart(i);
            // vertical services => should not occur in endcaps for now
            if (blueprint.isVertical()) {
                InactiveRing ring = mirrorRing(blueprint);
                if (blueprint.getFeederIndex() != -1) ring.setFeederIndex(2 * tracker.totalDiscs() - blueprint.getFeederIndex() - 1);
                is.addEndcapServicePart(ring);
            }
            // horizontal services => standard and only case for now
            else {
                InactiveTube tube = mirrorTube(blueprint);
                // calculate feeder index to match that of mirroring disc
                if (blueprint.getFeederIndex() != -1) tube.setFeederIndex(2 * tracker.totalDiscs() - blueprint.getFeederIndex() - 1);
                // calculate neighbour index to match that of a recently-reflected barrel service
                if (blueprint.getNeighbourType() == InactiveElement::barrel) {
                    tube.setNeighbourIndex(is.getBarrelServices().size() / 2 + blueprint.getNeighbourIndex());
                }
                // calculate neighbour index to match that of a recently-reflected endcap service
                else if (blueprint.getNeighbourType() == InactiveElement::endcap) {
                    tube.setNeighbourIndex(is.getEndcapServices().size() - i + blueprint.getNeighbourIndex());
                }
                // append the new volume to the list of endcap services
                is.addEndcapServicePart(tube);
            }
        }
        std::cout << "done." << std::endl << "Mirroring supports...";
        // number of support volumes that may need to be reflected
        half = is.getSupports().size();
        // supports loop
        for (unsigned int i = 0; i < half; i ++) {
            InactiveElement& blueprint = is.getSupportPart(i);
            // only supports that do not cross z=0 need to be reflected
            if (blueprint.getZOffset() > 0) {
                if (blueprint.isVertical()) {
                    InactiveRing ring = mirrorRing(blueprint);
                    // append the new volume to the list of supports
                    is.addSupportPart(ring);
                }
                else {
                    InactiveTube tube = mirrorTube(blueprint);
                    // append the new volume to the list of supports
                    is.addSupportPart(tube);
                }
            }
        }
        std::cout << "done." << std::endl;
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the service volumes in an UP configuration.
     * It first deals with the cut layers that are logically part of the barrel above but treated separately
     * by the analysis function. Then it works through the remaining barrels, then through the endcaps.
     * Feeder and neighbour relations within a single barrel or a single endcap are taken care of right
     * away, the interconnections between the two categories are done as a last step after everything
     * else is in place.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::servicesUp(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw; //zlength, zoffset, rinner, rwidth
        int k = 0, h; // global layer counter, helper variable
        // set the difference between number of barrels and number of endcaps
        if (tracker.nOfEndcaps() == 0) h = tracker.nOfBarrels() - 1;
        else h = tracker.nOfBarrels() - tracker.nOfEndcaps();
        // cut layers
        for (int i = 0; i < h; i ++) {
            zo = tracker.zOffsetBarrel(i) + epsilon;
            zl = tracker.zOffsetBarrel(h) - zo;
            for (int j = 0; j < tracker.nOfLayers(i); j++) {
                ri = tracker.innerRadiusLayer(k);
                rw = tracker.innerRadiusLayer(k + 1) - ri - epsilon;
                is = addBarrelServiceTube(is, zl, zo, ri, volume_width, false);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                is = addBarrelServiceRing(is, volume_width, zo + zl, ri, rw, false);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::barrel);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(is.getBarrelServices().size() - 2);
                if (i != 0) {
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 3);
                }
                k++;
            }
        }
        if ((k == (int)tracker.totalLayers()) && (tracker.nOfEndcaps() == 0) ) is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFinal(true);
        else {
            // regular inner barrels
            zl = volume_width;
            for (int i = h; i < tracker.nOfBarrels() - 1; i ++) {
                zo = tracker.zOffsetBarrel(i) + epsilon;
                for (int j = 0; j < tracker.nOfLayers(i); j++) {
                    ri = tracker.innerRadiusLayer(k);
                    rw = tracker.innerRadiusLayer(k + 1) - ri - volume_width - 2 * epsilon;
                    is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                    if (j == 0) {
                        if (i == h) {
                            if (is.getBarrelServices().size() > 1) {
                                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
                            }
                        }
                        else {
                            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                        }
                    }
                    else {
                        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
                    }
                    k++;
                }
            }
            // outmost barrel
            k = servicesOutmostBarrel(tracker, is, k);
            if ((tracker.nOfEndcaps() < 2) && (tracker.nOfBarrels() > 1)) {
                int sub = tracker.nOfLayers(tracker.nOfBarrels() - 1);
                is.getBarrelServicePart(is.getBarrelServices().size() - sub).setNeighbourType(InactiveElement::barrel);
                is.getBarrelServicePart(is.getBarrelServices().size() - sub).setNeighbourIndex(is.getBarrelServices().size() - sub - 1);
            }
            // endcaps
            k = 0;
            rw = volume_width;
            for (int i = 0; i < tracker.nOfEndcaps(); i++) {
                if (i < tracker.nOfEndcaps() - 1) {
                    int l = findBarrelInnerRadius(tracker.nOfBarrels() - tracker.nOfEndcaps() + 1, tracker);
                    ri = tracker.innerRadiusLayer(l) - 2 * rw - 2 * epsilon;
                }
                else ri = outer_radius - rw - epsilon;
                for (int j = 0; j < tracker.nOfDiscs(i); j++) {
                    if (k == 0) {
                        zl = tracker.zOffsetDisc(k) - tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps()) - volume_width - 2 * epsilon;
                        zo = tracker.zOffsetDisc(k) - zl;
                        if (tracker.totalDiscs() == 1) is = addEndcapServiceTube(is, zl, zo, ri, rw, true);
                        else is = addEndcapServiceTube(is, zl, zo, ri, rw, false);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(k));
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        if (i == tracker.nOfEndcaps() - 1) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                    }
                    else {
                        if (j == 0) {
                            double ro;
                            ro = tracker.innerRadiusLayer(findBarrelInnerRadius(tracker.nOfBarrels() - tracker.nOfEndcaps() + i, tracker));
                            ro = ro - 2 * rw - 2 * epsilon;
                            zo = tracker.zOffsetDisc(k - 1) + epsilon;
                            zl = tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps() + i) - zo + volume_width + epsilon;
                            if ((i == tracker.nOfEndcaps() - 1) && (j == tracker.nOfDiscs(i) - 1)) is = addEndcapServiceTube(is, zl, zo, ro, rw, true);
                            else is = addEndcapServiceTube(is, zl, zo, ro, rw, false);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
                            if (i == tracker.nOfEndcaps() - 1) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                            zo = tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps() + i) + volume_width + 2 * epsilon;
                            zl = tracker.zOffsetDisc(k) - zo;
                        }
                        else {
                            zl = tracker.zOffsetDisc(k) - tracker.zOffsetDisc(k - 1) - epsilon;
                            zo = tracker.zOffsetDisc(k) - zl;
                            
                        }
                        if ((i == tracker.nOfEndcaps() - 1) && (j == tracker.nOfDiscs(i) - 1)) is = addEndcapServiceTube(is, zl, zo, ri, rw, true);
                        else is = addEndcapServiceTube(is, zl, zo, ri, rw, false);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(k));
                        if (j > 0) {
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
                        }
                        else is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        if (i == tracker.nOfEndcaps() - 1) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                    }
                    k++;
                }
            }
            // join barrel and endcap services
            is = joinBarrelEndcapServices(tracker, is, 0, is.getBarrelServices().size(), 0, is.getEndcapServices().size(), tracker.totalDiscs());
        }
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the service volumes in an DOWN configuration.
     * It first deals with the cut layers that are logically part of the barrel above but treated separately
     * by the analysis function. Then it works through the remaining barrel-endcap pair. Feeder and
     * neighbour relations are taken care of right away since the structure is much simpler than that
     * of the UP configuration and allows for far less variation.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::servicesDown(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zrl, ztl, zbo, zeo, ri, rrw, rtw;
        int k = 0;
        // barrels
        zrl = volume_width;
        rtw = volume_width;
        if (findMaxBarrelZ(tracker) < (tracker.zOffsetDisc(tracker.totalDiscs() - 1) + tracker.lengthDisc(tracker.totalDiscs() - 1)))
            zeo = tracker.zOffsetDisc(tracker.totalDiscs() - 1) + tracker.lengthDisc(tracker.totalDiscs() - 1);
        else zeo = findMaxBarrelZ(tracker);
        // inner barrels
        for (int i = 0; i < tracker.nOfBarrels() - 1; i++) {
            zbo = tracker.zOffsetBarrel(i) + epsilon;
            ztl = zeo - zbo;
            for (int j = 0; j < tracker.nOfLayers(i); j++) {
                ri = tracker.innerRadiusLayer(k);
                rrw = tracker.innerRadiusLayer(k + 1) - ri - epsilon;
                is = addBarrelServiceTube(is, ztl, zbo, ri, rtw, true);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                k++;
            }
        }
        // outmost barrel
        k = servicesOutmostBarrel(tracker, is, k);
        // endcap
        ri = outer_radius - rtw - epsilon;
        for (int i = 0; i < tracker.totalDiscs(); i++) {
            if (i == 0) ztl = tracker.zOffsetDisc(i) - tracker.zOffsetBarrel(tracker.nOfBarrels() - 1) - volume_width - 2 * epsilon;
            else ztl = tracker.zOffsetDisc(i) - tracker.zOffsetDisc(i - 1) - epsilon;
            zeo = tracker.zOffsetDisc(i) - ztl;
            if (i == tracker.totalDiscs() - 1) is = addEndcapServiceTube(is, ztl, zeo, ri, rtw, true);
            else is = addEndcapServiceTube(is, ztl, zeo, ri, rtw, false);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(i));
            if (i == 0) {
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 1);
            }
            else {
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
            }
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
        }
        return is;
    }
    
    /**
     * This function provides a frame for the steps that are necessary to to build the supports on the z+ side.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsAll(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        // outer tube
        is = addSupportTube(is, 2 * max_length, 0.0 - max_length, outer_radius, volume_width);
        is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::o_sup);
        is.getSupportPart(is.getSupports().size() - 1).track(false);
        // barrels
        is = supportsRegularBarrels(tracker, is);
        // barrel tubes
        is = supportsBarrelTubes(tracker, is);
        // short barrels
        is = supportsShortBarrels(tracker, is);
        // endcaps
        is = supportsEndcaps(tracker, is);
        // rings from config file
        is = supportsUserDefined(tracker, is, geomfile);
        return is;
    }
    
// private
    /**
     * This convenience function bundles service creation and arrangement for the last barrel within a
     * tracker. It takes care of connecting feeders and neighbours to a certain extent but does not
     * connect anything to other active surface layers or discs outside the last barrel.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param layer The overall layer index of the first layer of the last barrel
     * @return One past the last overall layer index that occurs in the tracker, i.e. the total number of barrel layers
     */
    int Usher::servicesOutmostBarrel(TrackerIntRep& tracker, InactiveSurfaces& is, int layer) {
        int om_b;
        double zl, zo, ri, rw;
        om_b = tracker.nOfBarrels() - 1;
        zl = volume_width;
        // last barrel without last layer
        zo = tracker.zOffsetBarrel(om_b) + epsilon;
        for (int i = 0; i < tracker.nOfLayers(om_b) - 1; i++) {
            ri = tracker.innerRadiusLayer(layer);
            rw = tracker.innerRadiusLayer(layer + 1) - ri - epsilon;
            is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(layer));
            if (i != 0) {
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
            }
            layer++;
        }
        // last layer
        ri = tracker.innerRadiusLayer(layer);
        rw = outer_radius - ri - epsilon;
        if (tracker.nOfEndcaps() == 0) is = addBarrelServiceRing(is, zl, zo, ri, rw, true);
        else is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(layer));
        if (layer > 0) {
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
        }
        layer++;
        return layer;
    }
    
    /**
     * This is the core function that creates and arranges the support parts that hold up the regular barrels.
     * It places a disc at the end of each barrel and outside the service volumes that belong to that same barrel.
     * It also places an inner and an outer tube outside the first and the last layer of each barrel.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsRegularBarrels(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw;
        int k = 0;
        for (int i = 0; i < tracker.nOfBarrels(); i++) {
            if (tracker.nOfLayers(i) > 1) {
                zl = volume_width;
                zo = tracker.zOffsetBarrel(i) + volume_width + 2 * epsilon;
                ri = tracker.innerRadiusLayer(k);
                rw = tracker.innerRadiusLayer(k + tracker.nOfLayers(i) - 1) - ri;
                is = addSupportRing(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::b_sup);
            }
            k = k + tracker.nOfLayers(i);
        }
        return is;
    }
    
    InactiveSurfaces& Usher::supportsBarrelTubes(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double r, z, l;
        std::pair<int, int> aux, stst;
        std::pair<int, double> tmp;
        aux = findBarrelSupportParams(tracker, is.isUp());
        for (int i = 0; i < aux.second; i++) {
            tmp.first = i + 1;
            tmp.second = 0.0;
            stst = findSupportStartStop(tracker, tmp, aux, z, is.isUp());
            if ((stst.second > 0) && (stst.first < stst.second)) {
                r = tracker.innerRadiusLayer(stst.first) - volume_width - epsilon;
                if (aux.first == 0) z = tracker.zOffsetBarrel(i);
                else z = tracker.zOffsetBarrel(aux.first + i - 1);
                l = 2.0 * z;
                z = -z;
                is = addSupportTube(is, l, z, r, volume_width);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
                r = tracker.outerRadiusLayer(stst.second) + epsilon;
                is = addSupportTube(is, l, z, r, volume_width);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
                if (i == aux.second - 1) is.getSupportPart(is.getSupports().size() - 1).track(false);
            }
        }
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the extra support parts needed in short barrels.
     * It places a disc at the inner end of each of those, reaching in radius from the nearest long layer below
     * to the nearest long layer above.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsShortBarrels(TrackerIntRep& tracker, InactiveSurfaces& is) {
        bool regular;
        int start, stop;
        double r, w, z;
        std::list<std::pair<int, double> >::iterator next;
        std::list<std::pair<int, double> >::iterator iter = tracker.shortBarrelsList().begin();
        std::list<std::pair<int, double> >::iterator guard = tracker.shortBarrelsList().end();
        while (iter != guard) {
            next = iter;
            next++;
            z = iter->second - volume_width - epsilon;
            start = iter->first - 1;
            r = tracker.outerRadiusLayer(start) + epsilon;
            regular = (next == guard) && (iter->first < (int)tracker.totalLayers() - 1);
            regular = regular || ((next != guard) && (next->first > iter->first + 1));
            regular = regular || (iter->first == (int)tracker.totalLayers() - 1);
            if (!regular) {
                while ((next != guard) && (iter->first + 1 == next->first)) {
                    iter++;
                    next++;
                }
            }
            stop = iter->first + 1;
            if (stop == tracker.totalLayers()) w = outer_radius - r - epsilon;
            else w = tracker.innerRadiusLayer(stop) - r - epsilon;
            is = addSupportRing(is, volume_width, z, r, w);
            is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::b_sup);
            iter++;
        }
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the support parts that hold up the endcaps. It places
     * two tubes per endcap, one inside the innermost and one outside the outmost ring, outside the service volumes.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsEndcaps(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw;
        int k = 0;
        rw = volume_width;
        for (int i = 0; i < tracker.nOfEndcaps(); i++) {
            if (tracker.nOfDiscs(i) > 1) {
                ri = tracker.outerRadiusEndcap(i) + epsilon;
                zo = tracker.zOffsetDisc(k);
                zl = tracker.zOffsetDisc(k + tracker.nOfDiscs(i) - 1) - zo;
                is = addSupportTube(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::e_sup);
                if (i == tracker.nOfEndcaps() - 1) is.getSupportPart(is.getSupports().size() - 1).track(false);
                ri = tracker.innerRadiusEndcap(i) - rw - epsilon;
                is  = addSupportTube(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::e_sup);
            }
            k = k + tracker.nOfDiscs(i);
        }
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the user-defined support parts. It uses
     * a function in the <i>configParser</i> class to extract a collection of positions in z from
     * the geometry configuration file and goes on to place a series of support discs - one between
     * each layer pair - from the inner to the outer support tube across all barrels. Since these
     * discs are mirrored later along with everything else, only z positions equal or greater than 0
     * are considered in the config file.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file for the user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsUserDefined(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        std::pair<int, int> aux, stst;
        std::pair<int, double> tmp;
        double z;
        configParser persian;
        std::list<std::pair<int, double> >* extras = NULL;
        std::list<std::pair<int, double> >::iterator iter;
        extras = persian.parseSupportsFromFile(geomfile);
        if (extras) {
            aux = findBarrelSupportParams(tracker, is.isUp());
            for (iter = extras->begin(); iter != extras->end(); iter++) {
                if (iter->second < max_length) {
                    z = iter->second - volume_width / 2.0;
                    if (iter->first == 0) {
                        for (int i = 0; i < aux.second; i++) {
                            tmp.first = i + 1;
                            tmp.second = iter->second;
                            stst = findSupportStartStop(tracker, tmp, aux, z, is.isUp());
                            is = addBarrelSupportsUserDefined(is, tracker, stst.first, stst.second, z);
                        }
                    }
                    else {
                        if (iter->first <= aux.second) {
                            stst = findSupportStartStop(tracker, *iter, aux, z, is.isUp());
                            is = addBarrelSupportsUserDefined(is, tracker, stst.first, stst.second, z);
                        }
                    }
                }
            }
            delete extras;
        }
        return is;
    }
    
    /**
     * This function performs the vital task of connecting barrel and endcap services in neighbourhood
     * relations. Each barrel's last service volume is connected to its endcap's first, if it exists. By contrast,
     * each endcap's last service volume is connected to the first volume of the next barrel out, if such a
     * barrel exists.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param begin_b The index of the first barrel service that needs to be considered for connection to an endcap service
     * @param end_b One past the index of the last barrel service that needs to be considered for connection to an endcap service
     * @param begin_e The index of the first endcap service that needs to be considered for connection to a barrel service
     * @param end_e One past the index of the last endcap service that needs to be considered for connection to a barrel service
     * @param d_offset An index offset into the collection of endcap services; used when mirroring endcap services
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::joinBarrelEndcapServices(
    TrackerIntRep& tracker, InactiveSurfaces& is, int begin_b, int end_b, int begin_e, int end_e, int d_offset) {
        std::vector<std::pair<int, int> > first_last_barrel, first_last_endcap;
        int first_b = -1, last_b = -1, first_e = -1, last_e = -1;
        if (tracker.nOfEndcaps() > 0) {
            int h = tracker.nOfBarrels() - tracker.nOfEndcaps();
            for (int i = 0; i < tracker.nOfBarrels(); i++) {
                first_b = last_b + 1;
                last_b = last_b + tracker.nOfLayers(i);
                if (i >= h) first_last_barrel.push_back(std::pair<int, int>(first_b, last_b));
            }
            for (int i = 0; i < tracker.nOfEndcaps(); i++) {
                first_e = last_e + 1;
                last_e = last_e + tracker.nOfDiscs(i);
                first_last_endcap.push_back(std::pair<int, int>(first_e, last_e));
            }
            int b_out = begin_b, b_in = begin_b, e_out = begin_e, e_in = begin_e;
            for (unsigned int i = 0; i < first_last_endcap.size() - 1; i++) {
                b_out = b_in;
                while ((b_out < end_b) && (!is.getBarrelServicePart(b_out).isVertical()
                        || (is.getBarrelServicePart(b_out).getFeederType() != InactiveElement::tracker)
                        || (is.getBarrelServicePart(b_out).getFeederIndex() != first_last_barrel.at(i).second)))
                    b_out++;
                b_in = b_out;
                while ((b_in < end_b) && (!is.getBarrelServicePart(b_in).isVertical()
                        || (is.getBarrelServicePart(b_in).getFeederType() != InactiveElement::tracker)
                        || (is.getBarrelServicePart(b_in).getFeederIndex() != first_last_barrel.at(i + 1).first)))
                    b_in++;
                e_in = e_out;
                while ((e_in < end_e - 1) &&
                        (is.getEndcapServicePart(e_in).getFeederIndex() != (first_last_endcap.at(i).first + d_offset))) e_in++;
                e_out = e_in;
                while ((e_out < end_e - 2) &&
                        (is.getEndcapServicePart(e_out).getFeederIndex() != (first_last_endcap.at(i).second + d_offset))) e_out++;
                e_out++;
                if (!((b_out == end_b) || (b_in == end_b) || (e_in == end_e) || (e_out == end_e))) {
                    is.getEndcapServicePart(e_in).setNeighbourType(InactiveElement::barrel);
                    is.getEndcapServicePart(e_in).setNeighbourIndex(b_out);
                    is.getBarrelServicePart(b_in).setNeighbourType(InactiveElement::endcap);
                    is.getBarrelServicePart(b_in).setNeighbourIndex(e_out);
                }
            }
            b_out = b_in;
            while ((b_out < end_b) && (!is.getBarrelServicePart(b_out).isVertical()
                    || (is.getBarrelServicePart(b_out).getFeederType() != InactiveElement::tracker)
                    || (is.getBarrelServicePart(b_out).getFeederIndex() != first_last_barrel.back().second)))
                b_out++;
            e_in = e_out;
            while ((e_in < end_e) &&
                    (is.getEndcapServicePart(e_in).getFeederIndex() != (first_last_endcap.back().first + d_offset))) e_in++;
            if (!((b_out == end_b) || (e_in == end_e))) {
                is.getEndcapServicePart(e_in).setNeighbourType(InactiveElement::barrel);
                is.getEndcapServicePart(e_in).setNeighbourIndex(b_out);
            }
        }
        return is;
    }
    
    /**
     * This convenience function adds a ring of the given specifics to the collection of barrel services.
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param length The length in z
     * @param offset The position of the leftmost end in z
     * @param radius The inner radius
     * @param width The width in r
     * @param final A flag indicating if the volume will be neighbour to anything or not
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::addBarrelServiceRing(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final) {
        InactiveRing ir;
        ir.setZLength(length);
        ir.setZOffset(offset);
        ir.setInnerRadius(radius);
        ir.setRWidth(width);
        ir.setFinal(final);
        is.addBarrelServicePart(ir);
        return is;
    }
    
    /**
     * This convenience function adds a tube of the given specifics to the collection of barrel services.
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param length The length in z
     * @param offset The position of the leftmost end in z
     * @param radius The inner radius
     * @param width The width in r
     * @param final A flag indicating if the volume will be neighbour to anything or not
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::addBarrelServiceTube(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final) {
        InactiveTube it;
        it.setZLength(length);
        it.setZOffset(offset);
        it.setInnerRadius(radius);
        it.setRWidth(width);
        it.setFinal(final);
        is.addBarrelServicePart(it);
        return is;
    }
    
    /**
     * This convenience function adds a tube of the given specifics to the collection of endcap services.
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param length The length in z
     * @param offset The position of the leftmost end in z
     * @param radius The inner radius
     * @param width The width in r
     * @param final A flag indicating if the volume will be neighbour to anything or not
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::addEndcapServiceTube(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final) {
        InactiveTube it;
        it.setZLength(length);
        it.setZOffset(offset);
        it.setInnerRadius(radius);
        it.setRWidth(width);
        it.setFinal(final);
        is.addEndcapServicePart(it);
        return is;
    }
    
    /**
     * This convenience function adds a ring of the given specifics to the collection of supports.
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param length The length in z
     * @param offset The position of the leftmost end in z
     * @param radius The inner radius
     * @param width The width in r
     * @param final A flag indicating if the volume will be neighbour to anything or not
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::addSupportRing(InactiveSurfaces& is, double length, double offset, double radius, double width) {
        InactiveRing ir;
        ir.setZLength(length);
        ir.setZOffset(offset);
        ir.setInnerRadius(radius);
        ir.setRWidth(width);
        is.addSupportPart(ir);
        return is;
    }
    
    /**
     * This convenience function adds a tube of the given specifics to the collection of supports.
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param length The length in z
     * @param offset The position of the leftmost end in z
     * @param radius The inner radius
     * @param width The width in r
     * @param final A flag indicating if the volume will be neighbour to anything or not
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::addSupportTube(InactiveSurfaces& is, double length, double offset, double radius, double width) {
        InactiveTube it;
        it.setZLength(length);
        it.setZOffset(offset);
        it.setInnerRadius(radius);
        it.setRWidth(width);
        is.addSupportPart(it);
        return is;
    }
    
    InactiveSurfaces& Usher::addBarrelSupportsUserDefined(InactiveSurfaces& is, TrackerIntRep& tracker, int start, int stop, double z) {
        double r, w;
        for (int k = start; k < stop; k++) {
            r = tracker.outerRadiusLayer(k) + epsilon;
            w = tracker.innerRadiusLayer(k + 1) - r - epsilon;
            is = addSupportRing(is, volume_width, z, r, w);
            is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::u_sup);
        }
        return is;
    }
    
    /**
     * This convenience function creates a mirror image of a ring from a blueprint volume with respect to
     * the z origin. The feeder and neighbour indices may have to be adjusted further afterwards.
     * @param blueprint A reference to the original ring object
     * @return The mirror image constructed from the blueprint
     */
    InactiveRing Usher::mirrorRing(InactiveElement& blueprint) {
        InactiveRing ir;
        ir.setCategory(blueprint.getCategory());
        ir.setZOffset(0.0 - blueprint.getZOffset() - blueprint.getZLength());
        ir.setZLength(blueprint.getZLength());
        ir.setInnerRadius(blueprint.getInnerRadius());
        ir.setRWidth(blueprint.getRWidth());
        ir.setTotalMass(blueprint.getTotalMass());
        ir.setLocalMass(blueprint.getLocalMass());
        ir.setExitingMass(blueprint.getExitingMass());
        ir.setRadiationLength(blueprint.getRadiationLength());
        ir.setInteractionLength(blueprint.getInteractionLength());
        blueprint.copyMassVectors(ir);
        ir.setVertical(blueprint.isVertical());
        ir.setFeederType(blueprint.getFeederType());
        ir.setFeederIndex(blueprint.getFeederIndex());
        ir.setNeighbourType(blueprint.getNeighbourType());
        ir.setNeighbourIndex(blueprint.getNeighbourIndex());
        ir.setFinal(blueprint.isFinal());
        ir.track(blueprint.track());
        return ir;
    }
    
    /**
     * This convenience function creates a mirror image of a tube from a blueprint volume with respect to
     * the z origin. The feeder and neighbour indices may have to be adjusted further afterwards.
     * @param blueprint A reference to the original tube object
     * @return The mirror image constructed from the blueprint
     */
    InactiveTube Usher::mirrorTube(InactiveElement& blueprint) {
        InactiveTube it;
        it.setCategory(blueprint.getCategory());
        it.setZOffset(0.0 - blueprint.getZOffset() - blueprint.getZLength());
        it.setZLength(blueprint.getZLength());
        it.setInnerRadius(blueprint.getInnerRadius());
        it.setRWidth(blueprint.getRWidth());
        it.setTotalMass(blueprint.getTotalMass());
        it.setLocalMass(blueprint.getLocalMass());
        it.setExitingMass(blueprint.getExitingMass());
        it.setRadiationLength(blueprint.getRadiationLength());
        it.setInteractionLength(blueprint.getInteractionLength());
        blueprint.copyMassVectors(it);
        it.setVertical(blueprint.isVertical());
        it.setFeederType(blueprint.getFeederType());
        it.setFeederIndex(blueprint.getFeederIndex());
        it.setNeighbourType(blueprint.getNeighbourType());
        it.setNeighbourIndex(blueprint.getNeighbourIndex());
        it.setFinal(blueprint.isFinal());
        it.track(blueprint.track());
        return it;
    }
    
    /**
     * This convenience function calculates the inner radius of an entire barrel.
     * @param barrel The internal index of the barrel
     * @param tintrep A reference to the internal data collection that represents the measurements of the tracker object
     * @return The inner radius of the given barrel
     */
    int Usher::findBarrelInnerRadius(int barrel, TrackerIntRep& tintrep) {
        int l = 0;
        for (int i = 0; i < barrel; i ++) {
            l = l + tintrep.nOfLayers(i);
        }
        return l;
    }
    
    /**
     * This convenience function finds the maximum size in z of the tracker object.
     * @param tintrep A reference to the internal data collection that represents the measurements of the tracker object
     * @return The maxima extension in z of the active tracker volumes
     */
    double Usher::findMaxBarrelZ(TrackerIntRep& tintrep) {
        double z_max = -1.0;
        for (int i = 0; i < tintrep.nOfBarrels(); i++) {
            if (z_max < tintrep.zOffsetBarrel(i)) z_max = tintrep.zOffsetBarrel(i);
        }
        return z_max;
    }
    
    std::pair<int, int> Usher::findSupportStartStop(TrackerIntRep& tracker, std::pair<int, double> udef, std::pair<int, int> aux, double z, bool up) {
        std::pair<int, int> startstop(0, 0);
        if (up) {
            if (udef.first == 1) {
                int i = 0;
                startstop.first = 0;
                while (i < tracker.nOfBarrels()) {
                    if (tracker.zOffsetBarrel(i) < z) startstop.first = startstop.first + tracker.nOfLayers(i);
                    i++;
                }
                startstop.second = tracker.nOfLayers(aux.first);
                for (int i = 0; i < aux.first; i++) startstop.second = startstop.second + tracker.nOfLayers(i);
                if ((tracker.nOfEndcaps() > 0) && (tracker.nOfBarrels() == tracker.nOfEndcaps())) startstop.second--;
            }
            else if (udef.first == aux.second) {
                if (tracker.zOffsetBarrel(tracker.nOfBarrels() - 1) < z) {
                    startstop.first = 0;
                    startstop.second = 0;
                }
                else {
                    startstop.second = tracker.totalLayers();
                    startstop.first = startstop.second - tracker.nOfLayers(tracker.nOfBarrels() - 1);
                    startstop.second--;
                }
            }
            else {
                if (tracker.zOffsetBarrel(udef.first - 1) < z) {
                    startstop.first = 0;
                    startstop.second = 0;
                }
                else {
                    startstop.second = 0;
                    for (int i = 0; i < (aux.first + udef.first); i++) startstop.second = startstop.second + tracker.nOfLayers(i);
                    startstop.first = startstop.second - tracker.nOfLayers(aux.first + udef.first - 1);
                    startstop.second--;
                }
            }
        }
        else {
            if (udef.first == 1) {
                int i = 0;
                startstop.first = 0;
                while (i < tracker.nOfBarrels()) {
                    if (tracker.zOffsetBarrel(i) < z) startstop.first = startstop.first + tracker.nOfLayers(i);
                    i++;
                }
                startstop.second = tracker.totalLayers() - tracker.nOfLayers(aux.first) - 1;
            }
            else {
                startstop.first = tracker.totalLayers() - tracker.nOfLayers(aux.first);
                startstop.second = tracker.totalLayers() - 1;
            }
        }
        return startstop;
    }
    
    /**
     * Print the internal status of the active and inactive surfaces given to the usher to be arranged.
     * @param tintrep A reference to the internal data collection that represents the measurements of the tracker object
     * @param is A reference to the collection of inactive surfaces that belong to this tracker
     */
    void Usher::print(TrackerIntRep& tintrep, InactiveSurfaces& is, bool full_summary) {
        std::cout << std::endl << "Current state of the system is:" << std::endl;
        std::cout << "Printing tracker information..." << std::endl;
        tintrep.print();
        std::cout << "Printing inactive surfaces..." << std::endl;
        is.print(full_summary);
        std::cout << "Usher internal status printed." << std::endl;
    }
    
// nested class public
    /**
     * Get the number of barrels in a previously analysed tracker object. If no such object has been analysed
     * yet, the function returns the negative value of -1 to indicate this. The number reflects the logical number
     * of barrels relevant to inactive surface placement; it may be different from the conceptual number of
     * barrels defined in the geometry configuration file.
     * @return The number of barrels in the tracker as detected by the analysis function.
     */
    int Usher::TrackerIntRep::nOfBarrels() {
        if (post_analysis) return n_of_layers.size();
        else return -1;
    }
    
    /**
     * Get the number of endcaps in a previously analysed tracker object. If no such object has been analysed
     * yet, the function returns the negative value of -1 to indicate this. The number reflects the logical number
     * of endcaps relevant to inactive surface placement; it may be different from the conceptual number of
     * endcaps defined in the geometry configuration file. This is not the case at the moment but it may change
     * if necessary.
     * @return The number of endcaps in the tracker as detected by the analysis function.
     */
    int Usher::TrackerIntRep::nOfEndcaps() {
        if (post_analysis) return n_of_discs.size();
        else return -1;
    }
    
    /**
     * Get the number of layers in a given barrel of a previously analysed tracker object. If no such object has been
     * analysed yet or if the barrel index is out of range, the function returns the negative value of -1 to indicate this.
     * The number reflects the logical number of barrels relevant to inactive surface placement; it may be different
     * from the conceptual number of barrels defined in the geometry configuration file.
     * @param barrelindex The internal, logical index of the requested barrel
     * @return The number of layers in the requested logical barrel
     */
    int Usher::TrackerIntRep::nOfLayers(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < n_of_layers.size())) return n_of_layers.at(barrelindex);
        }
        return -1;
    }
    
    /**
     * Get the total number of layers across all barrels in a previously analysed tracker object. If no such object has
     * been analysed yet, the function returns the negative value of -1 to indicate this.
     * @return The total number of layers in the tracker
     */
    int Usher::TrackerIntRep::totalLayers() {
        if (post_analysis) {
            int sum = 0;
            for (int i = 0; i < (int)n_of_layers.size(); i++) sum = sum + n_of_layers.at(i);
            return sum;
        }
        return -1;
    }
    
    /**
     * Get the number of discs in a given endcap of a previously analysed tracker object. If no such object has been
     * analysed yet or if the endcap index is out of range, the function returns the negative value of -1 to indicate this.
     * The number reflects the logical number of endcaps relevant to inactive surface placement; it may be different
     * from the conceptual number of endcaps defined in the geometry configuration file. This is not the case at the
     * moment but may change if necessary.
     * @param barrelindex The internal, logical index of the requested endcap
     * @return The number of discs in the requested logical endcap
     */
    int Usher::TrackerIntRep::nOfDiscs(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < n_of_discs.size())) return n_of_discs.at(endcapindex);
        }
        return -1;
    }
    
    /**
     * Get the total number of discs across all endcaps on the z+ side in a previously analysed tracker object.
     * If no such object has been analysed yet, the function returns the negative value of -1 to indicate this.
     * @return The total number of discs on the z+ side of the tracker
     */
    int Usher::TrackerIntRep::totalDiscs() {
        if (post_analysis) {
            int sum = 0;
            for (int i = 0; i < (int)n_of_discs.size(); i++) sum = sum + n_of_discs.at(i);
            return sum;
        }
        return -1;
    }
    
    /**
     * Get the inner radius of a layer from the tracker object analysis data. If no analysis of a tracker object
     * has been performed yet or if the layer index is out of range, the function returns the negative value of
     * -1 to indicate this.
     * @param layerindex The layer index
     * @return The inner radius of the given layer
     */
    double Usher::TrackerIntRep::innerRadiusLayer(int layerindex) {
        if (post_analysis) {
            if ((layerindex >= 0) && ((unsigned int)layerindex < layers_io_radius.size())) return layers_io_radius.at(layerindex).first;
        }
        return -1;
    }
    
    /**
     * Get the outer radius of a layer from the tracker object analysis data. If no analysis of a tracker object
     * has been performed yet or if the layer index is out of range, the function returns the negative value of
     * -1 to indicate this.
     *@param layerindex The layer index
     * @return The inner radius of the given layer
     */
    double Usher::TrackerIntRep::outerRadiusLayer(int layerindex) {
        if (post_analysis) {
            if ((layerindex >= 0) && ((unsigned int)layerindex < layers_io_radius.size())) return layers_io_radius.at(layerindex).second;
        }
        return -1;
    }
    
    /**
     * Get the inner radius of an endcap from the tracker object analysis data. If no analysis of a tracker object
     * has been performed yet or if the endcap index is out of range, the function returns the negative value of
     * -1 to indicate this.
     * @param endcapindex The logical endcap index
     * @return The inner radius of the given endcap
     */
    double Usher::TrackerIntRep::innerRadiusEndcap(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < endcaps_io_radius.size())) return endcaps_io_radius.at(endcapindex).first;
        }
        return -1;
    }
    
    /**
     * Get the outer radius of an endcap from the tracker object analysis data. If no analysis of a tracker object
     * has been performed yet or if the endcap index is out of range, the function returns the negative value of
     * -1 to indicate this.
     * @param endcapindex The logical endcap index
     * @return The inner radius of the given endcap
     */
    double Usher::TrackerIntRep::outerRadiusEndcap(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < endcaps_io_radius.size())) return endcaps_io_radius.at(endcapindex).second;
        }
        return -1;
    }
    
    /**
     * Get the length of a barrel from the tracker object analysis data. If no analysis of a tracker object has been
     * performed yet or if the barrel index is out of range, the function returns the negative value of -1 to indicate
     * this.
     * @param barrelindex The logical barrel index
     * @return The total length of the given barrel
     */
    double Usher::TrackerIntRep::lengthBarrel(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < barrels_length_offset.size())) return barrels_length_offset.at(barrelindex).first;
        }
        return -1;
    }
    
    /**
     * Get the leftmost z of a barrel from the tracker object analysis data. If no analysis of a tracker object has been
     * performed yet or if the barrel index is out of range, the function returns the negative value of -1 to indicate
     * this.
     * @param barrelindex The logical barrel index
     * @return The leftmost z of the given barrel
     */
    double Usher::TrackerIntRep::zOffsetBarrel(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < barrels_length_offset.size())) return barrels_length_offset.at(barrelindex).second;
        }
        return -1;
    }
    
    /**
     * Get the width in z of a disc from the tracker object analysis data. If no analysis of a tracker object has been
     * performed yet or if the disc index is out of range, the function returns the negative value of -1 to indicate
     * this.
     * @param discindex The logical disc index
     * @return The total width in z of the given disc
     */
    double Usher::TrackerIntRep::lengthDisc(int discindex) {
        if (post_analysis) {
            if ((discindex >= 0) && ((unsigned int)discindex < discs_length_offset.size())) return discs_length_offset.at(discindex).first;
        }
        return -1;
    }
    
    /**
     * Get the leftmost z of a disc from the tracker object analysis data. If no analysis of a tracker object has been
     * performed yet or if the disc index is out of range, the function returns the negative value of -1 to indicate
     * this.
     * @param discindex The logical disc index
     * @return The leftmost z of the given disc
     */
    double Usher::TrackerIntRep::zOffsetDisc(int discindex) {
        if (post_analysis) {
            if ((discindex >= 0) && ((unsigned int)discindex < discs_length_offset.size())) return discs_length_offset.at(discindex).second;
        }
        return -1;
    }
    
    /**
     * Get the real index of a layer as opposed to the one recorded during analysis.
     * @param tintreplayer The layer as recorded internally during analysis
     * @return The position of the requested layer in the layer collection of the tracker object
     */
    int Usher::TrackerIntRep::realIndexLayer(int tintreplayer) {
        if (post_analysis) {
            if ((tintreplayer >= 0) && ((unsigned int)tintreplayer < real_index_layer.size())) return real_index_layer.at(tintreplayer);
        }
        return -1;
    }
    
    /**
     * Get the real index of a disc as opposed to the one recorded during analysis.
     * @param tintrepdisc The disc as recorded internally during analysis
     * @return The position of the requested disc in the layer collection of the tracker object
     */
    int Usher::TrackerIntRep::realIndexDisc(int tintrepdisc) {
        if (post_analysis) {
            if ((tintrepdisc >= 0) && ((unsigned int)tintrepdisc < real_index_disc.size())) return real_index_disc.at(tintrepdisc);
        }
        return -1;
    }
    
    /**
     * Get access to the entire list of short layers as recorded during analysis.
     * @return A reference to the internal list of pairs that contain the index and the leftmost z of the layer
     */
    std::list<std::pair<int, double> >& Usher::TrackerIntRep::shortBarrelsList() { return short_layers; }
    
    /**
     * This function provides a frame for those parts that analyse the geometrical properties of an existing tracker object.
     * The entire tracker is defined in terms of what is needed to build the inactive surfaces around it afterwards.
     * @param tracker A reference to the existing tracker object.
     * @return The polarity of the tracker: true if it is UP, false if it is DOWN
     */
    bool Usher::TrackerIntRep::analyze(Tracker& tracker) {
        bool up;
        n_of_layers = analyzeBarrels(*(tracker.getBarrelLayers()), layers_io_radius, barrels_length_offset, real_index_layer, short_layers);
        n_of_discs = analyzeEndcaps(*(tracker.getEndcapLayers()), endcaps_io_radius, discs_length_offset, real_index_disc);
        post_analysis = true;
        up = analyzePolarity();
        return up;
    }
    
    std::pair<int, int> Usher::findBarrelSupportParams(TrackerIntRep& tracker, bool up) {
        std::pair<int, int> res(0, 0);
        if (up) {
            if (tracker.nOfEndcaps() > 0) {
                if (tracker.nOfBarrels() == tracker.nOfEndcaps()) {
                    res.first = 0;
                    res.second = tracker.nOfBarrels();
                }
                else {
                    res.first = tracker.nOfBarrels() - tracker.nOfEndcaps() + 1;
                    res.second = tracker.nOfBarrels() - res.first + 1;
                }
            }
            else {
                res.first = tracker.nOfBarrels();
                res.second = 1;
            }
        }
        else {
            res.first = tracker.nOfBarrels() - 1;
            res.second = 2;
        }
        return res;
    }
    
    /**
     * Print a summary of the internal tracker representation.
     */
    void Usher::TrackerIntRep::print() {
        if (!post_analysis) std::cout << "Internal tracker representation has not been initialised yet." << std::endl;
        else {
            std::cout << std::endl << "Internal tracker representation consists of the following parameters:" << std::endl;
            std::cout << std::endl << "Number of barrels: " << nOfBarrels() << std::endl;
            for (int i = 0; i < nOfBarrels(); i++) {
                std::cout << "Barrel " << i << ": " << nOfLayers(i);
                std::cout << " layers of length " << lengthBarrel(i) << " reaching to z = " << zOffsetBarrel(i) << ".";
                std::cout << std::endl;
            }
            std::cout << std::endl << "Total number of layers: " << totalLayers() << std::endl;
            std::cout << "Layers at inner radius of:";
            for (int i = 0; i < totalLayers(); i++) std::cout << " " << innerRadiusLayer(i);
            std::cout << "." << std::endl;
            std::cout << std::endl << "Number of endcaps: " << nOfEndcaps() << std::endl;
            for (int i = 0; i < nOfEndcaps(); i++) {
                std::cout << "Endcap " << i << ": " << nOfDiscs(i);
                std::cout << " discs of radius " << outerRadiusEndcap(i) << ".";
                std::cout << std::endl;
            }
            std::cout << std::endl << "Total number of discs: " << totalDiscs() << std::endl;
            std::cout << "Discs at an offset of:";
            for (int i = 0; i < totalDiscs(); i++) std::cout << " " << zOffsetDisc(i);
            std::cout << "." << std::endl << std::endl;
            if (nOfEndcaps() < 1) std::cout << "There are no endcaps in this configuration." << std::endl;
            else std::cout << "Endcaps extend to z = " << (zOffsetDisc(totalDiscs() - 1) + lengthDisc(totalDiscs() - 1)) << "." << std::endl;
            std::cout << "Number of short layers: " << short_layers.size() << "." << std::endl;
            if (short_layers.size() > 0) {
                std::cout << "Short layers found at index";
                for (std::list<std::pair<int, double> >::iterator iter = short_layers.begin(); iter != short_layers.end(); iter++) std::cout << " " << iter->first;
                std::cout << "." << std::endl;
                std::cout << "Short layers start at ";
                for (std::list<std::pair<int, double> >::iterator iter = short_layers.begin(); iter != short_layers.end(); iter++) std::cout << iter->second << " ";
                std::cout << "." << std::endl;
            }
        }
    }
    
// nested class private
    /**
     * This core function analyses the geometrical properties of the barrels and regroups the values into a more easily accessible format.
     * @param barrel_layers A reference to the collection of layers in the existing tracker object; the source data
     * @param radius_list_io A reference to an empty vector of pairs that will contain the inner and outer radius of each layer
     * @param length_offset_list A reference to an empty vector of pairs that will contain the length and the leftmost point of each layer
     * @param real_index A reference to an empty vector that will contain the index of each layer within the collection of the tracker object
     * @param layers_short A reference to an empty list of pairs that will contain the internal index and the leftmost z of the short layers
     * @return A vector listing the number of layers per barrel
     */
    std::vector<int> Usher::TrackerIntRep::analyzeBarrels(std::vector<Layer*>& barrel_layers, std::vector<std::pair<double, double> >& radius_list_io,
            std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::list<std::pair<int, double> >& layers_short) {
        std::vector<int> layer_counters;
        if (!barrel_layers.empty() && !barrel_layers.at(0)->getModuleVector()->empty()) {
            int index = 0;
            std::vector<Layer*>::iterator current = barrel_layers.begin();
            double z_old, z_new, r_in, r_out, len;
            while ((current != barrel_layers.end()) && ((*current)->getMaxZ() < 0)) {
                current++;
                index++;
            }
            if (current != barrel_layers.end()) {
                z_old = (*current)->getMaxZ();
                len = z_old - (*current)->getMinZ();
                r_in = (*current)->getMinRho();
                r_out = (*current)->getMaxRho();
                layer_counters.push_back(1);
                radius_list_io.push_back(std::pair<double, double>(r_in, r_out));
                length_offset_list.push_back(std::pair<double, double>(len, z_old));
                real_index.push_back(index);
                current++;
                index++;
                while(current != barrel_layers.end()) {
                    if ((*current)->getMaxZ() > 0) {
                        r_in = (*current)->getMinRho();
                        r_out = (*current)->getMaxRho();
                        radius_list_io.push_back(std::pair<double, double>(r_in, r_out));
                        real_index.push_back(index);
                        if ((*current)->getMinZ() > 0) layers_short.push_back(std::pair<int, double>(radius_list_io.size() - 1, (*current)->getMinZ()));
                        z_new = (*current)->getMaxZ();
                        if (z_new == z_old) layer_counters.back()++;
                        else {
                            len = z_new - (*current)->getMinZ();
                            layer_counters.push_back(1);
                            length_offset_list.push_back(std::pair<double, double>(len, z_new));
                            z_old = z_new;
                        }
                    }
                    current++;
                    index++;
                }
            }
        }
        return layer_counters;
    }
    
    /**
     * This core function analyses the geometrical properties of the endcaps and regroups the values into a more easily accessible format.
     * @param endcap_layers A reference to the collection of discs in the existing tracker object; the source data
     * @param radius_list_io A reference to an empty vector of pairs that will contain the inner and outer radius of each disc
     * @param length_offset_list A reference to an empty vector of pairs that will contain the width in z and the leftmost point of each disc
     * @param real_index A reference to an empty vector that will contain the index of each disc within the collection of the tracker object
     * @return A vector listing the number of discs per endcap
     */
    std::vector<int> Usher::TrackerIntRep::analyzeEndcaps(std::vector<Layer*>& endcap_layers, std::vector<std::pair<double, double> >& radius_list_io,
            std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index) {
        std::vector<int> layer_counters;
        if (!endcap_layers.empty() && !(endcap_layers.at(0)->getModuleVector()->empty())) {
            int index = 0;
            std::vector<Layer*>::iterator current = endcap_layers.begin();
            double z_length, z_offset, r_min, r_old, r_new;
            while ((current != endcap_layers.end()) && ((*current)->getMinZ() < 0)) {
                current++;
                index++;
            }
            if (current != endcap_layers.end()) {
                z_offset = (*current)->getMinZ();
                z_length = (*current)->getMaxZ() - z_offset;
                r_min = (*current)->getMinRho();
                r_old = (*current)->getMaxRho();
                layer_counters.push_back(1);
                radius_list_io.push_back(std::pair<double, double>(r_min, r_old));
                length_offset_list.push_back(std::pair<double, double>(z_length, z_offset));
                real_index.push_back(index);
                current++;
                index++;
                while (current != endcap_layers.end()) {
                    z_offset = (*current)->getMinZ();
                    z_length = (*current)->getMaxZ() - z_offset;
                    r_new = (*current)->getMaxRho();
                    length_offset_list.push_back(std::pair<double, double>(z_length, z_offset));
                    real_index.push_back(index);
                    if (r_new == r_old) layer_counters.back()++;
                    else {
                        layer_counters.push_back(1);
                        r_min = (*current)->getMinRho();
                        radius_list_io.push_back(std::pair<double, double>(r_min, r_new));
                        r_old = r_new;
                    }
                    current++;
                    index++;
                }
            }
        }
        return layer_counters;
    }
    
    /**
     * This function attempts to decide if a given tracker geometry should be classified as UP or DOWN:
     * as long as each layer reaches equally far or farther in z as the previous one, it assumes UP. If it finds
     * one that is shorter but farther out in r, it decides on DOWN.
     * @return True if the given configuration is considered to be UP, false if it is considered to be DOWN
     */
    bool Usher::TrackerIntRep::analyzePolarity() {
        if (nOfBarrels() > 1) {
            for (int i = 0; i < nOfBarrels() - 1; i++) {
                if (zOffsetBarrel(i) > zOffsetBarrel(i + 1)) return false;
            }
        }
        return true;
    }
}
