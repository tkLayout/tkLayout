/**
 * @file Usher.cc
 * @brief This is the base class implementation for the algorithm that places the inactive elements around an existing <i>Tracker</i> object of active modules.
 */

#include <Usher.hh>
namespace insur {
    // public
    /**
     * This is the function framing the steps that are necessary to build up the inactive surfaces around a given tracker's
     * collection of active modules. The tracker object is analysed first, then the information is used to create the inactive
     * surfaces on the z+ side. In a last step, those surfaces are mirrored around the origin of z to create the complete
     * set of volumes.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @param printstatus A flag that turns the command line summary at the end of processing on or off
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrange(Tracker& tracker, InactiveSurfaces& is, const std::list<Support*>& supports, bool printstatus) {
        TrackerIntRep tintrep;
        is.setUp(tintrep.analyze(tracker));
        if (is.isUp()) is = arrangeUp(tintrep, is, geom_max_radius, supports);
        else is = arrangeDown(tintrep, is, geom_max_radius, supports);
        is = mirror(tintrep, is);
        if (printstatus) print(tintrep, is, false);
        return is;
    }
    
    /**
     * This is the function framing the steps that are necessary to build up the inactive surfaces around a given pixel
     * detector's collection of active modules. The pixel detector is analysed first, then the information is used to create the
     * inactive surfaces on the z+ side. In a last step, those surfaces are mirrored around the origin of z to create the
     * complete set of volumes
     * @param pixels A reference to an existing tracker object, representing a pixel detector, with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param printstatus A flag that turns the command line summary at the end of processing on or off
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrangePixels(Tracker& pixels, InactiveSurfaces& is, bool printstatus) {
        TrackerIntRep pintrep;
        startTaskClock("Arranging Pixel configuration");
        is.setUp(pintrep.analyze(pixels));
        if (!pixels.skipAllServices()) is = servicesUp(pintrep, is, geom_inner_strip_radius, true);
        if (!pixels.skipAllSupports()) is = supportsAll(pintrep, is, geom_inner_strip_radius, std::list<Support*>(), true);
        stopTaskClock();
        is = mirror(pintrep, is);
        if (printstatus) print(pintrep, is, false);
        return is;
    }
    
    // protected
    /**
     * This function provides a frame for the steps that are necessary to to build the inactive surfaces on the z+
     * side for the UP configuration.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrangeUp(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports) {
        startTaskClock("Arranging UP configuration");
        if (!tracker.skipAllServices()) is = servicesUp(tracker, is, r_outer, false);
        if (!tracker.skipAllSupports()) is = supportsAll(tracker, is, r_outer, supports, false);
        stopTaskClock();
        return is;
    }
    
    /**
     * This function provides a frame for the steps that are necessary to to build the inactive surfaces on the z+
     * side for the DOWN configuration.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::arrangeDown(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports) {
        startTaskClock("Arranging DOWN configuration");
        if (!tracker.skipAllServices()) is = servicesDown(tracker, is, r_outer, false);
        if (!tracker.skipAllSupports()) is = supportsAll(tracker, is, r_outer, supports, false);
        stopTaskClock();
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
        startTaskClock("Mirroring barrel services");
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
        stopTaskClock();
        startTaskClock("Mirroring endcap services");
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
        stopTaskClock();
        startTaskClock("Mirroring supports");
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
        stopTaskClock();
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
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param track_all A flag indicating whether all volumes should be considered as being inside the tracking volume; true if yes, false otherwise
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::servicesUp(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, bool track_all) {
        double zl, zo, ri, rw; //zlength, zoffset, rinner, rwidth
        int k = 0, h; // global layer counter, helper variable
        // set the difference between number of barrels and number of endcaps
        if (tracker.nOfEndcaps() == 0) h = tracker.nOfBarrels() - 1;
        else h = tracker.nOfBarrels() - tracker.nOfEndcaps();
        // cut layers
        for (int i = 0; i < h; i ++) {
            if (!tracker.barrelHasServices(i)) continue;
            zo = tracker.zOffsetBarrel(i) + geom_epsilon;
            zl = tracker.zOffsetBarrel(h) - zo;
            for (int j = 0; j < tracker.nOfLayers(i); j++) {
                // layer end to barrel end
                ri = tracker.innerRadiusLayer(k);
                rw = tracker.innerRadiusLayer(k + 1) - ri - geom_epsilon;
                is = addBarrelServiceTube(is, zl, zo, ri, geom_inactive_volume_width, false);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                is = addBarrelServiceRing(is,geom_inactive_volume_width, zo + zl, ri, rw, false);
                // lower layer to upper layer
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
            zl = geom_inactive_volume_width;
            // barrel loop
            for (int i = h; i < tracker.nOfBarrels() - 1; i ++) {
                if (!tracker.barrelHasServices(i)) continue;
                zo = tracker.zOffsetBarrel(i) + geom_epsilon;
                // layer loop
                for (int j = 0; j < tracker.nOfLayers(i); j++) {
                    ri = tracker.innerRadiusLayer(k);
                    rw = tracker.innerRadiusLayer(k + 1) - ri - geom_inactive_volume_width - 2 * geom_epsilon;
                    is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                    is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                    // first layer in  barrel
                    if (j == 0) {
                        // first regular barrel => connects to cut layer below if it exists
                        if (i == h) {
                            if (is.getBarrelServices().size() > 1) {
                                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
                            }
                        }
                        // other regular barrel => first layer connects to endcap service only
                        else {
                            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                        }
                    }
                    // other layer => always connects to barrel layer and barrel service
                    else {
                        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
                    }
                    k++;
                }
            }
            // outmost barrel
            k = servicesOutmostBarrel(tracker, is, r_outer, k);
            if ((tracker.nOfEndcaps() < 2) && (tracker.nOfBarrels() > 1)) {
                int sub = tracker.nOfLayers(tracker.nOfBarrels() - 1);
                is.getBarrelServicePart(is.getBarrelServices().size() - sub).setNeighbourType(InactiveElement::barrel);
                is.getBarrelServicePart(is.getBarrelServices().size() - sub).setNeighbourIndex(is.getBarrelServices().size() - sub - 1);
            }
            // endcaps
            k = 0;
            rw = geom_inactive_volume_width;
            // endcap loop
            for (int i = 0; i < tracker.nOfEndcaps(); i++) {
                if (!tracker.endcapHasServices(i)) continue;
                if (i < tracker.nOfEndcaps() - 1) {
                    int l = findBarrelInnerRadius(tracker.nOfBarrels() - tracker.nOfEndcaps() + 1, tracker);
                    ri = tracker.innerRadiusLayer(l) - 2 * rw - 2 * geom_epsilon;
                }
                else ri = r_outer - rw - geom_epsilon;
                // disc loop
                for (int j = 0; j < tracker.nOfDiscs(i); j++) {
                    // first disc in tracker
                    if (k == 0) {
                        zl = tracker.zOffsetDisc(k) - tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps()) - geom_inactive_volume_width - 2 * geom_epsilon;
                        zo = tracker.zOffsetDisc(k) - zl;
                        is = addEndcapServiceTube(is, zl, zo, ri, rw, false);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(k));
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        if ((!track_all) && (i == tracker.nOfEndcaps() - 1)) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                    }
                    // other disc in tracker
                    else {
                        // first disc in endcap
                        if (j == 0) {
                            double ro;
                            ro = tracker.innerRadiusLayer(findBarrelInnerRadius(tracker.nOfBarrels() - tracker.nOfEndcaps() + i, tracker));
                            ro = ro - 2 * rw - 2 * geom_epsilon;
                            zo = tracker.zOffsetDisc(k - 1) + geom_epsilon;
                            zl = tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps() + i) - zo + geom_inactive_volume_width + geom_epsilon;
                            is = addEndcapServiceTube(is, zl, zo, ro, rw, false);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
                            if ((!track_all) && (i == tracker.nOfEndcaps() - 1)) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                            zo = tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps() + i) + geom_inactive_volume_width + 2 * geom_epsilon;
                            zl = tracker.zOffsetDisc(k) - zo;
                        }
                        // other disc in endcap
                        else {
                            zl = tracker.zOffsetDisc(k) - tracker.zOffsetDisc(k - 1) - geom_epsilon;
                            zo = tracker.zOffsetDisc(k) - zl;
                            
                        }
                        // all discs except first in tracker
                        is = addEndcapServiceTube(is, zl, zo, ri, rw, false);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(k));
                        if (j > 0) {
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
                        }
                        else is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                        if ((!track_all) && (i == tracker.nOfEndcaps() - 1)) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
                    }
                    if ((i == tracker.nOfEndcaps() - 1) && (j == tracker.nOfDiscs(i) - 1)) {
                        is = addEndcapServiceTube(is, geom_max_length - zo - zl - geom_epsilon, zo + zl + geom_epsilon, ri, rw, true);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                        is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
                        if (!track_all) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
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
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param track_all A flag indicating whether all volumes should be considered as being inside the tracking volume; true if yes, false otherwise
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::servicesDown(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, bool track_all) {
      double ztl, zbo, zeo, ri, rtw;
      // double zrl;
      // double rrw;
        int k = 0;
        // barrels
        //zrl = volume_width;
        rtw = geom_inactive_volume_width;
        if (findMaxBarrelZ(tracker) < (tracker.zOffsetDisc(tracker.totalDiscs() - 1) + tracker.lengthDisc(tracker.totalDiscs() - 1)))
            zeo = tracker.zOffsetDisc(tracker.totalDiscs() - 1) + tracker.lengthDisc(tracker.totalDiscs() - 1);
        else zeo = findMaxBarrelZ(tracker);
        // inner barrels
        for (int i = 0; i < tracker.nOfBarrels() - 1; i++) {
            if (!tracker.barrelHasServices(i)) continue;
            zbo = tracker.zOffsetBarrel(i) + geom_epsilon;
            ztl = zeo - zbo;
            // layer loop
            for (int j = 0; j < tracker.nOfLayers(i); j++) {
                ri = tracker.innerRadiusLayer(k);
                //rrw = tracker.innerRadiusLayer(k + 1) - ri - epsilon;
                is = addBarrelServiceTube(is, ztl, zbo, ri, rtw, true);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(k));
                k++;
            }
        }
        // outmost barrel
        k = servicesOutmostBarrel(tracker, is, r_outer, k);
        // endcap (one by definition)
        ri = r_outer - rtw - geom_epsilon;
        // disc loop
        for (int i = 0; i < tracker.totalDiscs(); i++) {
            if (!tracker.endcapHasServices(tracker.endcapFromDisc(i))) continue;
            if (i == 0) ztl = tracker.zOffsetDisc(i) - tracker.zOffsetBarrel(tracker.nOfBarrels() - 1) - geom_inactive_volume_width - 2 * geom_epsilon;
            else ztl = tracker.zOffsetDisc(i) - tracker.zOffsetDisc(i - 1) - geom_epsilon;
            zeo = tracker.zOffsetDisc(i) - ztl;
            is = addEndcapServiceTube(is, ztl, zeo, ri, rtw, false);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederType(InactiveElement::tracker);
            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setFeederIndex(tracker.realIndexDisc(i));
            // first disc
            if (i == 0) {
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 1);
            }
            // other disc
            else {
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
            }
            if (i == tracker.totalDiscs() - 1) {
                is = addEndcapServiceTube(is, geom_max_length - zeo - ztl - geom_epsilon, zeo + ztl + geom_epsilon, ri, rtw, true);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
            }
            if (!track_all) is.getEndcapServicePart(is.getEndcapServices().size() - 1).track(false);
        }
        return is;
    }
    
    /**
     * This function provides a frame for the steps that are necessary to to build the supports on the z+ side.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param geomfile The name of the tracker's geometry configuration file in case there are user-defined supports
     * @param track_all A flag indicating whether all volumes should be considered as being inside the tracking volume; true if yes, false otherwise
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsAll(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports, bool track_all) {
        // outer tube
        is = addSupportTube(is, 2 * geom_max_length, 0.0 - geom_max_length, r_outer, geom_inactive_volume_width);
        is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::o_sup);
        if (!track_all) is.getSupportPart(is.getSupports().size() - 1).track(false);
        // barrels
        is = supportsRegularBarrels(tracker, is);
        // barrel tubes
        is = supportsBarrelTubes(tracker, is, track_all);
        // short barrels
        is = supportsShortBarrels(tracker, is, r_outer);
        // endcaps
        is = supportsEndcaps(tracker, is, track_all);
        // rings from config file
        is = supportsUserDefined(tracker, is, supports);
        return is;
    }
    
// private
    /**
     * This convenience function bundles service creation and arrangement for the last barrel within a
     * tracker. It takes care of connecting feeders and neighbours to a certain extent but does not
     * connect anything to other active surface layers or discs outside the last barrel.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param outer_r The outer radius of the enclosing tracker volume
     * @param layer The overall layer index of the first layer of the last barrel
     * @return One past the last overall layer index that occurs in the tracker, i.e. the total number of barrel layers
     */
    int Usher::servicesOutmostBarrel(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, int layer) {
        int om_b;
        double zl, zo, ri, rw;
        om_b = tracker.nOfBarrels() - 1;
        zl = geom_inactive_volume_width;
        // last barrel without last layer
        if (!tracker.barrelHasServices(om_b)) return layer;
        zo = tracker.zOffsetBarrel(om_b) + geom_epsilon;
        for (int i = 0; i < tracker.nOfLayers(om_b) - 1; i++) {
            ri = tracker.innerRadiusLayer(layer);
            rw = tracker.innerRadiusLayer(layer + 1) - ri - geom_epsilon;
            is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
            is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(layer));
            // layer is not first in barrel => neighbours need to be set now rather than later
            if (i != 0) {
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourType(InactiveElement::barrel);
                is.getBarrelServicePart(is.getBarrelServices().size() - 1).setNeighbourIndex(is.getBarrelServices().size() - 2);
            }
            layer++;
        }
        // last layer
        ri = tracker.innerRadiusLayer(layer);
        rw = r_outer - ri - geom_epsilon;
        if (tracker.nOfEndcaps() == 0) is = addBarrelServiceRing(is, zl, zo, ri, rw, true);
        else is = addBarrelServiceRing(is, zl, zo, ri, rw, false);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setCategory(MaterialProperties::b_ser);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederType(InactiveElement::tracker);
        is.getBarrelServicePart(is.getBarrelServices().size() - 1).setFeederIndex(tracker.realIndexLayer(layer));
        // neighbour relations if last layer is not the only one
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
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsRegularBarrels(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw;
        int k = 0;
        // barrel loop
        for (int i = 0; i < tracker.nOfBarrels(); i++) {
            // no supports modelled if barrel only has one layer (or less)
            if (tracker.nOfLayers(i) > 1) {
                zl = geom_inactive_volume_width;
                zo = findBarrelSupportRingZ(tracker, i); // tracker.zOffsetBarrel(i) + volume_width + 2 * epsilon;
                ri = tracker.innerRadiusLayer(k);
                rw = tracker.innerRadiusLayer(k + tracker.nOfLayers(i) - 1) - ri;
                is = addSupportRing(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::b_sup);
            }
            k = k + tracker.nOfLayers(i);
        }
        return is;
    }
    
    /**
     * This is the core function that creates and arranges the support tubes within and outside each regula barrel.
     * It places a tube below the minimal radius of the first and above the maximal radius of the last layer.
     * @param tracker A reference to the existing tracker object with the active modules in it
     * @param is A reference to the collection of inactive surfaces that needs to be built up
     * @param track_all A flag indicating whether all volumes should be considered as being inside the tracking volume; true if yes, false otherwise
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsBarrelTubes(TrackerIntRep& tracker, InactiveSurfaces& is, bool track_all) {
        double z = 0;
        double r, l;
        std::pair<int, int> aux, stst;
        std::pair<int, double> tmp;
        // find relevant barrel range taking into account cut layers
        aux = findBarrelSupportParams(tracker, is.isUp());
        // barrel loop
        for (int i = 0; i < aux.second; i++) {
            tmp.first = i + 1;
            tmp.second = 0.0;
            // find start and stop layer in barrel
            stst = findSupportStartStop(tracker, tmp, aux, z, is.isUp());
            // only continue if start and stop layer are meaningful
            if ((stst.second > 0) && (stst.first < stst.second)) {
                r = tracker.innerRadiusLayer(stst.first) - geom_inactive_volume_width - geom_epsilon;
                // correct for cut layers => logical barrels in internal representation but part of bigger picture in model
                if (aux.first == 0) z = tracker.zOffsetBarrel(i);
                else z = tracker.zOffsetBarrel(aux.first + i - 1);
                l = 2.0 * z;
                z = -z;
                is = addSupportTube(is, l, z, r, geom_inactive_volume_width);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
                r = tracker.outerRadiusLayer(stst.second) + geom_epsilon;
                is = addSupportTube(is, l, z, r, geom_inactive_volume_width);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
                // outmost barrel support tube is outside the tracking volume
                if ((!track_all) && (i == aux.second - 1)) is.getSupportPart(is.getSupports().size() - 1).track(false);
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
     * @param outer_r The outer radius of the enclosing tracker volume
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsShortBarrels(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer) {
        bool regular;
        int start, stop;
        double r, w, z;
        std::list<std::pair<int, double> >::iterator next;
        std::list<std::pair<int, double> >::iterator iter = tracker.shortBarrelsList().begin();
        std::list<std::pair<int, double> >::iterator guard = tracker.shortBarrelsList().end();
        // layer loop
        while (iter != guard) {
            next = iter;
            next++;
            z = iter->second - geom_inactive_volume_width - geom_epsilon;
            // find inner radius of support disc
            start = iter->first - 1;
            r = tracker.outerRadiusLayer(start) + geom_epsilon;
            // find width of support disc: check if layer above is long or not
            regular = (next == guard) && (iter->first < (int)tracker.totalLayers() - 1);
            regular = regular || ((next != guard) && (next->first > iter->first + 1));
            regular = regular || (iter->first == (int)tracker.totalLayers() - 1);
            // skip short layers above if they exist
            if (!regular) {
                while ((next != guard) && (iter->first + 1 == next->first)) {
                    iter++;
                    next++;
                }
            }
            stop = iter->first + 1;
            // set width of support disc
            if (stop == tracker.totalLayers()) w = r_outer - r - geom_epsilon;
            else w = tracker.innerRadiusLayer(stop) - r - geom_epsilon;
            // create volume
            is = addSupportRing(is, geom_inactive_volume_width, z, r, w);
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
     * @param track_all A flag indicating whether all volumes should be considered as being inside the tracking volume; true if yes, false otherwise
     * @return A reference to the modified collection of inactive surfaces
     */
    InactiveSurfaces& Usher::supportsEndcaps(TrackerIntRep& tracker, InactiveSurfaces& is, bool track_all) {
        double zl, zo, ri, rw;
        int k = 0;
        rw = geom_inactive_volume_width;
        // endcap loop
        for (int i = 0; i < tracker.nOfEndcaps(); i++) {
            // no supports if endcap has only one disc
            if (tracker.nOfDiscs(i) > 1) {
                ri = tracker.outerRadiusEndcap(i) + geom_epsilon;
                zo = tracker.zOffsetDisc(k);
                zl = tracker.zOffsetDisc(k + tracker.nOfDiscs(i) - 1) - zo;
                is = addSupportTube(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::e_sup);
                if ((!track_all) && (i == tracker.nOfEndcaps() - 1)) is.getSupportPart(is.getSupports().size() - 1).track(false);
                ri = tracker.innerRadiusEndcap(i) - rw - geom_epsilon;
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
    InactiveSurfaces& Usher::supportsUserDefined(TrackerIntRep& tracker, InactiveSurfaces& is, const std::list<Support*>& extras) {
        std::pair<int, int> aux, stst;
        std::pair<int, double> tmp;
        double z;
        // extract list of support positions
        // do nothing if parser returns a null pointer
        if (extras.size()) {
            // find relevant barrel range taking into account cut layers
            aux = findBarrelSupportParams(tracker, is.isUp());
            // support positions loop
            for (auto iter = extras.begin(); iter != extras.end(); iter++) {
                // skip support definition if position is outside tracker volume
                if ((*iter)->midZ() < geom_max_length) {
                    z = (*iter)->midZ() - geom_inactive_volume_width / 2.0;
                    // support definition applies across all layers
                    if ((*iter)->index() == 0) {
                        for (int i = 0; i < aux.second; i++) {
                            tmp.first = i + 1;
                            tmp.second = (*iter)->midZ();
                            stst = findSupportStartStop(tracker, tmp, aux, z, is.isUp());
                            is = addBarrelSupportsUserDefined(is, tracker, stst.first, stst.second, z);
                        }
                    }
                    // support is only relevant for barrel iter->first - 1
                    else {
                        // skip support definition if it applies to a barrel outside the relevant range
                        if ((*iter)->index() <= aux.second) {
                            stst = findSupportStartStop(tracker, (*iter)->toPair(), aux, z, is.isUp());
                            is = addBarrelSupportsUserDefined(is, tracker, stst.first, stst.second, z);
                        }
                    }
                }
            }
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
        // do nothing if there are no endcaps that need to be joined to the barrels
        if (tracker.nOfEndcaps() > 0) {
            int h = tracker.nOfBarrels() - tracker.nOfEndcaps();
            // find first and last layer index for each barrel
            for (int i = 0; i < tracker.nOfBarrels(); i++) {
                first_b = last_b + 1;
                last_b = last_b + tracker.nOfLayers(i);
                if (i >= h) first_last_barrel.push_back(std::pair<int, int>(first_b, last_b));
            }
            // find first and last disc index for each endcap
            for (int i = 0; i < tracker.nOfEndcaps(); i++) {
                first_e = last_e + 1;
                last_e = last_e + tracker.nOfDiscs(i);
                first_last_endcap.push_back(std::pair<int, int>(first_e, last_e));
            }
            int b_out = begin_b, b_in = begin_b, e_out = begin_e, e_in = begin_e;
            // connect endcaps to barrels
            for (unsigned int i = 0; i < first_last_endcap.size() - 1; i++) {
                b_out = b_in;
                /* TODO: to comment*/
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
            r = tracker.outerRadiusLayer(k) + geom_epsilon;
            w = tracker.innerRadiusLayer(k + 1) - r - geom_epsilon;
            is = addSupportRing(is, geom_inactive_volume_width, z, r, w);
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
    
    /**
     * This convenience function calculates the number of cut layers that are part of the innermost
     * barrel and the number of barrels as as defined by the geometry config file from a given internal
     * tracker summary.
     * @param tracker A reference to the existing tracker object
     * @up The polarity of the tracker: true if it is UP, false if it is DOWN
     * @return A pair of integers listing the calculated values
     */
    std::pair<int, int> Usher::findBarrelSupportParams(TrackerIntRep& tracker, bool up) {
        std::pair<int, int> res(0, 0);
        // UP configuration - complex case
        if (up) {
            // further investigations needed if there are endcaps
            if (tracker.nOfEndcaps() > 0) {
                // no cut layers, straightforward barrel definition
                if (tracker.nOfBarrels() == tracker.nOfEndcaps()) {
                    res.first = 0;
                    res.second = tracker.nOfBarrels();
                }
                // cut layers, complex barrel definition
                else {
                    res.first = tracker.nOfBarrels() - tracker.nOfEndcaps() + 1;
                    res.second = tracker.nOfBarrels() - res.first + 1;
                }
            }
            // straightforward case: barrel-only
            else {
                res.first = tracker.nOfBarrels();
                res.second = 1;
            }
        }
        // DOWN configuration - simple case
        else {
            res.first = tracker.nOfBarrels() - 1;
            res.second = 2;
        }
        return res;
    }
    
    /**
     * This convenience function finds the global indices of the first and last layer that are relevant
     * when placing supports for a barrel.
     * @param tracker A reference to the existing tracker object
     * @param udef A pair of parameters relating to the support volumes: the relevant barrel + 1 and the support position in z
     * @param aux The pair of values that is returned by a call to findBarrelSupportParams(); see there
     * @param z The extension in z of the associated barrel
     * @param up The polarity of the tracker: true if it is UP, false if it is DOWN
     * @return A pair of integers listing the calculated values
     */
    std::pair<int, int> Usher::findSupportStartStop(TrackerIntRep& tracker, std::pair<int, double> udef, std::pair<int, int> aux, double z, bool up) {
        std::pair<int, int> startstop(0, 0);
        // UP configuration
        if (up) {
            // only one barrel, or first barrel
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
            // last barrel
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
            // intermediate barrel
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
        // DOWN configuration
        else {
            // first barrel
            if (udef.first == 1) {
                int i = 0;
                startstop.first = 0;
                while (i < tracker.nOfBarrels()) {
                    if (tracker.zOffsetBarrel(i) < z) startstop.first = startstop.first + tracker.nOfLayers(i);
                    i++;
                }
                startstop.second = tracker.totalLayers() - tracker.nOfLayers(aux.first) - 1;
            }
            // second (and last) barrel
            else {
                startstop.first = tracker.totalLayers() - tracker.nOfLayers(aux.first);
                startstop.second = tracker.totalLayers() - 1;
            }
        }
        return startstop;
    }


    double Usher::findBarrelSupportRingZ(TrackerIntRep& tracker, int i) const { // returns the Z of the supports at the end of the barrel (so that they clear any service ring)
      double zDiff = i < tracker.nOfBarrels() - 1 ? tracker.zOffsetBarrel(i+1) - tracker.zOffsetBarrel(i) : 0.;
      if (fabs(zDiff) >= geom_z_threshold_service_zigzag) { // checking if services zigzag
        return tracker.zOffsetBarrel(i) + geom_inactive_volume_width + 2*geom_epsilon;
      } else {
        double maxZOffset = tracker.zOffsetBarrel(0);
        for (int i = 1; i < tracker.nOfBarrels(); i++) maxZOffset = MAX(maxZOffset, tracker.zOffsetBarrel(i));
        return maxZOffset + geom_inactive_volume_width + 2*geom_epsilon;
      } 
    }
    
    /**
     * Print the internal status of the active and inactive surfaces given to the usher to be arranged.
     * @param tintrep A reference to the internal data collection that represents the measurements of the tracker object
     * @param is A reference to the collection of inactive surfaces that belong to this tracker
     * @param full_summary A flag turning printing of detailed information about inactive surfaces on or off; mostly for debugging
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

    bool Usher::TrackerIntRep::barrelHasServices(int barrelindex) {
      return barrel_has_services.at(barrelindex);
    }

    bool Usher::TrackerIntRep::endcapHasServices(int endcapindex) {
      return endcap_has_services.at(endcapindex);
    }

    int Usher::TrackerIntRep::endcapFromDisc(int discindex) {
      int sum = 0;
      int i = 0;
      for (; i < n_of_discs.size(); i++) {
        sum += n_of_discs.at(i);
        if (sum > discindex) break;
      } 
      return i;
    }

    bool Usher::TrackerIntRep::skipAllServices() { return skipAllServices_; }
    bool Usher::TrackerIntRep::skipAllSupports() { return skipAllSupports_; }
    
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
        skipAllServices_ = tracker.skipAllServices();
        skipAllSupports_ = tracker.skipAllSupports();
        n_of_layers = analyzeBarrels(tracker, layers_io_radius, barrels_length_offset, real_index_layer, short_layers, barrel_has_services);
        n_of_discs = analyzeEndcaps(tracker, endcaps_io_radius, discs_length_offset, real_index_disc, endcap_has_services);
        post_analysis = true;
        up = tracker.servicesForcedUp() || analyzePolarity();
        return up;
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
    std::vector<int> Usher::TrackerIntRep::analyzeBarrels(Tracker& tracker, std::vector<std::pair<double, double> >& radius_list_io,
            std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::list<std::pair<int, double> >& layers_short, std::vector<bool>& barrel_has_services) {
        // do nothing if there are no barrel layers at all or if the first layer has no modules in it

        class BarrelVisitor : public ConstGeometryVisitor {
          std::vector<std::pair<double, double> >& radius_list_io_;
          std::vector<std::pair<double, double> >& length_offset_list_;
          std::vector<int>& real_index_;
          std::list<std::pair<int, double> >& layers_short_;
          std::vector<int> layer_counters_;
          std::vector<bool>& barrel_has_services_;
          
          double prevZ = -1.;
          int index = 0;
        public:
          BarrelVisitor(std::vector<std::pair<double, double> >& radius_list_io, 
                        std::vector<std::pair<double, double> >& length_offset_list, 
                        std::vector<int>& real_index, 
                        std::list<std::pair<int, double> >& layers_short,
                        std::vector<bool>& barrel_has_services) :
              radius_list_io_(radius_list_io),
              length_offset_list_(length_offset_list),
              real_index_(real_index),
              layers_short_(layers_short),
              barrel_has_services_(barrel_has_services) {}


          void visit(const Barrel& barrel) {
            int counter = 0;

            barrel_has_services_.push_back(!barrel.skipServices());

            for(const Layer& layer : barrel.layers()){
            //for(Barrel::Container::const_iterator layerIter = barrel.layers().begin(); layerIter != barrel.layers().end(); ++ layerIter){
              if (layer.maxZwithHybrids() > 0.0) {
                real_index_.push_back(index);
                radius_list_io_.push_back(std::make_pair(layer.minRwithHybrids(), layer.maxRwithHybrids()));

                counter ++;
              }
              index++;
            }
            layer_counters_.push_back(counter);
            length_offset_list_.push_back(std::make_pair(barrel.maxZwithHybrids() - barrel.minZwithHybrids(), barrel.maxZwithHybrids()));
            if (barrel.minZwithHybrids() > 0) layers_short_.push_back(std::pair<int, double>(radius_list_io_.size() - 1, barrel.minZwithHybrids()));
          }

/*
          void visit(const Layer& l) {
            if (l.maxZwithHybrids() > 0.) { // skip Z- mezzanine layers
              real_index_.push_back(index);
              radius_list_io_.push_back(std::make_pair(l.minRwithHybrids(), l.maxRwithHybrids()));
              if (prevZ == l.maxZwithHybrids()) layer_counters_.back()++;
              else {
                prevZ = l.maxZwithHybrids();
                double len = l.maxZwithHybrids() - l.minZwithHybrids();
                layer_counters_.push_back(1);
                length_offset_list_.push_back(std::make_pair(len, l.maxZwithHybrids()));
                if (l.minZwithHybrids() > 0) layers_short_.push_back(std::pair<int, double>(radius_list_io_.size() - 1, l.minZwithHybrids()));
              }
            }
            index++;
          }
*/

          const std::vector<int>& layer_counters() const { return layer_counters_; }
        };


        BarrelVisitor v(radius_list_io, length_offset_list, real_index, layers_short, barrel_has_services);
        tracker.accept(v);

        return v.layer_counters();
    }
    
    /**
     * This core function analyses the geometrical properties of the endcaps and regroups the values into a more easily accessible format.
     * @param endcap_layers A reference to the collection of discs in the existing tracker object; the source data
     * @param radius_list_io A reference to an empty vector of pairs that will contain the inner and outer radius of each disc
     * @param length_offset_list A reference to an empty vector of pairs that will contain the width in z and the leftmost point of each disc
     * @param real_index A reference to an empty vector that will contain the index of each disc within the collection of the tracker object
     * @return A vector listing the number of discs per endcap
     */
    std::vector<int> Usher::TrackerIntRep::analyzeEndcaps(Tracker& tracker, std::vector<std::pair<double, double> >& radius_list_io,
            std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::vector<bool>& endcap_has_services) {
        
        class EndcapVisitor : public ConstGeometryVisitor {
          std::vector<std::pair<double, double> >& radius_list_io_;
          std::vector<std::pair<double, double> >& length_offset_list_;
          std::vector<int>& real_index_;
          std::vector<int> layer_counters_;
          std::vector<bool>& endcap_has_services_;

          double prevR = -1.;
          int index = 0;
        public:
          EndcapVisitor(std::vector<std::pair<double, double> >& radius_list_io, 
                        std::vector<std::pair<double, double> >& length_offset_list, 
                        std::vector<int>& real_index,
                        std::vector<bool>& endcap_has_services) : 
              radius_list_io_(radius_list_io),
              length_offset_list_(length_offset_list),
              real_index_(real_index),
              endcap_has_services_(endcap_has_services) {}
        
          void visit(const Endcap& endcap) {
            int counter = 0;

            endcap_has_services_.push_back(!endcap.skipServices());

            for(const Disk& disk : endcap.disks()) {
              if (disk.maxZwithHybrids() > 0.0) {
                length_offset_list_.push_back(std::make_pair(disk.maxZwithHybrids() - disk.minZwithHybrids(), disk.minZwithHybrids()));
                real_index_.push_back(index);
                counter ++;
              }
              index ++;
            }

            layer_counters_.push_back(counter);
            radius_list_io_.push_back(std::make_pair(endcap.minRwithHybrids(), endcap.maxRwithHybrids()));
          }

          /*
          void visit(const Disk& d) {
            if (d.maxZwithHybrids() > 0.) {
              double len = d.maxZwithHybrids() - d.minZwithHybrids();
              length_offset_list_.push_back(std::make_pair(len, d.minZwithHybrids()));
              real_index_.push_back(index);
              if (prevR == d.maxRwithHybrids()) layer_counters_.back()++;
              else {
                prevR = d.maxRwithHybrids();
                layer_counters_.push_back(1);
                radius_list_io_.push_back(std::make_pair(d.minRwithHybrids(), d.maxRwithHybrids()));
              }
            }
            index++;
          }
           */

          const std::vector<int>& layer_counters() const { return layer_counters_; }
        };
        
        EndcapVisitor v(radius_list_io, length_offset_list, real_index, endcap_has_services);
        tracker.accept(v);

        return v.layer_counters();
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
