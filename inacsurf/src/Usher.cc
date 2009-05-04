/**
 * @file Usher.cc
 * @brief This is the base class implementation for the algorithm that places the inactive elements around an existing <i>Tracker</i> object of active modules.
 */

#include <Usher.h>
namespace insur {
    // public
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
    InactiveSurfaces& Usher::arrangeUp(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        std::cout << "Arranging UP configuration..." << std::endl;
        is = servicesUp(tracker, is);
        is = supportsAll(tracker, is, geomfile);
        std::cout << "Arrangement done." << std::endl;
        return is;
    }
    
    InactiveSurfaces& Usher::arrangeDown(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        std::cout << "Arranging DOWN configuration..." << std::endl;
        is = servicesDown(tracker, is);
        is = supportsAll(tracker, is, geomfile);
        std::cout << "Arrangement done." << std::endl;
        return is;
    }
    
    InactiveSurfaces& Usher::mirror(TrackerIntRep& tracker, InactiveSurfaces& is) {
        std::cout << "Mirroring barrel services...";
        unsigned int half = is.getBarrelServices().size();
        for (unsigned int i = 0; i < half; i++) {
            InactiveElement& blueprint = is.getBarrelServicePart(i);
            if (blueprint.isVertical()) {
                InactiveRing ring = mirrorRing(blueprint);
                if (ring.getFeederType() == InactiveElement::tracker) {
                    bool short_layer = false;
                    std::list<std::pair<int, double> >::const_iterator iter = tracker.shortBarrelsList().begin();
                    while (iter != tracker.shortBarrelsList().end()) {
                        if (blueprint.getFeederIndex() == tracker.realIndexLayer(iter->first)) {
                            short_layer = true;
                            break;
                        }
                        else iter++;
                    }
                    if (short_layer) ring.setFeederIndex(ring.getFeederIndex() - 1);
                }
                else if (ring.getFeederType() == InactiveElement::barrel) {
                    ring.setFeederIndex(is.getBarrelServices().size() - 1);
                }
                if (blueprint.getNeighbourType() == InactiveElement::barrel) {
                    ring.setNeighbourIndex(is.getBarrelServices().size() - i + blueprint.getNeighbourIndex());
                }
                else if (blueprint.getNeighbourType() == InactiveElement::endcap) {
                    ring.setNeighbourIndex(2 * tracker.totalDiscs() - blueprint.getNeighbourIndex() - 1);
                }
                is.addBarrelServicePart(ring);
            }
            else {
                InactiveTube tube = mirrorTube(blueprint);
                bool short_layer = false;
                std::list<std::pair<int, double> >::const_iterator iter = tracker.shortBarrelsList().begin();
                while (iter != tracker.shortBarrelsList().end()) {
                    if (blueprint.getFeederIndex() == iter->first) {
                        short_layer = true;
                        break;
                    }
                    else iter++;
                }
                if (short_layer) tube.setFeederIndex(tube.getFeederIndex() - 1);
                is.addBarrelServicePart(tube);
            }
        }
        std::cout << "done." << std::endl << "Mirroring endcap services...";
        half = is.getEndcapServices().size();
        for (unsigned int i = 0; i < half; i++) {
            InactiveElement& blueprint = is.getEndcapServicePart(i);
            if (blueprint.isVertical()) {
                InactiveRing ring = mirrorRing(blueprint);
                if (blueprint.getFeederIndex() != -1) ring.setFeederIndex(2 * tracker.totalDiscs() - blueprint.getFeederIndex() - 1);
                is.addEndcapServicePart(ring);
            }
            else {
                InactiveTube tube = mirrorTube(blueprint);
                if (blueprint.getFeederIndex() != -1) tube.setFeederIndex(2 * tracker.totalDiscs() - blueprint.getFeederIndex() - 1);
                if (blueprint.getNeighbourType() == InactiveElement::barrel) {
                    tube.setNeighbourIndex(is.getBarrelServices().size() / 2 + blueprint.getNeighbourIndex());
                }
                else if (blueprint.getNeighbourType() == InactiveElement::endcap) {
                    tube.setNeighbourIndex(is.getEndcapServices().size() - i + blueprint.getNeighbourIndex());
                }
                is.addEndcapServicePart(tube);
            }
        }
        std::cout << "done." << std::endl << "Mirroring supports...";
        half = is.getSupports().size();
        for (unsigned int i = 0; i < half; i ++) {
            InactiveElement& blueprint = is.getSupportPart(i);
            if (blueprint.getZOffset() > 0) {
                if (blueprint.isVertical()) {
                    InactiveRing ring = mirrorRing(blueprint);
                    is.addSupportPart(ring);
                }
                else {
                    InactiveTube tube = mirrorTube(blueprint);
                    is.addSupportPart(tube);
                }
            }
        }
        std::cout << "done." << std::endl;
        return is;
    }
    
    InactiveSurfaces& Usher::servicesUp(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw;
        int k = 0, h;
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
                    rw = tracker.innerRadiusLayer(k + 1) - ri - epsilon;
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
                    ri = tracker.innerRadiusLayer(l) - rw - epsilon;
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
                    }
                    else {
                        if (j == 0) {
                            double ro;
                            ro = tracker.innerRadiusLayer(findBarrelInnerRadius(tracker.nOfBarrels() - tracker.nOfEndcaps() + i, tracker));
                            ro = ro - rw - epsilon;
                            zo = tracker.zOffsetDisc(k - 1) + epsilon;
                            zl = tracker.zOffsetBarrel(tracker.nOfBarrels() - tracker.nOfEndcaps() + i) - zo + volume_width + epsilon;
                            if ((i == tracker.nOfEndcaps() - 1) && (j == tracker.nOfDiscs(i) - 1)) is = addEndcapServiceTube(is, zl, zo, ro, rw, true);
                            else is = addEndcapServiceTube(is, zl, zo, ro, rw, false);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setCategory(MaterialProperties::e_ser);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourType(InactiveElement::endcap);
                            is.getEndcapServicePart(is.getEndcapServices().size() - 1).setNeighbourIndex(is.getEndcapServices().size() - 2);
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
                    }
                    k++;
                }
            }
            // join barrel and endcap services
            is = joinBarrelEndcapServices(tracker, is, 0, is.getBarrelServices().size(), 0, is.getEndcapServices().size(), tracker.totalDiscs());
        }
        return is;
    }
    
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
        for (unsigned int i = 0; i < tracker.totalDiscs(); i++) {
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
        }
        return is;
    }
    
    InactiveSurfaces& Usher::supportsAll(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        // outer tube
        is = addSupportTube(is, max_length, 0.0 - max_length / 2.0, outer_radius, volume_width);
        is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
        // inner tube
        is = addSupportTube(is, max_length, 0.0 - max_length / 2.0, inner_radius - volume_width, volume_width);
        is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::t_sup);
        // barrel endpoints
        is = supportsRegularBarrels(tracker, is);
        // short barrels
        is = supportsShortBarrels(tracker, is);
        // endcaps
        is = supportsEndcaps(tracker, is);
        // rings from config file
        is = supportsUserDefined(tracker, is, geomfile);
        return is;
    }
    
// private
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
    
    InactiveSurfaces& Usher::supportsRegularBarrels(TrackerIntRep& tracker, InactiveSurfaces& is) {
        double zl, zo, ri, rw;
        int k = 0;
        zl = volume_width;
        for (int i = 0; i < tracker.nOfBarrels(); i++) {
            if (tracker.nOfLayers(i) > 1) {
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
    
    InactiveSurfaces& Usher::supportsShortBarrels(TrackerIntRep& tracker, InactiveSurfaces& is) {
        bool regular;
        unsigned int start, stop;
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
                ri = tracker.innerRadiusEndcap(i) - rw - epsilon;
                is  = addSupportTube(is, zl, zo, ri, rw);
                is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::e_sup);
            }
            k = k + tracker.nOfDiscs(i);
        }
        return is;
    }
    
    InactiveSurfaces& Usher::supportsUserDefined(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile) {
        int start, stop;
        double z, r, w;
        configParser persian;
        std::list<double>* extras = NULL;
        std::list<double>::iterator iter;
        extras = persian.parseSupportsFromFile(geomfile);
        if (extras) {
            for (iter = extras->begin(); iter != extras->end(); iter++) {
                if (*iter < max_length) {
                    start = 0;
                    stop = tracker.totalLayers() - 1;
                    z = *iter - volume_width / 2.0;
                    for (int i = 0; i < tracker.nOfBarrels(); i++) {
                        if (*iter > tracker.lengthBarrel(i) / 2.0) start = start + tracker.nOfLayers(i);
                        else break;
                    }
                    for (int i = tracker.nOfBarrels() - 1; i >= 0; i--) {
                        if (*iter > tracker.lengthBarrel(i) / 2.0) stop = stop - tracker.nOfLayers(i);
                        else break;
                    }
                    r = inner_radius + epsilon;
                    if (start > stop) w = tracker.innerRadiusLayer(tracker.totalLayers() - 1) - r - epsilon;
                    else {
                        if (stop < (int)tracker.totalLayers() - 1) {
                            r = tracker.outerRadiusLayer(stop) + epsilon;
                            w = tracker.innerRadiusLayer(tracker.totalLayers() - 1) - r - epsilon;
                            is = addSupportRing(is, volume_width, z, r, w);
                            is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::u_sup);
                            r = inner_radius + epsilon;
                        }
                        w = tracker.innerRadiusLayer(start) - r - epsilon;
                    }
                    is = addSupportRing(is, volume_width, z, r, w);
                    is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::u_sup);
                    for (int k = start; k < stop; k++) {
                        r = tracker.outerRadiusLayer(k) + epsilon;
                        w = tracker.innerRadiusLayer(k + 1) - r - epsilon;
                        is = addSupportRing(is, volume_width, z, r, w);
                        is.getSupportPart(is.getSupports().size() - 1).setCategory(MaterialProperties::u_sup);
                    }
                }
            }
            delete extras;
        }
        return is;
    }
    
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
    
    InactiveSurfaces& Usher::addSupportRing(InactiveSurfaces& is, double length, double offset, double radius, double width) {
        InactiveRing ir;
        ir.setZLength(length);
        ir.setZOffset(offset);
        ir.setInnerRadius(radius);
        ir.setRWidth(width);
        is.addSupportPart(ir);
        return is;
    }
    
    InactiveSurfaces& Usher::addSupportTube(InactiveSurfaces& is, double length, double offset, double radius, double width) {
        InactiveTube it;
        it.setZLength(length);
        it.setZOffset(offset);
        it.setInnerRadius(radius);
        it.setRWidth(width);
        is.addSupportPart(it);
        return is;
    }
    
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
        return ir;
    }
    
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
        return it;
    }
    
    int Usher::findBarrelInnerRadius(int barrel, TrackerIntRep& tintrep) {
        int l = 0;
        for (int i = 0; i < barrel; i ++) {
            l = l + tintrep.nOfLayers(i);
        }
        return l;
    }
    
    double Usher::findMaxBarrelZ(TrackerIntRep& tintrep) {
        double z_max = -1.0;
        for (int i = 0; i < tintrep.nOfBarrels(); i++) {
            if (z_max < tintrep.zOffsetBarrel(i)) z_max = tintrep.zOffsetBarrel(i);
        }
        return z_max;
    }
    
    void Usher::print(TrackerIntRep& tintrep, InactiveSurfaces& is, bool full_summary) {
        std::cout << std::endl << "Current state of the system is:" << std::endl;
        std::cout << "Printing tracker information..." << std::endl;
        tintrep.print();
        std::cout << "Printing inactive surfaces..." << std::endl;
        is.print(full_summary);
        std::cout << "Usher internal status printed." << std::endl;
    }
    
// nested class public
    int Usher::TrackerIntRep::nOfBarrels() {
        if (post_analysis) return n_of_layers.size();
        else return -1;
    }
    
    int Usher::TrackerIntRep::nOfEndcaps() {
        if (post_analysis) return n_of_discs.size();
        else return -1;
    }
    
    int Usher::TrackerIntRep::nOfLayers(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < n_of_layers.size())) return n_of_layers.at(barrelindex);
        }
        return -1;
    }
    
    unsigned int Usher::TrackerIntRep::totalLayers() {
        unsigned int sum = 0;
        for (unsigned int i = 0; i < n_of_layers.size(); i++) sum = sum + n_of_layers.at(i);
        return sum;
    }
    
    int Usher::TrackerIntRep::nOfDiscs(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < n_of_discs.size())) return n_of_discs.at(endcapindex);
        }
        return -1;
    }
    
    unsigned int Usher::TrackerIntRep::totalDiscs() {
        unsigned int sum = 0;
        for (unsigned int i = 0; i < n_of_discs.size(); i++) sum = sum + n_of_discs.at(i);
        return sum;
    }
    
    double Usher::TrackerIntRep::innerRadiusLayer(int layerindex) {
        if (post_analysis) {
            if ((layerindex >= 0) && ((unsigned int)layerindex < layers_io_radius.size())) return layers_io_radius.at(layerindex).first;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::outerRadiusLayer(int layerindex) {
        if (post_analysis) {
            if ((layerindex >= 0) && ((unsigned int)layerindex < layers_io_radius.size())) return layers_io_radius.at(layerindex).second;
        }
        return -1;
        return -1;
    }
    
    double Usher::TrackerIntRep::innerRadiusEndcap(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < endcaps_io_radius.size())) return endcaps_io_radius.at(endcapindex).first;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::outerRadiusEndcap(int endcapindex) {
        if (post_analysis) {
            if ((endcapindex >= 0) && ((unsigned int)endcapindex < endcaps_io_radius.size())) return endcaps_io_radius.at(endcapindex).second;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::lengthBarrel(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < barrels_length_offset.size())) return barrels_length_offset.at(barrelindex).first;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::zOffsetBarrel(int barrelindex) {
        if (post_analysis) {
            if ((barrelindex >= 0) && ((unsigned int)barrelindex < barrels_length_offset.size())) return barrels_length_offset.at(barrelindex).second;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::lengthDisc(int discindex) {
        if (post_analysis) {
            if ((discindex >= 0) && ((unsigned int)discindex < discs_length_offset.size())) return discs_length_offset.at(discindex).first;
        }
        return -1;
    }
    
    double Usher::TrackerIntRep::zOffsetDisc(int discindex) {
        if (post_analysis) {
            if ((discindex >= 0) && ((unsigned int)discindex < discs_length_offset.size())) return discs_length_offset.at(discindex).second;
        }
        return -1;
    }
    
    int Usher::TrackerIntRep::realIndexLayer(int tintreplayer) {
        if (post_analysis) {
            if ((tintreplayer >= 0) && ((unsigned int)tintreplayer < real_index_layer.size())) return real_index_layer.at(tintreplayer);
        }
        return -1;
    }
    
    int Usher::TrackerIntRep::realIndexDisc(int tintrepdisc) {
        if (post_analysis) {
            if ((tintrepdisc >= 0) && ((unsigned int)tintrepdisc < real_index_disc.size())) return real_index_disc.at(tintrepdisc);
        }
        return -1;
    }
    
    std::list<std::pair<int, double> >& Usher::TrackerIntRep::shortBarrelsList() { return short_layers; }
    
    bool Usher::TrackerIntRep::analyze(Tracker& tracker) {
        bool up;
        n_of_layers = analyzeBarrels(*(tracker.getBarrelLayers()), layers_io_radius, barrels_length_offset, real_index_layer, short_layers);
        n_of_discs = analyzeEndcaps(*(tracker.getEndcapLayers()), endcaps_io_radius, discs_length_offset, real_index_disc);
        post_analysis = true;
        up = analyzePolarity();
        return up;
    }
    
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
                for (unsigned int i = 0; i < totalLayers(); i++) std::cout << " " << innerRadiusLayer(i);
                std::cout << "." << std::endl;
                std::cout << std::endl << "Number of endcaps: " << nOfEndcaps() << std::endl;
                for (int i = 0; i < nOfEndcaps(); i++) {
                    std::cout << "Endcap " << i << ": " << nOfDiscs(i);
                    std::cout << " discs of radius " << outerRadiusEndcap(i) << ".";
                    std::cout << std::endl;
                }
                std::cout << std::endl << "Total number of discs: " << totalDiscs() << std::endl;
                std::cout << "Discs at an offset of:";
                for (unsigned int i = 0; i < totalDiscs(); i++) std::cout << " " << zOffsetDisc(i);
                std::cout << "." << std::endl << std::endl;
                std::cout << "Endcaps extend to z = " << (zOffsetDisc(totalDiscs() - 1) + lengthDisc(totalDiscs() - 1)) << "." << std::endl;
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
    
    bool Usher::TrackerIntRep::analyzePolarity() {
        if (nOfBarrels() > 1) {
            for (int i = 0; i < nOfBarrels() - 1; i++) {
                if (zOffsetBarrel(i) > zOffsetBarrel(i + 1)) return false;
            }
        }
        return true;
    }
}
