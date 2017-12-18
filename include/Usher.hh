//
// File:   Usher.h
// Author: ndemaio
//
// Created on October 28, 2008, 9:10 AM
//

/**
 * @file Usher.h
 * @brief This is the base class header file for the algorithm that places the inactive elements around an existing <i>Tracker</i> object of active modules
 */

#ifndef _USHER_H
#define	_USHER_H

#include <cmath>
#include <list>
#include <utility>
#include <Tracker.hh>
#include <Support.hh>
#include <InactiveSurfaces.hh>
#include <InactiveRing.hh>
#include <InactiveTube.hh>
#include <global_constants.hh>
#include <StopWatch.hh>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Barrel.hh"

namespace insur {
    /*
     * @class Usher
     * @brief The Usher class defines the core geometry creation algorithms for inactive surfaces.
     *
     * It analyses an existing collection of active modules, then goes on to decide where to place services and
     * support structures. It does this by building all the necessary parts on the z+ side and mirroring everything
     * that is not centered on the origin in z to the other side of the xy-plane. The two public functions build inactive
     * surfaces around a tracker and around a pixel detector, respectively. Their main difference is not in the
     * service and support parts themselves, but in some constraints regarding size and whether all parts are
     * inside the tracking volume or not.
     */
    class Usher {
    public:
        Usher() {}
        virtual ~Usher() {}
        virtual InactiveSurfaces& arrange(Tracker& tracker, InactiveSurfaces& is, const std::list<Support*>& supports, bool printstatus = false);
        virtual InactiveSurfaces& arrangePixels(Tracker& pixels, InactiveSurfaces& is, bool printstatus = false);
    protected:
        /**
         * @class TrackerIntRep
         * @brief This is a nested class that groups the geometry information extracted from the <i>Tracker</i>
         * collection of active modules.
         *
         * It analyses and calculates only those values that are relevant for the placement of the inactive surfaces.
         * This also includes some things not directly visible in the tracker information. The most important of these
         * is the polarity flag: for the algorithms in <i>Usher</i>, it is the main indicator of which placement
         * strategy is appropriate for the given module layout.
         */ 
        class TrackerIntRep {
        public:
            TrackerIntRep() { post_analysis = false; }
            virtual ~TrackerIntRep() {}
            int nOfBarrels();
            int nOfEndcaps();
            int nOfLayers(int barrelindex);
            int totalLayers();
            int nOfDiscs(int endcapindex);
            int totalDiscs();
            double innerRadiusLayer(int layerindex);
            double outerRadiusLayer(int layerindex);
            double innerRadiusEndcap(int endcapindex);
            double outerRadiusEndcap(int endcapindex);
            double lengthBarrel(int barrelindex);
            double zOffsetBarrel(int barrelindex);
            double lengthDisc(int discindex);
            double zOffsetDisc(int discindex);
            int realIndexLayer(int tintreplayer);
            int realIndexDisc(int tintrepdisc);
            bool barrelHasServices(int barrelindex);
            bool endcapHasServices(int endcapindex);
            int endcapFromDisc(int discindex);
            std::list<std::pair<int, double> >& shortBarrelsList();
            virtual bool analyze(Tracker& tracker);
            void print();
            bool skipAllServices();
            bool skipAllSupports();
        protected:
            std::vector<int> n_of_layers, n_of_discs, real_index_layer, real_index_disc;
            std::vector<bool> barrel_has_services, endcap_has_services;
            std::list<std::pair<int, double> > short_layers;
            std::vector<std::pair<double, double> > layers_io_radius, endcaps_io_radius, barrels_length_offset, discs_length_offset;
        private:
            bool post_analysis;
            bool skipAllServices_, skipAllSupports_;
            std::vector<int> analyzeBarrels(Tracker& tracker, std::vector<std::pair<double, double> >& radius_list_io,
                    std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::list<std::pair<int, double> >& layers_short, std::vector<bool>& barrel_has_services);
            std::vector<int> analyzeEndcaps(Tracker& tracker, std::vector<std::pair<double, double> >& radius_list_io,
                    std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::vector<bool>& endcap_has_services);
            bool analyzePolarity();
        };
        virtual InactiveSurfaces& arrangeUp(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports);
        virtual InactiveSurfaces& arrangeDown(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports);
        InactiveSurfaces& mirror(TrackerIntRep& tracker, InactiveSurfaces& is);
        virtual InactiveSurfaces& servicesUp(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, bool track_all);
        virtual InactiveSurfaces& servicesDown(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, bool track_all);
        virtual InactiveSurfaces& supportsAll(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, const std::list<Support*>& supports, bool track_all);
    private:
        int servicesOutmostBarrel(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer, int layer);
        InactiveSurfaces& supportsRegularBarrels(TrackerIntRep& tracker, InactiveSurfaces& is);
        InactiveSurfaces& supportsBarrelTubes(TrackerIntRep& tracker, InactiveSurfaces& is, bool track_all);
        InactiveSurfaces& supportsShortBarrels(TrackerIntRep& tracker, InactiveSurfaces& is, double r_outer);
        InactiveSurfaces& supportsEndcaps(TrackerIntRep& tracker, InactiveSurfaces& is, bool track_all);
        InactiveSurfaces& supportsUserDefined(TrackerIntRep& tracker, InactiveSurfaces& is, const std::list<Support*>& supports);
        InactiveSurfaces& joinBarrelEndcapServices(TrackerIntRep& tracker, InactiveSurfaces& is, int begin_b, int end_b, int begin_e, int end_e, int d_offset);
        InactiveSurfaces& addBarrelServiceRing(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final);
        InactiveSurfaces& addBarrelServiceTube(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final);
        InactiveSurfaces& addEndcapServiceTube(InactiveSurfaces& is, double length, double offset, double radius, double width, bool final);
        InactiveSurfaces& addSupportRing(InactiveSurfaces& is, double length, double offset, double radius, double width);
        InactiveSurfaces& addSupportTube(InactiveSurfaces& is, double length, double offset, double radius, double width);
        InactiveSurfaces& addBarrelSupportsUserDefined(InactiveSurfaces& is, TrackerIntRep& tracker, int start, int stop, double z);
        InactiveRing mirrorRing(InactiveElement& blueprint);
        InactiveTube mirrorTube(InactiveElement& blueprint);
        int findBarrelInnerRadius(int barrel, TrackerIntRep& tintrep);
        double findMaxBarrelZ(TrackerIntRep& tintrep);
        std::pair<int, int> findBarrelSupportParams(TrackerIntRep& tracker, bool up);
        std::pair<int, int> findSupportStartStop(TrackerIntRep& tracker, std::pair<int, double> udef, std::pair<int, int> aux, double z, bool up);
        double findBarrelSupportRingZ(TrackerIntRep& tracker, int i) const;
        void print(TrackerIntRep& tintrep, InactiveSurfaces& is, bool full_summary);
    };
}
#endif	/* _USHER_H */

