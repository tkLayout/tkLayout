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
#include <tracker.hh>
#include <configparser.hh>
#include <InactiveSurfaces.h>
#include <InactiveRing.h>
#include <InactiveTube.h>
#include <global_constants.h>
namespace insur {
    /*
     * @class Usher
     * @brief The Usher class defines the core geometry creation algorithm for inactive surfaces.
     *
     * It analyses an existing collection of active modules, then goes on to decide where to place services and
     * support structures. It does this by building all the necessary parts on the z+ side and mirroring everything
     * that is not centered on the origin in z to the other side of the xy-plane.
     */
    class Usher {
    public:
        Usher() {}
        virtual ~Usher() {}
        virtual InactiveSurfaces& arrange(Tracker& tracker, InactiveSurfaces& is, std::string geomfile, bool printstatus = false);
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
            std::list<std::pair<int, double> >& shortBarrelsList();
            virtual bool analyze(Tracker& tracker);
            void print();
        protected:
            std::vector<int> n_of_layers, n_of_discs, real_index_layer, real_index_disc;
            std::list<std::pair<int, double> > short_layers;
            std::vector<std::pair<double, double> > layers_io_radius, endcaps_io_radius, barrels_length_offset, discs_length_offset;
        private:
            bool post_analysis;
            std::vector<int> analyzeBarrels(std::vector<Layer*>& barrel_layers, std::vector<std::pair<double, double> >& radius_list_io,
                    std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index, std::list<std::pair<int, double> >& layers_short);
            std::vector<int> analyzeEndcaps(std::vector<Layer*>& endcap_layers, std::vector<std::pair<double, double> >& radius_list_io,
                    std::vector<std::pair<double, double> >& length_offset_list, std::vector<int>& real_index);
            bool analyzePolarity();
        };
        virtual InactiveSurfaces& arrangeUp(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile);
        virtual InactiveSurfaces& arrangeDown(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile);
        InactiveSurfaces& mirror(TrackerIntRep& tracker, InactiveSurfaces& is);
        virtual InactiveSurfaces& servicesUp(TrackerIntRep& tracker, InactiveSurfaces& is);
        virtual InactiveSurfaces& servicesDown(TrackerIntRep& tracker, InactiveSurfaces& is);
        virtual InactiveSurfaces& supportsAll(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile);
    private:
        int servicesOutmostBarrel(TrackerIntRep& tracker, InactiveSurfaces& is, int layer);
        InactiveSurfaces& supportsRegularBarrels(TrackerIntRep& tracker, InactiveSurfaces& is);
        InactiveSurfaces& supportsBarrelTubes(TrackerIntRep& tracker, InactiveSurfaces& is);
        InactiveSurfaces& supportsShortBarrels(TrackerIntRep& tracker, InactiveSurfaces& is);
        InactiveSurfaces& supportsEndcaps(TrackerIntRep& tracker, InactiveSurfaces& is);
        InactiveSurfaces& supportsUserDefined(TrackerIntRep& tracker, InactiveSurfaces& is, std::string geomfile);
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
        void print(TrackerIntRep& tintrep, InactiveSurfaces& is, bool full_summary);
    };
}
#endif	/* _USHER_H */

