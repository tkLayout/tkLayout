// 
// File:   TrackerActions.h
// Author: ndemaio
//
// Created on November 14, 2008, 10:09 AM
//

/**
 * @file TrackerActions.h
 * @brief This is the header file for the adapter class that encapsulates active volume creation
 */

#ifndef _TRACKERACTIONS_H
#define	_TRACKERACTIONS_H

#include <string>
#include <tracker.hh>
#include <configparser.hh>
namespace insur {
    /**
     * @class TrackerActions
     * @brief This class wraps the steps to create a <i>Tracker</i> object of active volumes into two simple functions.
     *
     * Their main difference is the number of configuration files: if both a geometry and a settings file are given, both
     * are passed through to the parser and the geometry is built according to the geometry file, then dressed according
     * to the module parameters in the settings file. If only the geometry file is supplied, the tracker is built with the default
     * parameters; no settings files or resource directories are created in the process.
     */
    class TrackerActions {
    public:
        Tracker* createActiveSurfaces(std::string geomfile, std::string settingsfile);
        Tracker* createActiveSurfaces(std::string geomfile);
    };
}
#endif	/* _TRACKERACTIONS_H */

