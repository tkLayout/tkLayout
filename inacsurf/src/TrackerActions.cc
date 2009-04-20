/**
 * @file TrackerActions.cc
 * @brief This is the implementation of the adapter class that encapsulates active volume creation
 */

#include <TrackerActions.h>
namespace insur {
    
    /**
     * This function encapsulates full geometry generation and module dressing, as speciefied in the two supplied configuration files.
     * @param geomfile The filename (absolute or relative) of the geometry configuration file
     * @param settingsfile The filename (absolute or relative) of the setting configuration file
     * @return A pointer to the newly created tracker object
     */
    Tracker* TrackerActions::createActiveSurfaces(std::string geomfile, std::string settingsfile) {
        Tracker* tracker;
        configParser persian;
        if ((tracker = persian.parseFile(geomfile))) persian.dressTracker(tracker, settingsfile);
        return tracker;
    }
    
    /**
     * This function encapsulates geometry generation only, as specified in the supplied configuration file.
     * The modules keep their initial default values for things like the number of strips/chips and module type.
     * @param geomfile The filename (absolute or relative) of the geometry configuration file
     * @return A pointer to the newly created tracker object
     */
    Tracker* TrackerActions::createActiveSurfaces(std::string geomfile) {
        configParser persian;
        return persian.parseFile(geomfile);
    }
}
