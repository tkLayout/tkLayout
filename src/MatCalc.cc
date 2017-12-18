/**
 * @file MatCalc.cc
 * @brief This class implements the algorithms that allocate materials to the various types of structures in the tracker
 */

#include <MatCalc.hh>
namespace insur {
    // public
    /**
     * This function returns the initialisation status.
     * @return True if the internal data structures have been initialised from a config file, false otherwise
     */
    bool MatCalc::initDone() { return init_done; }
    
    /**
     * This function sets the initialisation status.
     * @param yes True if the internal data structures contain information, false otherwise
     */
    void MatCalc::initDone(bool yes) { init_done = yes; }
    
    /**
     * This function resets the internal data structures.
     */
    void MatCalc::reset() {
        init_done = false;
    }
    
    /**
     * Getter for the internal global material table. The material table always exists but
     * it may be empty if the class has not been initialised.
     * @return A reference to the global material table
     */
    MaterialTable& MatCalc::getMaterialTable() { return mt; }
    
}
