/**
 * @file InactiveSurfaces.cc
 * @brief This is the implementation of the container class for inactive tracker elements
 */

#include <InactiveSurfaces.hh>
namespace insur {
    /*===== services =====*/
    /**
     * Add a single inactive element to the list of barrel services by copying it.
     * @param service The element that is appended to the list of barrel service parts
     */
    void InactiveSurfaces::addBarrelServicePart(InactiveElement service) {
        barrelservices.push_back(service);
    }
    
    /**
     * Access an individual element in the list of barrel services by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested barrel service part
     * @return A reference to the requested barrel service part
     */
    InactiveElement& InactiveSurfaces::getBarrelServicePart(int index) {
        return barrelservices.at(index);
    }
    
    /**
     * Remove a single barrel element identified by its index from the list. If the removed barrel service part
     * was the last on the list or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the barrel element that will be removed
     * @return An interator to the barrel element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeBarrelServicePart(int index) {
        if ((index >= 0) && ((unsigned int)index < barrelservices.size())) return barrelservices.erase(barrelservices.begin() + index);
        return barrelservices.end();
    }
    
    /**
     * Access the full list of barrel service parts at once.
     * @return A reference to the internal barrel service vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getBarrelServices() {
        return barrelservices;
    }
    
    /**
     * Add a single inactive element to the list of endcap services by copying it.
     * @param service The element that is appended to the list of endcap service parts
     */
    void InactiveSurfaces::addEndcapServicePart(InactiveElement service) {
        endcapservices.push_back(service);
    }
    
    /**
     * Access an individual element in the list of endcap services by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested endcap service part
     * @return A reference to the requested endcap service part
     */
    InactiveElement& InactiveSurfaces::getEndcapServicePart(int index) {
        return endcapservices.at(index);
    }
    
    /**
     * Remove a single endcap element identified by its index from the list. If the removed endcap service part
     * was the last on the list or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the endcap element that will be removed
     * @return An interator to the endcap element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeEndcapServicePart(int index) {
        if ((index >= 0) && ((unsigned int)index < endcapservices.size())) return endcapservices.erase(endcapservices.begin() + index);
        return endcapservices.end();
    }
    
    /**
     * Access the full list of endcap services at once.
     * @return A reference to the internal endcap services vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getEndcapServices() {
        return endcapservices;
    }
    /*===== supports =====*/
    /**
     * Add a single inactive element to the list of supports by copying it.
     * @param support The element that is appended to the list of support parts
     */
    void InactiveSurfaces::addSupportPart(InactiveElement support) {
        supports.push_back(support);
    }
    
    /**
     * Access an individual element in the list of supports by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested support part
     * @return A reference to the requested support part
     */
    InactiveElement& InactiveSurfaces::getSupportPart(int index) { // throws exception
        return supports.at(index);
    }
    
    /**
     *Remove a single element identified by its index from the list. If the removed support part was the last on the list
     * or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the element that will be removed
     * @return An interator to the element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeSupportPart(int index) {
        if ((index >= 0) && ((unsigned int)index < supports.size())) return supports.erase(supports.begin() + index);
        return supports.end();
    }
    
    /**
     * Access the full list of support parts at once.
     * @return A reference to the internal supports vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getSupports() { // may return empty vector
        return supports;
    }
    
    /*===== Flag and printing =====*/
    /**
     * Query the UP/DOWN flag.
     * @return True if the configuration is of type UP, false otherwise
     */
    bool InactiveSurfaces::isUp() { return is_up; }
    
    /**
     * Set the UP/DOWN flag.
     * @param up The new state of the flag
     */
    void InactiveSurfaces::setUp(bool up) { is_up = up; }
    
    /**
     * Print the contents of the collection.
     * @param full_summary A flag to switch verbose output on or off
     */
    void InactiveSurfaces::print(bool full_summary = true) {
        std::cout << "Number of barrel service elements: " << barrelservices.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < barrelservices.size(); i++) {
                std::cout << "Service element " << i << ":" << std::endl;
                barrelservices.at(i).print();
                std::cout << std::endl;
            }
        }
        std::cout << "Number of endcap service elements: " << endcapservices.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < endcapservices.size(); i++) {
                std::cout << "Service element " << i << ":" << std::endl;
                endcapservices.at(i).print();
                std::cout << std::endl;
            }
        }
        std::cout << "Number of support elements: " << supports.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < supports.size(); i++) {
                std::cout << "Support element " << i << ":" << std::endl;
                supports.at(i).print();
                std::cout << std::endl;
            }
        }
    }
}
