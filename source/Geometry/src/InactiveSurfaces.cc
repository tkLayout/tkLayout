/**
 * @file InactiveSurfaces.cc
 * @brief This is the implementation of the container class for inactive tracker elements
 */

#include <InactiveSurfaces.h>
#include <InactiveElement.h>

namespace insur {

  InactiveSurfaces::InactiveSurfaces() : m_isUp(false) {}
  InactiveSurfaces::~InactiveSurfaces() {}

    /*===== services =====*/
    /**
     * Add a single inactive element to the list of barrel services by copying it.
     * @param service The element that is appended to the list of barrel service parts
     */
    void InactiveSurfaces::addBarrelServicePart(InactiveElement service) {
        m_barrelServices.push_back(service);
    }
    
    /**
     * Access an individual element in the list of barrel services by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested barrel service part
     * @return A reference to the requested barrel service part
     */
    InactiveElement& InactiveSurfaces::getBarrelServicePart(int index) {
        return m_barrelServices.at(index);
    }
    
    /**
     * Remove a single barrel element identified by its index from the list. If the removed barrel service part
     * was the last on the list or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the barrel element that will be removed
     * @return An interator to the barrel element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeBarrelServicePart(int index) {
        if ((index >= 0) && ((unsigned int)index < m_barrelServices.size())) return m_barrelServices.erase(m_barrelServices.begin() + index);
        return m_barrelServices.end();
    }
    
    /**
     * Access the full list of barrel service parts at once.
     * @return A reference to the internal barrel service vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getBarrelServices() {
        return m_barrelServices;
    }
    
    /**
     * Add a single inactive element to the list of endcap services by copying it.
     * @param service The element that is appended to the list of endcap service parts
     */
    void InactiveSurfaces::addEndcapServicePart(InactiveElement service) {
        m_endcapServices.push_back(service);
    }
    
    /**
     * Access an individual element in the list of endcap services by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested endcap service part
     * @return A reference to the requested endcap service part
     */
    InactiveElement& InactiveSurfaces::getEndcapServicePart(int index) {
        return m_endcapServices.at(index);
    }
    
    /**
     * Remove a single endcap element identified by its index from the list. If the removed endcap service part
     * was the last on the list or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the endcap element that will be removed
     * @return An interator to the endcap element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeEndcapServicePart(int index) {
        if ((index >= 0) && ((unsigned int)index < m_endcapServices.size())) return m_endcapServices.erase(m_endcapServices.begin() + index);
        return m_endcapServices.end();
    }
    
    /**
     * Access the full list of endcap services at once.
     * @return A reference to the internal endcap services vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getEndcapServices() {
        return m_endcapServices;
    }
    /*===== supports =====*/
    /**
     * Add a single inactive element to the list of supports by copying it.
     * @param support The element that is appended to the list of support parts
     */
    void InactiveSurfaces::addSupportPart(InactiveElement support) {
        m_supports.push_back(support);
    }
    
    /**
     * Access an individual element in the list of supports by its index.
     * The internal vector will throw an exception if the index is out of range.
     * @param index The index of the requested support part
     * @return A reference to the requested support part
     */
    InactiveElement& InactiveSurfaces::getSupportPart(int index) { // throws exception
        return m_supports.at(index);
    }
    
    /**
     *Remove a single element identified by its index from the list. If the removed support part was the last on the list
     * or the given index is out of range, the returned iterator will point to <i>end()</i>.
     * @param index The index of the element that will be removed
     * @return An interator to the element immediately after the removed one
     */
    std::vector<InactiveElement>::iterator InactiveSurfaces::removeSupportPart(int index) {
        if ((index >= 0) && ((unsigned int)index < m_supports.size())) return m_supports.erase(m_supports.begin() + index);
        return m_supports.end();
    }
    
    /**
     * Access the full list of support parts at once.
     * @return A reference to the internal supports vector
     */
    std::vector<InactiveElement>& InactiveSurfaces::getSupports() { // may return empty vector
        return m_supports;
    }
    
    /*===== Flag and printing =====*/
    /**
     * Query the UP/DOWN flag.
     * @return True if the configuration is of type UP, false otherwise
     */
    bool InactiveSurfaces::isUp() { return m_isUp; }
    
    /**
     * Set the UP/DOWN flag.
     * @param up The new state of the flag
     */
    void InactiveSurfaces::setUp(bool up) { m_isUp = up; }
    
    /**
     * Print the contents of the collection.
     * @param full_summary A flag to switch verbose output on or off
     */
    void InactiveSurfaces::print(bool full_summary = true) {
        std::cout << "Number of barrel service elements: " << m_barrelServices.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < m_barrelServices.size(); i++) {
                std::cout << "Service element " << i << ":" << std::endl;
                m_barrelServices.at(i).print();
                std::cout << std::endl;
            }
        }
        std::cout << "Number of endcap service elements: " << m_endcapServices.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < m_endcapServices.size(); i++) {
                std::cout << "Service element " << i << ":" << std::endl;
                m_endcapServices.at(i).print();
                std::cout << std::endl;
            }
        }
        std::cout << "Number of support elements: " << m_supports.size() << std::endl;
        if (full_summary) {
            for (unsigned int i = 0; i < m_supports.size(); i++) {
                std::cout << "Support element " << i << ":" << std::endl;
                m_supports.at(i).print();
                std::cout << std::endl;
            }
        }
    }
}
