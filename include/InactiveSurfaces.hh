//
// File:   InactiveSurfaces.h
// Author: ndemaio
//
// Created on October 13, 2008, 4:55 PM
//

/**
 * @file InactiveSurfaces.h
 * @brief This is the header file for the container class for inactive tracker elements
 */


#ifndef _INACTIVESURFACES_H
#define	_INACTIVESURFACES_H

#include <map>
#include <string>
#include <vector>
#include <InactiveElement.hh>
namespace insur {
  /**
   * @class InactiveSurfaces
   * @brief This is the top-level container class  for the inactive surfaces.
   *
   * It contains lists of all subgroups of inactive volumes: the services list and the supporting parts list.
   * It provides access functions to them or their individual elements that typically copy a new element to
   * its place at the end of the vector or return a reference to a requested element. It also stores the type of
   * configuration (UP or DOWN) in a boolean flag. Some of the access functions to individual elements
   * may throw an exception if the requested index is out of range.
   */
  class InactiveSurfaces {
  public:
    InactiveSurfaces() {}
    virtual ~InactiveSurfaces() {}
    // services
    void addBarrelServicePart(InactiveElement service);
    InactiveElement& getBarrelServicePart(int index); // throws exception
    std::vector<InactiveElement>::iterator removeBarrelServicePart(int index);
    std::vector<InactiveElement>& getBarrelServices(); // may return empty vector
    void addEndcapServicePart(InactiveElement service);
    InactiveElement& getEndcapServicePart(int index); // throws exception
    std::vector<InactiveElement>::iterator removeEndcapServicePart(int index);
    std::vector<InactiveElement>& getEndcapServices(); // may return empty vector
    // supports
    void addSupportPart(InactiveElement support);
    InactiveElement& getSupportPart(int index); // throws exception
    std::vector<InactiveElement>::iterator removeSupportPart(int index);
    std::vector<InactiveElement>& getSupports(); // may return empty vector
    // layout flag
    bool isUp();
    void setUp(bool up);
    void print(bool full_summary);
  protected:
    //layout flag
    bool is_up;
    // element collections
    std::vector<InactiveElement> barrelservices, endcapservices, supports;
  private:

  };
}
#endif	/* _INACTIVESURFACES_H */

