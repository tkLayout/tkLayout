// 
// File:   MatCalcDummy.h
// Author: ndemaio
//
// Created on April 21, 2009, 2:43 PM
//

/**
 * @file MatCalcDummy.h
 * @brief This is the header file for a simplified material assignment class
 */

#ifndef _MATCALCDUMMY_H
#define	_MATCALCDUMMY_H

#include <string>
#include <MatCalc.hh>
namespace insur {
  /**
   * @class MatCalc
   * @brief This descendant of <i>MatCalc</i> is a highly simplified test class for assigning materials to volumes.
   *
   * Instead of parsing a config file and using its contents to calculate the material mix of the various volumes, this
   * class simply hands out a constant mass to every volume in one of the three categories <i>module</i>,
   * <i>service</i> or <i>support</i>. It uses its own, internal material table for this and never touches the global
   * material list that comes in a separate file either.
   */
  class MatCalcDummy : public MatCalc {
  public:
    MatCalcDummy();
    virtual ~MatCalcDummy() {}
    virtual bool calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps);
    virtual bool calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps);
    virtual bool calculateBarrelServiceMaterials(
      std::vector<std::vector<ModuleCap> >& barrelcaps,
      std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
    virtual bool calculateEndcapServiceMaterials(
      std::vector<std::vector<ModuleCap> >& endcapcaps,
      std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
    virtual bool calculateSupportMaterials(std::vector<InactiveElement>& supports);
  private:
    std::string ta, ts, tl;
    double ab, ae, sb, se, lb, le;
    void init();
  };
}
#endif	/* _MATCALCDUMMY_H */

