//
// File:   MaterialBudget.h
// Author: ndemaio
//
// Created on March 10, 2009, 6:40 PM
//

/**
 * @file MaterialBudget.h
 * @brief This is the header file for the material budget container class
 */

#ifndef _MATERIALBUDGET_H
#define	_MATERIALBUDGET_H

#include <vector>
#include "Tracker.hh"
#include <InactiveSurfaces.hh>
#include <ModuleCap.hh>
#include <MatCalc.hh>
namespace insur {
  /**
   * Errors that may occur during operations
   */
  static const std::string err_materials_supports = "Error calculating materials for supports.";
  static const std::string err_materials_bservices = "Error calculating materials for barrel services.";
  static const std::string err_materials_eservices = "Error calculating materials for endcap services.";
  static const std::string err_materials_bmodules = "Error calculating materials for barrel modules.";
  static const std::string err_materials_emodules = "Error calculating materials for endcap modules.";

  /**
   * @class MaterialBudget
   * @brief This class integrates information from a <i>Tracker</i> and an <i>InactiveSurface</i> instance with a collection of <i>ModuleCap</i> instances to provide the full material budget of a tracker.
   *
   * Its main function accepts an instance of a material calculator as an input parameter in order to use that calculator's
   * more specialised functions to assign mixtures of materials to the different categories of volumes in the tracker
   * geometry. Afterwards, each individual volume is ready to calculate its total mass, its radiation length and its interaction
   * length, which is also taken care of by the calculator class. A material budget that has been filled in this way can then
   * be passed on to be analysed by an instance of an <i>Analyzer</i> class.
   */
  class MaterialBudget {
  public:
    MaterialBudget(Tracker& tr, InactiveSurfaces& is);
    virtual ~MaterialBudget();
    Tracker& getTracker();
    InactiveSurfaces& getInactiveSurfaces();
    std::vector<std::vector<ModuleCap> >& getBarrelModuleCaps();
    std::vector<std::vector<ModuleCap> >& getEndcapModuleCaps();
    std::vector<InactiveElement> getAllServices();
    void print();
  protected:
    Tracker* tracker;
    InactiveSurfaces* inactive;
    std::vector<std::vector<ModuleCap> > capsbarrelmods, capsendmods;
    int onBoundary(std::vector<std::vector<ModuleCap> >& source, int layer); //throws exception
  private:
    MaterialBudget();
    MaterialBudget(const MaterialBudget& budget);
    MaterialBudget& operator=(const MaterialBudget& budget);
  };
}
#endif	/* _MATERIALBUDGET_H */

