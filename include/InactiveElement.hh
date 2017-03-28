//
// File:   InactiveElement.h
// Author: ndemaio
//
// Created on October 16, 2008, 5:11 PM
//

/**
 * @file InactiveElement.h
 * @brief This is the base class header file for a single inactive element
 */

#ifndef _INACTIVEELEMENT_H
#define	_INACTIVEELEMENT_H


#include "MaterialProperties.hh"
#include "global_constants.hh"
namespace insur {
  /**
   * @class InactiveElement
   * @brief This is the base class for the elements that make up the inactive surfaces.
   *
   * Since all inactive elements are simplified to tube shapes in the geometrical model, the geometry parameters have been
   * packed into the base class. All of these parameters apply to descendants of a different shape as well, though, because
   * they describe relations between the object and the origin, not between points within the object.
   */
  class InactiveElement : public MaterialProperties {
  public:
    /**
     * @enum InType A list of the various types of neighbour or feeder an element can have
     */
    enum InType { no_in, tracker, barrel, endcap };
    InactiveElement();
    virtual ~InactiveElement() {}
    virtual double getSurface() const;
    bool isVertical() const;
    void setVertical(bool vertical);
    bool isFinal();
    void setFinal(bool final);
    double getZOffset() const;
    void setZOffset(double zoffset);
    double getZLength() const;
    void setZLength(double zlength);
    double getInnerRadius() const;
    void setInnerRadius(double iradius);
    double getRWidth() const;
    void setRWidth(double width);
    virtual double getLength() const;
    double getVolume() const;
    int getFeederIndex();
    void setFeederIndex(int layer);
    InType getFeederType();
    void setFeederType(InType type);
    int getNeighbourIndex();
    InType getNeighbourType();
    void setNeighbourType(InType type);
    void setNeighbourIndex(int previous);
    void setTotalMass(double mass);
    void setLocalMass(double mass);
    void setRadiationLength(double rlength);
    void setInteractionLength(double ilength);
    std::pair<double, double> getEtaMinMax();
    virtual void print();
  protected:
    bool is_vertical, is_final;
    double z_offset, z_length, i_radius, w_radius;
    int feeder_index, neighbour_index;
    InType feeder_type, neighbour_type;
  };
}
#endif	/* _INACTIVEELEMENT_H */

