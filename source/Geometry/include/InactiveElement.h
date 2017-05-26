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
#ifndef INCLUDE_INACTIVEELEMENT_H_
#define INCLUDE_INACTIVEELEMENT_H_

#include "MaterialProperties.h"
#include <Math/Vector3Dfwd.h>

using ROOT::Math::XYZVector;

/**
* @class InactiveElement
* @brief This is the base class for the elements that make up the inactive surfaces (tubes & discs).
*
* Since all inactive elements are simplified to tube shapes (tubes or discs) in the geometrical model, the geometry parameters
* have been packed into this base class. All of these parameters apply to descendants of a different shape as well, though,
* they describe relations between the object and the origin, not between points within the object.
*/
class InactiveElement : public MaterialProperties {

 public:
    /**
     * @enum InType A list of the various types of neighbour or feeder an element can have
     */
    enum InType { no_in, tracker, barrel, endcap };

    //! Constructor: define starting z-pos, length in Z, inner radius, radial width & whether it's vertical (ring) or horizontal (disc)
    InactiveElement(double zOffset, double zLength, double rInner, double rWidth, bool isVertical);
    InactiveElement(InactiveElement& element);
    virtual ~InactiveElement() {}

    //! Check if track hit the inactive element -> if yes, return true with passed material & hit position vector
    bool checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, XYZVector& hitPos) const;

    // Getter methods
    virtual double getSurface() const { return m_surface; }
    virtual double getLength() const  { return m_isVertical ? m_rWidth : m_zLength; };
    virtual void   print();

    bool   isVertical() const         { return m_isVertical;          };
    bool   isHorizontal() const       { return !m_isVertical;         };
    double getZOffset() const         { return m_zOffset;             };
    double getZCentre() const         { return m_zOffset+m_zLength/2; };
    double getZLength() const         { return m_zLength;             };
    double getInnerRadius() const     { return m_rInner;              };
    double getRWidth() const          { return m_rWidth;              };
    double getOuterRadius() const     { return m_rOuter;              };
    double getVolume() const          { return m_volume;              };
    int    getFeederIndex() const     { return m_feederIndex;         };
    InType getFeederType() const      { return m_feederType;          };
    int    getNeighbourIndex() const  { return m_neighbourIndex;      };
    InType getNeighbourType() const   { return m_neighbourType;       };

    std::pair<double, double> getEtaMinMax() const;

    // Setter methods
    void setZOffset(double zOffset)          { m_zOffset       = zOffset; }; // Reset z offset -> no need to recalculate surface & volume (no effect)
    void setFeederIndex(int layer)           { m_feederIndex   = layer;   };
    void setFeederType(InType type)          { m_feederType    = type;    };
    void setNeighbourIndex(int previous)     { m_neighbourIndex= previous;};
    void setNeighbourType(InType type)       { m_neighbourType = type;    };

    void setTotalMass(double mass)           { total_mass = mass;        };
    void setLocalMass(double mass)           { local_mass = mass;        };
    void setExitingMass(double mass)         { exiting_mass = mass;      };
    void setRadiationLength(double rlength)  { r_length = rlength;       };
    void setInteractionLength(double ilength){ i_length = ilength;       };

  protected:

    bool   m_isVertical;//!< True if vertical, false if horizontal
    double m_zOffset;   //!< Starting Z position (end position = zOffset + zLength)
    double m_zLength;   //!< Length in Z direction
    double m_rInner;    //!< Inner radius
    double m_rWidth;    //!< Width in radial direction, rOuter = rInner+rWidth
    double m_rOuter;    //!< Outer radius
    double m_surface;   //!< Element avg. surface (dependent on whether vertical or horizontal object is defined)
    double m_volume;    //!< Element volume

    int    m_feederIndex, m_neighbourIndex;
    InType m_feederType , m_neighbourType;

}; // Class

#endif	/* INCLUDE_INACTIVEELEMENT_H_ */

