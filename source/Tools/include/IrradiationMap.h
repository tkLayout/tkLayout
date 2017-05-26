/**
 * @file IrradiationMap.h
 * @author Stefano Martina, Zbynek Drasal
 * @date 18/feb/2014
 */

#ifndef IRRADIATIONMAP_H_
#define IRRADIATIONMAP_H_

#include <string>
#include <fstream>
#include <utility>
#include <vector>

/**
 * @class IrradiationMap
 * @brief This class represent a single irradiation map defined either on (R,Z) or (Z,R) plane.
 * @details The map reads its configuration from the header (reading a new line starting by hash mark).
 * The maps are sortable with respect to the resolution using an operator <.
 * The map is either oriented in (R,Z) or (Z,R) plane - use calculateIrradiationRZ to access irradiation
 * value in the first case, calculateIrradiationZR in the latter case.
 * The map is read-in either in a format of histogram, where #bins = (Max-Min)/binWidth, values are
 * then centered with respect to bin center.
 * Or the map is read-in in a format of grid (mesh), where #bins = (Max-Min)/binWidth + 1, values are
 * then directly corresponding the input.
 */
class IrradiationMap {
public:
  /**
   * Constructor with feeding
   * @param irradiationMapFile is the path of the new file to feed
   */
  IrradiationMap(std::string irradiationMapFile);

  /**
   * Constructor without feeding
   */
  IrradiationMap();

  /**
   * Populate the map attributes reading the passed file
   * @param irradiationMapFile is the path of the raw file for the irradiation map to be read
   */
  void ingest(std::string irradiationMapFile);

  /**
   * Get the area of a bin of the map, identifies the resolution of the map
   * @return The area of the bin
   */
  double binArea() const;

  /**
   * Overriding of the operator < for allowing the sort of a vector of maps by their resolutions
   */
  bool operator < (const IrradiationMap& confrontedMap) const;

  /**
   * Test if a point is inside the area covered by the map
   * @param (z,rho) doubles that indicate a point in the plane ZxRho
   * @return True if the point is inside the map region, false otherwise
   */
  bool isInRegionZR(double zPos, double rPos) const;
  bool isInRegionZR(std::pair<double,double> coordinates) const {return isInRegionZR(coordinates.first , coordinates.second);};
  bool isInRegionRZ(std::pair<double,double> coordinates) const {return isInRegionZR(coordinates.second, coordinates.first);};

  /**
   * Test if a point in terms of bins is inside the number of bins covered by the map
   * @param (z,rho) ints that indicate a point in the plane ZxRho
   * @return True if the point is inside the map region, false otherwise
   */
  bool isInBinRegionZR(int& zPos, int& rPos) const;

  /**
   * Get the irradiation of the point if the map was saved in the plane ZxRho
   * @param (z,rho) doubles that indicate a point in the plane ZxRho
   * @return the value of the irradiation of the point, or 0 if the point is outside the map
   */
  double calculateIrradiationZR(double zPos, double rPos) const;

  /**
   * Get the irradiation of the point if the map was saved in the plane RhoxZ
   * @param (rho,z) doubles that indicate a point in the plane RhoxZ
   * @return the value of the irradiation of the point, or 0 if the point is outside the map
   */
  double calculateIrradiationRZ(double rPos, double zPos) const;
  double calculateIrradiationRZ(std::pair<double,double> coordinates) const {return calculateIrradiationRZ(coordinates.first, coordinates.second);};

  inline double getRMin()      const {return m_rMin;};
  inline double getRMax()      const {return m_rMax;};
  inline double getRBinWidth() const {return m_rBinWidth;};
  inline int    getRNBins()    const {return m_rBinNum;};
  inline double getZMin()      const {return m_zMin;};
  inline double getZMax()      const {return m_zMax;};
  inline double getZBinWidth() const {return m_zBinWidth;};
  inline int    getZNBins()    const {return m_zBinNum;};

  inline std::string getRUnit()    const {return m_rUnit;};
  inline std::string getZUnit()    const {return m_zUnit;};
  inline std::string getFluxUnit() const {return m_dataUnit;};

  inline bool isOfTypeHistogram() const {return m_typeHist;};
  inline bool isOfTypeMesh()      const {return m_typeMesh;};

private:
  const std::string comp_dataUnit    = "# Data unit: ";          /**< Prefix of the line of the header of the feeded file that precedes the value of flux unit*/
  const std::string comp_rUnit       = "# R unit: ";             /**< Prefix of the line of the header of the feeded file that precedes the value of rho unit*/
  const std::string comp_rMin        = "# R min: ";              /**< Prefix of the line of the header of the feeded file that precedes the value of min rho*/
  const std::string comp_rMax        = "# R max: ";              /**< Prefix of the line of the header of the feeded file that precedes the value of max rho*/
  const std::string comp_rBinWidth   = "# R bin width: ";        /**< Prefix of the line of the header of the feeded file that precedes the value of bin width in rho*/
  const std::string comp_rBinNum     = "# R number of bins: ";   /**< Prefix of the line of the header of the feeded file that precedes the value of the number of bins in rho*/
  const std::string comp_zUnit       = "# Z unit: ";             /**< Prefix of the line of the header of the feeded file that precedes the value of Z unit*/
  const std::string comp_zMin        = "# Z min: ";              /**< Prefix of the line of the header of the feeded file that precedes the value of min Z*/
  const std::string comp_zMax        = "# Z max: ";              /**< Prefix of the line of the header of the feeded file that precedes the value of max Z*/
  const std::string comp_zBinWidth   = "# Z bin width: ";        /**< Prefix of the line of the header of the feeded file that precedes the value of bin width in Z*/
  const std::string comp_zBinNum     = "# Z number of bins: ";   /**< Prefix of the line of the header of the feeded file that precedes the value of the number of bins in Z*/
  const std::string comp_norm        = "# normalization value: ";/**< Prefix of the line of the header of the feeded file that precedes the value of the normalization value in fb^-1*/
  const std::string comp_EscValue    = "/t";
  const std::string comp_EscLine     = "/n";
  const std::string comp_EscComment  = "#//;";

  std::string m_fileName;    /**< File name with irradiation data and predefined header*/
  std::string m_dataUnit;    /**< Flux unit*/
  double      m_dataFactor;  /**< Corresponding unit factor*/
  std::string m_rUnit;       /**< Rho unit*/
  double      m_rMin;        /**< Minimum radius range*/
  double      m_rMax;        /**< Maximum radius range*/
  double      m_rBinWidth;   /**< Defined bin width in radii*/
  long int    m_rBinNum;     /**< Defined number of bins in radii*/
  std::string m_zUnit;       /**< Z unit*/
  double      m_zMin;        /**< Minimum Z range*/
  double      m_zMax;        /**< Maximum Z range*/
  double      m_zBinWidth;   /**< Defined bin width in Z*/
  long int    m_zBinNum;     /**< Defined number of bins in Z*/
  double      m_norm;        /**< The value of normalization (per pp collision or per fb^-1?)*/
  bool        m_typeMesh;    /**< Irradiation map defined as a grid (mesh) - (zMax-zMin)/zBinWidth = zNumBins-1*/
  bool        m_typeHist;    /**< Irradiation map defined as a histogram - (zMax-zMin)/zBinWidth = zNumBins*/

  std::vector< std::vector<double> > m_irradiation;       /**< The matrix (rho * Z) or (Z*rho) -> see different get techniques based on the format. Matrix contains the irradiation values for each bin of the map*/
};


#endif /* IRRADIATIONMAP_H_ */
