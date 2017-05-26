/**
 * @file BFieldMap.h
 * @author Zbynek Drasal
 * @date 03/dec/2015
 */

#ifndef BFIELDMAP_H_
#define BFIELDMAP_H_

#include <string>
#include <fstream>
#include <utility>
#include <vector>
#include <cmath>

class TCanvas;

/**
 * @class BFieldMap
 * @brief This class represent a single B field (Bx, By, Bz) map defined in space by (X,Y,Z) coordinates.
 * @details The map reads its configuration from the header (reading a new line starting by hash mark).
 * The map is read-in either in a format of histogram, where #bins = (Max-Min)/binWidth, values are
 * then centered with respect to the bin center. Or the map is read-in in a format of grid (mesh), where
 * #bins = (Max-Min)/binWidth + 1, values are then directly corresponding the input.
 * Get the B field vector at given position using calculateBField method. B field is returned as a
 * std::vector<double> ([0]=X, [1]=Y, [2]=Z).
 */
class BFieldMap {
public:
  /**
   * Constructor with feeding
   * @param bfieldMapFile is the path of the new file to feed
   */
  BFieldMap(std::string bfieldMapFile);

  /**
   * Populate the map attributes reading the passed file
   * @param bfieldMapFile is the path of the raw file for the bfield map to be read
   */
  void ingest(std::string bfieldMapFile);

  /**
   * Test if a point is inside the area covered by the map
   * @param (x,y,z) doubles that indicate a point in the space
   * @return True if the point is inside the map region, false otherwise
   */
  bool isInRegionXYZ(double xPos, double yPos, double zPos) const;

  /**
   * Test if a point in terms of bins is inside the number of bins covered by the map
   * @param (x,y,z) ints that indicate a point in the space
   * @return True if the point is inside the map region, false otherwise
   */
  bool isInBinRegionXYZ(int& xPos, int& yPos, int& zPos) const;

  /**
   * Get the bfield corresponding to the point (X,Y,Z)
   * @param (x,y,z) doubles that indicate a point in the space
   * @return the value of the bfield of the point, or 0 if the point is outside the map
   */
  std::vector<double> calculateBField(double xPos, double yPos, double zPos) const;

  /**
   * Draw XZ projection of B field
   * @param TCanvas - new histogram will be created and drawn to an existing canvas
   * @return - if canvas doesn't exist False will be returned
   */
  bool drawXZBFieldProj(TCanvas& xzCanvas, std::string name, double minX, double maxX, double minZ, double maxZ);

  /**
   * Draw YZ projection of B field
   * @param TCanvas - new histogram will be created and drawn to an existing canvas
   * @return - if canvas doesn't exist False will be returned
   */
  bool drawYZBFieldProj(TCanvas& yzCanvas, std::string name, double minY, double maxY, double minZ, double maxZ);

  inline double getXMin()       const {return m_xMin;};
  inline double getXMax()       const {return m_xMax;};
  inline double getXBinWidth()  const {return m_xBinWidth;};
  inline int    getXNBins()     const {return m_xBinNum;};
  inline double getYMin()       const {return m_yMin;};
  inline double getYMax()       const {return m_yMax;};
  inline double getYBinWidth()  const {return m_yBinWidth;};
  inline int    getYNBins()     const {return m_yBinNum;};
  inline double getZMin()       const {return m_zMin;};
  inline double getZMax()       const {return m_zMax;};
  inline double getZBinWidth()  const {return m_zBinWidth;};
  inline int    getZNBins()     const {return m_zBinNum;};

  inline std::string getXUnit() const {return m_xUnit;};
  inline std::string getYUnit() const {return m_yUnit;};
  inline std::string getZUnit() const {return m_zUnit;};
  inline std::string getBUnit() const {return m_dataUnit;};

  inline bool isOfTypeHistogram() const {return m_typeHist;};
  inline bool isOfTypeMesh()      const {return m_typeMesh;};
  inline bool isOK()              const {return m_bFieldOK;}

private:
  const std::string c_fileDataUnit    = "# Data unit: ";          /**< Header prefix - value of flux unit*/
  const std::string c_fileXUnit       = "# X unit: ";             /**< Header prefix - value of x unit*/
  const std::string c_fileXMin        = "# X min: ";              /**< Header prefix - value of min x*/
  const std::string c_fileXMax        = "# X max: ";              /**< Header prefix - value of max x*/
  const std::string c_fileXBinWidth   = "# X bin width: ";        /**< Header prefix - value of bin width in x*/
  const std::string c_fileXBinNum     = "# X number of bins: ";   /**< Header prefix - value of the number of bins in x*/
  const std::string c_fileYUnit       = "# Y unit: ";             /**< Header prefix - value of y unit*/
  const std::string c_fileYMin        = "# Y min: ";              /**< Header prefix - value of min y*/
  const std::string c_fileYMax        = "# Y max: ";              /**< Header prefix - value of max y*/
  const std::string c_fileYBinWidth   = "# Y bin width: ";        /**< Header prefix - value of bin width in y*/
  const std::string c_fileYBinNum     = "# Y number of bins: ";   /**< Header prefix - value of the number of bins in y*/
  const std::string c_fileZUnit       = "# Z unit: ";             /**< Header prefix - value of Z unit*/
  const std::string c_fileZMin        = "# Z min: ";              /**< Header prefix - value of min Z*/
  const std::string c_fileZMax        = "# Z max: ";              /**< Header prefix - value of max Z*/
  const std::string c_fileZBinWidth   = "# Z bin width: ";        /**< Header prefix - value of bin width in Z*/
  const std::string c_fileZBinNum     = "# Z number of bins: ";   /**< Header prefix - value of the number of bins in Z*/
  const std::string c_fileEscValue    = "/t";                     /**< File - special character used between values*/
  const std::string c_fileEscLine     = "/n";                     /**< File - special character used at the end of line*/
  const std::string c_fileEscComment  = "#//;";                   /**< File - special character used for commet*/

  const double      c_arrowMaxLength  = 0.9;
  const double      c_arrowMinLength  = 0.4;
  const double      c_arrowMaxSize    = 0.012;
  const double      c_arrowMinSize    = 0.006;

  std::string m_fileName;    /**< File name with bfield data and predefined header*/
  std::string m_dataUnit;    /**< B field unit*/
  double      m_dataFactor;  /**< Corresponding unit factor*/
  std::string m_xUnit;       /**< x unit*/
  double      m_xMin;        /**< Minimum x range*/
  double      m_xMax;        /**< Maximum x range*/
  double      m_xBinWidth;   /**< Defined bin width in x*/
  long int    m_xBinNum;     /**< Defined number of bins in x*/
  std::string m_yUnit;       /**< y unit*/
  double      m_yMin;        /**< Minimum y range*/
  double      m_yMax;        /**< Maximum y range*/
  double      m_yBinWidth;   /**< Defined bin width in y*/
  long int    m_yBinNum;     /**< Defined number of bins in y*/
  std::string m_zUnit;       /**< Z unit*/
  double      m_zMin;        /**< Minimum Z range*/
  double      m_zMax;        /**< Maximum Z range*/
  double      m_zBinWidth;   /**< Defined bin width in Z*/
  long int    m_zBinNum;     /**< Defined number of bins in Z*/
  bool        m_typeMesh;    /**< B field map defined as a grid (mesh) - (xyzMax-xyzMin)/xyzBinWidth = xyzNumBins-1*/
  bool        m_typeHist;    /**< B field map defined as a histogram - (xyzMax-xyzMin)/xyzBinWidth = xyzNumBins*/

  bool        m_bFieldOK;    /**< Read B field correctly?*/
  std::vector<std::vector<std::vector< std::vector<double> >>> m_bField;       /**< The matrix X,Y,Z. Matrix contains the bfield vector for each bin of the map*/
};


#endif /* BFIELDMAP_H_ */
