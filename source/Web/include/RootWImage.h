/*
 * RootWImage.h
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWIMAGE_H_
#define INCLUDE_ROOTWIMAGE_H_

#include <string>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

#include "RootWItem.h"

// Forward declaration
class TCanvas;

// Typedefs
typedef std::string RootWImageSize;

/*
 * @class RootWImage
 * An html container for image - image input expected in ROOT TCanvas format. On Output the sw
 * provides one or more image(s) and their description. Only images supported by ROOT can be printed out,
 * for details see variable c_defaultAllowedExtensions.
 */
class RootWImage : public RootWItem {

 public:

  //! Default constructor - canvas is cloned and dynamically allocated
  RootWImage(const TCanvas& myCanvas);

  //! Constructor - clone canvas and define image width & height
  RootWImage(const TCanvas& myCanvas, int witdh, int height);

  //! Constructor - clone canvas and define image width & height, define relativeHtml directory wrt target directory
  RootWImage(const TCanvas& myCanvas, int witdh, int height, std::string relativeHtmlDirectory); // TODO: is this used for real?

  //! Default destructor
  ~RootWImage();

  // Setter methods
  void setComment(std::string newComment)                 {m_comment     = newComment;}
  void setName(std::string newName)                       {m_name        = newName; }
  void setZoomedSize(int witdh, int height)               {m_zoomedWidth = witdh; m_zoomedHeight = height;}
  void setRelativeHtmlDirectory(std::string newDirectory) {m_relativeHtmlDirectory = newDirectory;}
  void setTargetDirectory(std::string newDirectory)       {m_targetDirectory = newDirectory;}

  //! Image will be saved in default formats (setDefaultExtensions) + all formats defined by this method (needs
  //! to be a supported format, cf. c_defaultAllowedExtensions)
  bool addExtension(std::string newExt);

  // Getter methods
  std::string getName()      const {return m_name;}

  //! Return reference to the saved TCanvas to get access to plots etc. saved in the canvas
  const TCanvas& getCanvas() const;

  //! Method preceding dump to save files on disk and prepare images
  std::string saveFiles(int smallWidth, int smallHeight);
  std::string saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight);

  //! Dump method - printing image in the html format with description etc. (first call of saveFiles is needed)
  virtual std::ostream& dump(std::ostream& output);

 private:

  //! Private helper method cloning the given canvas and setting its properties
  void setCanvas(const TCanvas& myCanvas);

  RootWImageSize makeSizeCode(int sw, int sh, int lw, int lh);

  //! Private method - set default supported format: png (automatically) + those defined here
  void setDefaultExtensions();

  std::unique_ptr<TCanvas> m_canvas; //!< Image container internally implemented as ROOT TCanvas

  int            m_zoomedWidth;
  int            m_zoomedHeight;
  std::string    m_relativeHtmlDirectory;
  std::string    m_targetDirectory;
  RootWImageSize m_lastSize;
  std::string    m_comment;
  std::string    m_name;
  std::string    m_allowedExtensions; // Will be initialized in the constructor

  std::map<RootWImageSize, std::string> m_text;
  std::map<RootWImageSize, bool>        m_fileSaved;
  std::vector<std::string>              m_fileTypeV;

  static int                        s_imageCounter;
  static std::map<std::string, int> s_imageNameCounter;

  const double c_thumbCompression       = 2.;
  const char*  c_defaultAllowedExtensions = "|C|png|gif|svg|root|eps|pdf|ps|";
};



#endif /* INCLUDE_ROOTWIMAGE_H_ */
