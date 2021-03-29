//
// File:   tk2CMSSW.h
// Author: ndemaio
//
// Created on August 28, 2009, 5:23 PM
//

/**
 * @file tk2CMSSW.h
 * @brief This is the header file for the CMSSW XML translator class
 */

#ifndef _TK2CMSSW_H
#define _TK2CMSSW_H

#include <Extractor.hh>
#include <MainConfigHandler.hh>
#include <MaterialBudget.hh>
#include <MaterialTable.hh>
#include <XMLWriter.hh>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <iostream>
#include <pwd.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <time.h>
#include <tk2CMSSW_datatypes.hh>
#include <tk2CMSSW_strings.hh>
#include <unistd.h>

/**
 * A shorter alias for the filesystem library namespace
 */
namespace bfs = boost::filesystem;
namespace insur {
/**
 * @class tk2CMSSW
 * @brief This class is the main translator interface for generating XML output
 * for CMSSW from an existing material budget and table.
 *
 * It deals directly with setting up output paths and buffers for output files,
 * while it delegates analysis and XML formatting of the material budget to an
 * internal instance of <i>Extractor</i> and <i>XMLWriter</i>, respectively.
 * Existing files are not overwritten until analysis and XML formatting have
 * been completed successfully. Should any of them fail, any existing older
 * files are restored; preliminary results for the new files are discarded.
 */
class tk2CMSSW {
  mainConfigHandler &mainConfiguration;

public:
  tk2CMSSW(mainConfigHandler &mch) : mainConfiguration(mch) {}
  virtual ~tk2CMSSW() {}
  void translate(MaterialTable &mt, MaterialBudget &mb, XmlTags &trackerXmlTags,
                 std::string xmlDirectoryPath, std::string xmlOutputPath,
                 std::string xmlOutputName = "", bool wt = false);
  struct ConfigFile {
    std::string name, content;
  };
  void addConfigFile(const ConfigFile &file) { configFiles_.push_back(file); }

protected:
  CMSSWBundle trackerData;
  CMSSWBundle otstData;
  Extractor ex;
  XMLWriter wr;

private:
  std::vector<ConfigFile> configFiles_;
  void print(CMSSWBundle &data);
  void writeSimpleHeader(std::ostream &os, std::string &metadataFileName);
  void writeMetadata(std::ofstream &out);
  std::string currentDateTime(bool withTime) const;
  std::string fullUserName() const;
  const std::vector<ConfigFile> &getConfigFiles() const { return configFiles_; }
};
} // namespace insur
#endif /* _TK2CMSSW_H */
