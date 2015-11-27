/*
 * AnalyzerModule.h
 *
 * Base analyzis and visualization class -> inherit from this class whenever you want
 * to analyze data and visualize them in a user defined way. One has to implement the
 * content of two purly virtual methods: analyze() and visualize().
 *
 *  Created on: 27. 11. 2015
 *      Author: drasal
 */

#ifndef INCLUDE_ANALYZERMODULE_H_
#define INCLUDE_ANALYZERMODULE_H_

#include <vector>

class RootWSite;
class Tracker;
class SimParms;

class AnalyzerModule
{
 public:
  AnalyzerModule(std::vector<Tracker*> trackers);
  virtual ~AnalyzerModule();

  // Pure virtual analysis method -> analyzes data
  virtual bool analyze() = 0;
  // Pure virtual visualization method -> visualizes output
  virtual bool visualize(RootWSite& webSite, const SimParms* simParms) = 0;

 protected:
  // Vector of trackers -> const pointer, one can't and shouldn't change its content
  std::vector<const Tracker*> m_trackers;

};

#endif /* INCLUDE_ANALYZERMODULE_H_ */
