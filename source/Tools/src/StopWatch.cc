#include <StopWatch.h>

#include <iostream>

// Returns the instance (if already present) or creates one if needed
StopWatch& StopWatch::getInstance() {

  static StopWatch s_instance;
  return s_instance;
}

void StopWatch::setVerbosity(unsigned int newVerbosity, bool newPerformance) {

  m_verbosity  = newVerbosity;
  m_reportTime = newPerformance;
}

//
// Singleton constructor method -> must be private (singleton pattern)
//
StopWatch::StopWatch() {
  m_lastVerbosity = 0;
  m_verbosity     = 1000;
  m_reportTime    = true;
} 

//
// Destructor
//
StopWatch::~StopWatch() {
  std::cout << std::endl;
} 

void StopWatch::startCounter(std::string message) {

  m_startTimes.push_back(clock());
  if (m_startTimes.size()<=m_verbosity) {
    std::cout << std::endl;
    for (unsigned int i=1; i<m_startTimes.size(); ++i) std::cout << "  ";
    std::cout << message << " ... " << std::flush;
  }
}

double StopWatch::stopCounter() {
  double timeSeconds;
  if (m_startTimes.size()) {
    clock_t stopTime = clock();
    clock_t startTime = m_startTimes.back();
    m_startTimes.pop_back();
    timeSeconds = diffClock(stopTime, startTime);
    if (m_startTimes.size()<m_verbosity) {
      if (m_startTimes.size()<m_lastVerbosity) std::cout << std::endl;
      std::cout << "done" ;
      if (m_reportTime) std::cout << " [in " << timeSeconds << " s]";
      std::cout << std::flush;
      m_lastVerbosity=m_startTimes.size();
    }
  } else {
    timeSeconds = 0;
    logERROR("A timer stop was requested, not corresponding to any timer start");
  }
  return timeSeconds;
}

void StopWatch::addInfo(std::string message) {
  if (m_startTimes.size()<=m_verbosity) {
    std::cout << message << " " << std::flush;
  }
}

/* Returns the difference between start and stop time in s */
double StopWatch::diffClock(const clock_t& stopTime, const clock_t& startTime) {
  double diffticks=stopTime-startTime;
  double diffs=(diffticks)/CLOCKS_PER_SEC;
  return diffs;
}

