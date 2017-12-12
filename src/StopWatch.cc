#include <StopWatch.hh>

#include <iostream>

// Global static pointer used to ensure a single instance of the class
StopWatch* StopWatch::myInstance_ = NULL;

// Returns the instance (if already present) or creates one if needed
StopWatch* StopWatch::instance() {
  return myInstance_ ? myInstance_ : (myInstance_ = new StopWatch);
}

// Destroys the current instance
void StopWatch::destroy() {
  if (myInstance_) {
    delete myInstance_; 
    myInstance_ = NULL;
  }
}

void StopWatch::setVerbosity(unsigned int newVerbosity, bool newPerformance) {
  verbosity_ = newVerbosity;
  reportTime_ = newPerformance;
}

/* Object constructor */
StopWatch::StopWatch() {
  lastVerbosity_ = 0;
  verbosity_ = 1000;
  reportTime_ = true;
} 

/* Object destructor */
StopWatch::~StopWatch() {
  std::cout << std::endl;
} 

void StopWatch::startCounter(std::string message) {
  startTimes_.push_back(clock());
  if (startTimes_.size()<=verbosity_) {
    std::cout << std::endl;
    for (unsigned int i=1; i<startTimes_.size(); ++i) std::cout << "  ";
    std::cout << message << " ... " << std::flush;
  }
}

double StopWatch::stopCounter() {
  double timeSeconds;
  if (startTimes_.size()) {
    clock_t stopTime = clock();
    clock_t startTime = startTimes_.back();
    startTimes_.pop_back();
    timeSeconds = diffClock(stopTime, startTime);
    if (startTimes_.size()<verbosity_) {
      if (startTimes_.size()<lastVerbosity_) std::cout << std::endl;
      std::cout << "done" ;
      if (reportTime_) std::cout << " [in " << timeSeconds << " s]";
      std::cout << std::flush;
      lastVerbosity_=startTimes_.size();
    }
  } else {
    timeSeconds = 0;
    logERROR("A timer stop was requested, not corresponding to any timer start");
  }
  return timeSeconds;
}

void StopWatch::addInfo(std::string message) {
  if (startTimes_.size()<=verbosity_) {
    std::cout << message << " " << std::flush;
  }
}

/* Returns the difference between start and stop time in s */
double StopWatch::diffClock(const clock_t& stopTime, const clock_t& startTime) {
  double diffticks=stopTime-startTime;
  double diffs=(diffticks)/CLOCKS_PER_SEC;
  return diffs;
}

