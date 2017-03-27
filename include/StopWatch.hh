#ifndef StopWatch_h
#define StopWatch_h

#include <MessageLogger.hh>
#include <ctime>
#include <string>
#include <list>

#define startTaskClock(message) StopWatch::instance()->startCounter(message)
#define addTaskInfo(message) StopWatch::instance()->addInfo(message)
#define stopTaskClock() StopWatch::instance()->stopCounter()

/**
 * @class StopWatch
 * @brief This class measures the CPU time used between two moments
 *
 * It uses the ctime library to call time(), which is measuring the
 * CPU time used by the process so far. It works nicely as long as the
 * counter does not flip over (which should happen every 72 minutes at
 * most in a 32 bit CPU) The stop watch is started by the constructyor
 * (and by restart()) and it gives the elapsed CPU time in seconds
 * when asked
 */
class StopWatch {
 public:
  static StopWatch* instance();
  void startCounter(std::string message);
  double stopCounter();
  void setVerbosity(unsigned int newVerbosity, bool newPerformance);
  void addInfo(std::string message);
  static void destroy();
 private:
  StopWatch();
  ~StopWatch();
  static StopWatch* myInstance_;
  std::list<clock_t> startTimes_;
  double diffClock(const clock_t& stopTime, const clock_t& startTime);
  unsigned int verbosity_;
  unsigned int lastVerbosity_;
  bool reportTime_;
};

#endif
