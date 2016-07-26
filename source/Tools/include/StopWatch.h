#ifndef StopWatch_h
#define StopWatch_h

#include <MessageLogger.h>
#include <ctime>
#include <string>
#include <list>

#define startTaskClock(message) StopWatch::getInstance().startCounter(message)
#define addTaskInfo(message) StopWatch::getInstance().addInfo(message)
#define stopTaskClock() StopWatch::getInstance().stopCounter()

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

  //! StopWatch access method -> get instance of singleton class StopWatch
  static StopWatch& getInstance();

  //! Destructor
  ~StopWatch();

  //! Start the counter now with a given message
  void startCounter(std::string message);

  //! Stop the counter
  double stopCounter();

  void setVerbosity(unsigned int newVerbosity, bool newPerformance);

  void addInfo(std::string message);

 private:

  //! Singleton constructor method -> must be private (singleton pattern)
  StopWatch();

  //! Returns the difference between start and stop time in s
  double diffClock(const clock_t& stopTime, const clock_t& startTime);

  std::list<clock_t> m_startTimes;
  unsigned int       m_verbosity;
  unsigned int       m_lastVerbosity;
  bool               m_reportTime;
};

#endif
