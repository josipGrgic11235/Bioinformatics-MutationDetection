#include <chrono>

#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER

class ExecutionTimer
{
#if defined(__APPLE__) && defined(__MACH__)
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> start_value;
#else
  std::chrono::time_point<std::chrono::system_clock> start_value;
#endif

public:
  void start();

  void printExecutionTime();
};

#endif