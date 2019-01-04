#include <chrono>

#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER

class ExecutionTimer
{
  std::chrono::time_point<std::chrono::system_clock> start_value;

public:
  void start();

  void printExecutionTime();
};

#endif