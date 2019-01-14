// Created by Leon

#include <chrono>

#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER

/**
 * Class used to track execution time in milliseconds. 
 **/
class ExecutionTimer
{
#if defined(__APPLE__) && defined(__MACH__)
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> start_value;
#else
  std::chrono::time_point<std::chrono::system_clock> start_value;
#endif

public:
  /**
   * Start execution timer.
   **/
  void start();

  /**
   * Print current execution time.
   **/
  void printExecutionTime();
};

#endif