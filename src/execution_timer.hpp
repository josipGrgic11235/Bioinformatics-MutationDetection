#include <chrono>

class ExecutionTimer
{
    std::chrono::time_point<std::chrono::system_clock> start_value;

  public:
    void start();

    void printExecutionTime();
};