// Created by Leon

#include <execution_timer.hpp>
#include <iostream>

void ExecutionTimer::start()
{
    start_value = std::chrono::high_resolution_clock::now();
}

void ExecutionTimer::printExecutionTime()
{
    auto finish = std::chrono::high_resolution_clock::now();
    std::cout << "Execution time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start_value).count() / 1e3
              << "s";
}
