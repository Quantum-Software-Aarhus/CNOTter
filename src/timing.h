#include <chrono> // required for time stamps in logging
#include <iomanip> // required for setprecision

// here the program's starting time is measured
using namespace std::chrono;
const system_clock::time_point startTime = system_clock::now();
auto lifeTime = std::vector<system_clock::time_point>(omp_get_max_threads(),startTime);

uint64_t passedTime(system_clock::time_point stopwatch) {
    return duration_cast<seconds>(system_clock::now() - stopwatch).count();
}

uint64_t currentTime() {
    return duration_cast<seconds>(system_clock::now() - startTime).count();
}

void report(uint64_t level, uint64_t orbit) {
    std::cout   << std::setprecision(std::numeric_limits<double>::digits10)
                << "(" << currentTime() << "s) ("
                << level << " elts) (" << orbit << " orbits)" << std::endl;
}

void lifeBeat(int worker, uint64_t level, uint64_t orbit) {
    std::cout   << std::setprecision(std::numeric_limits<double>::digits10)
                << "...Worker " << worker
                << " (" << currentTime() << "s) ("
                << level << " elts) (" << orbit << " orbits)" << std::endl;
}
