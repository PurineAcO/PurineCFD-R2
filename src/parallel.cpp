#include "parallel.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <sched.h>
#include <set>
#include <string>
#include <unistd.h>
#include <utility>

namespace {

int physical_cores(){
    const int capacity = static_cast<int>(sysconf(_SC_NPROCESSORS_CONF));
    if(capacity <= 0){
        return 0;
    }
    cpu_set_t* mask = CPU_ALLOC(capacity);
    if(mask == nullptr){
        return 0;
    }
    const size_t bytes = CPU_ALLOC_SIZE(capacity);
    if(sched_getaffinity(0,bytes,mask) != 0){
        CPU_FREE(mask);
        return 0;
    }
    std::set<std::pair<int,int>> cores;
    for(int cpu=0;cpu<capacity;cpu++){
        if(!CPU_ISSET_S(cpu,bytes,mask)){
            continue;
        }
        const std::string base = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
        int socket = -1,core = -1;
        std::ifstream socket_file(base + "physical_package_id");
        std::ifstream core_file(base + "core_id");
        if(!(socket_file >> socket) || !(core_file >> core) || socket < 0 || core < 0){
            CPU_FREE(mask);
            return 0;
        }
        cores.emplace(socket,core);
    }
    CPU_FREE(mask);
    return static_cast<int>(cores.size());
}

}

const char* parallel::configure_threads(){
    omp_set_dynamic(0);
    if(const char* requested = getenv("OMP_NUM_THREADS"); requested && *requested){
        return "OMP_NUM_THREADS";
    }
    int cores = physical_cores();
    if(cores <= 0){
        cores = omp_get_num_procs();
    }
    omp_set_num_threads(std::min(cores,omp_get_thread_limit()));
    return "available-physical-cores";
}
