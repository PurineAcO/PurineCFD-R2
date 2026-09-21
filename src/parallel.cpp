#include "parallel.h"
#include <algorithm>
#include <cstdlib>
#include <omp.h>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fstream>
#include <sched.h>
#include <set>
#include <string>
#include <unistd.h>
#include <utility>
#endif

namespace {

#ifdef _WIN32

// 每个物理核心在 RelationProcessorCore 记录里出现一次, 与进程亲和掩码求交即可
// 排除已绑定外部的逻辑核, 语义与 Linux 下按 (socket, core) 去重一致。
int physical_cores(){
    DWORD_PTR process_mask = 0;
    DWORD_PTR system_mask = 0;
    if(!GetProcessAffinityMask(GetCurrentProcess(),&process_mask,&system_mask)){
        return 0;
    }
    DWORD bytes = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore,nullptr,&bytes);
    if(bytes == 0){
        return 0;
    }
    std::vector<char> buffer(bytes);
    if(!GetLogicalProcessorInformationEx(
           RelationProcessorCore,
           reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data()),
           &bytes)){
        return 0;
    }
    int cores = 0;
    for(DWORD offset=0;offset<bytes;){
        const auto* info = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(
            buffer.data() + offset);
        if(info->Relationship == RelationProcessorCore &&
           (info->Processor.GroupMask[0].Mask & process_mask) != 0){
            cores++;
        }
        offset += info->Size;
    }
    return cores;
}

#else

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

#endif

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
