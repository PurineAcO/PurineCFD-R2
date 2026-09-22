#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <string>
#include <sys/stat.h>

// 建立日志目录并重定向stdout
bool open_log(const char* path);
// 输出流场
bool dump_field(int step);

namespace {

inline bool make_dirs(const std::string& path){
    if(path.empty()){
        return true;
    }
    for(size_t i=1;i<path.size();i++){
        if(path[i] != '/'){
            continue;
        }
        std::string sub = path.substr(0,i);
        if(mkdir(sub.c_str(),0755) != 0 && errno != EEXIST){
            return false;
        }
    }
    return mkdir(path.c_str(),0755) == 0 || errno == EEXIST;
}

}

inline bool open_log(const char* path){
    std::string text(path);
    size_t slash = text.find_last_of('/');
    if(slash != std::string::npos && slash != 0 && !make_dirs(text.substr(0,slash))){
        fprintf(stderr,"Error: cannot create log directory for %s\n",path);
        return false;
    }
    if(!freopen(path,"w",stdout)){
        fprintf(stderr,"Error: cannot open log: %s\n",path);
        return false;
    }
    return true;
}

inline bool dump_field(int step){
    if(!make_dirs(cc::fieldpath)){
        fprintf(stderr,"Error: cannot create field directory: %s\n",cc::fieldpath.c_str());
        return false;
    }
    char tag[32];
    snprintf(tag,sizeof(tag),"step_%06d.dat",step);
    const std::string name = cc::fieldpath + "/" + tag;
    FILE* fp = fopen(name.c_str(),"w");
    if(!fp){
        fprintf(stderr,"Error: cannot create field file: %s\n",name.c_str());
        return false;
    }
    fprintf(fp,"TITLE=\"step %d\"\n",step);
    fprintf(fp,"VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\",\"nu_tilde\"\n");
    for(const cc::cell_class& cell : cc::CellList){
        fprintf(fp,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                cell.center.x,cell.center.y,cell.phy.rho,cell.phy.u,cell.phy.v,
                cell.phy.T,cell.phy.p,
                std::hypot(cell.phy.u,cell.phy.v)/cell.phy.a,cell.tur.miubl);
    }
    if(ferror(fp)){
        fclose(fp);
        fprintf(stderr,"Error: failed to write field file: %s\n",name.c_str());
        return false;
    }
    if(fclose(fp) != 0){
        fprintf(stderr,"Error: failed to close field file: %s\n",name.c_str());
        return false;
    }
    printf("[场输出] step %d  %s\n",step,name.c_str());
    return true;
}
