#include "io.h"
#include "classconfig.h"
#include "config.h"
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <string>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

bool is_separator(char c){
#ifdef _WIN32
    return c == '/' || c == '\\';
#else
    return c == '/';
#endif
}

bool make_dir(const std::string& path){
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(),0755) == 0 || errno == EEXIST;
#endif
}

bool make_dirs(const std::string& path){
    if(path.empty()){
        return true;
    }
    size_t start = 1;
#ifdef _WIN32
    if(path.size() >= 2 && path[1] == ':'){
        start = (path.size() >= 3 && is_separator(path[2])) ? 3 : 2;
    }
#endif
    for(size_t i=start;i<path.size();i++){
        if(!is_separator(path[i])){
            continue;
        }
        if(!make_dir(path.substr(0,i))){
            return false;
        }
    }
    return make_dir(path);
}

}

bool open_log(const char* path){
    std::string text(path);
#ifdef _WIN32
    const size_t slash = text.find_last_of("/\\");
#else
    const size_t slash = text.find_last_of('/');
#endif
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

bool dump_field(int step){
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
