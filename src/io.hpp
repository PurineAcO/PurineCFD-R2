#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include "udf.hpp"
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <nlohmann/json.hpp>
#include <set>
#include <string>
#include <system_error>
#include <vector>

// 建立日志目录并重定向stdout
bool open_log(const char* path);
// 输出流场
bool dump_field(int step);
// 读取config.json
bool config::load(const char* path);

static bool make_dirs(const std::string& path){
    if(path.empty()){
        return true;
    }
    std::error_code error;
    std::filesystem::create_directories(path,error);
    return !error;
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
                cell.phy.T,cell.otphy.p,
                std::hypot(cell.phy.u,cell.phy.v)/cell.otphy.a,cell.tur.miubl);
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

/*
Json 读取部分
*/

using Json = nlohmann::json;

static bool fail(const std::string& msg){
    fprintf(stderr,"Error: %s\n",msg.c_str());
    return false;
}

// 校验对象的键: 不在白名单内的报错, 必填键缺失也报错; optional 里的键可以缺席
static bool keys(const Json& object,std::initializer_list<const char*> allowed,
                 std::initializer_list<const char*> optional = {}){
    if(!object.is_object()){
        return fail("expected a configuration object");
    }
    const std::set<std::string> names(allowed.begin(),allowed.end());
    const std::set<std::string> extra(optional.begin(),optional.end());
    for(const auto& item : object.items()){
        if(!names.count(item.key()) && !extra.count(item.key())){
            return fail("unknown configuration key: " + item.key());
        }
    }
    for(const std::string& name : names){
        if(!object.contains(name)){
            return fail("missing configuration key: " + name);
        }
    }
    return true;
}

// 读取一个有限的数, 要求为正
static bool positive(const Json& object,const char* key,double& out){
    const auto item = object.find(key);
    if(item == object.end() || !item->is_number()){
        return fail(std::string(key) + " must be a number");
    }
    const double value = item->get<double>();
    if(!std::isfinite(value) || value <= 0.0){
        return fail(std::string(key) + " must be finite and positive");
    }
    out = value;
    return true;
}

// 读取一个有限的数, 符号不限
static bool number(const Json& object,const char* key,double& out){
    const auto item = object.find(key);
    if(item == object.end() || !item->is_number()){
        return fail(std::string(key) + " must be a number");
    }
    const double value = item->get<double>();
    if(!std::isfinite(value)){
        return fail(std::string(key) + " must be finite");
    }
    out = value;
    return true;
}

// 读取一个正整数, 按 32 位整数处理
static bool count(const Json& object,const char* key,int& out){
    double value = 0.0;
    if(!positive(object,key,value)){
        return false;
    }
    if(!object.find(key)->is_number_integer() || value > std::numeric_limits<int>::max()){
        return fail(std::string(key) + " must be a positive 32-bit integer");
    }
    out = static_cast<int>(value);
    return true;
}

// 读取路径字符串, 相对配置文件所在目录解析并做规范化
static bool path_value(const Json& object,const char* key,const std::filesystem::path& base,
                std::string& out,std::error_code& error){
    const auto item = object.find(key);
    if(item == object.end() || !item->is_string()){
        return fail(std::string(key) + " must be a path string");
    }
    const std::string text = item->get<std::string>();
    if(text.empty() || text.find('\0') != std::string::npos){
        return fail(std::string(key) + " is an invalid path");
    }
    std::filesystem::path full = std::filesystem::absolute(base/text,error);
    if(error){
        full = base/text;
    }
    out = full.lexically_normal().string();
    return true;
}

inline bool config::load(const char* path){
    // 打开配置文件
    std::ifstream stream(path);
    if(!stream){
        return fail("Cannot open configuration: " + std::string(path));
    }
    // 解析过程中逐层记录已出现的键, 用来发现重复键与过深嵌套
    std::vector<std::set<std::string>> object_keys;
    bool malformed = false;
    auto callback = [&](int depth,Json::parse_event_t event,Json& parsed){
        if(depth > 16){
            malformed = true;
        }
        if(event == Json::parse_event_t::object_start){
            object_keys.emplace_back();
        }else if(event == Json::parse_event_t::key){
            if(object_keys.empty() ||
               !object_keys.back().insert(parsed.get<std::string>()).second){
                malformed = true;
            }
        }else if(event == Json::parse_event_t::object_end && !object_keys.empty()){
            object_keys.pop_back();
        }
        return true;
    };
    const Json root = Json::parse(stream,callback,false);
    if(root.is_discarded()){
        return fail("Cannot parse configuration: " + std::string(path));
    }
    if(malformed){
        return fail("Duplicate configuration key or too deep nesting");
    }
    // 顶层只允许 io / solver / farfield 三节
    if(!keys(root,{"io","solver","farfield"})){
        return false;
    }
    const Json& io = root.at("io");
    const Json& solver = root.at("solver");
    const Json& far = root.at("farfield");
    // 各节的键白名单; io 的 structured 只在网格带结构化信息时出现
    if(!keys(io,{"mesh","log","field"},{"structured"}) ||
       !keys(solver,{"max_steps","cfl","dump_interval","convergence_interval"},{"threads"}) ||
       !keys(far,{"Ma","T","p","alpha"})){
        return false;
    }
    // 所有文件路径都相对配置文件所在目录解析
    std::error_code error;
    std::filesystem::path base = std::filesystem::path(path).parent_path();
    if(base.empty()){
        base = ".";
    }
    std::string mesh_path,log_path,field_path;
    if(!path_value(io,"mesh",base,mesh_path,error) ||
       !path_value(io,"log",base,log_path,error) ||
       !path_value(io,"field",base,field_path,error)){
        return false;
    }
    // structured 给出结构化邻接表路径并置起令牌; 缺省表示网格不带结构化信息
    if(io.contains("structured")){
        structer::ifstructer = true;
        if(!path_value(io,"structured",base,structer::adjacency,error)){
            return false;
        }
    }
    // 日志不得覆盖任何输入文件
    const std::string config_path =
        std::filesystem::absolute(path,error).lexically_normal().string();
    if(log_path == mesh_path || log_path == config_path){
        return fail("The log path must not overwrite an input file");
    }
    // 求解器参数
    int max_steps = 0,dump_interval = 0,conv_interval = 0;
    if(!count(solver,"max_steps",max_steps) ||
       !count(solver,"dump_interval",dump_interval) ||
       !count(solver,"convergence_interval",conv_interval)){
        return false;
    }
    if(conv_interval < 2){
        return fail("convergence_interval must be at least 2");
    }
    // threads 缺省表示按可用物理核心自动选择
    int threads = 0;
    if(solver.contains("threads") && !count(solver,"threads",threads)){
        return false;
    }
    // 来流条件
    double cfl = 0.0,ma = 0.0,T = 0.0,p = 0.0,alpha = 0.0;
    if(!positive(solver,"cfl",cfl) || !positive(far,"Ma",ma) ||
       !positive(far,"T",T) || !positive(far,"p",p) ||
       !number(far,"alpha",alpha)){
        return false;
    }
    if(ma >= 1.0){
        return fail("This solver requires a subsonic farfield: 0 < Ma < 1");
    }
    if(std::abs(alpha) >= 90.0){
        return fail("The farfield angle of attack must satisfy |alpha| < 90 degrees");
    }
    const double u_inf = ma*get_sonic_velocity(T);
    const double rho_inf = p/(cc::R*T);
    if(!std::isfinite(u_inf) || !std::isfinite(rho_inf) || rho_inf <= 0.0){
        return fail("Farfield values produce an invalid thermodynamic state");
    }
    // 全部校验通过后再一次性写入全局状态
    cc::meshpath = mesh_path;
    cc::testpath = log_path;
    cc::fieldpath = field_path;
    cc::max_step = max_steps;
    cc::threads = threads;
    fatime::CFL = cfl;
    dump_step = dump_interval;
    conv_step = conv_interval;
    const double rad = alpha*(std::acos(-1.0)/180.0);
    FAR_DEFINE.u = u_inf*std::cos(rad);
    FAR_DEFINE.v = u_inf*std::sin(rad);
    FAR_DEFINE.T = T;
    FAR_DEFINE.p = p;
    return true;
}
