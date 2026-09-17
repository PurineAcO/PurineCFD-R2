#include "config.h"
#include "physic.h"
#include "timarch.h"
#include "udf.h"
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

namespace {

using Json = nlohmann::json;

bool fail(const std::string& msg){
    fprintf(stderr,"Error: %s\n",msg.c_str());
    return false;
}

bool keys(const Json& object,std::initializer_list<const char*> allowed){
    if(!object.is_object()){
        return fail("expected a configuration object");
    }
    const std::set<std::string> names(allowed.begin(),allowed.end());
    for(const auto& item : object.items()){
        if(!names.count(item.key())){
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

bool positive(const Json& object,const char* key,double& out){
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

bool number(const Json& object,const char* key,double& out){
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

bool count(const Json& object,const char* key,int& out){
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

bool path_value(const Json& object,const char* key,const std::filesystem::path& base,
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

}

bool config::load(const char* path){
    std::ifstream stream(path);
    if(!stream){
        return fail("Cannot open configuration: " + std::string(path));
    }
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
    if(!keys(root,{"io","solver","farfield"})){
        return false;
    }
    const Json& io = root.at("io");
    const Json& solver = root.at("solver");
    const Json& far = root.at("farfield");
    if(!keys(io,{"mesh","log","field"}) ||
       !keys(solver,{"max_steps","cfl","dump_interval","convergence_interval"}) ||
       !keys(far,{"Ma","T","p","alpha"})){
        return false;
    }
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
    const std::string config_path =
        std::filesystem::absolute(path,error).lexically_normal().string();
    if(log_path == mesh_path || log_path == config_path){
        return fail("The log path must not overwrite an input file");
    }
    int max_steps = 0,dump_interval = 0,conv_interval = 0;
    if(!count(solver,"max_steps",max_steps) ||
       !count(solver,"dump_interval",dump_interval) ||
       !count(solver,"convergence_interval",conv_interval)){
        return false;
    }
    if(conv_interval < 2){
        return fail("convergence_interval must be at least 2");
    }
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
    cc::meshpath = mesh_path;
    cc::testpath = log_path;
    cc::fieldpath = field_path;
    cc::max_step = max_steps;
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
