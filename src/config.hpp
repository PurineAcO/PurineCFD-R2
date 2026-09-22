#pragma once

#include <string>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <nlohmann/json.hpp>
#include <set>
#include <system_error>
#include <vector>

namespace cc {

    // 程序基本参数
    inline int cell_num = 0;        // 网格数目
    inline int face_num = 0;        // 面数目
    inline std::string meshpath;    // 网格文件位置
    inline std::string testpath;    // 日志输出路径
    inline std::string fieldpath;   // 流场输出路径
    inline long long max_step = 0;  // 时间步数

    struct ivec2{int x = 0;int y = 0;ivec2() = default;ivec2(int x_,int y_):x(x_),y(y_){}};

    struct vec2{
        double x = 0.0, y = 0.0;
        vec2() = default;
        vec2(double x_,double y_):x(x_),y(y_){}
        vec2& operator+=(const vec2& o){ x += o.x; y += o.y; return *this; }
        vec2& operator-=(const vec2& o){ x -= o.x; y -= o.y; return *this; }
        vec2& operator*=(double s){ x *= s; y *= s; return *this; }
    };

    inline vec2 operator+(vec2 a, const vec2& b){ return a += b; }
    inline vec2 operator-(vec2 a, const vec2& b){ return a -= b; }
    inline vec2 operator*(vec2 a, double s){ return a *= s; }
    inline vec2 operator*(double s, vec2 a){ return a *= s; }
    inline double dot(const vec2& a, const vec2& b){ return a.x*b.x + a.y*b.y; }

    // 物理量矩阵
    struct physics{
        double rho = 0.0, u = 0.0, v = 0.0; // 密度, x/y方向速度
        double T = 0.0, a = 0.0;            // 温度, 声速
        double p = 0.0, e = 0.0;            // 压力, 单位质量总能量
        vec2 ugrad, vgrad, Tgrad;           // 速度与温度梯度
    };

    // 湍流变量矩阵
    struct turbulence{
        double miubl = 0.0;        // SA工作变量
        double miubl_former = 0.0; // 本步RK开始时的miubl
        double miubl_next = 0.0;   // 下一RK阶段的值,算完统一写回
        double sad = 0.0;          // 到最近壁面中点的距离
        vec2 miublgrad;            // miubl梯度
    };

    // 用于JST的人工耗散
    struct dissipation{
        double Y = 0.0;     // 激波捕捉因子
        double L[4] = {};   // 伪Laplace
        double Fd[4] = {};  // JST耗散
    };

    // 支持的边界条件
    inline constexpr short INTER = 0, WALL = 1, FAR = 4;

    // 常数
    inline constexpr double gamma = 1.4;    // 气体绝热常数
    inline constexpr double R = 287.05;     // 气体常数R
    inline constexpr double Cp = 1004.675;  // 定压热容
    inline constexpr double Cv = 717.645;   // 恒容热容
    inline constexpr double Pr = 0.71;      // 普朗特数

}

// 伪时间步参数
namespace fatime {
    inline double CFL = 1.0;
}

// 求解设置
namespace config {
    inline bool load(const char* path = "config.json"); // 读取config.json
    inline int dump_step = 1;   // 场输出间隔
    inline int conv_step = 1;   // 残差检查间隔
}

// 结构化网格参数
namespace structer{
    inline bool ifstructer = false; // 结构化网格令牌
    inline std::string adjacency;   // 结构化邻接表路径
    inline int S_MAX = 0;           // 环向单元数, 由邻接表表头给出
    inline int N_MAX = 0;           // 径向单元数, 由邻接表表头给出
    inline constexpr int HALO = 3;  // HALO网格层数
}

/*
Json 读取部分
*/

// 这两个模块依赖 cc 的基础类型, 必须在类型定义之后包含
#include "physic.hpp"
#include "udf.hpp"

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
       !keys(solver,{"max_steps","cfl","dump_interval","convergence_interval"}) ||
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
