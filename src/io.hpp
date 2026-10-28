#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
// 从流场文件读回初场, 并回传文件里记录的步号
bool load_field(const char* path, int& step);
// 读取config.json
bool config::load(const char* path);
// 写出壁面压力系数
void write_wall_cp(FILE* fp, int step, double t, cc::vecp refp, double refpressure);

static bool make_dirs(const std::string& path) {
    if (path.empty()) { return true; }
    std::error_code error;
    std::filesystem::create_directories(path, error);
    return !error;
}

inline bool open_log(const char* path) {
    std::string text(path);
    size_t slash = text.find_last_of('/');
    if (slash != std::string::npos && slash != 0 && !make_dirs(text.substr(0, slash))) {
        fprintf(stderr, "Error: cannot create log directory for %s\n", path);
        return false;
    }
    if (!freopen(path, "w", stdout)) {
        fprintf(stderr, "Error: cannot open log: %s\n", path);
        return false;
    }
    return true;
}

inline bool dump_field(int step) {
    if (!make_dirs(cc::fieldpath)) {
        fprintf(stderr, "Error: cannot create field directory: %s\n", cc::fieldpath.c_str());
        return false;
    }
    char tag[32];
    snprintf(tag, sizeof(tag), "step_%06d.dat", step);
    const std::string name = cc::fieldpath + "/" + tag;
    FILE* fp = fopen(name.c_str(), "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot create field file: %s\n", name.c_str());
        return false;
    }
    fprintf(fp, "TITLE=\"step %d\"\n", step);
    fprintf(fp, "VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\",\"nu_tilde\"\n");
    for (const cc::cell_class& cell : cc::CellList) {
        fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", cell.center.x, cell.center.y,
                cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.T, cell.otphy.p,
                std::hypot(cell.phy.u, cell.phy.v) / cell.otphy.a, cell.tur.miubl);
    }
    if (ferror(fp)) {
        fclose(fp);
        fprintf(stderr, "Error: failed to write field file: %s\n", name.c_str());
        return false;
    }
    if (fclose(fp) != 0) {
        fprintf(stderr, "Error: failed to close field file: %s\n", name.c_str());
        return false;
    }
    printf("[场输出] step %d  %s\n", step, name.c_str());
    return true;
}

/*
Json 读取部分
*/

using Json = nlohmann::json;

static bool fail(const std::string& msg) {
    fprintf(stderr, "Error: %s\n", msg.c_str());
    return false;
}

// 校验对象的键: 不在白名单内的报错, 必填键缺失也报错; optional 里的键可以缺席
static bool keys(const Json& object, std::initializer_list<const char*> allowed,
                 std::initializer_list<const char*> optional = {}) {
    if (!object.is_object()) { return fail("expected a configuration object"); }
    const std::set<std::string> names(allowed.begin(), allowed.end());
    const std::set<std::string> extra(optional.begin(), optional.end());
    for (const auto& item : object.items()) {
        if (!names.count(item.key()) && !extra.count(item.key())) {
            return fail("unknown configuration key: " + item.key());
        }
    }
    for (const std::string& name : names) {
        if (!object.contains(name)) { return fail("missing configuration key: " + name); }
    }
    return true;
}

// 读取一个有限的数, 要求为正
static bool positive(const Json& object, const char* key, double& out) {
    const auto item = object.find(key);
    if (item == object.end() || !item->is_number()) {
        return fail(std::string(key) + " must be a number");
    }
    const double value = item->get<double>();
    if (!std::isfinite(value) || value <= 0.0) {
        return fail(std::string(key) + " must be finite and positive");
    }
    out = value;
    return true;
}

// 读取一个有限的数, 符号不限
static bool number(const Json& object, const char* key, double& out) {
    const auto item = object.find(key);
    if (item == object.end() || !item->is_number()) {
        return fail(std::string(key) + " must be a number");
    }
    const double value = item->get<double>();
    if (!std::isfinite(value)) { return fail(std::string(key) + " must be finite"); }
    out = value;
    return true;
}

// 读取一个正整数, 按 32 位整数处理
static bool count(const Json& object, const char* key, int& out) {
    double value = 0.0;
    if (!positive(object, key, value)) { return false; }
    if (!object.find(key)->is_number_integer() || value > std::numeric_limits<int>::max()) {
        return fail(std::string(key) + " must be a positive 32-bit integer");
    }
    out = static_cast<int>(value);
    return true;
}

// 从dump_field写出的文件读回流场; x,y 用来确认这确实是本网格算出来的流场
inline bool load_field(const char* path, int& step) {
    std::ifstream stream(path);
    if (!stream) { return fail("Cannot open restart field: " + std::string(path)); }
    std::string line;
    step = 0;
    if (!std::getline(stream, line)) { return fail("restart field is empty"); }
    const size_t tag = line.find("step ");
    if (tag != std::string::npos) { step = (int)std::strtol(line.c_str() + tag + 5, nullptr, 10); }
    if (step < 0) { return fail("restart field has a negative step number"); }
    if (!std::getline(stream, line) || line.find("VARIABLES") == std::string::npos) {
        return fail("restart field is missing the VARIABLES header");
    }
    for (int i = 0; i < cc::cell_num; i++) {
        double x = 0.0, y = 0.0, rho = 0.0, u = 0.0, v = 0.0, T = 0.0, p = 0.0, ma = 0.0,
               nu = 0.0;
        if (!(stream >> x >> y >> rho >> u >> v >> T >> p >> ma >> nu)) {
            return fail("restart field has fewer cells than the mesh, stopped at cell " +
                        std::to_string(i + 1));
        }
        cc::cell_class& cell = cc::CellList[i];
        if (std::abs(x - cell.center.x) > 1e-6 * std::max(1.0, std::abs(cell.center.x)) ||
            std::abs(y - cell.center.y) > 1e-6 * std::max(1.0, std::abs(cell.center.y))) {
            return fail("restart field does not match the mesh, first at cell " +
                        std::to_string(i + 1));
        }
        if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(u) || !std::isfinite(v) ||
            !std::isfinite(T) || T <= 0.0 || !std::isfinite(nu) || nu < 0.0) {
            return fail("restart field has an invalid state at cell " + std::to_string(i + 1));
        }
        cell.phy = cc::vecp(rho, u, v, T);
        cell.form_otherphy();
        cell.form_conservative();
        cell.tur.miubl = nu;
    }
    double extra = 0.0;
    if (stream >> extra) { return fail("restart field has more cells than the mesh"); }
    for (cc::face_class& face : cc::FaceList) {
        if (face.type == cc::INTER) { face.face_physic_mid(); }
    }
    printf("流场读入 OK: %s  步号 %d  %d 个单元\n", path, step, cc::cell_num);
    return true;
}

// 读取路径字符串, 相对配置文件所在目录解析并做规范化
static bool path_value(const Json& object, const char* key, const std::filesystem::path& base,
                       std::string& out, std::error_code& error) {
    const auto item = object.find(key);
    if (item == object.end() || !item->is_string()) {
        return fail(std::string(key) + " must be a path string");
    }
    const std::string text = item->get<std::string>();
    if (text.empty() || text.find('\0') != std::string::npos) {
        return fail(std::string(key) + " is an invalid path");
    }
    std::filesystem::path full = std::filesystem::absolute(base / text, error);
    if (error) { full = base / text; }
    out = full.lexically_normal().string();
    return true;
}

inline bool config::load(const char* path) {
    // 打开配置文件
    std::ifstream stream(path);
    if (!stream) { return fail("Cannot open configuration: " + std::string(path)); }
    // 解析过程中逐层记录已出现的键, 用来发现重复键与过深嵌套
    std::vector<std::set<std::string>> object_keys;
    bool malformed = false;
    auto callback = [&](int depth, Json::parse_event_t event, Json& parsed) {
        if (depth > 16) { malformed = true; }
        if (event == Json::parse_event_t::object_start) {
            object_keys.emplace_back();
        } else if (event == Json::parse_event_t::key) {
            if (object_keys.empty() ||
                !object_keys.back().insert(parsed.get<std::string>()).second) {
                malformed = true;
            }
        } else if (event == Json::parse_event_t::object_end && !object_keys.empty()) {
            object_keys.pop_back();
        }
        return true;
    };
    const Json root = Json::parse(stream, callback, false);
    if (root.is_discarded()) { return fail("Cannot parse configuration: " + std::string(path)); }
    if (malformed) { return fail("Duplicate configuration key or too deep nesting"); }
    // 顶层只允许 io / solver / farfield 三节
    if (!keys(root, {"io", "solver", "farfield"})) { return false; }
    const Json& io = root.at("io");
    const Json& solver = root.at("solver");
    const Json& far = root.at("farfield");
    // 各节的键白名单; io 的 structured 只在网格带结构化信息时出现
    if (!keys(io, {"mesh", "log", "field"}, {"structured", "restart"}) ||
        !keys(solver, {"cfl", "dump_interval", "convergence_interval"},
              {"max_steps", "end_step", "threads", "urans", "miubl_inf_ratio"}) ||
        !keys(far, {"Ma", "T", "p", "alpha"})) {
        return false;
    }
    // 所有文件路径都相对配置文件所在目录解析
    std::error_code error;
    std::filesystem::path base = std::filesystem::path(path).parent_path();
    if (base.empty()) { base = "."; }
    std::string mesh_path, log_path, field_path;
    if (!path_value(io, "mesh", base, mesh_path, error) ||
        !path_value(io, "log", base, log_path, error) ||
        !path_value(io, "field", base, field_path, error)) {
        return false;
    }
    // structured 给出结构化邻接表路径并置起令牌; 缺省表示网格不带结构化信息
    if (io.contains("structured")) {
        structer::ifstructer = true;
        if (!path_value(io, "structured", base, structer::adjacency, error)) { return false; }
    }
    // restart 给出起始流场文件; 缺省表示从自由来流初始化
    if (io.contains("restart") &&
        !path_value(io, "restart", base, cc::restartpath, error)) {
        return false;
    }
    // 日志不得覆盖任何输入文件
    const std::string config_path =
        std::filesystem::absolute(path, error).lexically_normal().string();
    if (log_path == mesh_path || log_path == config_path) {
        return fail("The log path must not overwrite an input file");
    }
    // 求解器参数
    int dump_interval = 0, conv_interval = 0;
    if (!count(solver, "dump_interval", dump_interval) ||
        !count(solver, "convergence_interval", conv_interval)) {
        return false;
    }
    if (conv_interval < 2) { return fail("convergence_interval must be at least 2"); }
    // 步数两种给法: max_steps 是本趟步数(增量), end_step 是目标步号(绝对, 续算时最直观)
    long long max_steps = 0, end_step = 0;
    if (solver.contains("max_steps") && solver.contains("end_step")) {
        return fail("max_steps and end_step are mutually exclusive");
    }
    if (solver.contains("max_steps")) {
        int n = 0;
        if (!count(solver, "max_steps", n)) { return false; }
        max_steps = n;
    } else if (solver.contains("end_step")) {
        int n = 0;
        if (!count(solver, "end_step", n)) { return false; }
        end_step = n;
    } else {
        return fail("missing configuration key: max_steps or end_step");
    }
    // threads 缺省表示按可用物理核心自动选择
    int threads = 0;
    if (solver.contains("threads") && !count(solver, "threads", threads)) { return false; }
    // URANS参数, 缺省项保持config.hpp中的默认值
    double u_dt = 0.0, u_tol = urans::inner_tol, u_cfl = urans::steady_cfl;
    int u_inner = urans::inner, u_sweeps = urans::sweeps, u_steady = urans::steady_iters,
        u_wall = urans::wall_interval;
    if (solver.contains("urans")) {
        const Json& u = solver.at("urans");
        if (!keys(u, {"dt", "inner"},
                  {"inner_tol", "sweeps", "steady_iters", "steady_cfl", "wall_interval"}) ||
            !positive(u, "dt", u_dt) || !count(u, "inner", u_inner) ||
            (u.contains("inner_tol") && !positive(u, "inner_tol", u_tol)) ||
            (u.contains("sweeps") && !count(u, "sweeps", u_sweeps)) ||
            (u.contains("steady_iters") && !count(u, "steady_iters", u_steady)) ||
            (u.contains("steady_cfl") && !positive(u, "steady_cfl", u_cfl)) ||
            (u.contains("wall_interval") && !count(u, "wall_interval", u_wall))) {
            return false;
        }
    }
    // 来流条件
    double cfl = 0.0, ma = 0.0, T = 0.0, p = 0.0, alpha = 0.0;
    if (!positive(solver, "cfl", cfl) || !positive(far, "Ma", ma) || !positive(far, "T", T) ||
        !positive(far, "p", p) || !number(far, "alpha", alpha)) {
        return false;
    }
    double miubl_ratio = config::miubl_inf_ratio;
    if (solver.contains("miubl_inf_ratio") &&
        !positive(solver, "miubl_inf_ratio", miubl_ratio)) {
        return false;
    }
    if (ma >= 1.0) { return fail("This solver requires a subsonic farfield: 0 < Ma < 1"); }
    if (std::abs(alpha) >= 90.0) {
        return fail("The farfield angle of attack must satisfy |alpha| < 90 degrees");
    }
    const double u_inf = ma * get_sonic_velocity(T);
    const double rho_inf = p / (cc::R * T);
    if (!std::isfinite(u_inf) || !std::isfinite(rho_inf) || rho_inf <= 0.0) {
        return fail("Farfield values produce an invalid thermodynamic state");
    }

    // 全部校验通过后再一次性写入全局状态
    cc::meshpath = mesh_path;
    cc::testpath = log_path;
    cc::fieldpath = field_path;
    cc::max_step = max_steps;
    cc::end_step = end_step;
    cc::threads = threads;
    config::miubl_inf_ratio = miubl_ratio;
    fatime::CFL = cfl;
    dump_step = dump_interval;
    conv_step = conv_interval;
    cc::FAR_DEFINE.AOA = alpha * (std::acos(-1.0) / 180.0);
    cc::FAR_DEFINE.phy.u = u_inf * std::cos(cc::FAR_DEFINE.AOA);
    cc::FAR_DEFINE.phy.v = u_inf * std::sin(cc::FAR_DEFINE.AOA);
    cc::FAR_DEFINE.phy.T = T;
    cc::FAR_DEFINE.otphy.p = p;
    cc::FAR_DEFINE.phy.rho = p/cc::R/T;
    cc::FAR_DEFINE.otphy.form_otphy(cc::FAR_DEFINE.phy);
    cc::FAR_DEFINE.Ma = u_inf/cc::FAR_DEFINE.otphy.a;
    urans::dt = u_dt;
    urans::inner = u_inner;
    urans::inner_tol = u_tol;
    urans::sweeps = u_sweeps;
    urans::steady_iters = u_steady;
    urans::steady_cfl = u_cfl;
    urans::wall_interval = u_wall;
    return true;
}

inline void write_wall_cp(FILE* fp, int step, double t, cc::vecp refp, double refpressure) {
    const double q = dynamic_pressure(refp);
    fprintf(fp, "# step %d t %.8e\n", step, t);
    for (const cc::face_class* w : cc::WallFaces) {
        fprintf(fp, "%.6e %.6e %.6e\n", w->mid.x, w->mid.y, (w->otphy.p -refpressure) / q);
    }
    fflush(fp);
}

