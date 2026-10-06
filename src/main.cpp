#include "HALO.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "dualtime.hpp"
#include "dualtime_full.hpp"
#include "geometry.hpp"
#include "initialize.hpp"
#include "io.hpp"
#include "readmesh.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef PURINE_JAC_FULL
inline double flow_step_select(double dt,double cfl){ return dual_full::flow_step(dt,cfl); }
#else
inline double flow_step_select(double dt,double cfl){ return dual::flow_step(dt,cfl); }
#endif

// 来流动压
static double dynamic_pressure() {
    const double rho_inf = FAR_DEFINE.p / (cc::R * FAR_DEFINE.T);
    return 0.5 * rho_inf * (FAR_DEFINE.u * FAR_DEFINE.u + FAR_DEFINE.v * FAR_DEFINE.v);
}

// 壁面合力(压力+摩擦)的升阻力系数, 参考长度取1(圆柱直径/翼型弦长)
static void wall_forces(double& cl, double& cd) {
    const double q = dynamic_pressure();
    double fx = 0.0, fy = 0.0;
    for (cc::face_class* w : cc::WallFaces) {
        const cc::cell_class* c = cc::boundary_findcell(w);
        // s*nor 由流体指向物面
        const double s = cc::dot(w->nor, w->mid - c->center) > 0 ? 1.0 : -1.0;
        fx += s * (w->otphy.p * w->nor.x - w->visflux.x);
        fy += s * (w->otphy.p * w->nor.y - w->visflux.y);
    }
    const double a = std::atan2(FAR_DEFINE.v, FAR_DEFINE.u);
    cl = (-fx * std::sin(a) + fy * std::cos(a)) / q;
    cd = (fx * std::cos(a) + fy * std::sin(a)) / q;
}

static void write_wall_cp(FILE* fp, int step, double t) {
    const double q = dynamic_pressure();
    fprintf(fp, "# step %d t %.8e\n", step, t);
    for (const cc::face_class* w : cc::WallFaces) {
        fprintf(fp, "%.6e %.6e %.6e\n", w->mid.x, w->mid.y, (w->otphy.p - FAR_DEFINE.p) / q);
    }
    fflush(fp);
}

static bool field_valid() {
    for (const cc::cell_class& c : cc::CellList) {
        if (!std::isfinite(c.conser.c) || c.conser.c <= 0.0 || !std::isfinite(c.conser.e) ||
            !std::isfinite(c.tur.miubl)) {
            printf("Error: invalid state in cell #%d (%g,%g)\n", c.index, c.center.x, c.center.y);
            return false;
        }
    }
    return true;
}

// 初场叠加绕原点的高斯涡, 打破对称性使圆柱尽快起振
static void add_seed(double amp) {
    for (cc::cell_class& c : cc::CellList) {
        const double r = std::hypot(c.center.x, c.center.y);
        if (r < 0.5 || r > 3.0) continue;
        const double w = amp * std::exp(-r * r / 4.0);
        c.phy.u -= w * c.center.y;
        c.phy.v += w * c.center.x;
        c.form_otherphy();
        c.form_conservative();
    }
}

int main(int argc, char** argv) {
    const char* cfg = argc > 1 ? argv[1] : "config.json";
    if (!config::load(cfg)) return 1;
    if (urans::dt <= 0.0) {
        fprintf(stderr, "Error: solver.urans is required\n");
        return 1;
    }
    if (!make_dirs(cc::fieldpath)) {
        fprintf(stderr, "Error: cannot create field directory: %s\n", cc::fieldpath.c_str());
        return 1;
    }
    if (!open_log(cc::testpath.c_str())) return 1;
#ifdef _OPENMP
    // OMP_NUM_THREADS 优先, 其次 solver.threads
    const char* env_threads = getenv("OMP_NUM_THREADS");
    if ((!env_threads || !*env_threads) && cc::threads > 0) omp_set_num_threads(cc::threads);
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif
    if (!readmesh(cc::meshpath.c_str()) || !geometrymain()) return 1;
    std_initialize();
    allcell icell(i).form_conservative();
    // LSCB预处理必须在HALO之前: 虚网格没有几何
    allcell least_square_cell_based_preprocess(icell(i));
    HALO_structer_mesh();
    if (urans::seed > 0.0) add_seed(urans::seed * std::hypot(FAR_DEFINE.u, FAR_DEFINE.v));

    FILE* hist = fopen((cc::fieldpath + "/history.csv").c_str(), "w");
    FILE* wall =
        urans::wall_interval > 0 ? fopen((cc::fieldpath + "/wall_cp.dat").c_str(), "w") : nullptr;
    if (!hist || (urans::wall_interval > 0 && !wall)) {
        printf("Error: cannot create output files in %s\n", cc::fieldpath.c_str());
        return 1;
    }
    fprintf(hist, "phase,step,time,inner,res0,res,Cl,Cd\n");
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    double cl = 0.0, cd = 0.0;

    // 定常隐式迭代给URANS初场, CFL在前200步从2线性增大到steady_cfl
    for (int it = 1; it <= urans::steady_iters; it++) {
        const double cfl =
            std::min(urans::steady_cfl, 2.0 + (urans::steady_cfl - 2.0) * it / 200.0);
        const double r = flow_step_select(0.0, cfl);
        dual::sa_step(0.0, cfl);
        if (it % config::conv_step == 0 || it == urans::steady_iters) {
            wall_forces(cl, cd);
            printf("steady %6d  cfl %6.1f  res %.4e  Cl %.5f  Cd %.5f  [%.0fs]\n", it, cfl, r, cl,
                   cd, elapsed());
            fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", it, r, r, cl, cd);
            fflush(stdout);
            if (!field_valid()) return 1;
        }
    }
    if (urans::steady_iters > 0) {
        allcell icell(i).prim();
        if (!dump_field(0)) return 1;
    }

    // URANS: 每个物理步内做最多inner次隐式伪时间迭代
    double t = 0.0;
    for (int n = 1; n <= cc::max_step; n++) {
#pragma omp parallel for schedule(static)
        allcell {
            icell(i).lastconser = icell(i).conser;
            icell(i).tur.miubl_former = icell(i).tur.miubl;
        }
        double r0 = 0.0, r = 0.0;
        int k = 0;
        while (k < urans::inner) {
            r = flow_step_select(urans::dt, fatime::CFL);
            dual::sa_step(urans::dt, fatime::CFL);
            if (++k == 1)
                r0 = r;
            else if (r < urans::inner_tol * r0)
                break;
        }
        t += urans::dt;
        wall_forces(cl, cd);
        fprintf(hist, "U,%d,%.8e,%d,%.6e,%.6e,%.8f,%.8f\n", n, t, k, r0, r, cl, cd);
        fflush(hist);
        if (n % config::conv_step == 0 || n == 1) {
            printf("step %6d  t %.5e  inner %2d  res %.3e -> %.3e  Cl %.5f  Cd %.5f  [%.0fs]\n", n,
                   t, k, r0, r, cl, cd, elapsed());
            fflush(stdout);
            if (!field_valid()) return 1;
        }
        if (wall && n % urans::wall_interval == 0) write_wall_cp(wall, n, t);
        if (n % config::dump_step == 0) {
            allcell icell(i).prim();
            if (!dump_field(n)) return 1;
        }
    }
    fclose(hist);
    if (wall) fclose(wall);
    printf("Done in %.1fs\n", elapsed());
    return 0;
}
