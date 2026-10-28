#include "HALO.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "geometry.hpp"
#include "implict.hpp"
#include "initialize.hpp"
#include "io.hpp"
#include "readmesh.hpp"
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <string>

static volatile std::sig_atomic_t stop_requested = 0;
static void handle_stop(int) { stop_requested = 1; }

static std::vector<cc::vec4> snap_q;
static std::vector<double> snap_nu;
static double cfl_scale = 1.0;

static void save_state() {
    if ((int)snap_q.size() != cc::cell_num) {
        snap_q.resize(cc::cell_num);
        snap_nu.resize(cc::cell_num);
    }
    for (int i = 0; i < cc::cell_num; i++) {
        snap_q[i] = icell(i).conser;
        snap_nu[i] = icell(i).tur.miubl;
    }
}

static void restore_state() {
    for (int i = 0; i < cc::cell_num; i++) {
        icell(i).conser = snap_q[i];
        icell(i).tur.miubl = snap_nu[i];
        icell(i).prim();
        icell(i).form_otherphy();
    }
}

int main(int argc, char** argv) {

    // Prepare Period
    const char* cfg = argc > 1 ? argv[1] : "config.json";
    if (!config::load(cfg)) return 1;
    std::signal(SIGINT, handle_stop);
    std::signal(SIGTERM, handle_stop);
    if (urans::dt <= 0.0) {
        fprintf(stderr, "Error: solver.urans is required\n");
        return 1;
    }
    if (!make_dirs(cc::fieldpath)) return 1;
    if (!open_log(cc::testpath.c_str())) return 1;
#ifdef _OPENMP
    // OMP_NUM_THREADS 优先, 其次 solver.threads
    const char* env_threads = getenv("OMP_NUM_THREADS");
    if ((!env_threads || !*env_threads) && cc::threads > 0) omp_set_num_threads(cc::threads);
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif


    // Pre-Process
    if (!readmesh(cc::meshpath.c_str()) || !geometrymain()) return 1;
    std_initialize();
    allcell icell(i).form_conservative();
    allcell least_square_cell_based_preprocess(icell(i));
    HALO_structer_mesh();
    int restart_step = 0;
    if (!cc::restartpath.empty() && !load_field(cc::restartpath.c_str(), restart_step)) return 1;
    if (!cc::restartpath.empty()) {
        printf("续算: 从 step %d 起, 物理时间 t = %.6e s 起算\n", restart_step,
               restart_step * urans::dt);
        if (urans::steady_iters > 0) {
            printf("注意: steady_iters=%d, 读入的流场会先被定常迭代改掉\n", urans::steady_iters);
        }
    }
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

    // 保存流场
    auto finish = [&](int code, const char* why, int save_step) {
        if (save_step >= 0) {
            allcell icell(i).prim();
            if (!dump_field(save_step)) { code = 1; }
            printf("%s: 已保存 step %d 的流场到 %s\n", why, save_step, cc::fieldpath.c_str());
        }
        fclose(hist);
        if (wall) fclose(wall);
        if (code == 0) { printf("Done in %.1fs\n", elapsed()); }
        return code;
    };

    // RANS
    for (int it = 1; it <= urans::steady_iters;) {
        if (stop_requested) { return finish(0, "用户中断定常迭代", it - 1); }
        const double cfl =
            std::min(urans::steady_cfl, 2.0 + (urans::steady_cfl - 2.0) * it / 200.0) * cfl_scale;
        save_state();
        const double r = impl::flow_step(0.0, cfl);
        impl::turbulence_step(0.0, cfl);
        const int nb = bad_cells();
        if (nb > 0) {
            printf("steady %6d 发散: %d 个坏单元, 回退上一流场, CFL x%.4f -> x%.4f\n", it, nb,
                   cfl_scale, cfl_scale * 0.5);
            fflush(stdout);
            restore_state();
            cfl_scale *= 0.5;
            if (cfl_scale < 1e-3) {
                printf("CFL 已降到下限仍不收敛, 放弃\n");
                return finish(1, "定常迭代不收敛", it - 1);
            }
            continue;
        }
        if (cfl_scale < 1.0) cfl_scale = std::min(1.0, cfl_scale * 1.05);
        if (it % config::conv_step == 0 || it == urans::steady_iters) {
            wall_forces(cl,cd,cc::FAR_DEFINE.phy);
            printf("steady %6d  cfl %6.1f  res %.4e  Cl %.5f  Cd %.5f  clip %5d  [%.0fs]\n", it, cfl,
                   r, cl, cd, impl::n_clip, elapsed());
            fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", it, r, r, cl, cd);
            fflush(stdout);
        }
        ++it;
    }
    if (urans::steady_iters > 0) {
        allcell icell(i).prim();
        if (!dump_field(restart_step)) return 1;
    }

    // URANS
    const int step_0 = restart_step; // 续算时接着文件里的步号往下数
    const int last_step = cc::end_step > 0 ? (int)cc::end_step : step_0 + (int)cc::max_step;
    if (last_step <= step_0) {
        printf("Error: end_step %d 不大于起始步号 %d, 没有步要跑\n", last_step, step_0);
        fclose(hist);
        if (wall) fclose(wall);
        return 1;
    }
    printf("步号区间: %d -> %d, 共 %d 步\n", step_0 + 1, last_step, last_step - step_0);
    double t = (double)step_0 * urans::dt;
    for (int n = step_0 + 1; n <= last_step; n++) {
#pragma omp parallel for schedule(static)
        allcell {
            icell(i).lastconser = icell(i).conser;
            icell(i).tur.miubl_former = icell(i).tur.miubl;
        }
        const double cfl_nom = fatime::CFL;
        double r0 = 0.0, r = 0.0;
        int k = 0;
        for (;;) {
            save_state();
            fatime::CFL = cfl_nom * cfl_scale;
            r0 = 0.0;
            r = 0.0;
            k = 0;
            while (k < urans::inner) {
                r = impl::flow_step(urans::dt, fatime::CFL);
                impl::turbulence_step(urans::dt, fatime::CFL);
                if (++k == 1)
                    r0 = r;
                else if (r < urans::inner_tol * r0)
                    break;
            }
            const int nb = bad_cells();
            if (nb == 0) break;
            printf("step %6d 发散: %d 个坏单元, 回退上一流场, CFL x%.4f -> x%.4f\n", n, nb, cfl_scale,
                   cfl_scale * 0.5);
            fflush(stdout);
            restore_state();
            cfl_scale *= 0.5;
            if (cfl_scale < 1e-3) {
                printf("CFL 已降到下限仍不收敛, 放弃\n");
                return finish(1, "该步发散", n - 1);
            }
        }
        fatime::CFL = cfl_nom;
        if (cfl_scale < 1.0) cfl_scale = std::min(1.0, cfl_scale * 1.05);
        t += urans::dt;
        wall_forces(cl, cd,cc::FAR_DEFINE.phy);
        fprintf(hist, "U,%d,%.8e,%d,%.6e,%.6e,%.8f,%.8f\n", n, t, k, r0, r, cl, cd);
        fflush(hist);
        if (n % config::conv_step == 0 || n == step_0 + 1) {
            printf("step %6d  t %.5e  inner %2d  res %.3e -> %.3e  Cl %.5f  Cd %.5f  "
                   "clip %5d  [%.0fs]\n",
                   n, t, k, r0, r, cl, cd, impl::n_clip, elapsed());
            fflush(stdout);
        }
        if (wall && n % urans::wall_interval == 0) write_wall_cp(wall, n, t,cc::FAR_DEFINE.phy,cc::FAR_DEFINE.otphy.p);
        if (stop_requested) { return finish(0, "用户中断", n); }
        if (n % config::dump_step == 0) {
            allcell icell(i).prim();
            if (!dump_field(n)) return 1;
        }
    }
    return finish(0, nullptr, -1);
}
