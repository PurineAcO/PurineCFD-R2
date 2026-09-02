#include <cstdio>
#include <sys/stat.h>
#include "io.h"
#include "classconfig.h"
#include "config.h"

void dump_field(double t){
    static char dir[256];
    snprintf(dir, sizeof(dir), "%s", cc::fieldpath);
    mkdir(dir, 0755);
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/t_%06.2fms.dat", dir, t*1000.0);
    FILE* fp = fopen(fname, "w");
    if(!fp){ return; }
    fprintf(fp, "TITLE=\"t=%.6fs\"\n", t);
    fprintf(fp, "VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\"\n");
    for(cc::cell_class& cell : cc::CellList){
        fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                cell.center.x, cell.center.y, cell.phy.rho, cell.phy.u, cell.phy.v,
                cell.phy.T, cell.phy.p, cell.phy.u/cell.phy.a);
    }
    fclose(fp);
    printf("[场输出] t=%9.6fs  %s\n", t, fname);
}

void save_checkpoint(double t){
    static char dir[256];
    snprintf(dir, sizeof(dir), "%s", cc::fieldpath);
    mkdir(dir, 0755);
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/checkpoint.dat", dir);
    FILE* fp = fopen(fname, "w");
    if(!fp){ return; }
    fprintf(fp, "%.10e\n", t);
    for(cc::cell_class& cell : cc::CellList){
        fprintf(fp, "%.10e %.10e %.10e %.10e\n",
                cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.T);
    }
    fclose(fp);
    printf("[checkpoint] t=%9.6fs  %s\n", t, fname);
}

bool load_checkpoint(double& t){
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/checkpoint.dat", cc::fieldpath);
    FILE* fp = fopen(fname, "r");
    if(!fp){ return false; }
    if(fscanf(fp, "%lf", &t) != 1){ fclose(fp); return false; }
    for(cc::cell_class& cell : cc::CellList){
        if(fscanf(fp, "%lf %lf %lf %lf",
                  &cell.phy.rho, &cell.phy.u, &cell.phy.v, &cell.phy.T) != 4){
            fclose(fp); return false;
        }
        cell.form_physic();
    }
    fclose(fp);
    return true;
}
