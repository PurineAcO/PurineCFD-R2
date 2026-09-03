#include <cstdio>
#include <sys/stat.h>
#include "io.h"
#include "classconfig.h"
#include "config.h"

void dump_field(int step){
    static char dir[256];
    snprintf(dir, sizeof(dir), "%s", cc::fieldpath);
    mkdir(dir, 0755);
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/step_%06d.dat", dir, step);
    FILE* fp = fopen(fname, "w");
    if(!fp){ return; }
    fprintf(fp, "TITLE=\"step %d\"\n", step);
    fprintf(fp, "VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\"\n");
    for(cc::cell_class& cell : cc::CellList){
        fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                cell.center.x, cell.center.y, cell.phy.rho, cell.phy.u, cell.phy.v,
                cell.phy.T, cell.phy.p, cell.phy.u/cell.phy.a);
    }
    fclose(fp);
    printf("[场输出] step %d  %s\n", step, fname);
}

void save_checkpoint(int step){
    static char dir[256];
    snprintf(dir, sizeof(dir), "%s", cc::fieldpath);
    mkdir(dir, 0755);
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/checkpoint.dat", dir);
    FILE* fp = fopen(fname, "w");
    if(!fp){ return; }
    fprintf(fp, "%d\n", step);
    for(cc::cell_class& cell : cc::CellList){
        fprintf(fp, "%.10e %.10e %.10e %.10e\n",
                cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.T);
    }
    fclose(fp);
    printf("[checkpoint] step %d\n", step);
}

bool load_checkpoint(int& step){
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/checkpoint.dat", cc::fieldpath);
    FILE* fp = fopen(fname, "r");
    if(!fp){ return false; }
    if(fscanf(fp, "%d", &step) != 1){ fclose(fp); return false; }
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
