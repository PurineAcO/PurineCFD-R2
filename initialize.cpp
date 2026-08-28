#include "config.h"
#include "classconfig.h"
#include "physic.h"

void std_initialize(cc::cell_class &cell){
    cell.phy.T = cc::T;
    cell.phy.a = get_sonic_velocity(cell.phy.T);
    cell.phy.u = cell.phy.a * cc::Ma * deg2rad(cc::AOA);
    cell.phy.v = cell.phy.a * cc::Ma * deg2rad(cc::AOA);
    cell.phy.mu = SutherLand(cell.phy.T);
    cell.phy.rho = cc::Re * cell.phy.mu /(cell.phy.a * cc::Ma *cc::H);
    cell.phy.p = cell.phy.rho * cc::R * cell.phy.rho;
    cell.tur.miubl = cell.phy.mu * 0.1 /cell.phy.rho;
}