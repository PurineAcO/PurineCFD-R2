#pragma once
#include "mesh.h"

void prepare_convective_flux(cfd::Face& face);
void assemble_face_fluxes(cfd::Cell& cell);
