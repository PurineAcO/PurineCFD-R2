#define FAR_FIELD(p,Ma,T,AOA)
    extern void farfield_boundary_init(double p,double Ma,double T,double AOA);
#define Velocity_Inlet_Undersonic(u,v,T)
    extern void velocity_inlet_undersonic_boundary_init(double u,double v,double T);
#define Pressure_Outlet(p,T_stg)
    extern void pressure_outlet_boundary_init(double p,double T_stg);
#define OT_flowfield(eachstep)
    extern void output_flowfield(int eachstep);