#ifndef GLOBAL_H
#define GLOBAL_H

#include <QtGui>

// Note that in the Fortran DLL the chemokine constituent numbering starts at 1
#define MULTI -1
#define CFSE 0
#define OXYGEN 1
#define GLUCOSE 2
#define LACTATE 3
#define TRACER 4
#define DRUG_A_PARENT 5
#define DRUG_A_METAB_1 6
#define DRUG_A_METAB_2 7
#define DRUG_B_PARENT 8
#define DRUG_B_METAB_1 9
#define DRUG_B_METAB_2 10
#define GROWTH_RATE 11      // we pretend that this is a concentration
#define CELL_VOLUME 12
#define O2_BY_VOL 13
#define CYCLE_PHASE 14

#define PI_DEAD 15

// The intracellular (IC) dataIndex is the same as extracellular, from the tag determine which
#define IC_MULTI -1
#define IC_CFSE 0
#define IC_OXYGEN 1
#define IC_GLUCOSE 2
#define IC_LACTATE 3
#define IC_TRACER 4
#define IC_DRUG_A_PARENT 5
#define IC_DRUG_A_METAB_1 6
#define IC_DRUG_A_METAB_2 7
#define IC_DRUG_B_PARENT 8
#define IC_DRUG_B_METAB_1 9
#define IC_DRUG_B_METAB_2 10
#define IC_GROWTH_RATE 11      // we pretend that this is a concentration
#define IC_CELL_VOLUME 12
#define IC_O2_BY_VOL 13
#define IC_CYCLE_PHASE 14

#define DIST_NV 20
#define PI_NBINS 100

#define MAX_CELLS 200000
#define N_CELLINFO 7
//#define N_FACS_VARS 3

struct dist_set {
    bool used;
    double dv;
    double v0;
    double prob[DIST_NV];
};
typedef dist_set DIST_SET;

struct cell_data {
    int tag;
    double radius;
    double centre[3];
    int celltype;
    int status;
//    int highlight;
};
typedef cell_data CELL_DATA;

namespace Global
{
    extern QString GUI_build_version;
    extern QString DLL_build_version;

    extern int MAX_CHEMO;
    extern int N_EXTRA;
    extern int NX, NY, NZ;
    extern double DELTA_T;
    extern int NT_DISPLAY;
    extern double DELTA_X;
    extern double dfraction;
    extern int nt_vtk;
    extern int istep;
    extern bool leftb;

    extern int nvars_used;
    extern int nfieldvars_used;
    extern int GUI_to_DLL_index[32];
    extern int DLL_to_GUI_index[32];
    extern QString var_string[32];

    extern double *FACS_data;
    extern int nFACS_cells;
    extern int nFACS_dim;
    extern double FACS_vmin[32];
    extern double FACS_vmax[32];
    extern double FACS_vmin_log[32];
    extern double FACS_vmax_log[32];

    extern double *histo_data;
    extern double *histo_data_log;
    extern int nhisto_bins;
    extern int nhisto_dim;
    extern double histo_vmin[3*32];
    extern double histo_vmax[3*32];
    extern double histo_vmin_log[3*32];
    extern double histo_vmax_log[3*32];
    extern int histo_celltype;

//    extern int summaryData[100];
    extern double summaryData[100];
    extern int i_hypoxia_cutoff;
    extern int i_growth_cutoff;

    extern double concData[4000];
    extern double IC_concData[4000];
    extern int conc_nvars;
    extern int conc_nc_ex;
    extern int conc_nc_ic;
    extern double conc_dx_ex;
    extern double conc_dx_ic;
    extern QString casename;

    extern double volProb[100];
    extern int vol_nv;
    extern double vol_v0;
    extern double vol_dv;
    extern double oxyProb[100];
    extern int oxy_nv;
    extern double oxy_v0;
    extern double oxy_dv;

//    extern double distData[4000];
//    extern bool dist_used[20];
    extern int dist_nv;
    extern DIST_SET distParams[20];

    extern int cell_list[N_CELLINFO*MAX_CELLS];
    extern int ncell_list;

    extern bool showingVTK;
    extern bool recordingVTK;
    extern bool showingFACS;
    extern bool recordingFACS;
    extern bool showingField;
    extern bool recordingField;

    extern bool simulate_colony;
    extern double colony_days;
    extern double dist[40];
    extern double ddist;
    extern int ndist;

    extern double PI_fract[PI_NBINS];
    extern double PI_max_fract;
    extern double PI_max_fluor;

    extern bool volume_scaling;

} // namespace Global

#endif // GLOBAL_H
