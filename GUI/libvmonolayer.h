#ifndef LIBVMONOLAYER_H
#define LIBVMONOLAYER_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void get_dll_build_version(char *, int *);
void execute(int *,char *, int *,char *, int *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *, double *, int *, int *, int *, bool *);
//void get_summary(int *, int *, int *);
void get_summary(double *, int *, int *);
void get_nfacs(int *);
void get_facs(double *, double *, double *, double *, double *, bool);
void get_histo(int, double *, double *, double *, double *, double *, double *, bool);
void get_constituents(int *, int *, int *, char *, int *);
void make_colony_distribution(double *, double *, double *, int *, double *);

// For 3D display
//void get_scene(int *, int *);
// These are for profile plots
void get_concdata(int *, int *, double *, double *);
//void get_ic_concdata(int *, int *, double *, double *);
//void get_volprob(int *, double *, double *, double*);
//void get_oxyprob(int *, double *, double *, double *);
void get_pi_dist(int, double *, double *, double *);

void get_string(char **);

//
//
#ifdef __cplusplus
}
#endif

#endif // LIBVMONOLAYER_H
