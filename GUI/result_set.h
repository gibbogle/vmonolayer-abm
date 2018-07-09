#ifndef RESULT_SET_H
#define RESULT_SET_H
struct result_set {
	QString casename;
	int nsteps;
    double *tnow;      // simulation time (hours)
    double max_time;
	double *pData[16];
	double maxValue[16];
};
typedef result_set RESULT_SET;
#endif
