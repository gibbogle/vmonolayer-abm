#ifndef RESULT_SET_H
#define RESULT_SET_H
struct result_set {
	QString casename;
	int nsteps;
    double *tnow;      // simulation time (hours)
    double max_time;
    double *pData[maxGraphs];      // was 16, should be maxGraphs?
    double maxValue[maxGraphs];
};
typedef result_set RESULT_SET;
#endif
