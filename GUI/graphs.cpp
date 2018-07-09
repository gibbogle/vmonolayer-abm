#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
GRAPH_SET tsGraphSet[] = {

    {"nlive",
    "Live Cells",
    "No. of cells",
    "Number of live cells in the blob",
    1, false, 0, 1, 0, TS_TYPE},

    {"nviable",
    "Viable Cells",
    "No. of viable cells",
    "Number of viable cells in the blob",
    2, false, 0, 1, 0, TS_TYPE},

    {"nanoxiadead",
    "Anoxia-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by anoxia",
    3, false, 0, 1, 0, TS_TYPE},

    {"naglucosiadead",
    "Aglucosia-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by aglucosia",
    4, false, 0, 1, 0, TS_TYPE},

    {"ndrugAdead",
    "DrugA-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by drugA",
    5, false, 0, 1, 0, TS_TYPE},

    {"ndrugBdead",
    "DrugB-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by drugB",
    6, false, 0, 1, 0, TS_TYPE},

    {"nradiationdead",
    "Radiation-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by radiation",
    7, false, 0, 1, 0, TS_TYPE},

    {"nanoxiatagged",
    "Anoxia-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by anoxia",
    8, true, 0, 1, 0, TS_TYPE},

    {"naglucosiatagged",
    "Aglucosia-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by aglucosia",
    9, true, 0, 1, 0, TS_TYPE},

    {"ndrugAtagged",
    "DrugA-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by drugA",
    10, true, 0, 1, 0, TS_TYPE},

    {"ndrugBtagged",
    "DrugB-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by drugB",
    11, false, 0, 1, 0, TS_TYPE},

    {"nradiationtagged",
    "Radiation-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by radiation",
    12, false, 0, 1, 0, TS_TYPE},

    {"hypoxicfraction",
    "Hypoxic %",
    "%",
     "Percentage of cells with oxygen level below the specified threshold for hypoxia",
    13, true, 0, 0.1, 0, TS_TYPE},

    {"clonohypoxicfraction",
    "Clonogenic Hypoxic %",
    "%",
     "Percentage of clonogenic cells with oxygen level below the specified threshold for hypoxia",
    14, true, 0, 0.1, 0, TS_TYPE},

    {"growthfraction",
    "Slow-growth Fraction",
    "%",
     "Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits",
    15, false, 0, 0.1, 0, TS_TYPE},

    {"platingefficiency",
    "Plating Efficiency",
    "%",
     "Plating efficiency = 100*(number of viable cells)/(number of live cells)",
    16, true, 0, 0.1, 0, TS_TYPE},

    {"ECoxygen",
    "EC Oxygen",
    "Concentration",
    "Concentration of oxygen in the medium at the well bottom",
    17, true, 0, 0.001, 0, TS_TYPE},

    {"ECglucose",
    "EC Glucose",
    "Concentration",
     "Concentration of glucose in the medium at the well bottom",
    18, true, 0, 0.001, 0, TS_TYPE},

    {"EClactate",
    "EC Lactate",
    "Concentration",
     "Concentration of lactate in the medium at the well bottom",
    19, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugA",
    "EC Drug A",
    "Concentration",
     "Concentration of Drug A in the medium at the well bottom",
    20, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugAmet1",
    "EC Drug A metab1",
    "Concentration",
     "Concentration of Drug A metabolite 1 in the medium at the well bottom",
    21, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugAmet2",
    "EC Drug A metab2",
    "Concentration",
     "Concentration of Drug A metabolite 2 in the medium at the well bottom",
    22, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugB",
    "EC Drug B",
    "Concentration",
     "Concentration of Drug B in the medium at the well bottom",
    23, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugBmet1",
    "EC Drug B metab1",
    "Concentration",
     "Concentration of Drug B metabolite 1 in the medium at the well bottom",
    24, true, 0, 0.001, 0, TS_TYPE},

    {"ECdrugBmet2",
    "EC Drug B metab2",
    "Concentration",
     "Concentration of Drug B metabolite 2 in the medium at the well bottom",
    25, true, 0, 0.001, 0, TS_TYPE},

    {"ICoxygen",
    "IC Oxygen",
    "Concentration",
    "Intracellular concentration of oxygen",
    26, true, 0, 0.001, 0, TS_TYPE},

    {"ICglucose",
    "IC Glucose",
    "Concentration",
     "Intracellular concentration of glucose",
    27, true, 0, 0.001, 0, TS_TYPE},

    {"IClactate",
    "IC Lactate",
    "Concentration",
     "Intracellular concentration of lactate",
    28, true, 0, 0.001, 0, TS_TYPE},

    {"ICpyruvate",
    "IC Pyruvate",
    "Concentration",
     "Intracellular concentration of pyruvate",
    29, false, 0, 0.001, 0, TS_TYPE},

    {"ICdrugA",
    "IC Drug A",
    "Concentration",
     "Intracellular concentration of Drug A",
    30, true, 0, 0.001, 0, TS_TYPE},

    {"ICdrugAmet1",
    "IC Drug A metab1",
    "Concentration",
     "Intracellular concentration of Drug A metabolite 1",
    31, true, 0, 0.001, 0, TS_TYPE},

    {"ICdrugAmet2",
    "IC Drug A metab2",
    "Concentration",
     "Intracellular concentration of Drug A metabolite 2",
    32, true, 0, 0.001, 0, TS_TYPE},

    {"ICdrugB",
    "IC Drug B",
    "Concentration",
     "Intracellular concentration of Drug B",
    33, true, 0, 0.001, 0, TS_TYPE},

    {"ICdrugBmet1",
    "IC Drug B metab1",
    "Concentration",
     "Intracellular concentration of Drug B metabolite 1",
    34, true, 0, 0.001, 0, TS_TYPE},

    {"ICdrugBmet2",
    "IC Drug B metab2",
    "Concentration",
     "Intracellular concentration of Drug B metabolite 2",
    35, true, 0, 0.001, 0, TS_TYPE},

    {"Medoxygen",
    "Medium Oxygen",
    "Concentration",
    "Average concentration of oxygen in the medium",
    36, true, 0, 0.001, 0, TS_TYPE},

    {"Medglucose",
    "Medium Glucose",
    "Concentration",
     "Average concentration of glucose in the medium",
    37, true, 0, 0.001, 0, TS_TYPE},

    {"Medlactate",
    "Medium Lactate",
    "Concentration",
     "Average concentration of lactate in the medium",
    38, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugA",
    "Medium Drug A",
    "Concentration",
     "Average concentration of Drug A in the medium",
    39, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugAmet1",
    "Medium Drug A metab1",
    "Concentration",
     "Average concentration of Drug A metabolite 1 in the medium",
    40, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugAmet2",
    "Medium Drug A metab2",
    "Concentration",
     "Average concentration of Drug A metabolite 2 in the medium",
    41, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugB",
    "Medium Drug B",
    "Concentration",
     "Average concentration of Drug B in the medium",
    42, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugBmet1",
    "Medium Drug B metab1",
    "Concentration",
     "Average concentration of Drug B metabolite 1 in the medium",
    43, true, 0, 0.001, 0, TS_TYPE},

    {"MeddrugBmet2",
     "Medium Drug B metab2",
    "Concentration",
     "Average concentration of Drug B metabolite 2 in the medium",
    44, true, 0, 0.001, 0, TS_TYPE},

    {"doublingtime",
    "Ave Doubling time",
    "Hours",
    "Average doubling time",
    45, true, 0, 0.01, 0, TS_TYPE},

    {"Grate",
    "Glycolysis rate",
    "",
    "Normalised rate of glycolysis",
    46, false, 0, 0.001, 0, TS_TYPE},

    {"Prate",
    "Pyruvate utilisation rate",
    "",
    "Normalised rate of pyruvate utilisation",
    47, false, 0, 0.001, 0, TS_TYPE},

    {"Arate",
    "ATP production rate",
    "",
    "Normalised rate of ATP production",
    48, false, 0, 0.001, 0, TS_TYPE},

    {"Irate",
    "Intermediates production rate",
    "",
    "Normalised rate of production of anabolic intermediates",
    49, false, 0, 0.001, 0, TS_TYPE},

    {"f_G",
    "Glycolysis intermediates production factor",
    "",
    "Fraction of glucolysis rate going to production of anabolic intermediates",
    50, false, 0, 0.001, 0, TS_TYPE},

    {"f_P",
    "Pyruvate intermediates production factor",
    "",
    "Fraction of pyruvate utilisation rate going to production of anabolic intermediates",
    51, false, 0, 0.001, 0, TS_TYPE},

    {"HIF-1",
    "Normalised HIF-1 level",
    "",
    "Normalised HIF-1 level",
    52, false, 0, 0.001, 0, TS_TYPE},

    {"PDK1",
    "PDK1 factor level",
    "",
    "PDK1 factor level",
    53, false, 0, 0.001, 0, TS_TYPE},

    {"dividerate",
    "# divided/hour",
    "/hour",
    "# divided/hour",
    54, true, 0, 1, 0, TS_TYPE},


// Medium z-profiles

    {"MULTI",
    "Multi-constituent",
    "",
    "MULTI description",
    MULTI, true, 0, 1, 0, PROF_TYPE},

    {"Oxygen",
    "Oxygen Concentration",
    "",
    "Oxygen description",
    OXYGEN, false, 0, 1, 0, PROF_TYPE},

    {"Glucose",
    "Glucose Concentration",
    "",
    "Glucose description",
    GLUCOSE, false, 0, 1, 0, PROF_TYPE},

    {"Tracer",
    "Tracer Concentration",
    "",
    "Tracer description",
    TRACER, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A",
    "Drug A Concentration",
    "",
    "Drug_A description",
    DRUG_A_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A_metab1",
    "Drug A Metabolite 1 Concentration",
    "",
    "Drug_A_metab1 description",
    DRUG_A_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A_metab2",
    "Drug A Metabolite 2 Concentration",
    "",
    "Drug_A_metab2 description",
    DRUG_A_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B",
    "Drug B Concentration",
    "",
    "Drug_B description",
    DRUG_B_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab1",
    "Drug B Metabolite 1 Concentration",
    "",
    "Drug_B_metab1 description",
    DRUG_B_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab2",
    "Drug B Metabolite 2 Concentration",
    "",
    "Drug_B_metab2 description",
    DRUG_B_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"PI",
    "PI fixed fluorescence",
    "",
    "PI fluorescence description",
    PI_DEAD, false, 0, 1, 0, PROF_TYPE},

/*
// Intracellular profiles

    {"IC_MULTI",
    "IC Multi-constituent",
    "",
    "IC MULTI description",
    IC_MULTI, true, 0, 1, 0, PROF_TYPE},

    {"IC_Oxygen",
    "IC Oxygen Concentration",
    "",
    "IC Oxygen description",
    IC_OXYGEN, false, 0, 1, 0, PROF_TYPE},

    {"IC_Glucose",
    "IC Glucose Concentration",
    "",
    "IC Glucose description",
    IC_GLUCOSE, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A",
    "IC Drug A Concentration",
    "",
    "Drug_A description",
    IC_DRUG_A_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A_metab1",
    "IC Drug A Metabolite 1 Concentration",
    "",
    "Drug_A_metab1 description",
    IC_DRUG_A_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A_metab2",
    "IC Drug A Metabolite 2 Concentration",
    "",
    "Drug_A_metab2 description",
    IC_DRUG_A_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_B",
    "IC Drug B Concentration",
    "",
    "Drug_B description",
    IC_DRUG_B_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab1",
    "IC Drug B Metabolite 1 Concentration",
    "",
    "Drug_B_metab1 description",
    IC_DRUG_B_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab2",
    "IC Drug B Metabolite 2 Concentration",
    "",
    "Drug_B_metab2 description",
    IC_DRUG_B_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"IC_CFSE",
    "CFSE Concentration",
    "",
    "IC CFSE description",
    IC_CFSE, false, 0, 1, 0, PROF_TYPE},

    {"IC_growthrate",
    "Growth Rate",
    "",
    "Growth rate description",
    IC_GROWTH_RATE, false, 0, 1, 0, PROF_TYPE},

    {"IC_cellvolume",
    "Cell Volume",
    "",
    "Cell volume description",
    IC_CELL_VOLUME, false, 0, 1, 0, PROF_TYPE},

    {"IC_O2byvolume",
    "Cell O2xVolume",
    "",
    "Cell volume description",
    IC_O2_BY_VOL, false, 0, 1, 0, PROF_TYPE},
*/

};

    n_tsGraphs = sizeof(tsGraphSet)/sizeof(GRAPH_SET);
    tsGraphs = new GRAPH_SET[n_tsGraphs];
    for (int i=0; i<n_tsGraphs; i++) {
        tsGraphs[i] = tsGraphSet[i];
    }
    graphList = new GRAPH_SET[maxGraphs];
    nGraphs = maxGraphs;
}


GRAPH_SET Graphs::get_graph(int k)
{
	return graphList[k];
}

int Graphs::get_dataIndex(int k)
{
	return graphList[k].dataIndex;
}

QString Graphs::get_tag(int k)
{
	return graphList[k].tag;
}

QString Graphs::get_title(int k)
{
	return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
	return graphList[k].yAxisTitle;
}

QString Graphs::get_description(int k)
{
    return graphList[k].description;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

double Graphs::get_yscale(int k) {
    return graphList[k].yscale;
}

double Graphs::get_xscale(double xmax) {
    int n = 1;
    for (;;) {
        if (xmax <= n) break;
        n++;
    }
    return double(n);
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

int Graphs::get_type(int k) {
    return graphList[k].type;
}

bool Graphs::isTimeseries(int k)
{
    return (graphList[k].type == TS_TYPE);
}

bool Graphs::isProfile(int k)
{
    return (graphList[k].type == PROF_TYPE);
}

bool Graphs::isDistribution(int k)
{
    return (graphList[k].type == DIST_TYPE);
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList(int non_ts)
{
    int k = maxGraphs;
    int nts = 0;
    for (int i=0; i<n_tsGraphs; i++) {
        if (tsGraphs[i].active) {
            k--;
            graphList[k] = tsGraphs[i];
            nts++;
            if (nts == maxGraphs - non_ts) break;
        }
    }
    int ndummy = maxGraphs - nts - non_ts;
    for (k=0; k<ndummy; k++) {
        graphList[k].tag = "dummy";
        graphList[k].active = false;
        graphList[k].type = TS_TYPE;
        graphList[k].scaling = 1;
    }
    for (k=ndummy; k<ndummy + non_ts; k++) {
        graphList[k].tag = "non_ts";
        graphList[k].active = true;
        graphList[k].type = DIST_TYPE;  //????
        graphList[k].scaling = 1;
    }
    nGraphs = maxGraphs;

    char msg[128];
    sprintf(msg,"nGraphs: %d  non_ts: %d  nts: %d",nGraphs,non_ts,nts);
    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d scaling: %f",k,graphList[k].scaling);
//        LOG_MSG(msg);
//    }
}

