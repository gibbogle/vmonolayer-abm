#include <qstring.h>
#include "params.h"
#include "QDebug"

Params::Params()
{
    static infoStruct label_info[] = {
        {"PARENT_0", "Diffusion coefficient within the blob."},
        {"PARENT_1", "Diffusion coefficient in the medium."},
        {"PARENT_2", "Cell influx parameter Kin.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_3", "Cell efflux parameter Kout.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_4", "Half-life of the compound, used to calculate the decay rate.  This is the same in the cell and in the medium."},
        {"PARENT_CT1_0", "Kmet0 is the maximum rate of metabolism.  The actual rate is the product of drug concentration Cdrug, Kmet0 and a sigmoid function of O2 concentration C_O2, with parameters C2 and KO2:\n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_1", "C2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate.  When C_O2 = 0, the function = 1, when C_O2 >> KO2, the function = 1 - C2: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_2",  "KO2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_3", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
         {"PARENT_CT1_4", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
        {"PARENT_CT1_5", "Klesion is currently unused."},
        {"PARENT_CT1_6", "The O2 concentration in the kill experiment."},
        {"PARENT_CT1_7", "The drug concentration in the kill experiment."},
        {"PARENT_CT1_8", "The duration the kill experiment."},
        {"PARENT_CT1_9", "The kill fraction achieved in the kill experiment (1 - SF)."},
        {"PARENT_CT1_10", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_11", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_12", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_13",  "n_O2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
 metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
 This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_14", "The death probability of a drug-tagged cell at time of division."},
        {"PARENT_CT1_15", "The kill probability rate parameter.\n\
 E.g. for kill model 1, kill probability rate r = Kd*dM/dt = Kd*kmet*Cdrug, and for duration t, SF = exp(-rt)"},
        {"PARENT_CT1_16", "This box is ticked if the drug is cytotoxic and kill parameters are provided."},
        {"PARENT_CT1_17", "Using dMdt = Cdrug*(1 - C2 + C2*KO2/(KO2 + C_O2))*Kmet0, the kill probability Pkill in time dt for each model is: \n\
1. Kd*dMdt*dt  2. Kd*Cdrug*dMdt*dt  3. Kd*dMdt^2*dt  4. Kd*Cdrug*dt  5. Kd*Cdrug^2*dt"},
        {"PARENT_CT1_18", "This box is ticked if the drug sensitises the cells to radiation."},
    };

    PARAM_SET params[] = {

{"GUI_VERSION_NAME", 0, 0, 0,
 "GUI0.00",
 "GUI program version number."},

{"DLL_VERSION_NAME", 0, 0, 0,
 "DLL0.00",
 "DLL version number."},

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"USE_LOGNORMAL_DIST", 0, 0, 1,
"Use lognormal distribution",
"The divide time will be a random variate from a log-normal distribution. \n\
 Otherwise checkpoint times are exponentially distributed and base phase times are fixed"},

{"DIVIDE_TIME_1_MEDIAN", 18, 0, 0,
"Median (h)",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_1_SHAPE", 1.1, 0, 0,
"Shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE_TIME_2_MEDIAN", 18, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_2_SHAPE", 1.1, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 0, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 20.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[mins]"},

{"NT_CONC", 1, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver."},

{"VCELL_PL", 1.0, 0, 0,
"Cell volume",
"Nominal cell volume."},

{"WELL_AREA", 0.33, 0, 0,
"Well area",
"Cross-sectional area of the well."},

{"MEDIUM_VOLUME", 0.2, 0, 0,
"Medium volume",
"Volume of the medium in which the spheroid is growing."},

{"FULLY_MIXED", 0, 0, 1,
"Medium is fully mixed?",
"The medium is fully mixed"},

//{"UNSTIRRED_LAYER", 0.01, 0, 0,
//"Unstirred layer width",
//"Thickness of the unstirred layer around the spheroid (cm)."},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.3, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_THRESHOLD", 0.15, 0, 0,
"Tag threshold",
"A cell begins to experience starvation of oxygen (anoxia) or glucose (aglucosia) leading to cell death at the oxygen/glucose concentration given by this threshold value."},

{"ANOXIA_TAG_TIME", 3.0, 0, 0,
"Tag time limit",
"Length of time under anoxia (O2 < anoxia threshold) or aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of anoxia or aglucosia."},

{"ANOXIA_DEATH_TIME", 3.0, 0, 0,
"Death delay time",
"Time taken for a cell to die after it is tagged to die of anoxia or aglucosia."},

{"AGLUCOSIA_THRESHOLD", 0.15, 0, 0,
"Aglucosia threshold",
"A cell begins to experience aglucosia leading to cell death at the glucose concentration given by this threshold value."},

{"AGLUCOSIA_TAG_TIME", 3.0, 0, 0,
"Aglucosia time limit",
"Length of time under aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of aglucosia.\n\
[h]"},

{"AGLUCOSIA_DEATH_TIME", 3.0, 0, 0,
"Aglucosia death delay time",
"Time taken for a cell to die after it is tagged to die of aglucosia.\n\
[h]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NCELLTYPES", 2, 0, 0,
"Number of cell types",
"Maximum number of cell types in the spheroid.  The initial percentage of each type must be specified"},

{"CELLPERCENT_1", 100, 0, 100,
"Percentage of cell type 1",
"Percentage of cell type 1"},

{"CELLPERCENT_2", 0, 0, 100,
"Percentage of cell type 2",
"Percentage of cell type 2"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SHOW_PROGENY", 0, 0, 0,
 "Show descendants of cell #",
 "All the descendants of cell with the specified ID are highlighted.  (0 = no selection)"},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_GROWTH", 1, 0, 1,
"Oxygen growth?",
"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

{"OXYGEN_DEATH", 1, 0, 1,
"Anoxia death?",
"Oxygen controls death by anoxia"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Spheroid diffusion coeff",
 "Constituent diffusion coefficient in the spheroid"},

{"OXYGEN_MEDIUM_DIFF", 5.0e-5, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"OXYGEN_CELL_DIFF_IN", 600, 0, 0,
 "Cell influx parameter Kin",
 "Cell membrane diffusion constant Kin"},

{"OXYGEN_CELL_DIFF_OUT", 600, 0, 0,
 "Cell efflux parameter Kout",
 "Cell membrane diffusion constant Kout"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"OXYGEN_CONSUMPTION", 6.25e-17, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"OXYGEN_HILL_N", 1, 1, 2,
 "Hill function N",
 "Oxygen uptake rate Hill function N"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_GROWTH", 1, 0, 1,
"Glucose growth?",
"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

{"GLUCOSE_DEATH", 1, 0, 1,
"Aglucosia death?",
"Glucose controls death by aglucosia"},

{"GLUCOSE_DIFF_COEF", 3.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_MEDIUM_DIFF", 8.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUCOSE_CELL_DIFF_IN", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kin"},

{"GLUCOSE_CELL_DIFF_OUT", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kout"},

{"GLUCOSE_BDRY_CONC", 5.5, 0, 0,
 "Boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUCOSE_CONSUMPTION", 1.0e-16, 0, 0,
 "Max consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 400, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUCOSE_HILL_N", 1, 1, 2,
 "Hill function N",
 "Glucose uptake rate Hill function N"},

{"USE_LACTATE", 1, 0, 1,
"Use Lactate?",
"Lactate is simulated"},

//    {"LACTATE_GROWTH", 1, 0, 1,
//    "Lactate growth?",
//    "The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose, or something"},

    {"LACTATE_DIFF_COEF", 3.0e-7, 0, 0,
     "Spheroid diffusion coeff",
     "LACTATE diffusion coefficient"},

    {"LACTATE_MEDIUM_DIFF", 6.0e-6, 0, 0,
     "Medium diffusion coeff",
     "Constituent diffusion coefficient in the medium"},

    {"LACTATE_CELL_DIFF_IN", 400, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kin"},

    {"LACTATE_CELL_DIFF_OUT", 400, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kout"},

    {"LACTATE_BDRY_CONC", 0.4, 0, 0,
     "Boundary concentration",
     "LACTATE boundary concentration"},

    {"LACTATE_CONSUMPTION", 3.8e-17, 0, 0,
     "Max consumption rate",
     "LACTATE consumption rate"},

    {"LACTATE_MM_KM", 20, 0, 0,
     "Michaelis-Menten Km",
     "Michaelis-Menten Km (uM)"},

    {"LACTATE_HILL_N", 1, 1, 2,
     "Hill function N",
     "Lactate uptake rate Hill function N"},

{"USE_GLUTAMINE", 1, 0, 1,
"Use Glutamine?",
"Glutamine is simulated"},

{"GLUTAMINE_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "Glutamine diffusion coefficient"},

{"GLUTAMINE_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUTAMINE_CELL_DIFF_IN", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kin"},

{"GLUTAMINE_CELL_DIFF_OUT", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kout"},

{"GLUTAMINE_BDRY_CONC", 1.0, 0, 0,
 "Boundary concentration",
 "Glutamine boundary concentration"},

{"GLUTAMINE_CONSTANT", 1, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUTAMINE_CONSUMPTION", 1.0e-16, 0, 0,
 "Consumption rate",
 "Glutamine consumption rate"},

{"GLUTAMINE_MM_KM", 100, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUTAMINE_HILL_N", 1, 0, 2,
 "Hill function N",
 "Glutamine uptake rate Hill function N"},

    {"USE_OTHER", 1, 0, 1,
    "Use Othernutrient?",
    "Other nutrient is simulated"},

    {"OTHER_DIFF_COEF", 6.0e-7, 0, 0,
     "Spheroid diffusion coeff",
     "Other nutrient diffusion coefficient"},

    {"OTHER_MEDIUM_DIFF", 6.0e-6, 0, 0,
     "Medium diffusion coeff",
     "Constituent diffusion coefficient in the medium"},

    {"OTHER_CELL_DIFF_IN", 100, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kin"},

    {"OTHER_CELL_DIFF_OUT", 100, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kout"},

    {"OTHER_BDRY_CONC", 6.0, 0, 0,
     "Boundary concentration",
     "Other nutrient boundary concentration"},

    {"OTHER_CONSTANT", 1, 0, 1,
     "Constant concentration",
     "Extracellular concentration to be held constant everywhere at the specified boundary value"},

    {"OTHER_CONSUMPTION", 1.0e-16, 0, 0,
     "Consumption rate",
     "Other nutrient consumption rate"},

    {"OTHER_MM_KM", 100, 0, 0,
     "Michaelis-Menten Km",
     "Michaelis-Menten Km (uM)"},

    {"OTHER_HILL_N", 1, 0, 2,
     "Hill function N",
     "Other nutrient uptake rate Hill function N"},

//==========================
// Radiotherapy parameters
//==========================

{"RADIATION_ALPHA_H_1", 0.0738, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_1", 0.00725, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_1", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_1", 2.5, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_1", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_1", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"RADIATION_GROWTH_DELAY_FACTOR_1", 0.0, 0, 0,
"Growth delay factor",
"For a damaged cell, cell cycle is delayed for a number of hours given by this factor x radiation dose"},

{"RADIATION_GROWTH_DELAY_N_1", 0, 0, 0,
"Growth delay cycles",
"For a damaged cell, cell cycle delay persists for a number of cell cycles"},

{"RADIATION_ALPHA_H_2", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_2", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_2", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_2", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_2", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_2", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"RADIATION_GROWTH_DELAY_FACTOR_2", 0.0, 0, 0,
"Growth delay factor",
"For a damaged cell, cell cycle is delayed for a number of hours given by this factor x radiation dose"},

{"RADIATION_GROWTH_DELAY_N_2", 0, 0, 0,
"Growth delay cycles",
"For a damaged cell, cell cycle delay persists for a number of cell cycles"},

{"RADIATION_GROWTH_DELAY_ALL", 0, 0, 0,
"Delay growth of all cells",
"Cell cycle delay is also applied to cells that are not fated to die"},

    {"USE_CELL_CYCLE", 1,0,1,
     "Use cell cycle with G1, S, G2, M phases",
     "Cell cycle parameters determine the time spent in each phase.\n\
     In the case of G1 and G2, an exponentially distributed random delay is added"},

     {"USE_SYNCHRONISE", 0,0,1,
     "Synchronise cell cycles?",
     "Synchronise initial cell phases to start of M phase"},

     {"T_G1_1", 6, 0, 0,
     "G1 phase base time (h)",
     "Deterministic component of time spent in phase G1"},

     {"T_S_1", 8, 0, 0,
     "S phase base time (h)",
     "Time spent in phase S"},

     {"T_G2_1", 1, 0, 0,
     "G2 phase base time (h)",
     "Deterministic component of time spent in phase G2"},

     {"T_M_1", 0.5, 0, 0,
     "M phase base time (h)",
     "Time spent in phase M"},

     {"G1_MEAN_DELAY_1", 1.5, 0, 0,
     "G1 mean delay (h)",
     "Mean of the random component of time spent in phase G1 (exponentially distributed)"},

     {"S_MEAN_DELAY_1", 1.0, 0, 0,
     "S mean delay (h)",
     "Mean of the random component of time spent in phase S (exponentially distributed)"},

     {"G2_MEAN_DELAY_1", 1, 0, 0,
     "G2 mean delay (h)",
     "Mean of the random component of time spent in phase G2 (exponentially distributed)"},

     {"APOPTOSIS_RATE_1", 0.1, 0, 0,
     "Apoptosis rate/hr",
     "The rate of passage from state = DYING to state = DEAD, probability/hour"},

     {"ARREST_THRESHOLD_1", 0.1, 0, 0,
     "Arrest threshold",
     "Threshold level of fraction of maximum growth rate that causes arrest in S phase"},


     {"RMR_ETA_PL_1", 35, 0, 0,
     "PL lesion creation rate",
     "Coefficient of rate of creation of potentially lethal lesions: eta_PL"},

     {"RMR_ETA_IRL_1", 0.003, 0, 0,
     "IRL creation rate",
     "Coefficient of rate of creation of irrepairable lesions: eta_IRL"},

     {"RMR_KREP_BASE_1", 0.2, 0, 0,
     "Base PL lesion repair rate",
     "Base coefficient of rate of repair of potentially lethal lesions: Krepair_base (Curtis's epsilon_PL)\n\
     The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

     {"RMR_KREP_MAX_1", 2.0, 0, 0,
     "Max PL lesion repair rate",
     "Maximum coefficient of rate of repair of potentially lethal lesions: Krepair_max (Curtis's epsilon_PL)\n\
     The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

     {"RMR_KMIS_1", 0.00003, 0, 0,
     "Lesion misrepair rate",
     "Coefficient of rate of misrepair of PL lesions to Ch1 and Ch2 lesions: Kmisrepair"},

     {"RMR_F_MITOSIS_1", 4, 0, 0,
     "Mitosis factor",
     "Multiplying factor of misrepair during mitosis: mitosis_factor"},

     {"RMR_F_CH1_1", 0.5, 0, 0,
     "Ch1 fraction of misrepair",
     "Fraction of misrepaired lesions that are Ch1: fraction_Ch1"},

     {"RMR_PSURVIVE1_1", 0.45, 0, 0,
     "Ch1 prob of mitosis survival",
     "Probability that a cell with one Ch1 lesion survives mitosis: psurvive_Ch1"},

     {"RMR_PSURVIVE2_1", 0.05, 0, 0,
     "Ch2 prob of mitosis survival",
     "Probability that a cell with one Ch2 lesion survives mitosis: psurvive_Ch2"},

     {"RMR_ATCP_1", 50, 0, 0,
     "TCP Hill a",
     "TCP, the maximum time spent in the checkpoint, is a Hill function of # of L1 lesions, n. TCP(n) = b.n/(a + n)"},

      {"RMR_BTCP_1", 25, 0, 0,
      "TCP Hill b",
      "TCP, the maximum time spent in the checkpoint, is a Hill function of # of L1 lesions, n. TCP(n) = b.n/(a + n)"},

//     {"RMR_KCP_1", 0.13, 0, 0,
//     "Checkpoint time limit factor",
//     "Factor for computing maximum time spent in the checkpoint, as a function of # of L1 lesions: Kcp"},

      {"T_G1_2", 6, 0, 0,
      "G1 phase base time (h)",
      "Deterministic component of time spent in phase G1"},

      {"T_S_2", 8, 0, 0,
      "S phase base time (h)",
      "Time spent in phase S"},

      {"T_G2_2", 1, 0, 0,
      "G2 phase base time (h)",
      "Deterministic component of time spent in phase G2"},

      {"T_M_2", 0.5, 0, 0,
      "M phase base time (h)",
      "Time spent in phase M"},

      {"G1_MEAN_DELAY_2", 1.5, 0, 0,
      "G1 mean delay (h)",
      "Mean of the random component of time spent in phase G1 (exponentially distributed)"},

      {"S_MEAN_DELAY_2", 1.0, 0, 0,
      "S mean delay (h)",
      "Mean of the random component of time spent in phase S (exponentially distributed)"},

      {"G2_MEAN_DELAY_2", 1, 0, 0,
      "G2 mean delay (h)",
      "Mean of the random component of time spent in phase G2 (exponentially distributed)"},

      {"APOPTOSIS_RATE_2", 0.1, 0, 0,
      "Apoptosis rate/hr",
      "The rate of passage from state = DYING to state = DEAD, probability/hour"},

      {"ARREST_THRESHOLD_2", 0.1, 0, 0,
      "Arrest threshold",
      "Threshold level of fraction of maximum growth rate that causes arrest in S phase"},


      {"RMR_ETA_PL_2", 35, 0, 0,
      "PL lesion creation rate",
      "Coefficient of rate of creation of potentially lethal lesions: eta_PL"},

      {"RMR_ETA_IRL_2", 0.003, 0, 0,
      "IRL creation rate",
      "Coefficient of rate of creation of irrepairable lesions: eta_IRL"},

      {"RMR_KREP_BASE_2", 0.2, 0, 0,
      "Base PL lesion repair rate",
      "Base coefficient of rate of repair of potentially lethal lesions: Krepair_base (Curtis's epsilon_PL)\n\
      The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

      {"RMR_KREP_MAX_2", 2.0, 0, 0,
      "Max PL lesion repair rate",
      "Maximum coefficient of rate of repair of potentially lethal lesions: Krepair_max (Curtis's epsilon_PL)\n\
      The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

      {"RMR_KMIS_2", 0.00003, 0, 0,
      "Lesion misrepair rate",
      "Coefficient of rate of misrepair of PL lesions to Ch1 and Ch2 lesions: Kmisrepair"},

      {"RMR_F_MITOSIS_2", 4, 0, 0,
      "Mitosis factor",
      "Multiplying factor of misrepair during mitosis: mitosis_factor"},

      {"RMR_F_CH1_2", 0.5, 0, 0,
      "Ch1 fraction of misrepair",
      "Fraction of misrepaired lesions that are Ch1: fraction_Ch1"},

      {"RMR_PSURVIVE1_2", 0.45, 0, 0,
      "Ch1 prob of mitosis survival",
      "Probability that a cell with one Ch1 lesion survives mitosis: psurvive_Ch1"},

      {"RMR_PSURVIVE2_2", 0.05, 0, 0,
      "Ch2 prob of mitosis survival",
      "Probability that a cell with one Ch2 lesion survives mitosis: psurvive_Ch2"},

       {"RMR_ATCP_2", 50, 0, 0,
       "TCP Hill a",
       "TCP, the maximum time spent in the checkpoint, is a Hill function of # of L1 lesions, n. TCP(n) = b.n/(a + n)"},

       {"RMR_BTCP_2", 25, 0, 0,
       "TCP Hill b",
       "TCP, the maximum time spent in the checkpoint, is a Hill function of # of L1 lesions, n. TCP(n) = b.n/(a + n)"},


//      {"RMR_KCP_2", 0.13, 0, 0,
//      "Checkpoint time limit factor",
//      "Factor for computing maximum time spent in the checkpoint, as a function of # of L1 lesions: Kcp"},

// Metabolism parameters

      {"USE_METABOLISM", 1,0,1,
       "Use glucose metabolism",
       "If metabolism is turned on, cell cycle is used, and lactate is simulated.  If metabolism is turned off, lactate is not simulated"},

       {"F_GNORM_1", 0.5, 0, 0,
       "Normal fraction of glycolysis -> intermediates",
       "Fraction of glycolysis (r_G) going to make intermediates under conditions of full nutrition"},

       {"F_PNORM_1", 0.5, 0, 0,
       "Normal fraction of pyruvate -> intermediates",
        "Fraction of pyruvate (r_P) going to make intermediates under conditions of full nutrition"},

       {"F_GLNNORM_1", 0.5, 0, 0,
       "Normal fraction of glutamine -> intermediates",
        "Fraction of glutamine (r_Gn) going to make intermediates under conditions of full nutrition"},

      {"N_GA_1", 2, 0, 0,
      "ATP moles produced per glucose mole",
      "Number of ATP moles produced by the glycolysis of one glucose mole"},

      {"N_PA_1", 14, 0, 0,
      "ATP moles produced per pyruvate mole",
      "Number of ATP moles produced by the oxidation of one pyruvate mole"},

       {"N_GLNA_1", 8, 0, 0,
       "ATP moles produced per glutamine mole",
       "Number of ATP moles produced by the oxidation of one glutamine mole"},

       {"N_GP_1", 2, 0, 0,
       "Pyruvate moles produced per glucose mole",
       "Number of pyruvate moles produced by the glycolysis of one glucose mole (no intermediates produced)"},

      {"N_GI_1", 1, 0, 0,
      "Intermediate moles produced per glucose mole",
      "Number of moles of anabolic intermediates produced the glycolysis of one glucose mole"},

      {"N_PI_1", 1, 0, 0,
      "Intermediate moles produced per pyruvate mole",
      "Number of moles of anabolic intermediates produced the oxidation of one pyruvate mole"},

       {"N_GLNI_1", 1, 0, 0,
       "Intermediate moles produced per glutamine mole",
       "Number of moles of anabolic intermediates produced the oxidation of one glutamine mole"},

       {"N_ONI_1", 1.5, 0, 0,
       "Intermediate moles produced per othernutrient mole",
       "Number of moles of anabolic intermediates produced the oxidation of one other nutrient mole"},

      {"N_PO_1", 1, 0, 0,
      "Oxygen moles consumed per pyruvate mole",
      "Number of moles of oxygen consumed the oxidation of one pyruvate mole"},

       {"N_GLNO_1", 1, 0, 0,
       "Oxygen moles consumed per glutamine mole",
       "Number of moles of oxygen consumed the oxidation of one glutamine mole"},

      {"K_H1_1", 8, 0, 0,
      "Exponent N_H in Hss equation",
      "With x = C_O/C_O_max, steady-state HIF-1 level is Hss = (1-x)^N_H, where C_O = oxygen concentration"},

      {"K_H2_1", 0.00005, 0, 0,
      "HIF-1 rate constant K_H",
      "The rate of change of HIF-1 level H is given by: dH/dt = K_H*(Hss - H), where C_O = oxygen concentration, Hss is steady-state level for C_O"},

      {"K_HB_1", 3, 0, 0,
      "K_HB",
      "Glycolysis rate = K_HA*(1 + K_HB*H)*C_G^N/(C_G^N + Km^N)\n\
        where: H = HIF-1 level, C_G = glucose concentration, K_HA is maximum glucose consumption rate when H=0, Km and N are the glucose consumption Hill function parameters"},

      {"K_PDK_1", 4.63e-5, 0, 0,
      "K_PDK",
      "Representing PDK1 factor level by y in the range (PDKmin,1), the rate of change of y is: dy/dt = -K_PDK*(y - 1 + cH) where H = HIF-1 level, and c = 1 - PDKmin"},

       {"PDKMIN_1", 0.3, 0, 0,
       "PDKmin",
       "Minimum value of the PDK1 factor.\n\
       Representing PDK1 factor level by y in the range (PDKmin,1), the rate of change of y is: dy/dt = -K_PDK*(y - 1 + cH) where H = HIF-1 level, and c = 1 - PDKmin"},

       {"C_O2_NORM_1", 0.1, 0, 0,
       "Nominal normal IC O2 concentration",
       "Nominal normal IC O2 concentration, used to set normal metabolic rates for unconstrained growth"},

       {"C_G_NORM_1", 1.0, 0, 0,
       "Nominal normal IC glucose concentration",
       "Nominal normal IC glucose concentration, used to set normal metabolic rates for unconstrained growth"},

       {"C_L_NORM_1", 0.4, 0, 0,
       "Nominal normal IC lactate concentration",
       "Nominal normal IC lactate concentration, used to set normal metabolic rates for unconstrained growth"},

        {"C_GLN_NORM_1", 1.0, 0, 0,
        "Nominal normal IC glutamine concentration",
        "Nominal normal IC glutamine concentration, used to set normal metabolic rates for unconstrained growth"},

        {"C_ON_NORM_1", 6.0, 0, 0,
        "Nominal normal IC othernutrient concentration",
        "Nominal normal IC other nutrient concentration, used to set normal metabolic rates for unconstrained growth"},

      {"ATP_S_1", 0.3, 0, 0,
      "ATP production threshold for survival (fraction of peak)",
      "Cell death occurs when the ATP production rate falls below the fraction ATP_S of the maximum (no nutrient constraints) production rate"},

      {"ATP_G_1", 0.5, 0, 0,
       "ATP production threshold for growth (fraction of peak)",
       "Cell growth stops when the ATP production rate falls below the fraction ATP_G of the maximum (no nutrient constraints) production rate.\n\
        Intermediates production from glycolysis and pyruvate ramps down when ATP rate is below ATPramp*ATP_G, to 0 when ATP rate = ATP_G."},

      {"ATP_RAMP_1", 1.3, 0, 0,
        "Ramp factor for reducing r_G, r_P based on ATP",
       "Fractional rates of production of intermediates from glycolysis (r_G) and pyruvate (r_P) are reduced as ATP production rate falls.\n\
        If the normal ATP rate is r_An, and the current rate is r_A, then the fractional production rate f_A = r_A/r_An.\n\
        A multiplying factor on both r_G and r_P, f, is calculated as a ramp from 0 - 1 as f_A varies from ATP_G - ATPramp*ATP_G\n\
        i.e. the reduction starts when f_A = ATPramp*ATP_G (f = 1), and when f_A = ATP_G, f = 0.  f_G = f*f_Gn, f_P = f*f_Pn"},

       {"K_PL_1", 0.1, 0, 0,
       "Pyruvate -> lactate rate constant",
       "The forward rate constant K_PL of the pyruvate-lactate reaction, i.e. the rate constant for conversion of pyruvate to lactate"},

      {"K_LP_1", 0.1, 0, 0,
      "Lactate -> pyruvate rate constant",
      "The reverse rate constant K_LP of the pyruvate-lactate reaction, i.e. the rate constant for conversion of lactate to pyruvate"},

       {"PYRUVATE_MM_KM_1", 20, 0, 0,
       "Pyruvate Michaelis-Menten Km (uM)",
       "Pyruvate Michaelis-Menten Km (uM). fMM = C_P/(Km + C_P) is a multiple of pyruvate oxidation rate, ensuring that the rate -> 0 as C_P -> 0"},

       {"PASTEUR_MM_KM_1", 0.02, 0, 0,
       "Pasteur Effect Michaelis-Menten Km (normalised)",
       "The glycolysis rate is multiplied by cfactor, which is a function of the normalised total rate of oxygen consumption.\n\
        Normalised O2 rate = r_O2N = r_O2/r_O2max, cfactor = r_O2N/(Km + r_O2N)"},

       {"GLUTAMINE_BASERATE_1", 0, 0, 0,
       "Glutamine base consumption rate",
       "To account for a possible base rate of metabolism of glutamine.  The rate of consumption is made of a fixed base rate, r_base, + a rate that depends on glucose conc C_G.\n\
       With fbase = r_base/Vmax, and gfactor = fbase + (1 - fbase)**Km_GG/(Km_GG + C_G), where Km_GG = 0.2*Km_G,\n\
       the glutamine rate becomes: r_Gln = fPDK*Vmax*gfactor*C_Gln/(Km_Gln + C_Gln)"},

       {"F_N_1", 0.2, 0, 0,
       "Glutamine intermediates rate factor",
        "The factor f_N is the ratio of the rate of intermediates production from glutamine to the total intermediates production rate."},

       {"F_GP_SOLVER",1,0,0,
        "f_GP solver (1,2,3)",
        "Choices for f_GP solver: 1 = max ATP (tandem), 2 = max ATP (staged), 3 = survival (staged)"},

{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"HYPOXIA_THRESHOLD", 4.0, 0, 0,
"Hypoxia threshold",
"Hypoxia threshold"},

{"GROWTH_FRACTION_1", 0.25, 0, 0,
"Growth fraction threshold 1",
"Growth fraction threshold 1"},

{"GROWTH_FRACTION_2", 0.1, 0, 0,
"Growth fraction threshold 2",
"Growth fraction threshold 2"},

{"GROWTH_FRACTION_3", 0.01, 0, 0,
"Growth fraction threshold 3",
"Growth fraction threshold 3"},

{"DRUG_THRESHOLD", 1.0e-6, 0, 0,
 "Drug Threshold",
 "Threshold drug concentration - when all intracellular and extracellular concentrations fall below this level, the drug concentrations everywhere are set to zero"},

{"DRUG_LABEL_THRESHOLD", 0, 0, 0,
"Label Threshold",
"Threshold label-drug concentration - when a labelling drug (e.g. EDU) is used, this is the threshold for a cell to be considered as labelled"},

{"SPCRAD", 200.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

 {"SAVE_FACS_DATA",0,0,1,
  "Save FACS data",
  "Save data for FACS at a specified interval"},

 {"SAVE_FACS_DATA_FILE_NAME",0,0,0,
  "facs_data",
  "Base file name for saving FACS data"},

 {"SAVE_FACS_DATA_TSTART",0,0,0,
  "Start time",
  "Start time for saving FACS data"},

 {"SAVE_FACS_DATA_INTERVAL",0,0,0,
  "Interval",
  "Time interval for saving FACS data"},

 {"SAVE_FACS_DATA_NUMBER",1,0,0,
  "Number",
  "Number of times to save FACS data"},


// This is the end of the parameters that are actually read by the DLL
// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted,
// followed by the list of time-series and profile plots selected for this run.

{"DUMMY_HYPOXIA_THRESHOLD", 0, 0, 0,
"Hypoxia threshold",
"Select the intracellular O2 level below which the cell is counted as hypoxic"},

{"DUMMY_GROWTH_FRACTION", 0, 0, 0,
"Growth fraction",
"Select the threshold fraction of average growth rate (i.e. with no nutrient limits) used to count slow-growing cells"},

// Time-series plots
    {"nlive",                     1, 0,1,"","Number of live cells"},
    {"nviable",                   1, 0,1,"","Number of viable cells"},
    {"nonviable",                 1, 0,1,"","Total number of non-viable cells"},
    {"nATPdead",                  1, 0,1,"","Total number of cells that have been killed by low ATP"},
    {"ndrugAdead",                0, 0,1,"","Total number of cells that have been killed by drugA"},
    {"ndrugBdead",                0, 0,1,"","Total number of cells that have been killed by drugB"},
    {"nradiationdead",            0, 0,1,"","Total number of cells that have been killed by radiation"},
    {"ndead",                     1, 0,1,"","Total number of cellls that have died"},
    {"nATPtagged",                1, 0,1,"","Current number of cells tagged to die by lowATP"},
    {"ndrugAtagged",              0, 0,1,"","Current number of cells tagged to die by drugA"},
    {"ndrugBtagged",              0, 0,1,"","Current number of cells tagged to die by drugB"},
    {"nradiationtagged",          0, 0,1,"","Current number of cells tagged to die by radiation"},
    {"viablefraction",            1, 0,1,"","Fraction of cells that are viable"},
    {"hypoxicfraction",           1, 0,1,"","Fraction of cells with oxygen level below the specified threshold for hypoxia"},
    {"clonohypoxicfraction",      0, 0,1,"","Fraction of clonogenic cells with oxygen level below the specified threshold for hypoxia"},
    {"growthfraction",            0, 0,1,"","Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits"},
    {"nogrowfraction",            1, 0,1,"","Percentage of cells that are not growing (insufficient ATP rate for growth)"},
    {"clonofraction",             0, 0,1,"","Percentage of cells that are clonogenic (will give a colony >= 50)"},
    {"platingefficiency",         1, 0,1,"","Percentage of live cells that are viable"},
    {"ECoxygen",                  1, 0,1,"","EC concentration of oxygen in the medium (bottom)"},
    {"ECglucose",                 1, 0,1,"","EC concentration of glucose in the medium (bottom)"},
    {"EClactate",                 1, 0,1,"","EC concentration of lactate in the medium (bottom)"},
    {"ECglutamine",               1, 0,1,"","EC concentration of glutamine in the medium (bottom)"},
    {"ECother",                   0, 0,1,"","EC concentration of other nutrient in the medium (bottom)"},
    {"ECdrugA",                   0, 0,1,"","EC concentration of drug A in the medium (bottom)"},
    {"ECdrugAmet1",               0, 0,1,"","EC concentration of drug A metabolite 1 in the medium (bottom)"},
    {"ECdrugAmet2",               0, 0,1,"","EC concentration of drug A metabolite 2 in the medium (bottom)"},
    {"ECdrugB",                   0, 0,1,"","EC concentration of drug B in the medium (bottom)"},
    {"ECdrugBmet1",               0, 0,1,"","EC concentration of drug B metabolite 1 in the medium (bottom)"},
    {"ECdrugBmet2",               0, 0,1,"","EC concentration of drug B metabolite 2 in the medium (bottom)"},
    {"ICoxygen",                  1, 0,1,"","IC concentration of oxygen"},
    {"ICglucose",                 1, 0,1,"","IC concentration of glucose"},
    {"IClactate",                 1, 0,1,"","IC concentration of lactate"},
    {"ICglutamine",               1, 0,1,"","IC concentration of glutamine"},
    {"ICother",                   0, 0,1,"","IC concentration of other nutrient"},
    {"ICpyruvate",                1, 0,1,"","IC concentration of pyruvate"},
    {"ICATP",                     1, 0,1,"","IC concentration of ATP"},
    {"ICdrugA",                   0, 0,1,"","IC concentration of drug A"},
    {"ICdrugAmet1",               0, 0,1,"","IC concentration of drug A metabolite 1"},
    {"ICdrugAmet2",               0, 0,1,"","IC concentration of drug A metabolite 2"},
    {"ICdrugB",                   0, 0,1,"","IC concentration of drug B"},
    {"ICdrugBmet1",               0, 0,1,"","IC concentration of drug B metabolite 1"},
    {"ICdrugBmet2",               0, 0,1,"","IC concentration of drug B metabolite 2"},
    {"Medoxygen",                 0, 0,1,"","Average medium concentration of oxygen"},
    {"Medglucose",                0, 0,1,"","Average medium concentration of glucose"},
    {"Medlactate",                0, 0,1,"","Average medium concentration of lactate"},
    {"Medglutamine",              0, 0,1,"","Average medium concentration of glutamine"},
    {"Medother",                  0, 0,1,"","Average medium concentration of other nutrient"},
    {"MeddrugA",                  0, 0,1,"","Average medium concentration of drug A"},
    {"MeddrugAmet1",              0, 0,1,"","Average medium concentration of drug A metabolite 1"},
    {"MeddrugAmet2",              0, 0,1,"","Average medium concentration of drug A metabolite 2"},
    {"MeddrugB",                  0, 0,1,"","Average medium concentration of drug B"},
    {"MeddrugBmet1",              0, 0,1,"","Average medium concentration of drug B metabolite 1"},
    {"MeddrugBmet2",              0, 0,1,"","Average medium concentration of drug B metabolite 2"},
    {"doublingtime",              0, 0,1,"","Average doubling time"},
    {"Grate",                     1, 0,1,"","Normalised glycolysis rate"},
    {"Prate",                     0, 0,1,"","Normalised pyruvate utilisation rate"},
    {"Glnrate",                   0, 0,1,"","Normalised glutamine utilisation rate"},
    {"ONrate",                    0, 0,1,"","Normalised other nutrient utilisation rate"},
    {"Arate",                     1, 0,1,"","Normalised ATP production rate"},
    {"Irate",                     1, 0,1,"","Normalised rate of production of anabolic intermediates"},
    {"f_G",                       0, 0,1,"","f_G"},
    {"f_P",                       0, 0,1,"","f_P"},
    {"HIF-1",                     1, 0,1,"","HIF-1"},
    {"PDK1",                      1, 0,1,"","PDK1"},
    {"dividerate",                0, 0,1,"","# divided/hour"},
    {"G1_phase",                  0, 0,1,"","G1_phase"},
    {"G1_cp",                     0, 0,1,"","G1_cp"},
    {"S_phase",                   0, 0,1,"","S_phase"},
    {"G2_phase",                  0, 0,1,"","G2_phase"},
    {"G2_cp",                     0, 0,1,"","G2_cp"},
    {"M_phase",                   0, 0,1,"","M_phase"},
    {"S_phase_nonarrest",         0, 0,1,"","S_phase_nonarrest"},
    {"nmutations",                0, 0,1,"","nmutations"},


// Profile plots
    {"MULTI",                     1, 0,1,"","Selected constituent on a line through the blob centre"},
//    {"CFSE",                      0, 0,1,"","Extracellular CFSE concentration on a line through the blob centre"},
    {"Oxygen",                    0, 0,1,"","Extracellular oxygen concentration on a line through the blob centre"},
    {"Glucose",                   0, 0,1,"","Extracellular glucose concentration on a line through the blob centre"},
    {"Glutamine",                 0, 0,1,"","Extracellular glutamine concentration on a line through the blob centre"},
    {"Drug_A",                    0, 0,1,"","Extracellular drug A concentration on a line through the blob centre"},
    {"Drug_A_metab1",             0, 0,1,"","Extracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"Drug_A_metab2",             0, 0,1,"","Extracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"Drug_B",                    0, 0,1,"","Extracellular drug Bconcentration on a line through the blob centre"},
    {"Drug_B_metab1",             0, 0,1,"","Extracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"Drug_B_metab2",             0, 0,1,"","Extracellular drug B metabolite 2 concentration on a line through the blob centre"},
/*
    {"IC_MULTI",                  1, 0,1,"","Selected constituent on a line through the blob centre"},
    {"IC_Oxygen",                 0, 0,1,"","Intracellular oxygen concentration on a line through the blob centre"},
    {"IC_Glucose",                0, 0,1,"","Intracellular glucose concentration on a line through the blob centre"},
    {"IC_Glutamine",              0, 0,1,"","Intracellular glutamine concentration on a line through the blob centre"},
    {"IC_Drug_A",                 0, 0,1,"","Intracellular drug A concentration on a line through the blob centre"},
    {"IC_Drug_A_metab1",          0, 0,1,"","Intracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_A_metab2",          0, 0,1,"","Intracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"IC_Drug_B",                 0, 0,1,"","Intracellular drug Bconcentration on a line through the blob centre"},
    {"IC_Drug_B_metab1",          0, 0,1,"","Intracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_B_metab2",          0, 0,1,"","Intracellular drug B metabolite 2 concentration on a line through the blob centre"},
    {"IC_CFSE",                   0, 0,1,"","CFSE concentration on a line through the blob centre"},
    {"IC_growthrate",             0, 0,1,"","Cell growth rate on a line through the blob centre"},
    {"IC_cellvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
    {"IC_O2byvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
// Distribution plots
//    {"Oxygen",                    0, 0,1,"","Probability distribution of extracellular oxygen concentration"},
//    {"cellvolume",                0, 0,1,"","Probability distribution of cell volume fraction"}
*/

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}

    nInfolabel = sizeof(label_info)/sizeof(INFOSTRUCT);
    workingInfolabelList = new INFOSTRUCT[nInfolabel];
    for (int i=0; i<nInfolabel; i++) {
        workingInfolabelList[i] = label_info[i];
    }
    /*
    nInfocheckbox = sizeof(checkbox_info)/sizeof(INFOSTRUCT);
    workingInfocheckboxList = new INFOSTRUCT[nInfocheckbox];
    for (int i=0; i<nInfocheckbox; i++) {
        workingInfocheckboxList[i] = checkbox_info[i];
    }
    */
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}

void Params::get_labeltag(int i, QString *tag)
{
    *tag = workingInfolabelList[i].tag;
}

void Params::infoLabelInfo(QString tag, QString *info)
{
    for (int i=0; i<nInfolabel; i++) {
        if (tag == workingInfolabelList[i].tag) {
            *info = workingInfolabelList[i].info;
            return;
        } else {
            *info = "";
        }
    }
}

