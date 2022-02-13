## Estimating the cost-effectiveness of maternal vaccination and monoclonal antibodies for respiratory syncytial virus in Kenya and South Africa

These are files to perform cost-effectiveness analysis of public health interventions against RSV disease in Kenya and South Africa, using data provided by partners in these two countries. The intervention is either based on maternal vaccination (MV) or on monoclonal antibodies (mAB).

Data files are in the folder *custom_input*.  
**The main file to run the analysis to reproduce the results and figures in the manuscript is [run_cea_calc_parallel.R](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/blob/master/run_cea_calc_parallel.R).**

The code is based on the [previously developed R package McMarcel](https://zenodo.org/record/3663447), but it was modified to include:
1) user-provided incidence data
2) new efficacy data on vaccines and monoclonals, as well as an exponential decay model of the period of protection rather than an on-off model with a fixed duration
3) new cost estimates for both inpatient and outpatient care

Most of the modifications are in the file [get_burden_flexible_ari_sari](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/blob/master/functions/get_burden_flexible_ari_sari.R) which calculates the disease burden, costs and cost-effectiveness metrics following the interventions.

Incidence rates, efficacy estimates and some of the costs are generated from probability distributions that were inferred from the confidence intervals of the relevant data. Therefore the results are to some extent stochastic and will not be numerically identical to the figures reported in the paper, however they should be very close.

**Outputs** are in the folder [output/cea_plots](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots), the three folders corresponding to:  
1) the [default scenario](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/new_price_efficacy_KENdeaths_SAdeaths_CIs_SA_ILI_broader_effic_betafit): efficacy figures with published efficacy data, protection applies to a fixed period (3 months for MV and 5 months for mAb) then goes to zero.
2) [interim efficacy data for MV](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/new_price_efficacy_KENdeaths_SAdeaths_CIs_SA_ILI_broader_effic_betafit_interim): using unpublished data for MV (data presented by Pfizer at the [6th Resvinet conference](www.resvinet.org/6th-conference-2021.html))
3) [exponential waning model](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/new_price_efficacy_KENdeaths_SAdeaths_CIs_SA_ILI_broader_expwaning_effic_betafit): instead of an on-off model of protection, the efficacy wanes exponentially with a time-average equal to the reported mean efficacy over the period of protection


Contact: lshmk17 at lshtm.ac.uk
