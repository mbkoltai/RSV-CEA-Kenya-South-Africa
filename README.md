## Estimating the cost-effectiveness of maternal vaccination and monoclonal antibodies for respiratory syncytial virus (RSV) in Kenya and South Africa

This repository contains R scripts and data tables to perform a cost-effectiveness analysis (CEA) of public health interventions against RSV disease in Kenya and South Africa for children under the age of 5 years, accompanying [this manuscript](https://www.medrxiv.org/content/10.1101/2022.09.09.22279780v2). 

The simulated interventions are based either on maternal vaccination (MV) or on monoclonal antibodies (mAb).
We used probabilistic sensitivity analysis of a static cohort model, which is a modified version of a previous CEA by [Li 2020](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01537-6).

The file to run the analysis and reproduce the results and figures in the manuscript is [reprod_figs.R](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/blob/master/reprod_figs.R). Data files are in the folder *custom_input*.

The code is based on the [previously developed R package McMarcel](https://zenodo.org/record/3663447), which was modified to include:
1) user-provided incidence data (treated and untreated ARIs and SARIs and deaths)
2) new efficacy data on vaccines and monoclonals, as well as an exponential decay model of the period of protection rather than an on-off model with a fixed duration
3) country-specific cost estimates for both inpatient and outpatient care

Most of the modifications were made in the file [get_burden_flexible_ari_sari](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/blob/master/functions/get_burden_flexible_ari_sari.R) which calculates the disease burden, the costs and the cost-effectiveness metrics following the interventions.

Incidence rates, efficacy estimates and some of the costs are generated from probability distributions that were inferred from the confidence intervals of the relevant data. Therefore the results are to some extent stochastic and will not be numerically identical to the figures reported in the paper, but should be very close.

**Outputs** are in the folders:  
1) the [default scenario](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/SA_ILI_broader_coverage_ANC_effic_betafit_2022): coverage levels from [Baral 2020](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0237718), efficacy figures with the most recent efficacy data, protection lasting for a fixed period (3 months for MV and 5 months for mAb).
2) [BCG coverage scenario](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/SA_ILI_broader_coverage_BCG_effic_betafit_2022): same as previous scenario, but higher coverage level.
2) [exponential waning model](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots/SA_ILI_broader_expwaning_coverage_ANC_effic_betafit_2022): instead of an on-off model of protection, the efficacy wanes exponentially with a time-average equal to the reported mean efficacy over the period of protection. Coverage as in scenario #1.

For all **figures in the main text** the data tables can be found in the [CEA plots](https://github.com/mbkoltai/RSV-CEA-Kenya-South-Africa/tree/master/output/cea_plots) folder, named __figure_\*.csv__.

For questions or problems please post an issue or send an email ([to the authors](https://www.lshtm.ac.uk/aboutus/people/koltai.mihaly)).
