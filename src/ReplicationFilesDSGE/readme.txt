Before doing anything:
Must run split_join_clouds.jl

Default option is to use pre-saved estimation results to compute impulse responses/forecast/histograms and produce the figures. If you would like to change this, see instructions at bottom.

Replicating Figure 6.1 (IRFs with separate estimations half1 and half2):
 - Run specfiles_irfs/dsgevarifs_hrs_pi_ls_wpi_spec1002_10_2019Q4_1118.jl (this file computes IRFs in parallel and will require you to customize the procedure for adding workers to your personal computer/cluster (see line 106))
 - The IRF pdfs are saved in save_orig/output_data/m1002/ss10/forecast/figures (beginning with brookings_dsgevarif)

Replicating Figure 6.2 (histograms of selected parameters with separate estimations half1 and half2):
 - Run specfiles_parameters/analyze_separate_estimations_parameters_200129.jl
 - The histogram pdfs are saved in save_orig/output_data/m1002/ss10/estimate/figures

Replicating Figure 6.3 (IRFs with counterfactual policies with separate estimations half1 and half2):
 - Run specfiles_irfs/counter_fix_kappap_dsgevarirfs_hrs_pi_ls_wpi_spec1002_10_2019Q4_1118.jl
 - The IRF pdfs are saved in save_orig/output_data/m1002/ss10/forecast/figures (beginning with counter_fix_kappap)

Replicating Figure 6.4 (histograms of selected parameters with monetary policy regime-switching)
 - Run specfiles_parameters/analyze_regime_switching_parameters.jl
 - The histogram pdfs are saved in save_orig/output_data/m1002/ss24/estimate/figures

Replicating Figure 6.5 (IRFs with monetary policy regime-switching)
 - Run specfiles_irfs/dsgevarirfs_hrs_pic_ls_wpi_spec1002_24_2019Q4_1118.jl
 - The IRF pdfs are saved in save_orig/output_data/m1002/ss10/forecast/figures (names start with brookings_dsgevarirf)

Replicating Figure 6.6 (IRFs with Phillips curve slope regime-switching)
 - Run specfiles_irfs/dsgevarirfs_hrs_pic_ls_wpi_spec1002_22_2019Q4_1118.jl
 - The IRF pdfs are saved in save_orig/output_data/m1002/ss10/forecast/figures (names start with brook (names start with brookings_dsgevarirf)

Replicating Figure 7.1 (forecasts with alternative policies and with Phillips cruve slope regime-switching)
 - Run specfiles_estim/spec1002_forecasts.jl
 - The forecast pdfs are saved in save_orig/output_data/m1002/ss10/forecast/figures (can be found with ls *forecast*)

Replicating Figure C1 (histograms of selected parameters with separate estimations half1 and half2):
 - See Figure 6.2 (this script also replicates this figure)

Replicating Figure C2 (histograms of selected parameters with regime-switching monetary policy. PRE-ZLB):
- See figure 6.4 (that script also replicates this figure)

Replicating Figure C3 (IRFs with separate estimations half1 and half2):
 - See Figure 6.1 (this script also replicates this figure)

Replicating Figure C4 (histograms with regime-switching Phillips curve slope):
- See figure 6.4 (that script also replicates this figure)

set save_orig = false at the top of the spec
