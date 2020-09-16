using SMC
is_split = false

fp = "../../save_orig/output_data/m1002"

filenames = ["$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r1_preZLB=true_reg2start=900331_vint=191118.jld2",
             "$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r1_preZLB=true_reg2start=840331_vint=191118.jld2",
             "$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r1_preZLB=false_reg2start=900331_vint=191118.jld2",
             "$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r2_preZLB=true_reg2start=900331_vint=191118.jld2",
             "$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r2_preZLB=true_reg2start=840331_vint=191118.jld2",
             "$(fp)/ss10/estimate/raw/smc_cloud_npart=15000_period=r2_preZLB=false_reg2start=900331_vint=191118.jld2",
             "$(fp)/ss21/estimate/raw/smc_cloud_period=full_incZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss21/estimate/raw/smc_cloud_period=full_preZLB_reg2start=840331_reg=2_vint=191118.jld2",
             "$(fp)/ss21/estimate/raw/smc_cloud_period=full_preZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss22/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_incZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss22/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_preZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss22/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_preZLB_reg2start=840331_reg=2_vint=191118.jld2",
             "$(fp)/ss24/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_incZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss24/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_preZLB_reg2start=900331_reg=2_vint=191118.jld2",
             "$(fp)/ss24/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_preZLB_reg2start=840331_reg=2_vint=191118.jld2"]


for filename in filenames
    if is_split
        split_cloud(filename, 4)
    else
        join_cloud(filename, 4)
    end
end
