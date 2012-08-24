for E in pdf eps png; do

# Fig. 1
cp plots/sqr_acls_allexp_bydecay_logx.$E plots/PAS/mu95_subchannels_expected.$E -v
cp plots/sqr_acls_allexp_bydecay_zoom.$E plots/PAS/mu95_subchannels_expected_zoom.$E -v

# Fig. 2
cp plots/sqr_smacls_comb_logx.$E plots/PAS/CLs.$E -v
cp plots/sqr_smacls_comb_zoom.$E plots/PAS/CLs_zoom.$E -v


# Fig. 3
cp plots/sqr_acls_comb_logx.$E plots/PAS/mu95.$E -v
cp plots/sqr_acls_comb_zoom.$E plots/PAS/mu95_zoom.$E -v

# Fig. 4
cp plots/sqr_acls_all2_energy_logx.$E plots/PAS/mu95_7vs8.$E -v
cp plots/sqr_acls_all2_energy_zoom.$E plots/PAS/mu95_7vs8_zoom.$E -v

#Fig.5:
cp plots/sqr_pvala_all_energy.$E plots/PAS/pvalue_7vs8_zoom.$E -v
cp plots/sqr_pvala_all_bydecay.$E plots/PAS/pvalue_zoom.$E -v

#Fig.6:
#muhat_7vs8_zoom.pdf
#muhat_zoom.pdf

#Fig 7
#compatibility_11subcomb_7vs8.pdf
#compatibility_11subcomb.pdf

cp plots/sqr_mlzs_ccc_mH125.0.$E        plots/PAS/compatibility_11subcomb.$E
cp plots/sqr_mlzs_ccc_mH125.0_decay.$E  plots/PAS/compatibility_decay.$E
cp plots/sqr_mlzs_ccc_mH125.0_prod.$E   plots/PAS/compatibility_production.$E
cp plots/sqr_mlzs_ccc_mH125.0_energy.$E plots/PAS/compatibility_11subcomb_7vs8.$E
 

#Fig.8:
#compatibility_decay.pdf
#compatibility_production.pdf

#Fig.9:
cp plots/sqr_mass_scan_1d_all.$E plots/PAS/mH_2LLR_scan_gamgam_vs_zz4l.$E -v
cp plots/subchannels/sqr_mass_scan_1d_hires_comp.$E plots/PAS/mH_2LLR_scan_NoSyst_vs_WithSyst.$E -v

cp plots/sqr_mass_scan_2d_all_white.$E plots/PAS/mH_muhat_contours.$E -v
cp plots/subchannels/sqr_mass_scan_2d_hires.$E plots/PAS/mH_muhat_2LLR_scan.$E -v

#Fig.10:
#mH_muhat_contours.pdf
#mH_muhat_2LLR_scan.pdf

#Fig.11:
#cVcF.pdf

#Fig. 12
# MELA angles

#Fig.13:
cp plots/7/sqr_pvala_all_bydecay.$E plots/PAS/pvalue_7.$E -v
cp plots/8/sqr_pvala_all_bydecay.$E plots/PAS/pvalue_8.$E -v

#Fig.14:
cp plots/sqr_pvala_all_energy_hires.$E plots/PAS/pvalue_HighMassRes.$E -v
cp plots/sqr_pvala_all_energy_lowres.$E plots/PAS/pvalue_LowMassRes.$E -v


cp plots/subchannels/sqr_mass_scan_2d_ggH_hgg.$E   plots/PAS/appx_mH_muhat_2LLR_scan_gamgam_incl.$E
cp plots/subchannels/sqr_mass_scan_2d_hires.$E   plots/PAS/mH_muhat_2LLR_scan.$E
cp plots/subchannels/sqr_mass_scan_2d_hzz.$E   plots/PAS/appx_mH_muhat_2LLR_scan_zz4l.$E
cp plots/subchannels/sqr_mass_scan_2d_qqH_hgg.$E   plots/PAS/appx_mH_muhat_2LLR_scan_gamgam_vbf.$E
cp plots/subchannels/sqr_mass_bayes_2d_ggH_hgg.$E   plots/PAS/appx_mH_muhat_BayesPosterior_gamgam_incl.$E
cp plots/subchannels/sqr_mass_bayes_2d_hires.$E   plots/PAS/appx_mH_muhat_BayesPosterior.$E
cp plots/subchannels/sqr_mass_bayes_2d_hzz.$E   plots/PAS/appx_mH_muhat_BayesPosterior_zz4l.$E
cp plots/subchannels/sqr_mass_bayes_2d_qqH_hgg.$E   plots/PAS/appx_mH_muhat_BayesPosterior_gamgam_vbf.$E

cp plots/sqr_cvcf_scan_2d_comb.$E   plots/PAS/cVcF_2quadrants.$E
cp plots/sqr_cvcf_cut_scan_2d_comb.$E   plots/PAS/cVcF_1quadrant.$E
cp plots/sqr_cvcf_cut_bayes_2d_comb.$E   plots/PAS/appx_cVcF_1quadrant_BayesPosterior.$E
cp plots/sqr_cvcf_bayes_2d_comb.$E   plots/PAS/appx_cVcF_2quadrants_BayesPosterior.$E
done




