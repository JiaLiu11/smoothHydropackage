# generate the chi-square for glb and kln parameter search result

rm(list=ls())
# load libraries
suppressMessages(library(dplyr))


# experimental value
v2_ch_exp_value=0.0782 #from ATLAS_data/atlas_vn_data/meanNpart_meanV2.dat, pT>0.5GeV
v2_ch_exp_error=0.0019
v3_ch_exp_value=0.03165
v3_ch_exp_error=7.7e-04
pion_exp_value = 0.541850 # estimated from spectra
pion_exp_error = 0.017965
kaon_exp_value = 0.824950
kaon_exp_error = 0.027854
proton_exp_value=1.311448
proton_exp_error=0.043314
proton_pt2_exp_value = 2.085105
proton_pt2_exp_error = 0.070462
pik_ratio_exp_value = 6.69118
pik_ratio_exp_error = 0.67074
pip_ratio_exp_value = 21.6667
pip_ratio_exp_error = 2.29250
pil_ratio_exp_value = 26.7647
pil_ratio_exp_error = 3.63870
dn_ch_exp_value = 966 # do not use
dn_ch_exp_error = 37


chiSqGenerator<-function(infile, outfile, 
                         v2_exp, v2_exp_error,
                         v3_exp, v3_exp_error,
                         pion_exp,pion_exp_error,
                         kaon_exp,kaon_exp_error,
                         proton_exp, proton_exp_error,
                         proton_pt2_exp_value,proton_pt2_exp_error,
                         pik_ratio_exp_value ,pik_ratio_exp_error,
                         pip_ratio_exp_value ,pip_ratio_exp_error,
                         pil_ratio_exp_value ,pil_ratio_exp_error){
    # genearte the chi square for one type of run.
    # at each data point, we have 5 reference data: v2, v3, pion_meanPT, kaon_meanpT, proton_meanPT;
    # output format: taus eta/s Tdec chi_sq
    raw_data = read.table(infile, header=T)
    colnames(raw_data) = c('taus','eta_s', 'Tdec', 'visbulkNorm',
                           'v2_ch', 'v2_ch_err',
                           'v3_ch', 'v3_ch_err', 
                           'meanPT_pion', 'meanPT_pion_err',
                           'meanPT_kaon', 'meanPT_kaon_err', 
                           'meanPT_proton', 'meanPT_proton_err',
                           "proton_pt2_mean","proton_pt2_err",  
                           "pik_ratio_mean",  "pik_ratio_err",  
                           "pip_ratio_mean",  "pip_ratio_err",
                           "pil_ratio_mean", "pil_ratio_err")
    result.df = select(raw_data, taus, eta_s, Tdec, visbulkNorm)
    result.df$chiSq = (raw_data$v2_ch-v2_exp)^2/(v2_exp_error^2+raw_data$v2_ch_err^2)+
        (raw_data$v3_ch-v3_exp)^2/(v3_exp_error^2+raw_data$v3_ch_err^2)+
        (raw_data$meanPT_pion-pion_exp)^2/(pion_exp_error^2+raw_data$meanPT_pion_err^2)+
        (raw_data$meanPT_kaon-kaon_exp)^2/(kaon_exp_error^2+raw_data$meanPT_kaon_err^2)+
        (raw_data$meanPT_proton-proton_exp)^2/(proton_exp_error^2+raw_data$meanPT_proton_err^2)+
        (raw_data$proton_pt2_mean-proton_pt2_exp_value)^2/(proton_pt2_exp_error^2+raw_data$meanPT_proton_err^2)+
        (raw_data$pik_ratio_mean-pik_ratio_exp_value)^2/(pik_ratio_exp_error^2+raw_data$pik_ratio_err^2)+
        (raw_data$pip_ratio_mean-pip_ratio_exp_value)^2/(pip_ratio_exp_error^2+raw_data$pip_ratio_err^2)+
        (raw_data$pil_ratio_mean-pil_ratio_exp_value)^2/(pil_ratio_exp_error^2+raw_data$pil_ratio_err^2)
    raw_data$chiSq = result.df$chiSq
    raw_data = raw_data %>% arrange(chiSq)
    
    chiSqContribution = raw_data %>% 
        transmute(taus = taus, eta_s=eta_s, tdec = Tdec, visbulkNorm = visbulkNorm,
                  v2_chisq = (v2_ch-v2_exp)^2/(v2_exp_error^2+v2_ch_err^2),
                  v3_chisq = (v3_ch-v3_exp)^2/(v3_exp_error^2+v3_ch_err^2),
                  pion_chisq = (meanPT_pion-pion_exp)^2/(pion_exp_error^2+meanPT_pion_err^2),
                  kaon_chisq = (meanPT_kaon-kaon_exp)^2/(kaon_exp_error^2+meanPT_kaon_err^2),
                  proton_chisq = (meanPT_proton-proton_exp)^2/(proton_exp_error^2+meanPT_proton_err^2),
                  proton_pt2_chisq = (proton_pt2_mean-proton_pt2_exp_value)^2/(proton_pt2_exp_error^2+meanPT_proton_err^2),
                  pik_ratio_chisq = (pik_ratio_mean-pik_ratio_exp_value)^2/(pik_ratio_exp_error^2+pik_ratio_err^2),
                  pip_ratio_chisq = (pip_ratio_mean-pip_ratio_exp_value)^2/(pip_ratio_exp_error^2+pip_ratio_err^2),
                  pil_ratio_chisq = (pil_ratio_mean-pil_ratio_exp_value)^2/(pil_ratio_exp_error^2+pil_ratio_err^2),
                  chiSq = v2_chisq+v3_chisq+pion_chisq+kaon_chisq+proton_chisq+proton_pt2_chisq+pik_ratio_chisq+pip_ratio_chisq+pil_ratio_chisq) %>%
        arrange(chiSq)
    
    # add a comment to the first colname
    comment_symbol = "#" # for python
    colnames(result.df)[1] = paste(comment_symbol, colnames(result.df)[1],
                                   sep="")
    write.table(format(result.df, scientific = T), 
                outfile, quote=F, sep="\t",
                row.names=FALSE, col.names=TRUE)
    return(list("results"=result.df, 
                "data"=raw_data,
                "contribution"=chiSqContribution))
}

# process glb1
in_file = 'param_search_log_MCGlb_1_4params.dat'
out_file = 'hasfs_glb_tet_chiSq_4params.dat'
results_list = chiSqGenerator(in_file, out_file,
                              v2_ch_exp_value, v2_ch_exp_error,
                              v3_ch_exp_value, v3_ch_exp_error,
                              pion_exp_value, pion_exp_error,
                              kaon_exp_value, kaon_exp_error,
                              proton_exp_value, proton_exp_error,
                              proton_pt2_exp_value,proton_pt2_exp_error,
                              pik_ratio_exp_value ,pik_ratio_exp_error,
                              pip_ratio_exp_value ,pip_ratio_exp_error,
                              pil_ratio_exp_value ,pil_ratio_exp_error)
result_now = results_list$results
data_now = results_list$data
contribtuion_now = results_list$contribution

