# generate the chi-square for glb and kln parameter search result

rm(list=ls())
# load libraries
suppressMessages(library(dplyr))


# experimental value
v2_ch_exp_value=0.0782 #from ATLAS_data/atlas_vn_data/meanNpart_meanV2.dat, pT>0.5GeV
v2_ch_exp_error=0.0019
v3_ch_exp_value=0.03165
v3_ch_exp_error=7.7e-04
pion_exp_value = 0.517 # http://hepdata.cedar.ac.uk/view/ins1222333/d32;jsessionid=9g2su0gjgfvl
pion_exp_error = 0.017
kaon_exp_value = 0.871
kaon_exp_error = 0.027
proton_exp_value=1.311
proton_exp_error=0.034
dn_ch_exp_value = 966 # do not use
dn_ch_exp_error = 37


chiSqGenerator<-function(infile, outfile, 
                         v2_exp, v2_exp_error,
                         v3_exp, v3_exp_error,
                         pion_exp,pion_exp_error,
                         kaon_exp,kaon_exp_error,
                         proton_exp, proton_exp_error){
# genearte the chi square for one type of run.
# at each data point, we have 5 reference data: v2, v3, pion_meanPT, kaon_meanpT, proton_meanPT;
# output format: taus eta/s Tdec chi_sq
    raw_data = read.table(infile, header=T)
    colnames(raw_data) = c('taus','eta_s', 'Tdec', 'visbulkNorm',
                           'v2_ch', 'v2_ch_err',
                           'v3_ch', 'v3_ch_err', 
                           'meanPT_pion', 'meanPT_pion_err',
                           'meanPT_kaon', 'meanPT_kaon_err', 
                           'meanPT_proton', 'meanPT_proton_err')
    result.df = select(raw_data, taus, eta_s, Tdec, visbulkNorm)
    result.df$chiSq = (raw_data$v2_ch-v2_exp)^2/(v2_exp_error^2+raw_data$v2_ch_err^2)+
                      (raw_data$v3_ch-v3_exp)^2/(v3_exp_error^2+raw_data$v3_ch_err^2)+
                      (raw_data$meanPT_pion-pion_exp)^2/(pion_exp_error^2+raw_data$meanPT_pion_err^2)+
                      (raw_data$meanPT_kaon-kaon_exp)^2/(kaon_exp_error^2+raw_data$meanPT_kaon_err^2)+
                     (raw_data$meanPT_proton-proton_exp)^2/(proton_exp_error^2+raw_data$meanPT_proton_err^2)
    raw_data$chiSq = result.df$chiSq
    raw_data = raw_data %>% arrange(chiSq)
    raw_data$Tdec=NULL
    
    chiSqContribution = raw_data %>% 
        transmute(taus = taus, eta_s=eta_s, visbulkNorm = visbulkNorm,
                  v2_chisq = (v2_ch-v2_exp)^2/(v2_exp_error^2+v2_ch_err^2),
                  v3_chisq = (v3_ch-v3_exp)^2/(v3_exp_error^2+v3_ch_err^2),
                  pion_chisq = (meanPT_pion-pion_exp)^2/(pion_exp_error^2+meanPT_pion_err^2),
                  kaon_chisq = (meanPT_kaon-kaon_exp)^2/(kaon_exp_error^2+meanPT_kaon_err^2),
                  proton_chisq = (meanPT_proton-proton_exp)^2/(proton_exp_error^2+meanPT_proton_err^2),
                  chiSq = v2_chisq+v3_chisq+pion_chisq+kaon_chisq+proton_chisq) %>%
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



# # process glb0
# in_file = 'MCGlb_0/param_search_log_MCGlb_0.dat'
# out_file = 'nofs_glb_tet_chiSq.dat'
# results_list = chiSqGenerator(in_file, out_file,
#                         v2_ch_exp_value, v2_ch_exp_error,
#                         v3_ch_exp_value, v3_ch_exp_error,
#                         pion_exp_value, pion_exp_error,
#                         kaon_exp_value, kaon_exp_error,
#                         proton_exp_value, proton_exp_error)
# result_now = results_list$results
# data_now = results_list$data
# contribtuion_now = results_list$contribution

# process glb1
in_file = 'MCGlb_1/param_search_log_MCGlb_1.dat'
out_file = 'hasfs_glb_tet_chiSq.dat'
results_list = chiSqGenerator(in_file, out_file,
                              v2_ch_exp_value, v2_ch_exp_error,
                              v3_ch_exp_value, v3_ch_exp_error,
                              pion_exp_value, pion_exp_error,
                              kaon_exp_value, kaon_exp_error,
                              proton_exp_value, proton_exp_error)
result_now = results_list$results
data_now = results_list$data
contribtuion_now = results_list$contribution

# # process MCKLN0
# in_file = 'MCKLN_0/param_search_log_MCKLN_0.dat'
# out_file = 'nofs_kln_tet_chiSq.dat'
# result_now = chiSqGenerator(in_file, out_file,
#                             v2_ch_exp_value, v2_ch_exp_error,
#                             v3_ch_exp_value, v3_ch_exp_error,
#                             pion_exp_value, pion_exp_error,
#                             kaon_exp_value, kaon_exp_error,
#                             proton_exp_value, proton_exp_error)


# # process MCKLN1
# in_file = 'MCKLN_1/param_search_log_MCKLN_1.dat'
# out_file = 'hasfs_kln_tet_chiSq.dat'
# results_list = chiSqGenerator(in_file, out_file,
#                               v2_ch_exp_value, v2_ch_exp_error,
#                               v3_ch_exp_value, v3_ch_exp_error,
#                               pion_exp_value, pion_exp_error,
#                               kaon_exp_value, kaon_exp_error,
#                               proton_exp_value, proton_exp_error)
# result_now = results_list$results
# data_now = results_list$data
# contribution_now = results_list$contribution