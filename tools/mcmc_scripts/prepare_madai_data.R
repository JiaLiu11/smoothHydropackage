# generate model_output folders for MADAI tools 
# to train emulator
rm(list=ls())
# run options
model <- 'MCGlb_1'
madai_folder <- '/Users/Jia/code/madai/paramSearch'
parent_folder_name <- 'model_output'

# read in data
data_raw <- read.table(sprintf('%s/param_search_log_%s.dat', model, model),
                       header = T)
folders_total <- nrow(data_raw)
# redo colname
colnames(data_raw)<-c('tau_s', 'eta_over_s', 't_dec', 'bulk_norm',
                      'v2_ch', 'v2_ch_error', 'v3_ch', 'v3_ch_error',
                      'meanpt_pion', 'meanpt_pion_error',
                      'meanpt_kaon', 'meanpt_kaon_error',
                      'meanpt_proton', 'meanpt_proton_error')

# prepare folder structure
target_folder <- sprintf('%s/%s', madai_folder, model)
folder_names <- sprintf('run%04d', seq(1, folders_total))
if( file.exists(file.path(target_folder, parent_folder_name))==F){
    dir.create(file.path(target_folder, parent_folder_name))
}

# loop over folders
for(i in seq(1, folders_total)){
    ifolder_name <- folder_names[i]
    # create folder if not exists
    folder_now <- file.path(target_folder, 
                            parent_folder_name,
                            ifolder_name)
    if(file.exists(folder_now)==F){
        dir.create(folder_now)
    }
    # compose parameter data
    param_data_now <- data_raw[i, c(1,2,4)]
    write.table(t(param_data_now), 
                file=file.path(folder_now, 'parameters.dat'),
                col.names=F, quote=F)
    # compose results data
    sim_data_now <-data_raw[i, 5:14]
    sim_data_now <- matrix(sim_data_now, 
                           nrow=5, ncol=2,
                           byrow=T)
    row.names(sim_data_now)<-c('v2_ch', 'v3_ch', 'meanpt_pion', 
                               'meanpt_kaon', 'meanpt_proton')
    write.table(sim_data_now, 
                file=file.path(folder_now, 'results.dat'),
                col.names=F, quote=F)    
}
print(sprintf('All %d folders for %s model completes!', 
              folders_total, model))



