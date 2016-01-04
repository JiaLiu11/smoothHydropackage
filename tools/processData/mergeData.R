library(dplyr)
rm(list=ls())

# run options
mode = 'MCGlb_1'

# read in
v2_data <- read.table(paste(mode, 
                            sprintf('param_search_log_v2_%s.dat',tolower(mode)),sep='/'),
                      header = F, na.strings = 'nan')
v3_data <- read.table(paste(mode, 
                            sprintf('param_search_log_v3_%s.dat', tolower(mode)),sep='/'),
                      header = F,na.strings = 'nan')

# output
outputfile <- file.path(mode, 
                        paste0('param_search_log_',mode,'_test.dat'))

# join by 
joined <- left_join(v2_data, v3_data, by='V2')

# get new df
df_result <- data.frame(taus = joined$V1.x,
                        etas = joined$V2,
                        tdec = joined$V3.x,
                        VisBulkNorm = joined$V4.x,
                        v2_ch_mean = joined$V5.x,
                        v2_ch_err  = joined$V6.x,
                        v3_ch_mean = joined$V7.y,
                        v3_ch_err  = joined$V8.y,
                        pion_meanPT_mean = (joined$V9.x+joined$V9.y)/2,
                        pion_meanPT_err  = sqrt(joined$V10.x^2+joined$V10.y^2)/sqrt(2),
                        kaon_meanPT_mean = (joined$V11.x+joined$V11.y)/2,
                        kaon_meanPT_err  = sqrt(joined$V12.x^2+joined$V12.y^2)/sqrt(2),
                        proton_meanPT_mean = (joined$V13.x+joined$V13.y)/2,
                        proton_meanPT_err  = sqrt(joined$V14.x^2+joined$V14.y^2)/sqrt(2)
                        )

# output data
df_result_complete <- df_result[complete.cases(df_result),]
write.table(format(df_result_complete, scientific = T), 
            outputfile, quote=F, sep="\t",
            row.names=FALSE, col.names=TRUE)

