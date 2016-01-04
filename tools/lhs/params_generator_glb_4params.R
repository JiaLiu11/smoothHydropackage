################################################################
# Generate 3-parameter triplet (tau_s, eta_s, t_sw, zeta/s norm) 
# for MC-KLN or MC-Glb parameter search program

# Author: Jia Liu

################################################################
rm(list=ls())

library(lhs)
model <- 'glb'
saveFlag <- TRUE # save data?
visualizeFlag <- FALSE

# generate random numbers
set.seed(1) # make it reproducible
points_trial <- randomLHS(1024, 4)

########################################
# project it to MC-Glb or MC-KLN grids
########################################
if(model=='glb'){
    taus_bd  = c(0.1, 3.0)
    etas_bd  = c(0.0, 0.21)
    tsw_bd  = c(130, 170)
} else {
    taus_bd  = c(0.1, 2.0)
    etas_bd  = c(0.08, 0.28)
}
bNorm_bd = c(0,4)

# transform to required parameters
points_source <- points_trial
params.df = data.frame(
    tau_s = taus_bd[1]+(taus_bd[2]-taus_bd[1])*points_source[,1],
    eta_s = etas_bd[1]+(etas_bd[2]-etas_bd[1])*points_source[,2],
    t_sw = tsw_bd[1]+(tsw_bd[2]-tsw_bd[1])*points_source[,3],
    visbulknorm = bNorm_bd[1]+(bNorm_bd[2]-bNorm_bd[1])*points_source[,4])

# save data
if(saveFlag == TRUE){
    params.df.write = format(params.df, digits = 6) # make values for each variable distinct
    print("unique parameters:")
    print(sapply(params.df.write, function(x) length(unique(x))))
    # add comment symbol to the first line by hand
    colnames(params.df.write)[1] = paste("#", colnames(params.df.write)[1]) 
    write.table(params.df.write, file=paste0("params_list_",model,"_4params.dat"),
                quote=FALSE, sep="\t", row.names=FALSE)
    print(paste0('Table for ', model, ' has been dumped!'))    
}


########################################
# visualize the grid
########################################
if(visualizeFlag == TRUE){
    cmd <- sprintf('python %s_scatterMatrix_4params.py', model)
    system(cmd)
}



