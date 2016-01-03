########################################
# generate points
########################################
rm(list=ls())

library(lhs)
model <- 'glb'
saveFlag <- TRUE # save data?
visualizeFlag <- FALSE

# parameters in (tau_s, eta_s, zeta/s norm)
set.seed(1) # make it reproducible
points_trial <- randomLHS(32, 3)

# augment to 1024 events
points_augmented <- points_trial
for(i in seq(1,5)){
    set.seed(1) # make it reproducible
    points_augmented <- augmentLHS(points_augmented, 
                                   nrow(points_augmented))
}

########################################
# project it to MC-Glb or MC-KLN grids
########################################
if(model=='glb'){
    taus_bd  = c(0.1, 3.0)
    etas_bd  = c(0.0, 0.21)
} else {
    taus_bd  = c(0.1, 2.0)
    etas_bd  = c(0.08, 0.28)
}
bNorm_bd = c(0,3)

# transform to required parameters
points_source <- points_augmented
params.df = data.frame(
    tau_s = taus_bd[1]+(taus_bd[2]-taus_bd[1])*points_source[,1],
    eta_s = etas_bd[1]+(etas_bd[2]-etas_bd[1])*points_source[,2],
    t_dec = 155.0,
    visbulknorm = bNorm_bd[1]+(bNorm_bd[2]-bNorm_bd[1])*points_source[,3])

# save data
if(saveFlag == TRUE){
    params.df.write = format(params.df, digits = 6) # make values for each variable distinct
    print("unique parameters:")
    print(sapply(params.df.write, function(x) length(unique(x))))
    # add comment symbol to the first line by hand
    colnames(params.df.write)[1] = paste("#", colnames(params.df.write)[1]) 
    write.table(params.df.write, file=paste0("params_list_",model,".dat"),
                quote=FALSE, sep="\t", row.names=FALSE)
    print(paste0('Table for ', model, ' has been dumped!'))    
}


########################################
# visualize the grid
########################################
if(visualizeFlag == TRUE){
    cmd <- sprintf('python %s_scatterMatrix.py', model)
    system(cmd)
#     library(gridExtra)
#     library(ggplot2)
#     
#     # do the scatter+histogram plot
#     x_array = params.df$tau_s
#     y_array = params.df$eta_s
#     points_2D_total = nrow(params.df)
#     bin_edges=points_2D_total+1
#     x_breaks = seq(0.1, 2.0, length=40)
#     y_breaks = seq(0.0, 0.16, length=40)
#     x_name = "tau_s"
#     y_name = "eta_s"
#     
#     # histogram for x axis variable
#     hist_top <- ggplot()+geom_histogram(aes(x_array), breaks=x_breaks,
#                                         fill="red",colour="black")+
#         xlab(x_name)
#     # place holder block
#     empty <- ggplot()+geom_point(aes(1,1), colour="white") +
#         theme(                              
#             plot.background = element_blank(), 
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(), 
#             panel.background = element_blank(),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank(),
#             axis.text.x = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks = element_blank()
#         )
#     scatter <- ggplot()+geom_point(aes(x_array, y_array))+
#         xlab(x_name)+ylab(y_name)
#     hist_right <- ggplot()+geom_histogram(aes(y_array), breaks=y_breaks,
#                                           fill="red",colour="black")+
#         coord_flip()+xlab(y_name)
#     grid.arrange(hist_top, empty, scatter, hist_right, 
#                  ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}



