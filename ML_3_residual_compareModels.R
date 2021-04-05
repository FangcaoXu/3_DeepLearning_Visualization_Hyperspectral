#
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
# Author: Guido Cervone (cervone@psu.edu) and Fangcao Xu (xfangcao@psu.edu)
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

library(dplyr)
library(tibble)
library(gtools)
library(RColorBrewer)
library(plotly)
source("Z_global_variables.R")
source("Z_MODTRAN_functions.R")
source("Z_plot_functions.R")
source("Z_general_functions.R")


# Let's plot radiance as a function of wavelength, angles, and reflectivity
day <- 107
time <- c(2,6,10,14,18,22)  # here the time is the local time
reflec <- c(5, 10, 30, 50, 80, 100)


readfourComponentsbyDayTimeReflec <- function(predictfolder,originalfolder, day, time, reflec, down=FALSE){
  if(!down){
    predictfiles.grdw  <- unlist(lapply(reflec, function(z) lapply(time, function(y) lapply(day, function(x) paste(predictfolder, file.patterns["grnd"], sep="") %>% gsub("\\{day\\}", x, .)) %>%
                                                                     gsub("\\{hh\\}", y, .)) %>% gsub("\\{reflect\\}", z, .)))
  }else{
    predictfiles.grdw  <- unlist(lapply(reflec, function(z) lapply(time, function(y) lapply(day, function(x) paste(predictfolder, file.patterns["down"], sep="") %>% gsub("\\{day\\}", x, .)) %>%
                                                                     gsub("\\{hh\\}", y, .)) %>% gsub("\\{reflect\\}", z, .)))
  }
  predictfiles.up1 <- unlist(lapply(reflec, function(z) lapply(time, function(y) lapply(day, function(x) paste(predictfolder, file.patterns["up1"], sep="") %>% gsub("\\{day\\}", x, .)) %>%
                                                          gsub("\\{hh\\}", y, .)) %>% gsub("\\{reflect\\}", z, .)))
  predictfiles.up2 <- unlist(lapply(reflec, function(z) lapply(time, function(y) lapply(day, function(x) paste(predictfolder, file.patterns["up2"], sep="") %>% gsub("\\{day\\}", x, .)) %>%
                                                          gsub("\\{hh\\}", y, .)) %>% gsub("\\{reflect\\}", z, .)))
  originalfiles.grdw <- find.originalfiles(predictfiles.grwd, originalfolder) # ground-reflected or downwelling 
  originalfiles.up1 <- find.originalfiles(predictfiles.up1, originalfolder)
  originalfiles.up2<- find.originalfiles(predictfiles.up2, originalfolder)
  totalInputfiles <- find.originaltotalInput(predictfiles.up1, originalfolder)
  
  predorigfiles.grdw.read <- read.predorigmatrix(predictfiles.grwd, originalfiles.grdw)
  predorigfiles.up1.read <- read.predorigmatrix(predictfiles.up1, originalfiles.up1)
  predorigfiles.up2.read <- read.predorigmatrix(predictfiles.up2, originalfiles.up2)
  totalInputfiles.read <- read.originalmatrix(totalInputfiles)
  return(list(predorigfiles.grdw.read, predorigfiles.up1.read, predorigfiles.up2.read, totalInputfiles.read))
}

min.max.all <- function(predictedorigin){
  grwd.diff <-  mapply(function(x,y) diff.matrices(x, y), predictedorigin[[1]][[1]], predictedorigin[[1]][[2]], SIMPLIFY = FALSE)
  up1.diff <- mapply(function(x,y) diff.matrices(x, y), predictedorigin[[2]][[1]], predictedorigin[[2]][[2]], SIMPLIFY = FALSE)
  up2.diff <- mapply(function(x,y) diff.matrices(x, y), predictedorigin[[3]][[1]], predictedorigin[[3]][[2]], SIMPLIFY = FALSE)
  return(list(grwd.diff, up1.diff, up2.diff))
}

june30th.gr <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedCSV/MidLatitude", "~/disk10TB/DARPA/MatrixCSV/MidLatitude", day, time, reflec, down=FALSE) # returned grnd, up1, up2, total
june30th.dw <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedCSV/MidLatitude","~/disk10TB/DARPA/MatrixCSV/MidLatitude", day, time, reflec, down=TRUE) # returned down, up1, up2, total
aug5th <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedThree/Aug5th", "~/disk10TB/DARPA/MatrixCSV/MidLatitude", day, time, reflec, down=FALSE) # returned grnd, up1, up2, total
aug9th <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedThree/Aug9th", "~/disk10TB/DARPA/MatrixCSV/MidLatitude",day, time, reflec, down=TRUE) # returned down, up1, up2, total
aug10th.gr <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedThree/Aug10th/groundreflected", "~/disk10TB/DARPA/MatrixCSV/MidLatitude", day, time, reflec, down=FALSE) # returned grnd, up1, up2, total
aug10th.dw <- readfourComponentsbyDayTimeReflec("~/disk10TB/DARPA/predictedThree/Aug10th/downwelling", "~/disk10TB/DARPA/MatrixCSV/MidLatitude", day, time, reflec, down=TRUE) # returned down, up1, up2, total
#june30th[[1]]

########################################################### Statistical Tables to Check Different Models Residuals #################################################################################
# Aug 5th, ground reflected, up1, up2
model.changes5th <- as.data.frame(t(mapply(function(x,y) c(max(abs(x)), max(abs(y)), min(abs(x)), min(abs(y))), unlist(min.max.all(june30th.gr), recursive = FALSE), 
                                     unlist(min.max.all(aug5th), recursive = FALSE))), stringsAsFactors = FALSE)
model.changes5th <-rownames_to_column(model.changes5th)
colnames(model.changes5th) <- c("filename", "June30thModel_max", "Aug5thModel_max","June30thModel_min","Aug5thModel__min")
write.csv(model.changes5th, "/amethyst/s0/fbx5002/geolab_storage_V3/data/DARPA/Plots/Others/Changes5th.csv")
# Aug 9th, downwelling, up1, up2
model.changes9th <- as.data.frame(t(mapply(function(x,y) c(max(abs(x)), max(abs(y)), min(abs(x)), min(abs(y))), unlist(min.max.all(june30th.dw), recursive = FALSE), 
                                          unlist(min.max.all(aug9th), recursive = FALSE))), stringsAsFactors = FALSE)
model.changes9th <-rownames_to_column(model.changes9th)
colnames(model.changes9th) <- c("filename", "June30thModel_max", "Aug9thModel_max","June30thModel_min","Aug9thModel__min")
write.csv(model.changes9th, "/amethyst/s0/fbx5002/geolab_storage_V3/data/DARPA/Plots/Others/Changes9th.csv")
# Aug 10th, ground reflected, up1, up2
model.changes10th.gr <- as.data.frame(t(mapply(function(x,y) c(max(abs(x)), max(abs(y)), min(abs(x)), min(abs(y))), unlist(min.max.all(june30th.gr), recursive = FALSE), 
                                           unlist(min.max.all(aug10th.gr), recursive = FALSE))), stringsAsFactors = FALSE)
model.changes10th.gr <-rownames_to_column(model.changes10th.gr)
colnames(model.changes10th.gr) <- c("filename", "June30thModel_max", "Aug10thModel_max","June30thModel_min","Aug10thModel__min")
write.csv(model.changes10th.gr, "/amethyst/s0/fbx5002/geolab_storage_V3/data/DARPA/Plots/Others/Changes10th_gr.csv")
# Aug 10th, downwelling, up1, up2
model.changes10th.dw <- as.data.frame(t(mapply(function(x,y) c(max(abs(x)), max(abs(y)), min(abs(x)), min(abs(y))), unlist(min.max.all(june30th.dw), recursive = FALSE), 
                                           unlist(min.max.all(aug10th.dw), recursive = FALSE))), stringsAsFactors = FALSE)
model.changes10th.dw <-rownames_to_column(model.changes10th.dw)
colnames(model.changes10th.dw) <- c("filename", "June30thModel_max", "Aug10thModel_max","June30thModel_min","Aug10thModel__min")
write.csv(model.changes10th.dw, "/amethyst/s0/fbx5002/geolab_storage_V3/data/DARPA/Plots/Others/Changes10th_dw.csv")



########################################################################################## Comparison Plots #########################################################################################
diff.minmax <- minmax.crosstables(c(unlist(min.max.all(june30th.gr), recursive = FALSE),unlist(min.max.all(june30th.dw), recursive = FALSE),
                                     unlist(min.max.all(aug5th), recursive = FALSE), unlist(min.max.all(aug9th), recursive = FALSE),
                                     unlist(min.max.all(aug10th.gr), recursive = FALSE),unlist(min.max.all(aug10th.dw), recursive = FALSE)))
########## june 30th
# ground reflected component
comparisonMatrixPlot(june30th.gr[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/GroundReflected/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# downwelling component
comparisonMatrixPlot(june30th.dw[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/Downwelling/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal emission component
comparisonMatrixPlot(june30th.gr[[2]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(june30th.gr[[3]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = diff.minmax)
########## aug5th
# ground reflected component
comparisonMatrixPlot(aug5th[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug5th/GroundReflected/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal emission component
comparisonMatrixPlot(aug5th[[2]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug5th/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(aug5th[[3]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug5th/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = diff.minmax)
########## aug9th
# downwelling component
comparisonMatrixPlot(aug9th[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug9th/Downwelling/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal emission component
comparisonMatrixPlot(aug9th[[2]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug9th/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(aug9th[[3]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug9th/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = diff.minmax)
########## aug10th
# downwelling component
comparisonMatrixPlot(aug10th.dw[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/downwelling/Downwelling/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal emission component
comparisonMatrixPlot(aug10th.dw[[2]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/downwelling/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(aug10th.dw[[3]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/downwelling/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = diff.minmax)
# ground reflected component
comparisonMatrixPlot(aug10th.gr[[1]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/groundreflected/GroundReflected/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal emission component
comparisonMatrixPlot(aug10th.gr[[2]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/groundreflected/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(aug10th.gr[[3]], "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/groundreflected/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = diff.minmax)

# total radiance component
mapply(function(x,y) MODTRAN.matplot(as.matrix(x), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/TotalRadiance/", y, ".jpg",sep = "")), june30th.dw[[4]], names(june30th.dw[[4]]))
mapply(function(x,y) MODTRAN.matplot(as.matrix(x), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug5th/TotalRadiance/", y, ".jpg",sep = "")), aug5th[[4]], names(aug5th[[4]]))
mapply(function(x,y) MODTRAN.matplot(as.matrix(x), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug9th/TotalRadiance/", y, ".jpg",sep = "")), aug9th[[4]], names(aug9th[[4]]))
mapply(function(x,y) MODTRAN.matplot(as.matrix(x), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Aug10th/TotalRadiance/", y, ".jpg",sep = "")), aug10th.dw[[4]], names(aug10th.dw[[4]]))


### read and plot the predicted transmission for day 107 with autumn winter model 
trans.ss.byaw <- list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude", pattern="*.csv", full.names = TRUE)
originalfile.trans.ss.byaw <- find.originalfiles(trans.ss.byaw, "~/disk10TB/DARPA/MatrixCSV/MidLatitude")
trans.ss.byaw.read <- rowcol.renames(read.csv(trans.ss.byaw, header = FALSE, stringsAsFactors = FALSE), wv, an)
originalfile.trans.ss.byaw.read <- rowcol.renames(read.csv(originalfile.trans.ss.byaw, skip = 3, header = FALSE, stringsAsFactors = FALSE)[,-1][,angles.index], wv, an)
zlim.predorig <- minmax.crosstables(list(trans.ss.byaw.read,originalfile.trans.ss.byaw.read))
MODTRAN.matplot(as.matrix(trans.ss.byaw.read),zlim=zlim.predorig, filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/useAWPredictSS/Radiance_107_14_10_trans_predict.jpg")
MODTRAN.matplot(as.matrix(originalfile.trans.ss.byaw.read),zlim=zlim.predorig, filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/useAWPredictSS/Radiance_107_14_10_trans_original.jpg")
MODTRAN.matplot(as.matrix(diff.matrices(trans.ss.byaw.read, originalfile.trans.ss.byaw.read)),filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/useAWPredictSS/Radiance_107_14_10_trans_residual.jpg", residual="T")

