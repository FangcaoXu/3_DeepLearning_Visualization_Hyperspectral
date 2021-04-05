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
# # error sources: 1) ML training and predicting model (predicted up1, up2, down, transmissivity); 
# # 2) RTE simplification; 3) provided true reflectivity which are discrete values

source("Z_global_variables.R")
source("Z_MODTRAN_functions.R")
source("Z_plot_functions.R")
source("Z_general_functions.R")


########################################################## Codes start here ########################################################## 
# read lambertian surface targets for the midlatitude
predictfiles.trans <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/Stage1_7.5_12um/MidLatitude/TransmissionCSV",full.names=TRUE))
predictfiles.down <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/Stage1_7.5_12um/MidLatitude/DownwellingCSV",full.names=TRUE))
predictfiles.up1 <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/Stage1_7.5_12um/MidLatitude/PathThermalEmissionCSV",full.names=TRUE))
predictfiles.up2 <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/Stage1_7.5_12um/MidLatitude/PathThermalScatteringCSV",full.names=TRUE))
# original surface emission and total radiance files
originalfiles.surf <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/MidLatitude/SurfaceEmissionCSV",full.names=TRUE))
totalInputfiles <- mixedsort(list.files('~/disk10TB/DARPA/MatrixCSV/MidLatitude/Stage1_7.5_12um/TotalRadianceCSV',full.names=TRUE))
# predicted files stored in disk10TB are not complete and same. #15530 total files 
vars.trans <-return.index(predictfiles.trans, 2:4, names = c("Day","Time","Reflectivity"))
vars.down<-return.index(predictfiles.down, 2:4, names = c("Day","Time","Reflectivity"))
vars.up1 <-return.index(predictfiles.up1, 2:4, names = c("Day","Time","Reflectivity"))
vars.up2 <-return.index(predictfiles.up2, 2:4, names = c("Day","Time","Reflectivity"))
vars.surf<-return.index(originalfiles.surf, 2:4, names = c("Day","Time","Reflectivity"))
vars.total <- return.index(totalInputfiles, 2:4, names = c("Day","Time","Reflectivity"))
# Let's plot radiance as a function of wavelength, angles, and reflectivity
day <- 107
time <- 14  # here the time is the local time
# prepare the sliced input data based on the assigned variables
valid.index <- which(vars.trans$Day == day & vars.trans$Time == time)
input.trans <- read.predictedmatrix(predictfiles.trans[valid.index])
# down
valid.index <- which(vars.down$Day == day & vars.down$Time == time)
input.down <- read.predictedmatrix(predictfiles.down[valid.index])
# up1
valid.index <- which(vars.up1$Day == day & vars.up1$Time == time)
input.up1 <- read.predictedmatrix(predictfiles.up1[valid.index])
# up2
valid.index <- which(vars.up2$Day == day & vars.up2$Time == time)
input.up2 <- read.predictedmatrix(predictfiles.up2[valid.index])
# surf
valid.index <- which(vars.surf$Day == day & vars.surf$Time == time)
input.surf <- read.originalmatrix(originalfiles.surf[valid.index])
# total radiance data
valid.index <- which(vars.total$Day == day & vars.total$Time == time)
input.total <- read.originalmatrix(totalInputfiles[valid.index])
# Lambertian reflectivity retrieval
for(i in 2:length(input.total)){ # input.total[[1]] is for reflectivity 0
  retrieved_reflectivity <- as.matrix(emissivity.retrieval(input.total[[i]], input.up1[[i]], input.up2[[i]], input.down[[i-1]], input.trans[[i]]))
  plot_retrievedreflectivity_boxplot(retrieved_reflectivity, reflectivity.default[i]/100, 
                                     paste("~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/","Lb_retrievedReflect_", reflectivity.default[i], ".jpg", sep = ""), 
                                     ylim2=c(0, 0.1))
}


# Polyethylene reflectivity retrieval
# predicted trans, down, up1, up2
polyethylene.trans<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_trans.csv")[[1]]
polyethylene.down<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_down.csv")[[1]]
polyethylene.up1<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_up1.csv")[[1]]
polyethylene.up2<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_up2.csv")[[1]]
# ground-truth simulated surface emission and total radiance
originalpolyethylene.surf <- read.originalmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_surf.csv")[[1]]
polyethylene.total <- read.originalmatrix('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_total.csv')[[1]]
# retrieve the reflectivity of the polyethylene and plot it against ground truth values
retrieved_reflectivity <- as.matrix(emissivity.retrieval(polyethylene.total, polyethylene.up1, polyethylene.up2, polyethylene.down, polyethylene.trans))
polyethylene_reflectivity <- read.csv('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/polyethylene_info/polyethtlene_reflectivity.csv', row.names = 1)
# plot_retrievedreflectivity(retrieved_reflectivity, polyethylene_reflectivity, "~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/Polyethylene_retrievedReflect.jpg")
# if returnRMSE, then store them in RMSE.differentmaterials
plot_retrievedreflectivity_boxplot(retrieved_reflectivity, polyethylene_reflectivity[,2], '~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/Polyethylene_retrievedReflectBoxPlot.jpg',
                                   wv=polyethylene_reflectivity[,1], ylim1=c(0, 0.3), ylim2=c(0, 0.1), ablineRMSE = 0.08)

# Acetone reflectivity retrieval
# predicted trans, down, up1, up2
acetone.trans<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/predictedCSV/Radiance_107_14_trans.csv")[[1]]
acetone.down<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/predictedCSV/Radiance_107_14_down.csv")[[1]]
acetone.up1<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/predictedCSV/Radiance_107_14_up1.csv")[[1]]
acetone.up2<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/predictedCSV/Radiance_107_14_up2.csv")[[1]]
acetone.total <- read.originalmatrix('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/Radiance_107_14_total.csv')[[1]]
# retrieve the reflectivity of the polyethylene and plot it against ground truth values
retrieved_reflectivity <- emissivity.retrieval(acetone.total, acetone.up1, acetone.up2, acetone.down, acetone.trans)
acetone_reflectivity <- read.csv('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Acetone/acetone_info/acetone_absorptance.csv', row.names = 1)
acetone_reflectivity$reflec <- 1-acetone_reflectivity$absorp
acetone_reflectivity <- acetone_reflectivity[, 3]
plot_retrievedreflectivity_boxplot(retrieved_reflectivity, acetone_reflectivity, '~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/Acetone_retrievedReflectBoxPlot.jpg')

# Ammonia reflectivity retrieval
# predicted trans, down, up1, up2
ammonia.trans<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/predictedCSV/Radiance_107_14_trans.csv")[[1]] #unnamed dataframe
ammonia.down<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/predictedCSV/Radiance_107_14_down.csv")[[1]]
ammonia.up1<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/predictedCSV/Radiance_107_14_up1.csv")[[1]]
ammonia.up2<- read.predictedmatrix("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/predictedCSV/Radiance_107_14_up2.csv")[[1]]
ammonia.total <- read.originalmatrix('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/Radiance_107_14_total.csv')[[1]]
# retrieve the reflectivity of the polyethylene and plot it against ground truth values
retrieved_reflectivity <- as.matrix(emissivity.retrieval(ammonia.total, ammonia.up1, ammonia.up2, ammonia.down, ammonia.trans))
# read ammonia reflectivity retreived from the NIST book
ammonia_reflectivity <- read.csv('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Ammonia/ammonia_info/ammonia_absorptance.csv', row.names = 1)
ammonia_reflectivity$reflec <- 1-ammonia_reflectivity$absorp
ammonia_reflectivity <- ammonia_reflectivity[, 3]
plot_retrievedreflectivity_boxplot(retrieved_reflectivity, ammonia_reflectivity, '~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/Ammonia_retrievedReflectBoxPlot.jpg')


# Create empty vector to store the RMSE of retrieved reflectivity for different materials
RMSE.differentmaterials <- vector('numeric')
# exoscan reflectivity retrieval which contains subfolders of different materialsz
originalfolders.exoscan <- list.files('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/exoscan/originalCSV', full.names = TRUE)
predictfolders.exoscan <- list.files('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/exoscan/predictedCSV', full.names = TRUE)
reflectivities.exoscan <- list.files('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/exoscan/exoscan_info', pattern="*.csv", full.names = TRUE)
for(i in 1:length(originalfolders.exoscan)){
  # predictCSV each forlder has 6 predicted components for a specific material
  # originalCSV each folder has 7 components including the total radiance for a specific material
  originalfolder.exoscan <- originalfolders.exoscan[i]
  predictfolder.exoscan <- predictfolders.exoscan[i]
  reflectivity <- read.csv(reflectivities.exoscan[i], row.names = 1)[,2] # two columns, whose first column is wv while second column is reflectivity
  if(max(reflectivity[,2]) > 0.5){
    ylim=c(0, 1.1)
  }else if(max(reflectivity[,2]) > 0.3){
    ylim=c(0, 0.5)
  }else{
    ylim=c(0, 0.3)
  }
  # read predicted or original files in the folder
  predict.trans<- read.predictedmatrix(paste(predictfolder.exoscan, "/Radiance_107_14_trans.csv", sep=""))[[1]] #unnamed dataframe
  predict.down<- read.predictedmatrix(paste(predictfolder.exoscan, "/Radiance_107_14_down.csv", sep=""))[[1]]
  predict.up1<- read.predictedmatrix(paste(predictfolder.exoscan, "/Radiance_107_14_up1.csv", sep=""))[[1]]
  predict.up2<- read.predictedmatrix(paste(predictfolder.exoscan, "/Radiance_107_14_up2.csv", sep=""))[[1]]
  original.total <- read.originalmatrix(paste(originalfolder.exoscan, "/Radiance_107_14_total.csv", sep=""))[[1]]
  # retrieve the reflectivity
  retrieved_reflectivity <- as.matrix(emissivity.retrieval(original.total, predict.up1, predict.up2, predict.down, predict.trans))
  # plot
  RMSE.differentmaterials <-cbind(RMSE.differentmaterials, plot_retrievedreflectivity_boxplot(retrieved_reflectivity, reflectivity, paste('~/geolab_storage_V3/data/DARPA/Plots/retrievedReflectance/', basename(originalfolder.exoscan),'_retrievedReflectBoxPlot.jpg', sep=""),
                                     ylim1=ylim, ylim2=c(0, 0.1),returnRMSE = TRUE))
}

# random error added
for(i in 1:length(reflectivities.exoscan)){
  reflectivity <- read.csv(reflectivities.exoscan[i], row.names = 1)
  plot_randomerror_reflectivity_boxplot(reflectivity, paste('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/exoscan/exoscan_info/', 
                                              sub('\\.csv$', '', basename(reflectivities.exoscan))[i], '_randomError.jpg', sep=""))
}
# write the average spectral RMSE into CSV
colnames(RMSE.differentmaterials) <- c("Polyethylene", basename(originalfolders.exoscan))
avgRMSE.differentmaterials <- colMeans(RMSE.differentmaterials)
write.csv(avgRMSE.differentmaterials, "averageRMSE_differentMaterials.csv")


############################################################################################################################################
############################################################ Blue Heron Extracted Data #####################################################
############################################################################################################################################
############################################################################################################################################
day = 108; time=14
trainingfolders <- '~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHSimulatedMatrix'
training.trans.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_trans.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
training.down.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_down.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
training.up1.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_up1.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
training.up2.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_up2.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
training.surf.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_surf.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
training.total.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_total.csv", sep='') %>% gsub("\\{day\\}", day, .) %>% gsub("\\{time\\}", time, .),  trainingfolders)
# read
training.trans <- read.originalmatrix(training.trans.files, wv.BH.default, an.BH.default)
training.down <- read.originalmatrix(training.down.files, wv.BH.default, an.BH.default)
training.up1 <- read.originalmatrix(training.up1.files, wv.BH.default, an.BH.default)
training.up2 <- read.originalmatrix(training.up2.files, wv.BH.default, an.BH.default)
training.surf <- read.originalmatrix(training.surf.files, wv.BH.default, an.BH.default)
training.total <- read.originalmatrix(training.total.files, wv.BH.default, an.BH.default)
# BH read
BHExtracted.total<-  read.csv("~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted.csv", row.names = 1, header = TRUE,check.names=FALSE, stringsAsFactors = FALSE)
wv=as.numeric(rownames(BHExtracted.total))
an=as.numeric(colnames(BHExtracted.total))
BHExtracted.down<- rowcol.renames(read.csv("~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/Radiance_down.csv", header = FALSE, check.names=FALSE, stringsAsFactors = FALSE), wv, an)
BHExtracted.trans<- rowcol.renames(read.csv("~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork//Radiance_trans.csv", header = FALSE, check.names=FALSE, stringsAsFactors = FALSE), wv, an)
BHExtracted.up1<- rowcol.renames(read.csv("~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/Radiance_up1.csv", header = FALSE, check.names=FALSE, stringsAsFactors = FALSE), wv, an)
BHExtracted.up2 <- rowcol.renames(read.csv("~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/Radiance_up2.csv", header = FALSE, check.names=FALSE, stringsAsFactors = FALSE), wv, an)
# get zlim for the matrix plot
training.minmax.total <- minmax.crosstables(training.total)
minmax.total <- c(min(training.minmax.total[1], min(BHExtracted.total)), max(training.minmax.total[2], max(BHExtracted.total)))
training.minmax.trans <- minmax.crosstables(training.trans)
minmax.trans <- c(min(training.minmax.trans[1], min(BHExtracted.trans)), max(training.minmax.trans[2], max(BHExtracted.trans)))
training.minmax.down <- minmax.crosstables(training.down)
minmax.down <- c(min(training.minmax.down[1], min(BHExtracted.down)), max(training.minmax.down[2], max(BHExtracted.down)))
training.minmax.up1 <- minmax.crosstables(training.up1)
minmax.up1 <- c(min(training.minmax.up1[1], min(BHExtracted.up1)), max(training.minmax.up1[2], max(BHExtracted.up1)))
training.minmax.up2 <- minmax.crosstables(training.up2)
minmax.up2 <- c(min(training.minmax.up2[1], min(BHExtracted.up2)), max(training.minmax.up2[2], max(BHExtracted.up2)))
# surf, which is not predicted
minmax.surf <- minmax.crosstables(training.surf)

# generate matrix plot
outputfolder= "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/SimulatedPlots/"
for(i in 1:length(training.total)){
  MODTRAN.matplot(as.matrix(training.total[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                  zlim=minmax.total, filename = paste(outputfolder, "TotalRadiance/", names(training.total)[i], ".jpg", sep="") )
  MODTRAN.matplot(as.matrix(training.trans[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), zname = "Transmission", 
                  zlim=c(0,1), filename = paste(outputfolder,"Transmission/", names(training.trans)[i], ".jpg", sep=""))
  MODTRAN.matplot(as.matrix(training.down[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                  zlim=minmax.down, filename = paste(outputfolder, "Downwelling/", names(training.down)[i], ".jpg", sep="") )
  MODTRAN.matplot(as.matrix(training.up1[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                  zlim=minmax.up1, filename = paste(outputfolder, "PathThermalEmission/", names(training.up1)[i], ".jpg", sep="") )
  MODTRAN.matplot(as.matrix(training.up2[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                  zlim=minmax.up2, filename = paste(outputfolder, "PathThermalScattering/", names(training.up2)[i], ".jpg", sep="") )
  MODTRAN.matplot(as.matrix(training.surf[[i]]), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                  zlim=minmax.surf, filename = paste(outputfolder, "SurfaceEmission/", names(training.surf)[i], ".jpg", sep="") )
}

MODTRAN.matplot(as.matrix(BHExtracted.total), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), 
                zlim=minmax.total, filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted_total.jpg")
MODTRAN.matplot(as.matrix(BHExtracted.trans), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5), zname = "Transmission", zlim=c(0,1),
                filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted_trans.jpg")
MODTRAN.matplot(as.matrix(BHExtracted.down), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5),
                zlim=minmax.down, filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted_down.jpg")
MODTRAN.matplot(as.matrix(BHExtracted.up1), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5),
                zlim=minmax.up1, filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted_up1.jpg")
MODTRAN.matplot(as.matrix(BHExtracted.up2), pretty.wv = seq(7.5, 13.5, 0.5), pretty.an = seq(30, 60, 5),
                zlim=minmax.up2, filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHExtracted_up2.jpg")
  

BH.reflectivity <- emissivity.retrieval(BHExtracted.total, BHExtracted.up1, BHExtracted.up2, BHExtracted.down, BHExtracted.trans)
plot_retrievedreflectivity_boxplot(BH.reflectivity, filename = "~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHExtracted/1_BHnetwork/BHRetrievedReflectivity.jpg", wv=wv, an=an, pretty.wv = seq(7.5, 13.5, 0.5)) 

