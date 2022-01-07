# R code to create year-by-year tables from original "snapshot" input tables
# Written by S. Gilbert, last updated 1-1-2022



# Load R libraries needed to run analyses ---------------------------------

rm(list = ls())				# clear r's memory of any stored objects 							  
library(popbio)				# popbio has a lot of good tools for population dynamics
library(reshape2)			# this lets us reshape matrices and arrays (change the dimensions)
library(abind)				# this lets us bind arrays (like matrices) into 3-D stacks, like pages of a book


# Define directories where data are stored --------------------------------
# This requires input data for habitat variables such as deer K and road length/density
# if these change over time, need to load that too
# make sure to set your working directory!

Dir <- "./data/"	                                            # Where on the computer is the main folder with input tables?
DirRoad <- paste(Dir, "roads", sep="")												# Where in main folder are the .csvs for roads for each scenario?  
DirK <- paste(Dir, "deer-k", sep="")													# Where in main folder are the .csvs for deer K "


# Make scenarios for the future ---------------------------------

roads.Base = read.table(paste(DirRoad, "Road_No_Change.csv", sep="/"), sep=",", header=T)
roads.A = read.table(paste(DirRoad, "Road_Decom_Plan.csv", sep="/"), sep=",", header=T)
roads.B <- roads.C <- roads.D <- roads.A
roads.E = read.table(paste(DirRoad, "Road_Construction.csv", sep="/"), sep=",", header=T)
roads.mid.decom = read.table(paste(DirRoad, "Road_Decom_Mid.csv", sep="/"), sep=",", header=T)
roads.max.decom = read.table(paste(DirRoad, "Road_Decom_Max.csv", sep="/"), sep=",", header=T)

K.Base = read.table(paste(DirK, "K_Scen_Baseline.csv", sep="/"), sep=",", header=T)
K.Nochange = read.table(paste(DirK, "K_Scen_Nochange.csv", sep="/"), sep=",", header=T)
K.A <- K.Base
K.B = read.table(paste(DirK, "K_Scen_B.csv", sep="/"), sep=",", header=T)
K.C = read.table(paste(DirK, "K_Scen_C.csv", sep="/"), sep=",", header=T)
K.D = read.table(paste(DirK, "K_Scen_D.csv", sep="/"), sep=",", header=T)
K.E = read.table(paste(DirK, "K_Scen_E.csv", sep="/"), sep=",", header=T)

# Create a function to fill in data columns for missing yrs via extrapolation ---------------------------------

fill.in.yrs <- function(data, years.row, years.col){				# tell it which data, which row contains yrs (should be 1), which cols contain data
dat.row = 2:nrow(data)
years = as.numeric(data[years.row,years.col])
dat = data[dat.row, years.col]

for(i in 2:length(years)){
	year.seq = years[i-1]:years[i]
	run = years[i]-years[i-1]
	rise = dat[,i] - dat[,(i-1)]
	slope = rise/run
		for(k in 2:(length(year.seq)-1)){
			if(k==2){val.new = dat[,i-1] + slope}else{val.new=val.new+slope}
			if(k==2){fill.yrs=val.new}else{fill.yrs=cbind(fill.yrs, val.new)}
			}
	filled = data.frame(dat[,i-1], fill.yrs)
	if(i==2){filled.dat = filled}else{filled.dat = cbind(filled.dat, filled)}
}
filled.dat = cbind(filled.dat, dat[,length(years)])
colnames(filled.dat)<- as.character(as.numeric(years[1]):as.numeric(years[length(years)]))

return(filled.dat)
}

# fill in missing years columns for roads and deer K ---------------------------------

fill.roads.Base <- fill.in.yrs(roads.Base, 1, 2:5)								# fill in year values for roads.Base
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.Base), sep=""))
fill.roads.Base <- cbind(seq(1:nrow(fill.roads.Base)), fill.roads.Base)
colnames(fill.roads.Base) <- name

fill.roads.A <- fill.in.yrs(roads.A, 1, 2:5)
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.A), sep=""))				# fill in year values for roads.A
fill.roads.A <- cbind(seq(1:nrow(fill.roads.A)), fill.roads.A)
colnames(fill.roads.A) <- name

fill.roads.B <- fill.in.yrs(roads.B, 1, 2:5)
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.B), sep=""))				# fill in year values for roads.B
fill.roads.B <- cbind(seq(1:nrow(fill.roads.B)), fill.roads.B)
colnames(fill.roads.B) <- name

fill.roads.C <- fill.in.yrs(roads.C, 1, 2:5)									# fill in year values for roads.C
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.C), sep=""))
fill.roads.C <- cbind(seq(1:nrow(fill.roads.C)), fill.roads.C)
colnames(fill.roads.C) <- name

fill.roads.D <- fill.in.yrs(roads.D, 1, 2:5)									# fill in year values for roads.D
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.D), sep=""))
fill.roads.D <- cbind(seq(1:nrow(fill.roads.D)), fill.roads.D)
colnames(fill.roads.D) <- name

fill.roads.E <- fill.in.yrs(roads.E, 1, 2:5)									# fill in year values for roads.E
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.E), sep=""))
fill.roads.E <- cbind(seq(1:nrow(fill.roads.E)), fill.roads.E)
colnames(fill.roads.E) <- name

fill.roads.mid.decom <- fill.in.yrs(roads.mid.decom, 1, 2:5)									# fill in year values for roads.mid.decom
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.mid.decom), sep=""))
fill.roads.mid.decom <- cbind(seq(1:nrow(fill.roads.mid.decom)), fill.roads.mid.decom)
colnames(fill.roads.mid.decom) <- name

fill.roads.max.decom <- fill.in.yrs(roads.max.decom, 1, 2:5)									# fill in year values for roads.max.decom
name <- c("Pack_ID", paste("Road_", colnames(fill.roads.max.decom), sep=""))
fill.roads.max.decom <- cbind(seq(1:nrow(fill.roads.max.decom)), fill.roads.max.decom)
colnames(fill.roads.max.decom) <- name

fill.K.Base <- fill.in.yrs(K.Base, 1, 2:5)										# fill in year values for K.Base
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.Base), sep=""))
fill.K.Base <- cbind(seq(1:nrow(fill.K.Base)), fill.K.Base)
colnames(fill.K.Base) <- name

fill.K.Nochange <- fill.in.yrs(K.Nochange, 1, 2:5)										# fill in year values for K.Nochange
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.Nochange), sep=""))
fill.K.Nochange <- cbind(seq(1:nrow(fill.K.Nochange)), fill.K.Nochange)
colnames(fill.K.Nochange) <- name

fill.K.A <- fill.in.yrs(K.A, 1, 2:5)											# fill in year values for K.A
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.A), sep=""))
fill.K.A <- cbind(seq(1:nrow(fill.K.A)), fill.K.A)
colnames(fill.K.A) <- name

fill.K.B <- fill.in.yrs(K.B, 1, 2:5)											# fill in year values for K.B
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.B), sep=""))
fill.K.B <- cbind(seq(1:nrow(fill.K.B)), fill.K.B)
colnames(fill.K.B) <- name

fill.K.C <- fill.in.yrs(K.C, 1, 2:5)											# fill in year values for K.C
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.C), sep=""))
fill.K.C <- cbind(seq(1:nrow(fill.K.C)), fill.K.C)
colnames(fill.K.C) <- name

fill.K.D <- fill.in.yrs(K.D, 1, 2:5)											# fill in year values for K.D
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.D), sep=""))
fill.K.D <- cbind(seq(1:nrow(fill.K.D)), fill.K.D)
colnames(fill.K.D) <- name

fill.K.E <- fill.in.yrs(K.E, 1, 2:5)											# fill in year values for K.E
name <- c("Pack_ID", paste("Deer_", colnames(fill.K.E), sep=""))
fill.K.E <- cbind(seq(1:nrow(fill.K.E)), fill.K.E)
colnames(fill.K.E) <- name

# now write output, processed data, to file---------------------------------


write.csv(fill.roads.Base, paste("./data-processed/roads", "Data.roads.base.csv", sep="/"), row.names=F)
write.csv(fill.roads.A, paste("./data-processed/roads", "Data.roads.A.csv", sep="/"), row.names=F)
write.csv(fill.roads.B, paste("./data-processed/roads", "Data.roads.B.csv", sep="/"), row.names=F)
write.csv(fill.roads.C, paste("./data-processed/roads", "Data.roads.C.csv", sep="/"), row.names=F)
write.csv(fill.roads.D, paste("./data-processed/roads", "Data.roads.D.csv", sep="/"), row.names=F)
write.csv(fill.roads.E, paste("./data-processed/roads", "Data.roads.E.csv", sep="/"), row.names=F)
write.csv(fill.roads.mid.decom, paste("./data-processed/roads", "Data.roads.mid.decom.csv", sep="/"), row.names=F)
write.csv(fill.roads.max.decom, paste("./data-processed/roads", "Data.roads.max.decom.csv", sep="/"), row.names=F)

write.csv(fill.K.Base, paste("./data-processed/deer-k", "Data.K.Base.csv", sep="/"), row.names=F)
write.csv(fill.K.Nochange, paste("./data-processed/deer-k", "Data.K.Nochange.csv", sep="/"), row.names=F)
write.csv(fill.K.A, paste("./data-processed/deer-k", "Data.K.A.csv", sep="/"), row.names=F)
write.csv(fill.K.B, paste("./data-processed/deer-k", "Data.K.B.csv", sep="/"), row.names=F)
write.csv(fill.K.C, paste("./data-processed/deer-k", "Data.K.C.csv", sep="/"), row.names=F)
write.csv(fill.K.D, paste("./data-processed/deer-k", "Data.K.D.csv", sep="/"), row.names=F)
write.csv(fill.K.E, paste("./data-processed/deer-k", "Data.K.E.csv", sep="/"), row.names=F)





