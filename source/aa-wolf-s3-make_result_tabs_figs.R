### R code to create figures and calculate summary statistics for the Alexander Archipelago Wolf Project
### By Sophie Gilbert 
### Last updated 1-1-2022



#### Libraries to be loaded (if already installed), otherwise, need to install----


rm(list = ls())				# clear r's memory of any stored objects (for instance, if we've been manually messing around with names in r). 
							# This is always a good idea at the start, FYI
library(ggplot2)			# for plotting
library(gridExtra)			# for plotting
library(fields)
library(colorRamps)



#### First define where to get the data------

Dir.scen =  "./results/"															# where within g drive is the scenario data?
Dir.out = "./results/tabs/"



#### Create a function to calculate probability of quasi-extinction, when different thresholds are chosen----

extinction.prob <- function(data, q.thresh){
for(i in 1:ncol(data)){
	if(length(data[which(data[,i]<=q.thresh),i])>0){q.count=1}else{q.count=0} # does pop cross thresh. in iteration i?
	if(length(data[which(data[,i]<=1),i])>0){e.count=1}else{e.count=0}		  # does pop <=1 in iteration i? 
	p.count=(length(data[which(data[,i]<=q.thresh),i]))/nrow(data)			  # what percentage of years below thresh?
	if(i==1){q.tot = q.count}else{q.tot=c(q.tot, q.count)}
	if(i==1){e.tot = e.count}else{e.tot=c(e.tot, e.count)}
	if(i==1){p.tot =p.count}else{p.tot = c(p.tot, p.count)}
}
q.total = sum(q.tot)/ncol(data)
e.total = sum(e.tot)/ncol(data)
p.total = sum(p.tot)/ncol(data)
q.lcl = q.total -1.96*(sqrt(q.total/(1-q.total)/ncol(data)))
q.ucl = q.total +1.96*(sqrt(q.total/(1-q.total)/ncol(data)))
e.lcl = e.total -1.96*(sqrt(e.total/(1-e.total)/ncol(data)))
e.ucl = e.total +1.96*(sqrt(e.total/(1-e.total)/ncol(data)))
p.lcl = quantile(p.tot, probs=0.025)
p.ucl = quantile(p.tot, probs=0.975)


q.output = data.frame(q.threshold.n = q.thresh,
						q.prob = q.total, q.lcl = q.lcl, q.ucl=q.ucl,
						e.prob=e.total, e.lcl = e.lcl, e.ucl = e.ucl,
						perc.yrs = p.total, perc.lcl = p.lcl, perc.ucl=p.ucl)
return(q.output)
}

## A function to calculate percentage change by the end of timesteps

pct.change.end <-function(data){
	Lcl = min(data$lcl[nrow(data)], data$ucl[nrow(data)])
	Ucl = max(data$lcl[nrow(data)], data$ucl[nrow(data)])
	pct.out <- data.frame(Med.pct=(data$median[nrow(data)]-data$median[1])/data$median[1], Lcl.pct=(Lcl-data$median[1])/data$median[1], Ucl.pct=(Ucl-data$median[1])/data$median[1])
	return(pct.out)
}



##### Plot of empirical change in wolf abundance
dat.emp = read.csv(paste(Dir, "Pow_trend.csv", sep=""))

par(las=1)
plot(dat.emp$year, dat.emp$popN.est, pch=19, col="black", 
					xlab="Year", ylab="Fall wolf population size",
					ylim=c(0, 600))
arrows(dat.emp$year, dat.emp$popN.est, dat.emp$year, dat.emp$lcl, length=0.05, angle=90)					
arrows(dat.emp$year, dat.emp$popN.est, dat.emp$year, dat.emp$ucl, length=0.05, angle=90)



############## Create scenarios results ###########################


## Read in the data from scenario results----

wolf.base.s <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
wolf.A.s <- read.table(paste(Dir.scen, "scenarios", "ScenA", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.B.s <- read.table(paste(Dir.scen, "scenarios", "ScenB", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.C.s <- read.table(paste(Dir.scen, "scenarios", "ScenC", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.D.s <- read.table(paste(Dir.scen, "scenarios", "ScenD", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.E.s <- read.table(paste(Dir.scen, "scenarios", "ScenE", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.base <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
wolf.A <- read.table(paste(Dir.scen, "scenarios", "ScenA", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.B <- read.table(paste(Dir.scen, "scenarios", "ScenB", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.C <- read.table(paste(Dir.scen, "scenarios", "ScenC", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.D <- read.table(paste(Dir.scen, "scenarios", "ScenD", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.E <- read.table(paste(Dir.scen, "scenarios", "ScenE", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

deer.base <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
deer.A <- read.table(paste(Dir.scen, "scenarios", "ScenA", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.B <- read.table(paste(Dir.scen, "scenarios", "ScenB", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.C <- read.table(paste(Dir.scen, "scenarios", "ScenC", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.D <- read.table(paste(Dir.scen, "scenarios", "ScenD", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.E <- read.table(paste(Dir.scen, "scenarios", "ScenE", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

packs.base <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
packs.A <- read.table(paste(Dir.scen, "scenarios", "ScenA", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
packs.B <- read.table(paste(Dir.scen, "scenarios", "ScenB", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
packs.C <- read.table(paste(Dir.scen, "scenarios", "ScenC", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
packs.D <- read.table(paste(Dir.scen, "scenarios", "ScenD", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
packs.E <- read.table(paste(Dir.scen, "scenarios", "ScenE", "extinct.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.base.monte <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
wolf.A.monte <- read.table(paste(Dir.scen, "scenarios", "ScenA", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.B.monte <- read.table(paste(Dir.scen, "scenarios", "ScenB", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.C.monte <- read.table(paste(Dir.scen, "scenarios", "ScenC", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.D.monte <- read.table(paste(Dir.scen, "scenarios", "ScenD", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.E.monte <- read.table(paste(Dir.scen, "scenarios", "ScenE", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

hunt.base <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
hunt.A <- read.table(paste(Dir.scen, "scenarios", "ScenA", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.B <- read.table(paste(Dir.scen, "scenarios", "ScenB", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.C <- read.table(paste(Dir.scen, "scenarios", "ScenC", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.D <- read.table(paste(Dir.scen, "scenarios", "ScenD", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.E <- read.table(paste(Dir.scen, "scenarios", "ScenE", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

browse.base <- read.table(paste(Dir.scen, "scenarios", "ScenBase", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
browse.A <- read.table(paste(Dir.scen, "scenarios", "ScenA", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.B <- read.table(paste(Dir.scen, "scenarios", "ScenB", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.C <- read.table(paste(Dir.scen, "scenarios", "ScenC", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.D <- read.table(paste(Dir.scen, "scenarios", "ScenD", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.E <- read.table(paste(Dir.scen, "scenarios", "ScenE", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)



### Calculate percent change for deer, wolves, hunting (browsing % change makes no sense)----

d.p.base <- pct.change.end(deer.base)
d.p.A <- pct.change.end(deer.A)
d.p.B <- pct.change.end(deer.B)
d.p.C <- pct.change.end(deer.C)
d.p.D <- pct.change.end(deer.D)
d.p.E <- pct.change.end(deer.E)

w.p.base <- pct.change.end(wolf.base)
w.p.A <- pct.change.end(wolf.A)
w.p.B <- pct.change.end(wolf.B)
w.p.C <- pct.change.end(wolf.C)
w.p.D <- pct.change.end(wolf.D)
w.p.E <- pct.change.end(wolf.E)

w.p.base.s <- pct.change.end(wolf.base.s)
w.p.A.s <- pct.change.end(wolf.A.s)
w.p.B.s <- pct.change.end(wolf.B.s)
w.p.C.s <- pct.change.end(wolf.C.s)
w.p.D.s <- pct.change.end(wolf.D.s)
w.p.E.s <- pct.change.end(wolf.E.s)

h.p.base <- pct.change.end(hunt.base)
h.p.A <- pct.change.end(hunt.A)
h.p.B <- pct.change.end(hunt.B)
h.p.C <- pct.change.end(hunt.C)
h.p.D <- pct.change.end(hunt.D)
h.p.E <- pct.change.end(hunt.E)


### Make scenario results tables, deer, wolves, hunting, browsing, pack extinction ----

d.scen.tab <- data.frame(Scen = c("base", "A", "B", "C", "D", "E"), 
	Med.pct = c(d.p.base$Med.pct, d.p.A$Med.pct, d.p.B$Med.pct, d.p.C$Med.pct, d.p.D$Med.pct, d.p.E$Med.pct),
	Lcl.pct  = c(min(d.p.base$Lcl.pct, d.p.base$Ucl.pct), min(d.p.A$Lcl.pct, d.p.A$Ucl.pct), min(d.p.B$Lcl.pct, d.p.B$Ucl.pct), min(d.p.C$Lcl.pct, d.p.C$Ucl.pct), min(d.p.D$Lcl.pct, d.p.D$Ucl.pct), min(d.p.E$Lcl.pct, d.p.E$Ucl.pct)),
	Ucl.pct  = c(max(d.p.base$Lcl.pct, d.p.base$Ucl.pct), max(d.p.A$Lcl.pct, d.p.A$Ucl.pct), max(d.p.B$Lcl.pct, d.p.B$Ucl.pct), max(d.p.C$Lcl.pct, d.p.C$Ucl.pct), max(d.p.D$Lcl.pct, d.p.D$Ucl.pct), max(d.p.E$Lcl.pct, d.p.E$Ucl.pct)),
	Med = c(deer.base$med[31], deer.A$med[31], deer.B$med[31], deer.C$med[31], deer.D$med[31], deer.E$med[31]),
	Lcl = c(deer.base$lcl[31], deer.A$lcl[31], deer.B$lcl[31], deer.C$lcl[31], deer.D$lcl[31], deer.E$lcl[31]),
	Ucl = c(deer.base$ucl[31], deer.A$ucl[31], deer.B$ucl[31], deer.C$ucl[31], deer.D$ucl[31], deer.E$ucl[31]))
	
	
w.scen.tab <- data.frame(Scen = c("base", "A", "B", "C", "D", "E"), 
	Med.pct = c(w.p.base$Med.pct, w.p.A$Med.pct, w.p.B$Med.pct, w.p.C$Med.pct, w.p.D$Med.pct, w.p.E$Med.pct),
	Lcl.pct  = c(min(w.p.base$Lcl.pct, w.p.base$Ucl.pct), min(w.p.A$Lcl.pct, w.p.A$Ucl.pct), min(w.p.B$Lcl.pct, w.p.B$Ucl.pct), min(w.p.C$Lcl.pct, w.p.C$Ucl.pct), min(w.p.D$Lcl.pct, w.p.D$Ucl.pct), min(w.p.E$Lcl.pct, w.p.E$Ucl.pct)),
	Ucl.pct  = c(max(w.p.base$Lcl.pct, w.p.base$Ucl.pct), max(w.p.A$Lcl.pct, w.p.A$Ucl.pct), max(w.p.B$Lcl.pct, w.p.B$Ucl.pct), max(w.p.C$Lcl.pct, w.p.C$Ucl.pct), max(w.p.D$Lcl.pct, w.p.D$Ucl.pct), max(w.p.E$Lcl.pct, w.p.E$Ucl.pct)),
	Med = c(wolf.base$med[31], wolf.A$med[31], wolf.B$med[31], wolf.C$med[31], wolf.D$med[31], wolf.E$med[31]),
	Lcl = c(wolf.base$lcl[31], wolf.A$lcl[31], wolf.B$lcl[31], wolf.C$lcl[31], wolf.D$lcl[31], wolf.E$lcl[31]),
	Ucl = c(wolf.base$ucl[31], wolf.A$ucl[31], wolf.B$ucl[31], wolf.C$ucl[31], wolf.D$ucl[31], wolf.E$ucl[31]))
	
w.scen.tab.s <- data.frame(Scen = c("base", "A", "B", "C", "D", "E"), 
	Med.pct = c(w.p.base.s$Med.pct, w.p.A.s$Med.pct, w.p.B.s$Med.pct, w.p.C.s$Med.pct, w.p.D.s$Med.pct, w.p.E.s$Med.pct),
	Lcl.pct  = c(min(w.p.base.s$Lcl.pct, w.p.base.s$Ucl.pct), min(w.p.A.s$Lcl.pct, w.p.A.s$Ucl.pct), min(w.p.B.s$Lcl.pct, w.p.B.s$Ucl.pct), min(w.p.C.s$Lcl.pct, w.p.C.s$Ucl.pct), min(w.p.D.s$Lcl.pct, w.p.D.s$Ucl.pct), min(w.p.E.s$Lcl.pct, w.p.E.s$Ucl.pct)),
	Ucl.pct  = c(max(w.p.base.s$Lcl.pct, w.p.base.s$Ucl.pct), max(w.p.A.s$Lcl.pct, w.p.A.s$Ucl.pct), max(w.p.B.s$Lcl.pct, w.p.B.s$Ucl.pct), max(w.p.C.s$Lcl.pct, w.p.C.s$Ucl.pct), max(w.p.D.s$Lcl.pct, w.p.D.s$Ucl.pct), max(w.p.E.s$Lcl.pct, w.p.E.s$Ucl.pct)),
	Med = c(wolf.base.s$med[31], wolf.A.s$med[31], wolf.B.s$med[31], wolf.C.s$med[31], wolf.D.s$med[31], wolf.E.s$med[31]),
	Lcl = c(wolf.base.s$lcl[31], wolf.A.s$lcl[31], wolf.B.s$lcl[31], wolf.C.s$lcl[31], wolf.D.s$lcl[31], wolf.E.s$lcl[31]),
	Ucl = c(wolf.base.s$ucl[31], wolf.A.s$ucl[31], wolf.B.s$ucl[31], wolf.C.s$ucl[31], wolf.D.s$ucl[31], wolf.E.s$ucl[31]))

h.scen.tab <- data.frame(Scen = c("base", "A", "B", "C", "D", "E"), 
	Med.pct = c(h.p.base$Med.pct, h.p.A$Med.pct, h.p.B$Med.pct, h.p.C$Med.pct, h.p.D$Med.pct, h.p.E$Med.pct),
	Lcl.pct  = c(min(h.p.base$Lcl.pct, h.p.base$Ucl.pct), min(h.p.A$Lcl.pct, h.p.A$Ucl.pct), min(h.p.B$Lcl.pct, h.p.B$Ucl.pct), min(h.p.C$Lcl.pct, h.p.C$Ucl.pct), min(h.p.D$Lcl.pct, h.p.D$Ucl.pct), min(h.p.E$Lcl.pct, h.p.E$Ucl.pct)),
	Ucl.pct  = c(max(h.p.base$Lcl.pct, h.p.base$Ucl.pct), max(h.p.A$Lcl.pct, h.p.A$Ucl.pct), max(h.p.B$Lcl.pct, h.p.B$Ucl.pct), max(h.p.C$Lcl.pct, h.p.C$Ucl.pct), max(h.p.D$Lcl.pct, h.p.D$Ucl.pct), max(h.p.E$Lcl.pct, w.p.E$Ucl.pct)),
	Med = c(hunt.base$med[31], hunt.A$med[31], hunt.B$med[31], hunt.C$med[31], hunt.D$med[31], hunt.E$med[31]),
	Lcl = c(hunt.base$lcl[31], hunt.A$lcl[31], hunt.B$lcl[31], hunt.C$lcl[31], hunt.D$lcl[31], hunt.E$lcl[31]),
	Ucl = c(hunt.base$ucl[31], hunt.A$ucl[31], hunt.B$ucl[31], hunt.C$ucl[31], hunt.D$ucl[31], hunt.E$ucl[31]))


p.scen.tab <- data.frame(Scen =  c("base", "A", "B", "C", "D", "E"), 
	Med = c(packs.base$med[31], packs.A$med[31], packs.B$med[31], packs.C$med[31], packs.D$med[31], packs.E$med[31]),
	Lcl = c(packs.base$lcl[31], packs.A$lcl[31], packs.B$lcl[31], packs.C$lcl[31], packs.D$lcl[31], packs.E$lcl[31]),
	Ucl = c(packs.base$ucl[31], packs.A$ucl[31], packs.B$ucl[31], packs.C$ucl[31], packs.D$ucl[31], packs.E$ucl[31]))


b.scen.tab <- data.frame(Scen = c("base", "A", "B", "C", "D", "E"), 

	Med = c(browse.base$med[31], browse.A$med[31], browse.B$med[31], browse.C$med[31], browse.D$med[31], browse.E$med[31]),
	Lcl = c(browse.base$lcl[31], browse.A$lcl[31], browse.B$lcl[31], browse.C$lcl[31], browse.D$lcl[31], browse.E$lcl[31]),
	Ucl = c(browse.base$ucl[31], browse.A$ucl[31], browse.B$ucl[31], browse.C$ucl[31], browse.D$ucl[31], browse.E$ucl[31]))



##### calculate extinction probabilities for a range of quasi-extinction thresholds

q.range = seq(0, 150, 1)

for(i in 1:length(q.range)){
	q.th = q.range[i]
	e.base = extinction.prob(wolf.base.monte, q.thresh= q.th)
	e.A = extinction.prob(wolf.A.monte, q.thresh= q.th)
	e.B = extinction.prob(wolf.B.monte, q.thresh= q.th)
	e.C = extinction.prob(wolf.C.monte, q.thresh= q.th)
	e.D = extinction.prob(wolf.D.monte, q.thresh= q.th)
	e.E = extinction.prob(wolf.E.monte, q.thresh= q.th)
	
	if(i==1){
		ext.base = e.base;
		ext.A = e.A;
		ext.B = e.B;
		ext.C = e.C;
		ext.D = e.D;
		ext.E = e.E
			}else{
				ext.base = rbind(ext.base, e.base);
				ext.A = rbind(ext.A, e.A);
				ext.B = rbind(ext.B, e.B);
				ext.C = rbind(ext.C, e.C);
				ext.D = rbind(ext.D, e.D);
				ext.E = rbind(ext.E, e.E)}

}








############### Create sensitivity analysis results ##################

## Read in the data from sensitivity results----

# wolves, spring

wolf.clim.low.s <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
wolf.clim.mid.s <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.clim.max.s <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.diet.9.5.s <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.15.s <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet15","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.20.5.s <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.26.s <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.nohunt.s <- read.table(paste(Dir.scen, "sensitivities", "SenHunt", "wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.hunt.s <- wolf.B.s

wolf.K.base.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.maxharvest.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.contharvest.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.midharvest.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.nochange.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.transition.s <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KTransition","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.roads.maxbuild.s <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.maxdecom.s <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.middecom.s <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.nochange.s <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.plandecom.s <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadPlandecom","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.trap.0.s <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.20.s <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap20","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.30.s <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.none.s <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.100.s <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.gone.s <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","wolf.spring.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)


# wolves, fall

wolf.clim.low <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
wolf.clim.mid <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.clim.max <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.diet.9.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.15 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet15","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.20.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.26 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.nohunt <- read.table(paste(Dir.scen, "sensitivities", "SenHunt", "wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.hunt <- wolf.B

wolf.K.base <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.maxharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.contharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.midharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.transition <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KTransition","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.roads.maxbuild <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.maxdecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.middecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.plandecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadPlandecom","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.trap.0 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.20 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap20","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.30 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.none <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.100 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.gone <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","wolf.fall.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)




# deer

deer.clim.low <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
deer.clim.mid <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.clim.max <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

deer.diet.9.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.diet.15 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet15","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.diet.20.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.diet.26 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

deer.nohunt <- read.table(paste(Dir.scen, "sensitivities", "SenHunt", "deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.hunt <- deer.B

deer.K.base <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.K.maxharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.K.contharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.K.midharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.K.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.K.transition <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KTransition","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

deer.roads.maxbuild <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.roads.maxdecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.roads.middecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.roads.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.roads.plandecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadPlandecom","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

deer.trap.0 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.trap.20 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap20","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.trap.30 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.trap.none <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.trap.100 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
deer.gone <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","deer.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

# Hunting
hunt.clim.low <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
hunt.clim.mid <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.clim.max <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

hunt.diet.9.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.diet.15 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet15","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.diet.20.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.diet.26 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

hunt.K.base <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.K.maxharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.K.contharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.K.midharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.K.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.K.transition <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KTransition","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

hunt.roads.maxbuild <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.roads.maxdecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.roads.middecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.roads.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.roads.plandecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadPlandecom","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

hunt.trap.0 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.trap.20 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap20","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.trap.30 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.trap.none <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.trap.100 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
hunt.gone <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","hunt.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)


# Browsing
browse.clim.low <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	# read in scenarios results
browse.clim.mid <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.clim.max <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

browse.diet.9.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.diet.15 <- browse.B
browse.diet.20.5 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.diet.26 <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

browse.nohunt <- read.table(paste(Dir.scen, "sensitivities", "SenHunt", "browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.hunt <- browse.B

browse.K.base <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.K.maxharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.K.contharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.K.midharvest <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.K.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.K.transition <- browse.B

browse.roads.maxbuild <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.roads.maxdecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.roads.middecom <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.roads.nochange <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.roads.plandecom <- browse.B

browse.trap.0 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.trap.20 <- browse.B
browse.trap.30 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.trap.none <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.trap.100 <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
browse.gone <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","browse.pack.perc.mn.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

# wolf extinction stuff

wolf.clim.low.monte <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowLow", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)	
wolf.clim.mid.monte <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMid", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.clim.max.monte <- read.table(paste(Dir.scen, "sensitivities", "SenClimate", "SnowMax","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.diet.9.5.monte <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet9.5","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.15.monte <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet15","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.20.5.monte <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet20.5","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.diet.26.monte <- read.table(paste(Dir.scen, "sensitivities", "SenDiet", "Diet26","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.nohunt.monte <- read.table(paste(Dir.scen, "sensitivities", "SenHunt", "wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.hunt.monte <- wolf.B.monte

wolf.K.base.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KBase","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.maxharvest.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMaxharvest","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.contharvest.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KContharvest","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.midharvest.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KMidharvest","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.nochange.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KNochange","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.K.transition.monte <- read.table(paste(Dir.scen, "sensitivities", "SenK", "KTransition","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.roads.maxbuild.monte <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxbuild","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.maxdecom.monte <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMaxdecom","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.middecom.monte <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadMiddecom","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.nochange.monte <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadNochange","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.roads.plandecom.monte <- read.table(paste(Dir.scen, "sensitivities", "SenRoads", "RoadPlandecom","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

wolf.trap.0.monte <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap0","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.20.monte <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap20","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.30.monte <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap30","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.none.monte <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "TrapNone","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.trap.100.monte <- read.table(paste(Dir.scen, "sensitivities", "SenTrap", "Trap100","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)
wolf.gone.monte <- read.table(paste(Dir.scen, "sensitivities", "SenWolf", "WolfGone","wolf.spring.monte.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)


### calculate percent change for deer, wolves, hunting (browsing % change makes no sense)----

d.p.clim.low <- pct.change.end(deer.clim.low)
d.p.clim.mid <- pct.change.end(deer.clim.mid)
d.p.clim.max <- pct.change.end(deer.clim.max)

d.p.diet.9.5 <- pct.change.end(deer.diet.9.5)
d.p.diet.15 <- pct.change.end(deer.diet.15)
d.p.diet.20.5 <- pct.change.end(deer.diet.20.5)
d.p.diet.26 <- pct.change.end(deer.diet.26)

d.p.nohunt <- pct.change.end(deer.nohunt)
d.p.hunt <- pct.change.end(deer.B)

d.p.K.base <- pct.change.end(deer.K.base)
d.p.K.maxharvest <- pct.change.end(deer.K.maxharvest)
d.p.K.midharvest <- pct.change.end(deer.K.midharvest)
d.p.K.contharvest <- pct.change.end(deer.K.contharvest)
d.p.K.nochange <- pct.change.end(deer.K.nochange)
d.p.K.transition <- pct.change.end(deer.K.transition)

d.p.roads.maxdecom <- pct.change.end(deer.roads.maxdecom)
d.p.roads.plandecom <- pct.change.end(deer.roads.plandecom)
d.p.roads.middecom <- pct.change.end(deer.roads.middecom)
d.p.roads.nochange <- pct.change.end(deer.roads.nochange)
d.p.roads.maxbuild <- pct.change.end(deer.roads.maxbuild)

d.p.trap.0 <- pct.change.end(deer.trap.0)
d.p.trap.20 <- pct.change.end(deer.trap.20)
d.p.trap.30 <- pct.change.end(deer.trap.30)
d.p.trap.none <- pct.change.end(deer.trap.none)
d.p.trap.100 <- pct.change.end(deer.trap.100)
d.p.gone <- pct.change.end(deer.gone)

w.p.clim.low <- pct.change.end(wolf.clim.low)
w.p.clim.mid <- pct.change.end(wolf.clim.mid)
w.p.clim.max <- pct.change.end(wolf.clim.max)

w.p.diet.9.5 <- pct.change.end(wolf.diet.9.5)
w.p.diet.15 <- pct.change.end(wolf.diet.15)
w.p.diet.20.5 <- pct.change.end(wolf.diet.20.5)
w.p.diet.26 <- pct.change.end(wolf.diet.26)

w.p.nohunt <- pct.change.end(wolf.nohunt)
w.p.hunt <- pct.change.end(wolf.B)

w.p.K.base <- pct.change.end(wolf.K.base)
w.p.K.maxharvest <- pct.change.end(wolf.K.maxharvest)
w.p.K.contharvest <- pct.change.end(wolf.K.contharvest)
w.p.K.midharvest <- pct.change.end(wolf.K.midharvest)
w.p.K.nochange <- pct.change.end(wolf.K.nochange)
w.p.K.transition <- pct.change.end(wolf.K.transition)

w.p.roads.maxdecom <- pct.change.end(wolf.roads.maxdecom)
w.p.roads.plandecom <- pct.change.end(wolf.roads.plandecom)
w.p.roads.middecom <- pct.change.end(wolf.roads.middecom)
w.p.roads.nochange <- pct.change.end(wolf.roads.nochange)
w.p.roads.maxbuild <- pct.change.end(wolf.roads.maxbuild)

w.p.trap.0 <- pct.change.end(wolf.trap.0)
w.p.trap.20 <- pct.change.end(wolf.trap.20)
w.p.trap.30 <- pct.change.end(wolf.trap.30)
w.p.trap.none <- pct.change.end(wolf.trap.none)
w.p.trap.100 <- pct.change.end(wolf.trap.100)
w.p.gone <- pct.change.end(wolf.gone)

w.p.clim.low.s <- pct.change.end(wolf.clim.low.s)
w.p.clim.mid.s <- pct.change.end(wolf.clim.mid.s)
w.p.clim.max.s <- pct.change.end(wolf.clim.max.s)

w.p.diet.9.5.s <- pct.change.end(wolf.diet.9.5.s)
w.p.diet.15.s <- pct.change.end(wolf.diet.15.s)
w.p.diet.20.5.s <- pct.change.end(wolf.diet.20.5.s)
w.p.diet.26.s <- pct.change.end(wolf.diet.26.s)

w.p.nohunt.s <- pct.change.end(wolf.nohunt.s)
w.p.hunt.s <- pct.change.end(wolf.B.s)

w.p.K.base.s <- pct.change.end(wolf.K.base.s)
w.p.K.maxharvest.s <- pct.change.end(wolf.K.maxharvest.s)
w.p.K.contharvest.s <- pct.change.end(wolf.K.contharvest.s)
w.p.K.midharvest.s <- pct.change.end(wolf.K.midharvest.s)
w.p.K.nochange.s <- pct.change.end(wolf.K.nochange.s)
w.p.K.transition.s <- pct.change.end(wolf.K.transition.s)

w.p.roads.maxdecom.s <- pct.change.end(wolf.roads.maxdecom.s)
w.p.roads.plandecom.s <- pct.change.end(wolf.roads.plandecom.s)
w.p.roads.middecom.s <- pct.change.end(wolf.roads.middecom.s)
w.p.roads.nochange.s <- pct.change.end(wolf.roads.nochange.s)
w.p.roads.maxbuild.s <- pct.change.end(wolf.roads.maxbuild.s)

w.p.trap.0.s <- pct.change.end(wolf.trap.0.s)
w.p.trap.20.s <- pct.change.end(wolf.trap.20.s)
w.p.trap.30.s <- pct.change.end(wolf.trap.30.s)
w.p.trap.none.s <- pct.change.end(wolf.trap.none.s)
w.p.trap.100.s <- pct.change.end(wolf.trap.100.s)
w.p.gone.s <- pct.change.end(wolf.gone.s)

h.p.clim.low <- pct.change.end(hunt.clim.low)
h.p.clim.mid <- pct.change.end(hunt.clim.mid)
h.p.clim.max <- pct.change.end(hunt.clim.max)

h.p.diet.9.5 <- pct.change.end(hunt.diet.9.5)
h.p.diet.15 <- pct.change.end(hunt.diet.15)
h.p.diet.20.5 <- pct.change.end(hunt.diet.20.5)
h.p.diet.26 <- pct.change.end(hunt.diet.26)

h.p.K.base <- pct.change.end(hunt.K.base)
h.p.K.maxharvest <- pct.change.end(hunt.K.maxharvest)
h.p.K.contharvest <- pct.change.end(hunt.K.contharvest)
h.p.K.midharvest <- pct.change.end(hunt.K.midharvest)
h.p.K.nochange <- pct.change.end(hunt.K.nochange)
h.p.K.transition <- pct.change.end(hunt.K.transition)

h.p.roads.maxdecom <- pct.change.end(hunt.roads.maxdecom)
h.p.roads.plandecom <- pct.change.end(hunt.roads.plandecom)
h.p.roads.middecom <- pct.change.end(hunt.roads.middecom)
h.p.roads.nochange <- pct.change.end(hunt.roads.nochange)
h.p.roads.maxbuild <- pct.change.end(hunt.roads.maxbuild)

h.p.trap.0 <- pct.change.end(hunt.trap.0)
h.p.trap.20 <- pct.change.end(hunt.trap.20)
h.p.trap.30 <- pct.change.end(hunt.trap.30)
h.p.trap.none <- pct.change.end(hunt.trap.none)
h.p.trap.100 <- pct.change.end(hunt.trap.100)
h.p.gone <- pct.change.end(hunt.gone)

### create results table, deer sensitivity to parms ----

d.climate.sens.tab <- data.frame(Sen = c("clim.low", "clim.mid" ,"clim.max"), 
	Med.pct = c(d.p.clim.low$Med.pct, d.p.clim.mid$Med.pct, d.p.clim.max$Med.pct),
	Lcl.pct  = c(d.p.clim.low$Lcl.pct, d.p.clim.mid$Lcl.pct, d.p.clim.max$Lcl.pct),
	Ucl.pct  = c(d.p.clim.low$Ucl.pct, d.p.clim.mid$Ucl.pct, d.p.clim.max$Ucl.pct), 
	Med = c(deer.clim.low$med[31], deer.clim.mid$med[31], deer.clim.max$med[31]),
	Lcl = c(deer.clim.low$lcl[31], deer.clim.mid$lcl[31], deer.clim.max$lcl[31]),
	Ucl = c(deer.clim.low$ucl[31], deer.clim.mid$ucl[31], deer.clim.max$ucl[31]))

d.diet.sens.tab <- data.frame(Sen = c("diet.9.5", "diet.15" ,"diet.20.5", "diet.26"), 
	Med.pct = c(d.p.diet.9.5$Med.pct, d.p.diet.15$Med.pct, d.p.diet.20.5$Med.pct, d.p.diet.26$Med.pct),
	Lcl.pct  = c(d.p.diet.9.5$Lcl.pct, d.p.diet.15$Lcl.pct, d.p.diet.20.5$Lcl.pct, d.p.diet.26$Lcl.pct),
	Ucl.pct  = c(d.p.diet.9.5$Ucl.pct, d.p.diet.15$Ucl.pct, d.p.diet.20.5$Ucl.pct, d.p.diet.26$Ucl.pct), 
	Med = c(deer.diet.9.5$med[31], deer.diet.15$med[31], deer.diet.20.5$med[31], deer.diet.26$med[31]),
	Lcl = c(deer.diet.9.5$lcl[31], deer.diet.15$lcl[31], deer.diet.20.5$lcl[31], deer.diet.26$lcl[31]),
	Ucl = c(deer.diet.9.5$ucl[31], deer.diet.15$ucl[31], deer.diet.20.5$ucl[31], deer.diet.26$ucl[31]))

d.hunt.sens.tab <- data.frame(Sen = c("No.deer.hunt", "Deer.hunt"), 
	Med.pct = c(d.p.nohunt$Med.pct, d.p.hunt$Med.pct),
	Lcl.pct  = c(d.p.nohunt$Lcl.pct, d.p.hunt$Lcl.pct),
	Ucl.pct  = c(d.p.nohunt$Ucl.pct, d.p.hunt$Ucl.pct), 
	Med = c(deer.nohunt$med[31], deer.hunt$med[31]),
	Lcl = c(deer.nohunt$lcl[31], deer.hunt$lcl[31]),
	Ucl = c(deer.nohunt$ucl[31], deer.hunt$ucl[31]))

d.K.sens.tab <- data.frame(Sen = c("K.nochange", "K.base" ,"K.transition", "K.contharvest", "K.midharvest", "K.maxharvest"), 
	Med.pct = c(d.p.K.nochange$Med.pct, d.p.K.base$Med.pct, d.p.K.transition$Med.pct, d.p.K.contharvest$Med.pct, d.p.K.midharvest$Med.pct, d.p.K.maxharvest$Med.pct),
	Lcl.pct  = c(d.p.K.nochange$Lcl.pct, d.p.K.base$Lcl.pct, d.p.K.transition$Lcl.pct, d.p.K.contharvest$Lcl.pct, d.p.K.midharvest$Lcl.pct, d.p.K.maxharvest$Lcl.pct),
	Ucl.pct  = c(d.p.K.nochange$Ucl.pct, d.p.K.base$Ucl.pct, d.p.K.transition$Ucl.pct, d.p.K.contharvest$Ucl.pct, d.p.K.midharvest$Ucl.pct, d.p.K.maxharvest$Ucl.pct), 
	Med = c(deer.K.nochange$med[31], deer.K.base$med[31], deer.K.transition$med[31], deer.K.contharvest$med[31], deer.K.midharvest$med[31], deer.K.maxharvest$med[31]),
	Lcl = c(deer.K.nochange$lcl[31], deer.K.base$lcl[31], deer.K.transition$lcl[31], deer.K.contharvest$lcl[31], deer.K.midharvest$lcl[31], deer.K.maxharvest$lcl[31]),
	Ucl = c(deer.K.nochange$ucl[31], deer.K.base$ucl[31], deer.K.transition$ucl[31], deer.K.contharvest$ucl[31], deer.K.midharvest$ucl[31], deer.K.maxharvest$ucl[31]))

d.roads.sens.tab <- data.frame(Sen = c("roads.maxdecom", "roads.middecom" ,"roads.plandecom", "roads.nochange", "roads.maxbuild"), 
	Med.pct = c(d.p.roads.maxdecom$Med.pct, d.p.roads.middecom$Med.pct, d.p.roads.plandecom$Med.pct, d.p.roads.nochange$Med.pct, d.p.roads.maxbuild$Med.pct),
	Lcl.pct  = c(d.p.roads.maxdecom$Lcl.pct, d.p.roads.middecom$Lcl.pct, d.p.roads.plandecom$Lcl.pct, d.p.roads.nochange$Lcl.pct, d.p.roads.maxbuild$Lcl.pct),
	Ucl.pct  = c(d.p.roads.maxdecom$Ucl.pct, d.p.roads.middecom$Ucl.pct, d.p.roads.plandecom$Ucl.pct, d.p.roads.nochange$Ucl.pct, d.p.roads.maxbuild$Ucl.pct), 
	Med = c(deer.roads.maxdecom$med[31], deer.roads.middecom$med[31], deer.roads.plandecom$med[31], deer.roads.nochange$med[31], deer.roads.maxbuild$med[31]),
	Lcl = c(deer.roads.maxdecom$lcl[31], deer.roads.middecom$lcl[31], deer.roads.plandecom$lcl[31], deer.roads.nochange$lcl[31], deer.roads.maxbuild$lcl[31]),
	Ucl = c(deer.roads.maxdecom$ucl[31], deer.roads.middecom$ucl[31], deer.roads.plandecom$ucl[31], deer.roads.nochange$ucl[31], deer.roads.maxbuild$ucl[31]))

d.trap.sens.tab <- data.frame(Sen = c("TrapNone", "Trap0" ,"Trap20", "Trap30", "Trap100", "Gone"), 
	Med.pct = c(d.p.trap.none$Med.pct, d.p.trap.0$Med.pct, d.p.trap.20$Med.pct, d.p.trap.30$Med.pct, d.p.trap.100$Med.pct, d.p.gone$Med.pct),
	Lcl.pct  = c(d.p.trap.none$Lcl.pct, d.p.trap.0$Lcl.pct, d.p.trap.20$Lcl.pct, d.p.trap.30$Lcl.pct, d.p.trap.100$Lcl.pct, d.p.gone$Lcl.pct),
	Ucl.pct  = c(d.p.trap.none$Ucl.pct, d.p.trap.0$Ucl.pct, d.p.trap.20$Ucl.pct, d.p.trap.30$Ucl.pct, d.p.trap.100$Ucl.pct, d.p.gone$Ucl.pct), 
	Med = c(deer.trap.none$med[31], deer.trap.0$med[31], deer.trap.20$med[31], deer.trap.30$med[31], deer.trap.100$med[31], deer.gone$med[31]),
	Lcl = c(deer.trap.none$lcl[31], deer.trap.0$lcl[31], deer.trap.20$lcl[31], deer.trap.30$lcl[31], deer.trap.100$lcl[31], deer.gone$lcl[31]),
	Ucl = c(deer.trap.none$ucl[31], deer.trap.0$ucl[31], deer.trap.20$ucl[31], deer.trap.30$ucl[31], deer.trap.100$ucl[31], deer.gone$ucl[31]))





### create results table, wolf sensitivity to parms, fall ----

w.climate.sens.tab <- data.frame(Sen = c("clim.low", "clim.mid" ,"clim.max"), 
	Med.pct = c(w.p.clim.low$Med.pct, w.p.clim.mid$Med.pct, w.p.clim.max$Med.pct),
	Lcl.pct  = c(w.p.clim.low$Lcl.pct, w.p.clim.mid$Lcl.pct, w.p.clim.max$Lcl.pct),
	Ucl.pct  = c(w.p.clim.low$Ucl.pct, w.p.clim.mid$Ucl.pct, w.p.clim.max$Ucl.pct), 
	Med = c(wolf.clim.low$med[31], wolf.clim.mid$med[31], wolf.clim.max$med[31]),
	Lcl = c(wolf.clim.low$lcl[31], wolf.clim.mid$lcl[31], wolf.clim.max$lcl[31]),
	Ucl = c(wolf.clim.low$ucl[31], wolf.clim.mid$ucl[31], wolf.clim.max$ucl[31]))

w.diet.sens.tab <- data.frame(Sen = c("diet.9.5", "diet.15" ,"diet.20.5", "diet.26"), 
	Med.pct = c(w.p.diet.9.5$Med.pct, w.p.diet.15$Med.pct, w.p.diet.20.5$Med.pct, w.p.diet.26$Med.pct),
	Lcl.pct  = c(w.p.diet.9.5$Lcl.pct, w.p.diet.15$Lcl.pct, w.p.diet.20.5$Lcl.pct, w.p.diet.26$Lcl.pct),
	Ucl.pct  = c(w.p.diet.9.5$Ucl.pct, w.p.diet.15$Ucl.pct, w.p.diet.20.5$Ucl.pct, w.p.diet.26$Ucl.pct), 
	Med = c(wolf.diet.9.5$med[31], wolf.diet.15$med[31], wolf.diet.20.5$med[31], wolf.diet.26$med[31]),
	Lcl = c(wolf.diet.9.5$lcl[31], wolf.diet.15$lcl[31], wolf.diet.20.5$lcl[31], wolf.diet.26$lcl[31]),
	Ucl = c(wolf.diet.9.5$ucl[31], wolf.diet.15$ucl[31], wolf.diet.20.5$ucl[31], wolf.diet.26$ucl[31]))

w.hunt.sens.tab <- data.frame(Sen = c("No.wolf.hunt", "wolf.hunt"), 
	Med.pct = c(w.p.nohunt$Med.pct, w.p.hunt$Med.pct),
	Lcl.pct  = c(w.p.nohunt$Lcl.pct, w.p.hunt$Lcl.pct),
	Ucl.pct  = c(w.p.nohunt$Ucl.pct, w.p.hunt$Ucl.pct), 
	Med = c(wolf.nohunt$med[31], wolf.hunt$med[31]),
	Lcl = c(wolf.nohunt$lcl[31], wolf.hunt$lcl[31]),
	Ucl = c(wolf.nohunt$ucl[31], wolf.hunt$ucl[31]))

w.K.sens.tab <- data.frame(Sen = c("K.nochange", "K.base" ,"K.transition", "K.contharvest", "K.midharvest", "K.maxharvest"), 
	Med.pct = c(w.p.K.nochange$Med.pct, w.p.K.base$Med.pct, w.p.K.transition$Med.pct, w.p.K.contharvest$Med.pct, w.p.K.midharvest$Med.pct, w.p.K.maxharvest$Med.pct),
	Lcl.pct  = c(w.p.K.nochange$Lcl.pct, w.p.K.base$Lcl.pct, w.p.K.transition$Lcl.pct, w.p.K.contharvest$Lcl.pct, w.p.K.midharvest$Lcl.pct, w.p.K.maxharvest$Lcl.pct),
	Ucl.pct  = c(w.p.K.nochange$Ucl.pct, w.p.K.base$Ucl.pct, w.p.K.transition$Ucl.pct, w.p.K.contharvest$Ucl.pct, w.p.K.midharvest$Ucl.pct, w.p.K.maxharvest$Ucl.pct), 
	Med = c(wolf.K.nochange$med[31], wolf.K.base$med[31], wolf.K.transition$med[31], wolf.K.contharvest$med[31], wolf.K.midharvest$med[31], wolf.K.maxharvest$med[31]),
	Lcl = c(wolf.K.nochange$lcl[31], wolf.K.base$lcl[31], wolf.K.transition$lcl[31], wolf.K.contharvest$lcl[31], wolf.K.midharvest$lcl[31], wolf.K.maxharvest$lcl[31]),
	Ucl = c(wolf.K.nochange$ucl[31], wolf.K.base$ucl[31], wolf.K.transition$ucl[31], wolf.K.contharvest$ucl[31], wolf.K.midharvest$ucl[31], wolf.K.maxharvest$ucl[31]))

w.roads.sens.tab <- data.frame(Sen = c("roads.maxdecom", "roads.middecom" ,"roads.plandecom", "roads.nochange", "roads.maxbuild"), 
	Med.pct = c(w.p.roads.maxdecom$Med.pct, w.p.roads.middecom$Med.pct, w.p.roads.plandecom$Med.pct, w.p.roads.nochange$Med.pct, w.p.roads.maxbuild$Med.pct),
	Lcl.pct  = c(w.p.roads.maxdecom$Lcl.pct, w.p.roads.middecom$Lcl.pct, w.p.roads.plandecom$Lcl.pct, w.p.roads.nochange$Lcl.pct, w.p.roads.maxbuild$Lcl.pct),
	Ucl.pct  = c(w.p.roads.maxdecom$Ucl.pct, w.p.roads.middecom$Ucl.pct, w.p.roads.plandecom$Ucl.pct, w.p.roads.nochange$Ucl.pct, w.p.roads.maxbuild$Ucl.pct), 
	Med = c(wolf.roads.maxdecom$med[31], wolf.roads.middecom$med[31], wolf.roads.plandecom$med[31], wolf.roads.nochange$med[31], wolf.roads.maxbuild$med[31]),
	Lcl = c(wolf.roads.maxdecom$lcl[31], wolf.roads.middecom$lcl[31], wolf.roads.plandecom$lcl[31], wolf.roads.nochange$lcl[31], wolf.roads.maxbuild$lcl[31]),
	Ucl = c(wolf.roads.maxdecom$ucl[31], wolf.roads.middecom$ucl[31], wolf.roads.plandecom$ucl[31], wolf.roads.nochange$ucl[31], wolf.roads.maxbuild$ucl[31]))


w.trap.sens.tab <- data.frame(Sen = c("TrapNone", "Trap0" ,"Trap20", "Trap30", "Trap100"), 
	Med.pct = c(w.p.trap.none$Med.pct, w.p.trap.0$Med.pct, w.p.trap.20$Med.pct, w.p.trap.30$Med.pct, w.p.trap.100$Med.pct),
	Lcl.pct  = c(w.p.trap.none$Lcl.pct, w.p.trap.0$Lcl.pct, w.p.trap.20$Lcl.pct, w.p.trap.30$Lcl.pct, w.p.trap.100$Lcl.pct),
	Ucl.pct  = c(w.p.trap.none$Ucl.pct, w.p.trap.0$Ucl.pct, w.p.trap.20$Ucl.pct, w.p.trap.30$Ucl.pct, w.p.trap.100$Ucl.pct), 
	Med = c(wolf.trap.none$med[31], wolf.trap.0$med[31], wolf.trap.20$med[31], wolf.trap.30$med[31], wolf.trap.100$med[31]),
	Lcl = c(wolf.trap.none$lcl[31], wolf.trap.0$lcl[31], wolf.trap.20$lcl[31], wolf.trap.30$lcl[31], wolf.trap.100$lcl[31]),
	Ucl = c(wolf.trap.none$ucl[31], wolf.trap.0$ucl[31], wolf.trap.20$ucl[31], wolf.trap.30$ucl[31], wolf.trap.100$ucl[31]))




### create results table, wolf sensitivity to parms, spring----

w.climate.sens.tab.s <- data.frame(Sen = c("clim.low", "clim.mid" ,"clim.max"), 
	Med.pct = c(w.p.clim.low.s$Med.pct, w.p.clim.mid.s$Med.pct, w.p.clim.max.s$Med.pct),
	Lcl.pct  = c(w.p.clim.low.s$Lcl.pct, w.p.clim.mid.s$Lcl.pct, w.p.clim.max.s$Lcl.pct),
	Ucl.pct  = c(w.p.clim.low.s$Ucl.pct, w.p.clim.mid.s$Ucl.pct, w.p.clim.max.s$Ucl.pct), 
	Med = c(wolf.clim.low.s$med[31], wolf.clim.mid.s$med[31], wolf.clim.max.s$med[31]),
	Lcl = c(wolf.clim.low.s$lcl[31], wolf.clim.mid.s$lcl[31], wolf.clim.max.s$lcl[31]),
	Ucl = c(wolf.clim.low.s$ucl[31], wolf.clim.mid.s$ucl[31], wolf.clim.max.s$ucl[31]))

w.diet.sens.tab.s <- data.frame(Sen = c("diet.9.5", "diet.15" ,"diet.20.5", "diet.26"), 
	Med.pct = c(w.p.diet.9.5.s$Med.pct, w.p.diet.15.s$Med.pct, w.p.diet.20.5.s$Med.pct, w.p.diet.26.s$Med.pct),
	Lcl.pct  = c(w.p.diet.9.5.s$Lcl.pct, w.p.diet.15.s$Lcl.pct, w.p.diet.20.5.s$Lcl.pct, w.p.diet.26.s$Lcl.pct),
	Ucl.pct  = c(w.p.diet.9.5.s$Ucl.pct, w.p.diet.15.s$Ucl.pct, w.p.diet.20.5.s$Ucl.pct, w.p.diet.26.s$Ucl.pct), 
	Med = c(wolf.diet.9.5.s$med[31], wolf.diet.15.s$med[31], wolf.diet.20.5.s$med[31], wolf.diet.26.s$med[31]),
	Lcl = c(wolf.diet.9.5.s$lcl[31], wolf.diet.15.s$lcl[31], wolf.diet.20.5.s$lcl[31], wolf.diet.26.s$lcl[31]),
	Ucl = c(wolf.diet.9.5.s$ucl[31], wolf.diet.15.s$ucl[31], wolf.diet.20.5.s$ucl[31], wolf.diet.26.s$ucl[31]))

w.hunt.sens.tab.s <- data.frame(Sen = c("No.wolf.hunt", "wolf.hunt"), 
	Med.pct = c(w.p.nohunt.s$Med.pct, w.p.hunt.s$Med.pct),
	Lcl.pct  = c(w.p.nohunt.s$Lcl.pct, w.p.hunt.s$Lcl.pct),
	Ucl.pct  = c(w.p.nohunt.s$Ucl.pct, w.p.hunt.s$Ucl.pct), 
	Med = c(wolf.nohunt.s$med[31], wolf.hunt.s$med[31]),
	Lcl = c(wolf.nohunt.s$lcl[31], wolf.hunt.s$lcl[31]),
	Ucl = c(wolf.nohunt.s$ucl[31], wolf.hunt.s$ucl[31]))

w.K.sens.tab.s <- data.frame(Sen = c("K.nochange", "K.base" ,"K.transition", "K.contharvest", "K.midharvest", "K.maxharvest"), 
	Med.pct = c(w.p.K.nochange.s$Med.pct, w.p.K.base.s$Med.pct, w.p.K.transition.s$Med.pct, w.p.K.contharvest.s$Med.pct, w.p.K.midharvest.s$Med.pct, w.p.K.maxharvest.s$Med.pct),
	Lcl.pct  = c(w.p.K.nochange.s$Lcl.pct, w.p.K.base.s$Lcl.pct, w.p.K.transition.s$Lcl.pct, w.p.K.contharvest.s$Lcl.pct, w.p.K.midharvest.s$Lcl.pct, w.p.K.maxharvest.s$Lcl.pct),
	Ucl.pct  = c(w.p.K.nochange.s$Ucl.pct, w.p.K.base.s$Ucl.pct, w.p.K.transition.s$Ucl.pct, w.p.K.contharvest.s$Ucl.pct, w.p.K.midharvest.s$Ucl.pct, w.p.K.maxharvest.s$Ucl.pct), 
	Med = c(wolf.K.nochange.s$med[31], wolf.K.base.s$med[31], wolf.K.transition.s$med[31], wolf.K.contharvest.s$med[31], wolf.K.midharvest.s$med[31], wolf.K.maxharvest.s$med[31]),
	Lcl = c(wolf.K.nochange.s$lcl[31], wolf.K.base.s$lcl[31], wolf.K.transition.s$lcl[31], wolf.K.contharvest.s$lcl[31], wolf.K.midharvest.s$lcl[31], wolf.K.maxharvest.s$lcl[31]),
	Ucl = c(wolf.K.nochange.s$ucl[31], wolf.K.base.s$ucl[31], wolf.K.transition.s$ucl[31], wolf.K.contharvest.s$ucl[31], wolf.K.midharvest.s$ucl[31], wolf.K.maxharvest.s$ucl[31]))
	
w.roads.sens.tab.s <- data.frame(Sen = c("roads.maxdecom", "roads.middecom" ,"roads.plandecom", "roads.nochange", "roads.maxbuild"), 
	Med.pct = c(w.p.roads.maxdecom.s$Med.pct, w.p.roads.middecom.s$Med.pct, w.p.roads.plandecom.s$Med.pct, w.p.roads.nochange.s$Med.pct, w.p.roads.maxbuild.s$Med.pct),
	Lcl.pct  = c(w.p.roads.maxdecom.s$Lcl.pct, w.p.roads.middecom.s$Lcl.pct, w.p.roads.plandecom.s$Lcl.pct, w.p.roads.nochange.s$Lcl.pct, w.p.roads.maxbuild.s$Lcl.pct),
	Ucl.pct  = c(w.p.roads.maxdecom.s$Ucl.pct, w.p.roads.middecom.s$Ucl.pct, w.p.roads.plandecom.s$Ucl.pct, w.p.roads.nochange.s$Ucl.pct, w.p.roads.maxbuild.s$Ucl.pct), 
	Med = c(wolf.roads.maxdecom.s$med[31], wolf.roads.middecom.s$med[31], wolf.roads.plandecom.s$med[31], wolf.roads.nochange.s$med[31], wolf.roads.maxbuild.s$med[31]),
	Lcl = c(wolf.roads.maxdecom.s$lcl[31], wolf.roads.middecom.s$lcl[31], wolf.roads.plandecom.s$lcl[31], wolf.roads.nochange.s$lcl[31], wolf.roads.maxbuild.s$lcl[31]),
	Ucl = c(wolf.roads.maxdecom.s$ucl[31], wolf.roads.middecom.s$ucl[31], wolf.roads.plandecom.s$ucl[31], wolf.roads.nochange.s$ucl[31], wolf.roads.maxbuild.s$ucl[31]))

w.trap.sens.tab.s <- data.frame(Sen = c("TrapNone", "Trap0" ,"Trap20", "Trap30", "Trap100"), 
	Med.pct = c(w.p.trap.none.s$Med.pct, w.p.trap.0.s$Med.pct, w.p.trap.20.s$Med.pct, w.p.trap.30.s$Med.pct, w.p.trap.100.s$Med.pct),
	Lcl.pct  = c(w.p.trap.none.s$Lcl.pct, w.p.trap.0.s$Lcl.pct, w.p.trap.20.s$Lcl.pct, w.p.trap.30.s$Lcl.pct, w.p.trap.100.s$Lcl.pct),
	Ucl.pct  = c(w.p.trap.none.s$Ucl.pct, w.p.trap.0.s$Ucl.pct, w.p.trap.20.s$Ucl.pct, w.p.trap.30.s$Ucl.pct, w.p.trap.100.s$Ucl.pct), 
	Med = c(wolf.trap.none.s$med[31], wolf.trap.0.s$med[31], wolf.trap.20.s$med[31], wolf.trap.30.s$med[31], wolf.trap.100.s$med[31]),
	Lcl = c(wolf.trap.none.s$lcl[31], wolf.trap.0.s$lcl[31], wolf.trap.20.s$lcl[31], wolf.trap.30.s$lcl[31], wolf.trap.100.s$lcl[31]),
	Ucl = c(wolf.trap.none.s$ucl[31], wolf.trap.0.s$ucl[31], wolf.trap.20.s$ucl[31], wolf.trap.30.s$ucl[31], wolf.trap.30.s$ucl[31]))



### create results table, hunt sensitivity to parms ----

h.climate.sens.tab <- data.frame(Sen = c("clim.low", "clim.mid" ,"clim.max"), 
	Med.pct = c(h.p.clim.low$Med.pct, h.p.clim.mid$Med.pct, h.p.clim.max$Med.pct),
	Lcl.pct  = c(h.p.clim.low$Lcl.pct, h.p.clim.mid$Lcl.pct, h.p.clim.max$Lcl.pct),
	Ucl.pct  = c(h.p.clim.low$Ucl.pct, h.p.clim.mid$Ucl.pct, h.p.clim.max$Ucl.pct), 
	Med = c(hunt.clim.low$med[31], hunt.clim.mid$med[31], hunt.clim.max$med[31]),
	Lcl = c(hunt.clim.low$lcl[31], hunt.clim.mid$lcl[31], hunt.clim.max$lcl[31]),
	Ucl = c(hunt.clim.low$ucl[31], hunt.clim.mid$ucl[31], hunt.clim.max$ucl[31]))

h.diet.sens.tab <- data.frame(Sen = c("diet.9.5", "diet.15" ,"diet.20.5", "diet.26"), 
	Med.pct = c(h.p.diet.9.5$Med.pct, h.p.diet.15$Med.pct, h.p.diet.20.5$Med.pct, h.p.diet.26$Med.pct),
	Lcl.pct  = c(h.p.diet.9.5$Lcl.pct, h.p.diet.15$Lcl.pct, h.p.diet.20.5$Lcl.pct, h.p.diet.26$Lcl.pct),
	Ucl.pct  = c(h.p.diet.9.5$Ucl.pct, h.p.diet.15$Ucl.pct, h.p.diet.20.5$Ucl.pct, h.p.diet.26$Ucl.pct), 
	Med = c(hunt.diet.9.5$med[31], hunt.diet.15$med[31], hunt.diet.20.5$med[31], hunt.diet.26$med[31]),
	Lcl = c(hunt.diet.9.5$lcl[31], hunt.diet.15$lcl[31], hunt.diet.20.5$lcl[31], hunt.diet.26$lcl[31]),
	Ucl = c(hunt.diet.9.5$ucl[31], hunt.diet.15$ucl[31], hunt.diet.20.5$ucl[31], hunt.diet.26$ucl[31]))


h.K.sens.tab <- data.frame(Sen = c("K.nochange", "K.base" ,"K.transition", "K.contharvest", "K.midharvest", "K.maxharvest"), 
	Med.pct = c(h.p.K.nochange$Med.pct, h.p.K.base$Med.pct, h.p.K.transition$Med.pct, h.p.K.contharvest$Med.pct, h.p.K.midharvest$Med.pct, h.p.K.maxharvest$Med.pct),
	Lcl.pct  = c(h.p.K.nochange$Lcl.pct, h.p.K.base$Lcl.pct, h.p.K.transition$Lcl.pct, h.p.K.contharvest$Lcl.pct, h.p.K.midharvest$Lcl.pct, h.p.K.maxharvest$Lcl.pct),
	Ucl.pct  = c(h.p.K.nochange$Ucl.pct, h.p.K.base$Ucl.pct, h.p.K.transition$Ucl.pct, h.p.K.contharvest$Ucl.pct, h.p.K.midharvest$Ucl.pct, h.p.K.maxharvest$Ucl.pct), 
	Med = c(hunt.K.nochange$med[31], hunt.K.base$med[31], hunt.K.transition$med[31], hunt.K.contharvest$med[31], hunt.K.midharvest$med[31], hunt.K.maxharvest$med[31]),
	Lcl = c(hunt.K.nochange$lcl[31], hunt.K.base$lcl[31], hunt.K.transition$lcl[31], hunt.K.contharvest$lcl[31], hunt.K.midharvest$lcl[31], hunt.K.maxharvest$lcl[31]),
	Ucl = c(hunt.K.nochange$ucl[31], hunt.K.base$ucl[31], hunt.K.transition$ucl[31], hunt.K.contharvest$ucl[31], hunt.K.midharvest$ucl[31], hunt.K.maxharvest$ucl[31]))

h.roads.sens.tab <- data.frame(Sen = c("roads.maxdecom", "roads.middecom" ,"roads.plandecom", "roads.nochange", "roads.maxbuild"), 
	Med.pct = c(h.p.roads.maxdecom$Med.pct, h.p.roads.middecom$Med.pct, h.p.roads.plandecom$Med.pct, h.p.roads.nochange$Med.pct, h.p.roads.maxbuild$Med.pct),
	Lcl.pct  = c(h.p.roads.maxdecom$Lcl.pct, h.p.roads.middecom$Lcl.pct, h.p.roads.plandecom$Lcl.pct, h.p.roads.nochange$Lcl.pct, h.p.roads.maxbuild$Lcl.pct),
	Ucl.pct  = c(h.p.roads.maxdecom$Ucl.pct, h.p.roads.middecom$Ucl.pct, h.p.roads.plandecom$Ucl.pct, h.p.roads.nochange$Ucl.pct, h.p.roads.maxbuild$Ucl.pct), 
	Med = c(hunt.roads.maxdecom$med[31], hunt.roads.middecom$med[31], hunt.roads.plandecom$med[31], hunt.roads.nochange$med[31], hunt.roads.maxbuild$med[31]),
	Lcl = c(hunt.roads.maxdecom$lcl[31], hunt.roads.middecom$lcl[31], hunt.roads.plandecom$lcl[31], hunt.roads.nochange$lcl[31], hunt.roads.maxbuild$lcl[31]),
	Ucl = c(hunt.roads.maxdecom$ucl[31], hunt.roads.middecom$ucl[31], hunt.roads.plandecom$ucl[31], hunt.roads.nochange$ucl[31], hunt.roads.maxbuild$ucl[31]))


h.trap.sens.tab <- data.frame(Sen = c("TrapNone", "Trap0" ,"Trap20", "Trap30", "Trap100", "Gone"), 
	Med.pct = c(h.p.trap.none$Med.pct, h.p.trap.0$Med.pct, h.p.trap.20$Med.pct, h.p.trap.30$Med.pct, h.p.trap.100$Med.pct, h.p.gone$Med.pct),
	Lcl.pct  = c(h.p.trap.none$Lcl.pct, h.p.trap.0$Lcl.pct, h.p.trap.20$Lcl.pct, h.p.trap.30$Lcl.pct, h.p.trap.100$Lcl.pct, h.p.gone$Lcl.pct),
	Ucl.pct  = c(h.p.trap.none$Ucl.pct, h.p.trap.0$Ucl.pct, h.p.trap.20$Ucl.pct, h.p.trap.30$Ucl.pct, h.p.trap.100$Ucl.pct, h.p.gone$Ucl.pct), 
	Med = c(hunt.trap.none$med[31], hunt.trap.0$med[31], hunt.trap.20$med[31], hunt.trap.30$med[31], hunt.trap.100$med[31], hunt.gone$med[31]),
	Lcl = c(hunt.trap.none$lcl[31], hunt.trap.0$lcl[31], hunt.trap.20$lcl[31], hunt.trap.30$lcl[31], hunt.trap.100$lcl[31], hunt.gone$lcl[31]),
	Ucl = c(hunt.trap.none$ucl[31], hunt.trap.0$ucl[31], hunt.trap.20$ucl[31], hunt.trap.30$ucl[31], hunt.trap.100$ucl[31], hunt.gone$ucl[31]))



### create results table, browse sensitivity to parms #####

b.climate.sens.tab <- data.frame(Sen = c("clim.low", "clim.mid" ,"clim.max"),  
	Med = c(browse.clim.low$med[31], browse.clim.mid$med[31], browse.clim.max$med[31]),
	Lcl = c(browse.clim.low$lcl[31], browse.clim.mid$lcl[31], browse.clim.max$lcl[31]),
	Ucl = c(browse.clim.low$ucl[31], browse.clim.mid$ucl[31], browse.clim.max$ucl[31]))

b.diet.sens.tab <- data.frame(Sen = c("diet.9.5", "diet.15" ,"diet.20.5", "diet.26"), 
	Med = c(browse.diet.9.5$med[31], browse.diet.15$med[31], browse.diet.20.5$med[31], browse.diet.26$med[31]),
	Lcl = c(browse.diet.9.5$lcl[31], browse.diet.15$lcl[31], browse.diet.20.5$lcl[31], browse.diet.26$lcl[31]),
	Ucl = c(browse.diet.9.5$ucl[31], browse.diet.15$ucl[31], browse.diet.20.5$ucl[31], browse.diet.26$ucl[31]))

b.hunt.sens.tab <- data.frame(Sen = c("No.hunt", "hunt"), 
	Med = c(browse.nohunt$med[31], browse.hunt$med[31]),
	Lcl = c(browse.nohunt$lcl[31], browse.hunt$lcl[31]),
	Ucl = c(browse.nohunt$ucl[31], browse.hunt$ucl[31]))

b.K.sens.tab <- data.frame(Sen = c("K.nochange", "K.base" ,"K.transition", "K.contharvest", "K.midharvest", "K.maxharvest"), 
	Med = c(browse.K.nochange$med[31], browse.K.base$med[31], browse.K.transition$med[31], browse.K.contharvest$med[31], browse.K.midharvest$med[31], browse.K.maxharvest$med[31]),
	Lcl = c(browse.K.nochange$lcl[31], browse.K.base$lcl[31], browse.K.transition$lcl[31], browse.K.contharvest$lcl[31], browse.K.midharvest$lcl[31], browse.K.maxharvest$lcl[31]),
	Ucl = c(browse.K.nochange$ucl[31], browse.K.base$ucl[31], browse.K.transition$ucl[31], browse.K.contharvest$ucl[31], browse.K.midharvest$ucl[31], browse.K.maxharvest$ucl[31]))

b.roads.sens.tab <- data.frame(Sen = c("roads.maxdecom", "roads.middecom" ,"roads.plandecom", "roads.nochange", "roads.maxbuild"), 
		Med = c(browse.roads.maxdecom$med[31], browse.roads.middecom$med[31], browse.roads.plandecom$med[31], browse.roads.nochange$med[31], browse.roads.maxbuild$med[31]),
	Lcl = c(browse.roads.maxdecom$lcl[31], browse.roads.middecom$lcl[31], browse.roads.plandecom$lcl[31], browse.roads.nochange$lcl[31], browse.roads.maxbuild$lcl[31]),
	Ucl = c(browse.roads.maxdecom$ucl[31], browse.roads.middecom$ucl[31], browse.roads.plandecom$ucl[31], browse.roads.nochange$ucl[31], browse.roads.maxbuild$ucl[31]))


b.trap.sens.tab <- data.frame(Sen = c("TrapNone", "Trap0" ,"Trap20", "Trap30", "Trap100", "Gone"), 
	Med = c(browse.trap.none$med[31], browse.trap.0$med[31], browse.trap.20$med[31], browse.trap.30$med[31], browse.trap.100$med[31], browse.gone$med[31]),
	Lcl = c(browse.trap.none$lcl[31], browse.trap.0$lcl[31], browse.trap.20$lcl[31], browse.trap.30$lcl[31], browse.trap.100$lcl[31], browse.gone$lcl[31]),
	Ucl = c(browse.trap.none$ucl[31], browse.trap.0$ucl[31], browse.trap.20$ucl[31], browse.trap.30$ucl[31], browse.trap.100$ucl[31], browse.gone$ucl[31]))




### calculate extinction probabilities for a range of quasi-extinction thresholds----

q.range = seq(0, 60, 1)

for(i in 1:length(q.range)){
	q.th = q.range[i]
	
	e.clim.low = extinction.prob(wolf.clim.low.monte, q.thresh= q.th)
	e.clim.mid = extinction.prob(wolf.clim.mid.monte, q.thresh= q.th)
	e.clim.max = extinction.prob(wolf.clim.max.monte, q.thresh= q.th)
	
	e.diet.9.5 = extinction.prob(wolf.diet.9.5.monte, q.thresh= q.th)
	e.diet.15  = extinction.prob(wolf.diet.15.monte, q.thresh= q.th)
	e.diet.20.5 =extinction.prob(wolf.diet.20.5.monte, q.thresh= q.th)
	e.diet.26  = extinction.prob(wolf.diet.26.monte, q.thresh= q.th)
	
	e.nohunt = extinction.prob(wolf.nohunt.monte, q.thresh= q.th)
	e.hunt   = extinction.prob(wolf.hunt.monte, q.thresh= q.th)
	
	e.K.nochange = extinction.prob(wolf.K.nochange.monte, q.thresh= q.th)
	e.K.base    = extinction.prob(wolf.K.base.monte, q.thresh= q.th)
	e.K.transition = extinction.prob(wolf.K.transition.monte, q.thresh= q.th)
	e.K.contharvest = extinction.prob(wolf.K.contharvest.monte, q.thresh= q.th)
	e.K.midharvest = extinction.prob(wolf.K.midharvest.monte, q.thresh= q.th)
	e.K.maxharvest = extinction.prob(wolf.K.maxharvest.monte, q.thresh= q.th)
	
	e.roads.maxdecom = extinction.prob(wolf.roads.maxdecom.monte, q.thresh= q.th)
	e.roads.middecom = extinction.prob(wolf.roads.middecom.monte, q.thresh= q.th)
	e.roads.plandecom = extinction.prob(wolf.roads.plandecom.monte, q.thresh= q.th)
	e.roads.nochange = extinction.prob(wolf.roads.nochange.monte, q.thresh= q.th)
	e.roads.maxbuild = extinction.prob(wolf.roads.maxbuild.monte, q.thresh= q.th)
	
	e.trap.none = extinction.prob(wolf.trap.none.monte, q.thresh= q.th)
	e.trap.0 = extinction.prob(wolf.trap.0.monte, q.thresh= q.th)
	e.trap.20 = extinction.prob(wolf.trap.20.monte, q.thresh= q.th)
	e.trap.30 = extinction.prob(wolf.trap.30.monte, q.thresh= q.th)
	e.trap.100 = extinction.prob(wolf.trap.100.monte, q.thresh= q.th)
	
	if(i==1){		
		ext.clim.low = e.clim.low;
		ext.clim.mid = e.clim.mid;
		ext.clim.max = e.clim.max;
	
		ext.diet.9.5 = e.diet.9.5;
		ext.diet.15 = e.diet.15;
		ext.diet.20.5 = e.diet.20.5;
		ext.diet.26  = e.diet.26;
	
		ext.nohunt = e.nohunt;
		ext.hunt   = e.hunt;
	
		ext.K.nochange = e.K.nochange;
		ext.K.base    = e.K.base;
		ext.K.transition = e.K.transition;
		ext.K.contharvest = e.K.contharvest;
		ext.K.midharvest = e.K.midharvest;
		ext.K.maxharvest = e.K.maxharvest;
	
		ext.roads.maxdecom = e.roads.maxdecom;
		ext.roads.middecom = e.roads.middecom;
		ext.roads.plandecom = e.roads.plandecom;
		ext.roads.nochange = e.roads.nochange;
		ext.roads.maxbuild = e.roads.maxbuild;
	
		ext.trap.none = e.trap.none;
		ext.trap.0 = e.trap.0;
		ext.trap.20 = e.trap.20;
		ext.trap.30 = e.trap.30;
		ext.trap.100 = e.trap.100

			}else{
				ext.clim.low = rbind(ext.clim.low, e.clim.low);
				ext.clim.mid = rbind(ext.clim.mid, e.clim.mid);
				ext.clim.max = rbind(ext.clim.max, e.clim.max);
	
				ext.diet.9.5 = rbind(ext.diet.9.5, e.diet.9.5);
				ext.diet.15 = rbind(ext.diet.15, e.diet.15);
				ext.diet.20.5 = rbind(ext.diet.20.5, e.diet.20.5);
				exn.diet.26  = rbind(ext.diet.26, e.diet.26);
	
				ext.nohunt = rbind(ext.nohunt, e.nohunt);
				ext.hunt   = rbind(ext.hunt, e.hunt);
	
				ext.K.nochange = rbind(ext.K.nochange, e.K.nochange);
				ext.K.base    = rbind(ext.K.base, e.K.base);
				ext.K.transition = rbind(ext.K.transition, e.K.transition);
				ext.K.contharvest = rbind(ext.K.contharvest, e.K.contharvest);
				ext.K.midharvest = rbind(ext.K.midharvest, e.K.midharvest);
				ext.K.maxharvest = rbind(ext.K.maxharvest, e.K.maxharvest);
	
				ext.roads.maxdecom = rbind(ext.roads.maxdecom, e.roads.maxdecom);
				ext.roads.middecom = rbind(ext.roads.middecom, e.roads.middecom);
				ext.roads.plandecom = rbind(ext.roads.plandecom, e.roads.plandecom);
				ext.roads.nochange = rbind(ext.roads.nochange, e.roads.nochange);
				ext.roads.maxbuild = rbind(ext.roads.maxbuild, e.roads.maxbuild);
	
				ext.trap.none = rbind(ext.trap.none, e.trap.none);
				ext.trap.0 = rbind(ext.trap.0, e.trap.0);
				ext.trap.20 = rbind(ext.trap.20, e.trap.20);
				ext.trap.30 = rbind(ext.trap.30, e.trap.30);
				ext.trap.100 = rbind(ext.trap.100, e.trap.100)
				}
}





############### Write the results to file #############################

## scenarios----

write.table(w.scen.tab, paste(Dir.out, "wolf.scen.tab.fall.csv", sep='/'), sep=",", row.names=F)
write.table(w.scen.tab.s, paste(Dir.out, "wolf.scen.tab.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.scen.tab, paste(Dir.out, "deer.scen.tab.fall.csv", sep='/'), sep=",", row.names=F)
write.table(p.scen.tab, paste(Dir.out, "packs.scen.tab.csv", sep='/'), sep=",", row.names=F)
write.table(h.scen.tab, paste(Dir.out, "huntdeer.scen.tab.csv", sep='/'), sep=",", row.names=F)
write.table(b.scen.tab, paste(Dir.out, "browse.scen.tab.csv", sep='/'), sep=",", row.names=F)

write.table(ext.base, paste(Dir.out, "wolf.extinct.base.csv", sep='/'), sep=",", row.names=F)
write.table(ext.A, paste(Dir.out, "wolf.extinct.A.csv", sep='/'), sep=",", row.names=F)
write.table(ext.B, paste(Dir.out, "wolf.extinct.B.csv", sep='/'), sep=",", row.names=F)
write.table(ext.C, paste(Dir.out, "wolf.extinct.C.csv", sep='/'), sep=",", row.names=F)
write.table(ext.D, paste(Dir.out, "wolf.extinct.D.csv", sep='/'), sep=",", row.names=F)
write.table(ext.E, paste(Dir.out, "wolf.extinct.E.csv", sep='/'), sep=",", row.names=F)


## sensitivities----

write.table(w.climate.sens.tab, paste(Dir.out, "wolf.climate.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.climate.sens.tab.s, paste(Dir.out, "wolf.climate.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.climate.sens.tab, paste(Dir.out, "deer.climate.sens.csv", sep='/'), sep=",", row.names=F)
write.table(h.climate.sens.tab, paste(Dir.out, "hunt.climate.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.climate.sens.tab, paste(Dir.out, "browse.climate.sens.csv", sep='/'), sep=",", row.names=F)

write.table(w.diet.sens.tab, paste(Dir.out, "wolf.diet.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.diet.sens.tab.s, paste(Dir.out, "wolf.diet.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.diet.sens.tab, paste(Dir.out, "deer.diet.sens.csv", sep='/'), sep=",", row.names=F)
write.table(h.diet.sens.tab, paste(Dir.out, "hunt.diet.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.diet.sens.tab, paste(Dir.out, "browse.diet.sens.csv", sep='/'), sep=",", row.names=F)

write.table(w.hunt.sens.tab, paste(Dir.out, "wolf.hunt.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.hunt.sens.tab.s, paste(Dir.out, "wolf.hunt.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.hunt.sens.tab, paste(Dir.out, "deer.hunt.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.hunt.sens.tab, paste(Dir.out, "browse.hunt.sens.csv", sep='/'), sep=",", row.names=F)

write.table(w.K.sens.tab, paste(Dir.out, "wolf.K.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.K.sens.tab.s, paste(Dir.out, "wolf.K.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.K.sens.tab, paste(Dir.out, "deer.K.sens.csv", sep='/'), sep=",", row.names=F)
write.table(h.K.sens.tab, paste(Dir.out, "hunt.K.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.K.sens.tab, paste(Dir.out, "browse.K.sens.csv", sep='/'), sep=",", row.names=F)

write.table(w.roads.sens.tab, paste(Dir.out, "wolf.roads.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.roads.sens.tab.s, paste(Dir.out, "wolf.roads.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.roads.sens.tab, paste(Dir.out, "deer.roads.sens.csv", sep='/'), sep=",", row.names=F)
write.table(h.roads.sens.tab, paste(Dir.out, "hunt.roads.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.roads.sens.tab, paste(Dir.out, "browse.roads.sens.csv", sep='/'), sep=",", row.names=F)

write.table(w.trap.sens.tab, paste(Dir.out, "wolf.trap.sens.csv", sep='/'), sep=",", row.names=F)
write.table(w.trap.sens.tab.s, paste(Dir.out, "wolf.trap.sens.spring.csv", sep='/'), sep=",", row.names=F)
write.table(d.trap.sens.tab, paste(Dir.out, "deer.trap.sens.csv", sep='/'), sep=",", row.names=F)
write.table(h.trap.sens.tab, paste(Dir.out, "hunt.trap.sens.csv", sep='/'), sep=",", row.names=F)
write.table(b.trap.sens.tab, paste(Dir.out, "browse.trap.sens.csv", sep='/'), sep=",", row.names=F)

write.table(ext.clim.low, paste(Dir.out, "wolf.extinct.clim.low.csv", sep='/'), sep=",", row.names=F)
write.table(ext.clim.mid, paste(Dir.out, "wolf.extinct.clim.mid.csv", sep='/'), sep=",", row.names=F)
write.table(ext.clim.max, paste(Dir.out, "wolf.extinct.clim.max.csv", sep='/'), sep=",", row.names=F)

write.table(ext.diet.9.5, paste(Dir.out, "wolf.extinct.diet.9.5.csv", sep='/'), sep=",", row.names=F)
write.table(ext.diet.15, paste(Dir.out, "wolf.extinct.diet.15.csv", sep='/'), sep=",", row.names=F)
write.table(ext.diet.20.5, paste(Dir.out, "wolf.extinct.diet.20.5.csv", sep='/'), sep=",", row.names=F)
write.table(ext.diet.26, paste(Dir.out, "wolf.extinct.diet.26.csv", sep='/'), sep=",", row.names=F)

write.table(ext.hunt, paste(Dir.out, "wolf.extinct.hunt.csv", sep='/'), sep=",", row.names=F)
write.table(ext.nohunt, paste(Dir.out, "wolf.extinct.nohunt.csv", sep='/'), sep=",", row.names=F)

write.table(ext.K.nochange, paste(Dir.out, "wolf.extinct.K.nochange.csv", sep='/'), sep=",", row.names=F)
write.table(ext.K.base, paste(Dir.out, "wolf.extinct.K.base.csv", sep='/'), sep=",", row.names=F)
write.table(ext.K.transition, paste(Dir.out, "wolf.extinct.K.transition.csv", sep='/'), sep=",", row.names=F)
write.table(ext.K.midharvest, paste(Dir.out, "wolf.extinct.K.midharvest.csv", sep='/'), sep=",", row.names=F)
write.table(ext.K.maxharvest, paste(Dir.out, "wolf.extinct.K.maxharvest.csv", sep='/'), sep=",", row.names=F)

write.table(ext.roads.maxdecom, paste(Dir.out, "wolf.extinct.roads.maxdecom.csv", sep='/'), sep=",", row.names=F)
write.table(ext.roads.middecom, paste(Dir.out, "wolf.extinct.roads.middecom.csv", sep='/'), sep=",", row.names=F)
write.table(ext.roads.plandecom, paste(Dir.out, "wolf.extinct.roads.plandecom.csv", sep='/'), sep=",", row.names=F)
write.table(ext.roads.nochange, paste(Dir.out, "wolf.extinct.roads.nochange.csv", sep='/'), sep=",", row.names=F)
write.table(ext.roads.maxbuild, paste(Dir.out, "wolf.extinct.roads.maxbuild.csv", sep='/'), sep=",", row.names=F)

write.table(ext.trap.none, paste(Dir.out, "wolf.extinct.trap.none.csv", sep='/'), sep=",", row.names=F)
write.table(ext.trap.0, paste(Dir.out, "wolf.extinct.trap.0.csv", sep='/'), sep=",", row.names=F)
write.table(ext.trap.20, paste(Dir.out, "wolf.extinct.trap.20.csv", sep='/'), sep=",", row.names=F)
write.table(ext.trap.30, paste(Dir.out, "wolf.extinct.trap.30.csv", sep='/'), sep=",", row.names=F)
write.table(ext.trap.100, paste(Dir.out, "wolf.extinct.trap.100.csv", sep='/'), sep=",", row.names=F)






############### Finally, let's make some figures! ###################


## Scenarios plots----

### Figure 3. Result probability distributions, etc., at end of 30 yrs----

name.arg = c("No New Actions", "Scenario A", "Scenario B", "Scenario C", "Scenario D", "Scenario E")
col1 = "gray20"
col2 = "black"
col3 = "chartreuse4"
col4 = "royalblue"
col5 = "darkorchid3"
col6 = "red3"
col7="orange3"
col8 = "darkred4"
col9 = "salmon2"
col10 = "olivegreen2"


#quartz(width=7, height=4.75)
pdf("./results/figs/figure-3.pdf", width = 7.4, height = 5)

par(mfrow=c(2,3), mar=c(4,4.1,4,2.9)+.1)

lwd=1
lty=1
ylim=c(0,400)
ylim2 = c(0,1)


### Wolf pop through time

plot(wolf.base.s$year, wolf.base.s$mean, type="l", lwd=lwd, col="transparent", bg ="gray90", xlim=c(2015,2045), ylim=ylim, xlab="Year", ylab="Wolf abudance (median)")
lines(wolf.A.s$year, wolf.A.s$median, lwd=lwd, col=col3)
lines(wolf.B.s$year, wolf.B.s$median, lwd=lwd, col=col4)
lines(wolf.C.s$year, wolf.C.s$median, lwd=lwd, col=col5)
lines(wolf.D.s$year, wolf.D.s$median, lwd=lwd, col=col6)
lines(wolf.E.s$year, wolf.E.s$median, lwd=lwd, col=col7)
lines(wolf.base.s$year, wolf.base.s$median, lwd=lwd, col=col2)
mtext("a", at = 2003, line = 2)

### Wolf abundance at end

plot(density(as.numeric(wolf.base.monte[30,])), lwd=lwd, xlim=c(0,525), ylim= c(0, 0.015), xlab="Wolf abundance in 2045", ylab= "Probability density", main = "")
lines(density(as.numeric(wolf.A.monte[30,])), lwd=lwd, col = col3)
lines(density(as.numeric(wolf.B.monte[30,])), lwd=lwd, col = col4)
lines(density(as.numeric(wolf.C.monte[30,])), lwd=lwd, col = col5)
lines(density(as.numeric(wolf.D.monte[30,])), lwd=lwd, col = col6)
lines(density(as.numeric(wolf.E.monte[30,])), lwd=lwd, col = col7)
abline(v = 119, col="grey30", lty=2)
mtext("b", at = -100, line = 2)

### percentage of vacant packs

plot(packs.base$year, packs.base$median, type="l", col="transparent", bg ="gray90", xlim=c(2015,2045), ylim=ylim2, xlab="Year", ylab="% Vacant packs")

lines(packs.A$year, packs.A$median, lwd=lwd, col=col3)
lines(packs.B$year, packs.B$median, lwd=lwd, col=col4)
lines(packs.C$year, packs.C$median, lwd=lwd, col=col5)
lines(packs.D$year, packs.D$median, lwd=lwd, col=col6)
lines(packs.E$year, packs.E$median, lwd=lwd, col=col7)
lines(packs.base$year, packs.base$median, lwd=lwd, col=col2)
mtext("c", at = 2003, line = 2)

### Deer pop through time

ylim=c(30,60)
plot(deer.base$year, deer.base$mean/1000, type="l", lwd=lwd, col="transparent", bg ="gray90", xlim=c(2015,2045), ylim=ylim, xlab="Year", ylab="Deer abudance (median), 1000's")
lines(deer.A$year, deer.A$median/1000, lwd=lwd, col=col3)
lines(deer.B$year, deer.B$median/1000, lwd=lwd, col=col4)
lines(deer.C$year, deer.C$median/1000, lwd=lwd, col=col5)
lines(deer.D$year, deer.D$median/1000, lwd=lwd, col=col6)
lines(deer.E$year, deer.E$median/1000, lwd=lwd, col=col7)
lines(deer.base$year, deer.base$median/1000, lwd=lwd, col=col2)
mtext("d", at = 2003, line = 2)

## Add quasi-extinction probability, up to N = 30, and Densit/probability distributions of scenarios, yr 30

plot(ext.base$q.threshold.n, ext.base$perc.yrs, lwd=lwd, type="l", ylim=c(0,1), xlab = "Quasi-extinction threshold (N)", ylab = "% years below threshold")
lines(ext.A$q.threshold.n, ext.A$perc.yrs, type="l", lwd=lwd, col = col3)
lines(ext.B$q.threshold.n, ext.B$perc.yrs, type="l", lwd=lwd, col = col4)
lines(ext.C$q.threshold.n, ext.C$perc.yrs, type="l", lwd=lwd, col = col5)
lines(ext.D$q.threshold.n, ext.D$perc.yrs, type="l", lwd=lwd, col = col6)
lines(ext.E$q.threshold.n, ext.E$perc.yrs, type="l", lwd=lwd, col = col7)
abline(v = 119, col = "grey30", lty = 2)
mtext("e", at = -20, line = 2)

plot(wolf.base$year, wolf.base$mean, col="transparent", xaxt="null", yaxt="null", xlab="", ylab="", bty="n")

legend(2015, 200, legend=name.arg, box.col="transparent", lwd=c(rep(lwd,6),1), lty=c(rep(lty, 6), 4), col=c(col2, col3, col4, col5, col6, col7, "dodgerblue"))

dev.off()





### Figure 4. Barplots of deer and wolf change in 2045, side by side ----

#quartz(width=3.5, height=6)
pdf("./results/figs/figure-4.pdf", width = 3.5, height = 6)
par(mfrow=c(2, 1), oma= c(1,1,0,0))


wolf.col = "grey70"
deer.col = "saddlebrown"
#quartz(width=6, height=4)
par(mar=c(4,4,2,0)+.1)

b <- barplot(rbind(w.scen.tab$Med.pct*100, d.scen.tab$Med.pct*100), beside=T, col=c(wolf.col, deer.col), ylim=c(-100, 250), ylab="", las=1)
axis(1, labels=F, tick=F, cex = 0.8)
#text(x=colMeans(b), y=-115, srt= 60, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.scen.tab$Med.pct*100, d.scen.tab$Med.pct*100), b, rbind(w.scen.tab$Lcl.pct*100, d.scen.tab$Lcl.pct*100), angle=90, length=0.05)
arrows(b, rbind(w.scen.tab$Med.pct*100, d.scen.tab$Med.pct*100), b, rbind(w.scen.tab$Ucl.pct*100, d.scen.tab$Ucl.pct*100), angle=90, length=0.05)
legend(6, 250, legend=c("Wolf abundance", "Deer abundance"), fill=c(wolf.col, deer.col), box.col="transparent", cex=0.8)
mtext("a", at = -1, line = 1)

### now barplots of hunting and browsing change in 2045, side by side

hunt.col =  "salmon"
browse.col = "darkolivegreen"
#quartz(width=6, height=4)
par(mar=c(7,4,0,0)+.1)

b <- barplot(rbind(b.scen.tab$Med*100, h.scen.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=c(-100, 100), ylab="", las=1)
axis(1, labels=F, tick=F, cex = 0.8)
text(x=colMeans(b), y=-115, srt= 60, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(b.scen.tab$Med*100, h.scen.tab$Med.pct*100), b, rbind(b.scen.tab$Lcl*100, h.scen.tab$Lcl.pct*100), angle=90, length=0.05)
arrows(b, rbind(b.scen.tab$Med*100, h.scen.tab$Med.pct*100), b, rbind(b.scen.tab$Ucl*100, h.scen.tab$Ucl.pct*100), angle=90, length=0.05)
legend(1, 100, legend=c("Browse impacted wolf ranges", "Change in deer harvested"), fill=c(browse.col, hunt.col), box.col="transparent", cex=0.8)
mtext("b", at = -1, line = 1)
mtext("Percent (%)", side=2, line=0, outer=TRUE)

dev.off()

### Figure S1. Sensitivities ----

Tab.pert = read.table(paste(Dir.scen, "sensitivities","Sensitivity perturbation amounts.csv", sep="/"), sep=",", header=T, stringsAsFactors=F)

## First, abundance of wolf and deer sensitivities

#quartz(width=7.4, height=5)
pdf("./results/figs/figure-s1.pdf", width = 7.4, height = 5)

par(mar=c(6.75,2,3.5,1)+.1, mfrow=c(2,3), oma= c(0,3,0,0))
max.group=6
arrow.length = 0.02
ylim=c(-180, 450)
y = -150

name.arg = c("Restore", "No future harvest", "*Transition SG", "Continued OG", "Increased OG", "Max OG")
b <- barplot(rbind(w.K.sens.tab.s$Med.pct*100, d.K.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.K.sens.tab.s$Med.pct*100, d.K.sens.tab$Med.pct*100), b, rbind(w.K.sens.tab.s$Lcl.pct*100, d.K.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(w.K.sens.tab.s$Med.pct*100, d.K.sens.tab$Med.pct*100), b, rbind(w.K.sens.tab.s$Ucl.pct*100, d.K.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("a", at = -2, line = 2)
legend(4, 400, legend=c("Wolf abundance", "Deer abundance"), fill=c(wolf.col, "saddlebrown"), box.col="transparent")


name.arg = c("Max decom", "Mid decom","*Planned decom", "No change", "Construction")
b <- barplot(rbind(w.roads.sens.tab.s$Med.pct*100, d.roads.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.roads.sens.tab.s$Med.pct*100, d.roads.sens.tab$Med.pct*100), b, rbind(w.roads.sens.tab.s$Lcl.pct*100, d.roads.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(w.roads.sens.tab.s$Med.pct*100, d.roads.sens.tab$Med.pct*100), b, rbind(w.roads.sens.tab.s$Ucl.pct*100, d.roads.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("b", at = -2, line = 2)

name.arg = c("Low freq.", "*Ave. freq.", "High freq.")
b <- barplot(rbind(w.climate.sens.tab.s$Med.pct*100, d.climate.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.climate.sens.tab.s$Med.pct*100, d.climate.sens.tab$Med.pct*100), b, rbind(w.climate.sens.tab.s$Lcl.pct*100, d.climate.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(w.climate.sens.tab.s$Med.pct*100, d.climate.sens.tab$Med.pct*100), b, rbind(w.climate.sens.tab.s$Ucl.pct*100, d.climate.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("c", at = -2, line = 2)

name.arg = c("9.5 deer/yr", "*15 deer/yr", "20.5 deer/yr", "26 deer/yr")
b <- barplot(rbind(w.diet.sens.tab.s$Med.pct*100, d.diet.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.diet.sens.tab.s$Med.pct*100, d.diet.sens.tab$Med.pct*100), b, rbind(w.diet.sens.tab.s$Lcl.pct*100, d.diet.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(w.diet.sens.tab.s$Med.pct*100, d.diet.sens.tab$Med.pct*100), b, rbind(w.diet.sens.tab.s$Ucl.pct*100, d.diet.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("d", at = -2, line = 2)

name.arg = c("No harvest", "0% legal", "*20% legal cap", "30% legal cap", "No limit", "Extinct")
b <- barplot(rbind(c(w.trap.sens.tab.s$Med.pct*100, -100), d.trap.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y=y, srt= 45, adj = 1, labels= name.arg[1:7], xpd=T)
arrows(b, rbind(c(w.trap.sens.tab.s$Med.pct*100, NA), d.trap.sens.tab$Med.pct*100), b, rbind(c(w.trap.sens.tab.s$Lcl.pct*100, NA), d.trap.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(c(w.trap.sens.tab.s$Med.pct*100, NA), d.trap.sens.tab$Med.pct*100), b, rbind(c(w.trap.sens.tab.s$Ucl.pct*100, NA), d.trap.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("e", at = -2, line = 2)

name.arg = c("Hunt no deer", "*Hunt regular")
b <- barplot(rbind(w.hunt.sens.tab.s$Med.pct*100, d.hunt.sens.tab$Med.pct*100), beside=T, col=c(wolf.col, "saddlebrown"), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1)
axis(1, labels=F, tick=F, las=1)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(w.hunt.sens.tab.s$Med.pct*100, d.hunt.sens.tab$Med.pct*100), b, rbind(w.hunt.sens.tab.s$Lcl.pct*100, d.hunt.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(w.hunt.sens.tab.s$Med.pct*100, d.hunt.sens.tab$Med.pct*100), b, rbind(w.hunt.sens.tab.s$Ucl.pct*100, d.hunt.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("f", at = -2, line = 2)

mtext("Change in abundance (%)", side=2, line=1, outer=TRUE)
dev.off()

### Figure S2. Ecosystem service (browse and hunting) sensitivities ----

#quartz(width=7.4, height=5)
pdf("./results/figs/figure-s2.pdf", width = 7.4, height = 5)

par(mar=c(6.75,3,3.5,0)+.1, mfrow=c(2,3), oma= c(0,3,0,0))

#par(mar=c(8,2,3.5,1)+.1, mfrow=c(2,3), oma= c(0,5,0,0))
max.group=6
arrow.length = 0.02
ylim=c(-100, 100)
y=-100

name.arg = c("Restore", "No future harvest", "*Transition SG", "Continued OG", "Increased OG", "Max OG")
b <- barplot(rbind(b.K.sens.tab$Med*100, h.K.sens.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(b.K.sens.tab$Med*100, h.K.sens.tab$Med.pct*100), b, rbind(b.K.sens.tab$Lcl*100, h.K.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(b.K.sens.tab$Med*100, h.K.sens.tab$Med.pct*100), b, rbind(b.K.sens.tab$Ucl*100, h.K.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("a", at = -2, line = 2)
legend(1, 100, legend=c("Browse impacted pack ranges", "Change in deer harvested"), fill=c(browse.col, hunt.col), box.col="transparent")


name.arg = c("Max decom", "Mid decom","*Planned decom", "No change", "Construction")
b <- barplot(rbind(b.roads.sens.tab$Med*100, h.roads.sens.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(b.roads.sens.tab$Med*100, h.roads.sens.tab$Med.pct*100), b, rbind(b.roads.sens.tab$Lcl*100, h.roads.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(b.roads.sens.tab$Med*100, h.roads.sens.tab$Med.pct*100), b, rbind(b.roads.sens.tab$Ucl*100, h.roads.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("b", at = -2, line = 2)

name.arg = c("Low freq.", "*Ave. freq.", "High freq.")
b <- barplot(rbind(b.climate.sens.tab$Med*100, h.climate.sens.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(b.climate.sens.tab$Med*100, h.climate.sens.tab$Med.pct*100), b, rbind(b.climate.sens.tab$Lcl*100, h.climate.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(b.climate.sens.tab$Med*100, h.climate.sens.tab$Med.pct*100), b, rbind(b.climate.sens.tab$Ucl*100, h.climate.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("c", at = -2, line = 2)


name.arg = c("9.5 deer/yr", "*15 deer/yr", "20.5 deer/yr", "26 deer/yr")
b <- barplot(rbind(b.diet.sens.tab$Med*100, h.diet.sens.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:6], xpd=T)
arrows(b, rbind(b.diet.sens.tab$Med*100, h.diet.sens.tab$Med.pct*100), b, rbind(b.diet.sens.tab$Lcl*100, h.diet.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(b.diet.sens.tab$Med*100, h.diet.sens.tab$Med.pct*100), b, rbind(b.diet.sens.tab$Ucl*100, h.diet.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("d", at = -2, line = 2)


name.arg = c("No harvest", "0% legal", "*20% legal cap", "30% legal cap", "No limit", "Extinct")
b <- barplot(rbind(b.trap.sens.tab$Med*100, h.trap.sens.tab$Med.pct*100), beside=T, col=c(browse.col, hunt.col), ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=colMeans(b), y = y, srt= 45, adj = 1, labels= name.arg[1:7], xpd=T)
arrows(b, rbind(b.trap.sens.tab$Med*100, h.trap.sens.tab$Med.pct*100), b, rbind(b.trap.sens.tab$Lcl*100, h.trap.sens.tab$Lcl.pct*100), angle=90, length=arrow.length)
arrows(b, rbind(b.trap.sens.tab$Med*100, h.trap.sens.tab$Med.pct*100), b, rbind(b.trap.sens.tab$Ucl*100, h.trap.sens.tab$Ucl.pct*100), angle=90, length=arrow.length)
mtext("e", at = -2, line = 2)


name.arg = c("","Hunt no deer", "*Hunt regular")
b <- barplot(c(NA, b.hunt.sens.tab$Med*100), space = c(0.5, 0.5), beside=T, col=browse.col, ylim=ylim, ylab="", xlim=c(0, max.group*3), width=1, las = 1)
axis(1, labels=F, tick=F)
text(x=b, y = y, srt= 45, adj = 1, labels= name.arg, xpd=T)
arrows(b, c(NA, b.hunt.sens.tab$Med*100), b, c(NA, b.hunt.sens.tab$Lcl*100), angle=90, length=arrow.length)
arrows(b, c(NA, b.hunt.sens.tab$Med*100), b, c(NA, b.hunt.sens.tab$Ucl*100), angle=90, length=arrow.length)
mtext("f", at = -2, line = 2)

mtext("Proportion (%)", side=2, line=1, outer=TRUE)



dev.off()


