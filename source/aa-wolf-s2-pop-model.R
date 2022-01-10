##############################################################################################
#########################################################################################
# This code runs a predator-prey-habitat model, based on wolf-deer dyanmics in Southeast Alaska
# The code was produced by Sophie Gilbert (U Idaho), 1-7-2022
# GitHub repository for code, data, etc., at https://github.com/sophielgilbert/alexander-archipelago-wolf

# title         : aa-wolf-s2-pop-model.r
# purpose       : Predicts changes in wolf and deer abundance and ecosystem services from 2015- 2045
# data inputs	  : Tables of road density and deer carrying capacity for different scenarios developed in stakeholder meeting (see Gilbert et al. 2015)
# author        : Sophie Gilbert
# input         : Binary 16-bit signed integer big endian byte-order
# output        : Prediction tables for Monte Carlo iterations of model runs for each scenario
# additional	  : Seperate R code to build figures from model outputs is available at: https://github.com/sophielgilbert/alexander-archipelago-wolf

# publication 1	: Gilbert, S., M. Lindberg, T. Haynes, M. Kissling, D. Albert. 2015. 
#					Future population trends and drivers of change for Alexander Archipelago wolves on and near Prince of Wales Island, Alaska. 
#					Final Report to the U.S. Fish & Wildlife Service, Anchorage, AK.
#					Available at: https://www.fws.gov/r7/fisheries/endangered/pdf/aa_wolf/Gilbert%20et%20al%202015%20Population%20Model_final.pdf

# publication 2	: Gilbert, S., M. Lindberg, T. Haynes, M. Kissling, D. Albert, D. Person. In Press. 
#					Future population trends and drivers of change for Alexander Archipelago wolves on and near Prince of Wales Island, Alaska. 
#					In preparation for Frontiers in Ecology and Evolution, special edition on wolf management
#					Available at: INSERT URL HERE ONCE PUBLISHED
 
################################################################################

########### 1. Load R libraries needed to run analyses ###########

rm(list = ls())																				# clear r's memory of any stored objects 							  
library(popbio)																				# popbio has a lot of good tools for population dynamics
library(reshape2)																			# this lets us reshape matrices and arrays (change the dimensions)
library(abind)																				# this lets us bind arrays (like matrices) into 3-D stacks, like book pgs



########### 2. Define directories where data are stored ###########

# This requires input data for habitat variables such as deer K and road length/density
# if these change over time, need to load that too

DirRoad <- "./data-processed/roads/"														# Where in main folder are the .csvs for roads for each scenario?  
DirK <- "./data-processed/deer-k/"														  # Where in main folder are the .csvs for deer K "
DirOut <- "./results/"											 	                  # Specify where you want to store results. Be careful here!


########### 3. Read in conditions for scenarios and sensitivities ###########

clim.tab <- c(0.07, 0.08, 0.10)																		# probs. of severe winter under clim. scenarios
trap.tab <- c(0, 0.2, 0.3, 1.0)																		# caps on reported harvest *NOTE 9/1/21: now allowing elimination of legal harvest
packs.tab <- read.table("./data/wolf-inputs/pack-characteristics.csv", sep=",", header=T)			# pack characteristics

roads.Base = read.table(paste(DirRoad, "Data.roads.base.csv", sep="/"), sep=",", header=T)			# roads for each pack and year, scenario base
roads.A = read.table(paste(DirRoad, "Data.roads.A.csv", sep="/"), sep=",", header=T)				# roads for each pack and year, scenario A
roads.B <- roads.C <- roads.D <- roads.A															# roads for each pack and year, scenarios B, C, D = A
roads.E = read.table(paste(DirRoad, "Data.roads.E.csv", sep="/"), sep=",", header=T)				# roads for each pack and year, scenario E
roads.mid.decom = read.table(paste(DirRoad, "Data.roads.mid.decom.csv", sep="/"), sep=",", header=T) # roads for each pack and year, scenario E
roads.max.decom = read.table(paste(DirRoad, "Data.roads.max.decom.csv", sep="/"), sep=",", header=T) # roads for each pack and year, scenario E

K.Base = read.table(paste(DirK, "Data.K.base.csv", sep="/"), sep=",", header=T)						# K for each pack and year, scenario base
K.Nochange = read.table(paste(DirK, "Data.K.Nochange.csv", sep="/"), sep=",", header=T)				# K for each pack and year, scenario base
K.A <- K.Base																						                              # K for each pack and year, scenario A = base
K.B = read.table(paste(DirK, "Data.K.B.csv", sep="/"), sep=",", header=T)							# K for each pack and year, scenario B
K.C = read.table(paste(DirK, "Data.K.C.csv", sep="/"), sep=",", header=T)							# K for each pack and year, scenario C
K.D = read.table(paste(DirK, "Data.K.D.csv", sep="/"), sep=",", header=T)							# K for each pack and year, scenario D
K.E = read.table(paste(DirK, "Data.K.E.csv", sep="/"), sep=",", header=T)							# K for each pack and year, scenario E



########### 4. Make scenarios ###########

Scen.names = c("roads", "K", "climate", "trap", "wolf.harv", "deer.harv", "wolf.diet")				# Name the parts of each scenario's list elements

Scen.Base 	= list(roads.Base, K.Base, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)					# Baseline scenario (i.e., no change except sucession)
Scen.A 		= list(roads.A, K.A, clim.tab[1], trap.tab[1], TRUE, TRUE, 15)							# Scenario A
Scen.B 		= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)							# Scenario B
Scen.C 		= list(roads.C, K.C, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)							# Scenario C
Scen.D 		= list(roads.D, K.D, clim.tab[3], trap.tab[3], TRUE, TRUE, 15)							# Scenario D
Scen.E 		= list(roads.E, K.E, clim.tab[3], trap.tab[3], TRUE, TRUE, 15)							# Scenario E
names(Scen.Base) <-names(Scen.A) <-names(Scen.B) <-names(Scen.C) <-names(Scen.D) <-names(Scen.E) <-Scen.names



########### 5. Make sensitivity analyses combinations ###########

Sensitivity.deer.harvest 	= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, FALSE, 15)			# Sensitivity to deer harvest being eliminated
Sensitivity.roads.noch 		= list(roads.Base, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)		# Sens., no change to roads
Sensitivity.roads.mid.decom 	= list(roads.mid.decom, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)	# Sens., increased decommission
Sensitivity.roads.max.decom = list(roads.max.decom, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)	# Sens., max decomission
Sensitivity.roads.max.build = list(roads.E, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., increased building

Sensitivity.K.baseline 		= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., only succession
Sensitivity.K.nochange 		= list(roads.B, K.Nochange, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)	# Sens., no change in K
Sensitivity.K.cont.harvest	= list(roads.B, K.C, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., continued old growth harvest
Sensitivity.K.mid.harvest 	= list(roads.B, K.D, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., mid old growth harvest
Sensitivity.K.max.harvest	= list(roads.B, K.E, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., max old growth harvest

Sensitivity.climate.low 	= list(roads.B, K.B, clim.tab[1], trap.tab[2], TRUE, TRUE, 15)			# Sens., low snow
Sensitivity.climate.mid 	= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 15)			# Sens., mid snow
Sensitivity.climate.high	 = list(roads.B, K.B, clim.tab[3], trap.tab[2], TRUE, TRUE, 15)			# Sens., max snow

Sensitivity.diet.9.5		= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 9.5)			# Sens., 25% deer in diet
Sensitivity.diet.20.5		= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 20.5)		# Sens., 50% deer in  diet
Sensitivity.diet.26			= list(roads.B, K.B, clim.tab[2], trap.tab[2], TRUE, TRUE, 26)			# Sens., 75% deer in diet

Sensitivity.wolftrap.none	= list(roads.B, K.B, clim.tab[2], trap.tab[1], FALSE, TRUE, 15)			# Sens., no wolf harvest at all (no poaching, not really possible)
Sensitivity.wolftrap.0		= list(roads.B, K.B, clim.tab[2], trap.tab[1], TRUE, TRUE, 15)			# Sens., 0 reported wolf harvest (poaching only)
Sensitivity.wolftrap.30		= list(roads.B, K.B, clim.tab[2], trap.tab[3], TRUE, TRUE, 15)			# Sens., 30% reported wolf harvest (plus poaching)
Sensitivity.wolftrap.100	= list(roads.B, K.B, clim.tab[2], trap.tab[4], TRUE, TRUE, 15)			# Sens., up to 100% wolf harvest possible (no regs)

names(Sensitivity.deer.harvest)  <- Scen.names
names(Sensitivity.roads.noch) <- names(Sensitivity.roads.mid.decom) <- names(Sensitivity.roads.max.decom) <-names(Sensitivity.roads.max.build) <- Scen.names
names(Sensitivity.K.baseline) <-names(Sensitivity.K.nochange) <- names(Sensitivity.K.cont.harvest) <- names(Sensitivity.K.mid.harvest) <- names(Sensitivity.K.max.harvest) <- Scen.names
names(Sensitivity.climate.low) <-names(Sensitivity.climate.mid) <- names(Sensitivity.climate.high) <- Scen.names
names(Sensitivity.diet.9.5) <-names(Sensitivity.diet.20.5) <- names(Sensitivity.diet.26) <- Scen.names
names(Sensitivity.wolftrap.none) <-names(Sensitivity.wolftrap.0) <- names(Sensitivity.wolftrap.30) <- names(Sensitivity.wolftrap.100)<- Scen.names

Sensitivity.wolfgone <- Scen.B

Scen.to.run = list(Scen.Base, Scen.A, Scen.B, Scen.C, Scen.D, Scen.E)
Scen.to.run.folds = c("ScenBase", "ScenA", "ScenB", "ScenC", "ScenD", "ScenE")

Sens.to.run = list(	Sensitivity.deer.harvest,
					Sensitivity.roads.noch,
					Sensitivity.roads.mid.decom,
					Sensitivity.roads.max.decom,
					Sensitivity.roads.max.build,
					Sensitivity.K.baseline,
					Sensitivity.K.nochange,
					Sensitivity.K.cont.harvest,
					Sensitivity.K.mid.harvest,
					Sensitivity.K.max.harvest,
					Sensitivity.climate.low,
					Sensitivity.climate.mid, 		
					Sensitivity.climate.high,	 
					Sensitivity.diet.9.5,		
					Sensitivity.diet.20.5,		
					Sensitivity.diet.26,
					Sensitivity.wolftrap.none,	
					Sensitivity.wolftrap.0,
					Sensitivity.wolftrap.30,
					Sensitivity.wolftrap.100,
					Sensitivity.wolfgone)
					
Sens.to.run.names <- c("Sensitivity.deer.harvest",
						"Sensitivity.roads.noch",
						"Sensitivity.roads.mid.decom",
						"Sensitivity.roads.max.decom",
						"Sensitivity.roads.max.build",
						"Sensitivity.K.baseline",
						"Sensitivity.K.nochange",
						"Sensitivity.K.cont.harvest",
						"Sensitivity.K.mid.harvest",
						"Sensitivity.K.max.harvest",
						"Sensitivity.climate.low",
						"Sensitivity.climate.mid", 		
						"Sensitivity.climate.high",	 
						"Sensitivity.diet.9.5",		
						"Sensitivity.diet.20.5",		
						"Sensitivity.diet.26",
						"Sensitivity.wolftrap.none",	
						"Sensitivity.wolftrap.0",
						"Sensitivity.wolftrap.30",
						"Sensitivity.wolftrap.100",
						"Sensitivity.wolfgone")	

Sens.to.run.folds <- c("SenHunt",
						"SenRoads/RoadNochange",
						"SenRoads/RoadMiddecom",
						"SenRoads/RoadMaxdecom",
						"SenRoads/RoadMaxbuild",
						"SenK/KBase",
						"SenK/KNochange",
						"SenK/KContharvest",
						"SenK/KMidharvest",
						"SenK/KMaxharvest",
						"SenClimate/SnowLow",
						"SenClimate/SnowMid",
						"SenClimate/SnowMax",
						"SenDiet/Diet9.5",
						"SenDiet/Diet20.5",
						"SenDiet/Diet26",
						"SenTrap/TrapNone",
						"SenTrap/Trap0",
						"SenTrap/Trap30",
						"SenTrap/Trap100",
						"SenWolf/WolfGone")
						


########### Running the models #################################################

########### 6. Choose which individual scenario or sensitivity analysis to run #########


Scenario = Sensitivity.wolfgone					# chose scenario to run
Kind = "Sensitivities"									# Lets us choose where to store results
Which.Scen = paste(Kind, "SenWolf/WolfGone", sep="/")	# for storing results; see end of code


#Scenario = Scen.D										# chose scenario to run
#Kind = "Scenarios"										# Lets us choose where to store results
#Which.Scen = paste(Kind, "ScenD", sep="/")				# for storing results; see end of code






##########  7. Or we can loop through a list of all the scenarios or sensitivities, running each in turn ####
# Rather than adding yet another loop to run through scenarios and sensitivities collectively
# The user must manually set Kind = ... (below) to choose whether to run all scenarios OR all sensitivities
# Outputs are saved automatically to the /results/scenarios or sensitivities  folders

#Kind = "Scenarios"                             # NEED TO MANUALLY CHOOSE HERE, run both lines
Kind = "Sensitivities"                          # NEED TO MANUALLY CHOOSE HERE, run both lines

if(Kind =="Scenarios"){													# chose the correct inputs (scenarios or sensitivities) ond output folders
	To.run = Scen.to.run;
	To.run.folds = Scen.to.run.folds}else{
	To.run = Sens.to.run;
	To.run.folds = Sens.to.run.folds}




for(f in 1:length(To.run)){
	Scenario = To.run[[f]]
	Which.Scen =paste(Kind, To.run.folds[[f]], sep="/")
	print(Which.Scen)



########### Basic parameters of wolves and deer population ###########

### Overall parameters

monte.length 	= 1000													# chose number of Monte Carlo simulations to run
tmax 			= 30												      	# number of years										
npack 			= nrow(packs.tab)										# number of packs		
years 			= seq(from=2015, to = 2015 + tmax + 1, by = 1)			# sequence of years to simulate through
winter.sev 		= Scenario$climate										# frequency of severe winter, from scenario


### Deer parameters

rmax 				= 0.6												# maximum reproductive rate for deer
theta 				= 2													# time lag in yrs for density dependence in deer population
bear.fawn 			= 0.46												# mean fawn prodation rate by bears
bear.fawn.sd			= 0.023											# SE of fawn predation rate by bears
bear.dd	 			= 0.5												# maximum amount of bear predation on fawns that can be additive
bear.adult 			= 0.03												# mean adult deer predation rate by bears
bear.adult.se 		= 0.0015											# SE adult deer predation rate by bears

if(Scenario$deer.harv ==TRUE){											# if deer hunting allowed at all...
h = 0.012}else{ h = 0}													# baseline hunting rate of deer, from Person 1997

start.deer.min 		= 0.5												# choose "floor" for deer start relative to carrynig capacity... 0-1 = possible range
start.deer.max		= 1													# choose ceiling


### Wolf parameters


mean.density.2014 	= 9.9												# wolf population survey from 2014 (ADF&G 2015), wolves/1000km2
sd.density.2014 	= 3.0												# sd of wolf density	
density 			= mean.density.2014									# density from 2014 estimates ADFG	
most.recent.harvest = round(7*31/18)									# reported harvest season (2014/15), times unreported scalar
trap.cap 			= Scenario$trap										# cap (%) on reported wolf harvest, from scenario
pack.areas 			= packs.tab$land_sqkm								# pack areas, in km sq.
pack.ocean 			= packs.tab$avg_shore_dist							# pack mean distance of all shoreline to nearest town									
pack.roads 			= Scenario$roads[,2:length(years)]					# roads per pack per year, from scenario
pack.k				= Scenario$K[,2:length(years)]						# deer K per pack per year, from scenario	
C 					= Scenario$wolf.diet								# predation rate, deer/wolf/year
C.sd 				= 4													# SD of predation rate
litter 				= 4.1												# average litter size
litter.max 			= 11												# maximum litter size
pack.mean 			= 6													# mean wolf pack size
pack.max 			= 18												# maximum reasonable pack size for wolves
wolf.lag 			= 2													# how many years should wolves lag deer
cmort 				= 0.5												# mean chronic mortality rate at high densities
cmort.sd 			= 0.3												# SD of chronic mortality rate at high densities
disp 				= 0.5												# probability of dispersal at high densities
disp.sd 			= 0.3												# sd of dispersal prob
surv.disp 			= 0.34												# mean survival of dispersers
surv.disp.sd 		= 0.3												# sd survival of dispersers
road.den.in 		= packs.tab$roads_2015_km/packs.tab$land_sqkm		# road density at start
qextinct			= 10												# quasi-extinction threshold

if(Scenario$wolf.harv==TRUE){											# if wolf harvest allowed at all...
unreported.scalar = 31/18;												# tot. wolves killed by humans/wolves killed legally, Person & Russell 2008 p. 1545 
trap.prop.rep.mort = (0.23)/(0.23 + 0.19 + 0.04); 						# proportion of mort. from legal trapping, for disperser pool. Person & Russell 2008 p. 1545
trap.prop.mort = (0.23 + 0.19)/(0.23 + 0.19 + 0.04)}else{				# proportion of mort. that is from trapping, both legal and illegal, for disperser pool.
	unreported.scalar = 0;
	trap.prop.rep.mort = 0;
	trap.prop.mort = 0}

Tab.5.regression 	= TRUE												# If TRUE, will use the Tab 5 single regression. If false, Tab 6 regressions, Person & Russell 2008

main.t = "density= ADFG 2014, \n pack.mean = 6, \n start.deer.min = 0.5, \n unreported.scalar = 1.72, \n Tab 5 trap rergession"


### Predation services parameters
# 30 deer/km2, sustained for 50 years, results on almost complete loss of understory, including cedar (Jean-Louise Martin work, and personnal communication, Haida Gwaii)
# Howevever, deer K in our model is based off of USFS deer HCI, which may simply not allow deer to reach this density
# And deer impacts to browse are likely to be more about proximity to K rather than absolute density
# So we will invoke browse impacts of deer are at 90% of K or above
deer.high.K = 	0.90			# If deer reach 90% of K, induce browsing effects		
deer.browse.time = 20			# sustain high density for 20 years to lose cedar, and to delay conifer regenaration
deer.low.K = 0.10	# Deer must go down to 10% of K (3/35, 3-4 deer per km^2 was the threshold on Haida Gwaii, 35/km was max reported 
deer.conifer.time = 20			# If deer are above threshold density, conifer regeneration is delayed ~20 yrs
								# Not sure how to parameterize soil carbon loss yet, still reading and thinking...



# SELECT ALL CODE BELOW THIS POINT AND RUN IT TO SEE THE MONTE CARLO RESULTS GRAPHED 

########### Monte Carlo iterations of model for chosen scenario ###########

																
pdf(paste(DirOut,Which.Scen,"wolf_deer_plot.pdf", sep="/"), height=6, width=6)
	for(m in 1:monte.length){														####################loop through monte carlo interations ###################

		######## make starting wolf population

				
		density.pop = rnorm(10, mean=density, sd = sd.density.2014)						# generate 10 wolf densities per 1000km1, but we need only one
		density.pop = density.pop[which(density.pop>0)]									# wolf density can't be below zero, rule these out
		density.pop = density.pop[1]													# take first density of remaining, non-zero densities
		pop.area = sum(pack.areas)/1000													# total population (study) area, converted to 1000 sq. km units
		pop.sum = round(pop.area*density.pop	)										# total wolf pop in study area at start
		pop.sum = pop.sum - most.recent.harvest											# reduce starting pop size (fall est.) by winter harvest
		if(pop.sum<0){pop.sum=0}														# can't have a starting population size of less than zero

		P.t.all = rnorm(length(pack.areas), mean=pack.mean)								# generate starting wolf sizes for each pack based on mean pack size
		P.t.all = round(P.t.all)

		P.t.all = data.frame(packs.tab, nstart = P.t.all, road.rank=rank(road.den.in))	# store these starting pack sizes
		p.rand = sample(1:nrow(P.t.all), nrow(P.t.all))									# random order for stocking packs, if desired
		P.t.all = P.t.all[order(p.rand),]												        # stock packs randomely
		#P.t.all = P.t.all[order(P.t.all$road.rank),]									  # stock packs in order of low to high road density
		
		for(p in 1:nrow(P.t.all)){														# loop to assign wolves to the population
			if(p==1){accumulate = P.t.all$nstart[p]}else{								# lowest-roaded pack gets all its wolves
			accumulate=c(accumulate, sum(P.t.all$nstart[1:p]))}							# count how many total wolves being assigned down the list
			}

		P.t.all$accumulate = accumulate													# cumulative wolves assigned across packs
		if(max(accumulate)>=pop.sum){too.much = min(which(accumulate>=pop.sum))}else{
		too.much = 31}																	# which pack pushes the total pop above random start pop?
		if(too.much <=1){too.much = 2}													# we need to fill at least 1 packs worth w/wolves
		pack.marg = round(pop.sum -accumulate[too.much-1])  							# this pack will get a few wolves, just not all starting wolves
		if(pack.marg<1){pack.marg = 0}
		P.t.all$nstart[too.much:nrow(P.t.all)] = 0										# no wolves above threshold, pack by pack
		P.t.all$accumulate[too.much:nrow(P.t.all)] = 0									# no wolves above threshold, in cumulative pop
		P.t.all$nstart[too.much] = pack.marg											# give marginal wolves to marginal pack		
		if(pack.marg > 0) {P.t.all$accumulate[too.much] = 								# add marginal wolves to cumulative pop
			P.t.all$accumulate[too.much - 1] + pack.marg	}								
		P.t.all <- P.t.all[order(P.t.all$Pack_ID),]	
		
		if(Which.Scen==paste(Kind, "SenWolf/WolfGone", sep="/")){                 # If running a "no wolf" simulation, start w/1 wolf
		  P.t.all$nstart = c(1, rep(0, (nrow(P.t.all)-1)))}
				
		for(t in 1:(length(years)-1)){														###################### loop through time, 1:tmax ################
			
			winter = rbinom(1,1, prob= winter.sev)											# generate winter severity for each year, binomial distribution
					
				for(i in 1:npack){																###################### loop through packs each pack, 1:npack, in year t	
					
					K.t = pack.k[i,t]															# deer K for time t, pack i
					road.t = pack.roads[i,t]													# roads for time t, pack i
					pack.area = pack.areas[i]															
			
				##### initial wolf population size, from previous time-step
					if(t==1){P.t = P.t.all$nstart[i]}										# get randomly-generated pack starts sizes, from above				
					if(t==1 & P.t < 0){P.t = 0}												# truncate the number in pack i in year t to >= 0
					if(t==1 & P.t > 18){P.t = 18}											# truncate the number in pack i in year t to <= 18
					if(t>1){p.prev = data.frame(t(packs.time[i,,t-1]))}						# access the stored data for this pack from t-1
					if(t>1){P.t = (as.numeric(p.prev$fall.wolf) 							# after year 1, wolf pack pop at start of interval t = ...
					- as.numeric(p.prev$wolf.harv.tot) 										# fall walf pack pop in (t-1) minus wolf pack harvest (t-1)...
					- as.numeric(p.prev$dispersers) 										# minus wolf pack dispersal (t-1)...
					- as.numeric(p.prev$cmort) 												# minus wolf pack chronic mort (t-1)...
					+ as.numeric(p.prev$immigrants))}										# plus wolf pack immigrants (t-1)							
					if(P.t <0){P.t = 0}														# no negative wolf numbers allowed
					P.t = round(P.t, digits=0)												# no decimal wolf numbers allowed

				##### Initial deer K and pop size, from previous time-step
					start.prop = runif(1, start.deer.min, start.deer.max)					# year 1 deer proximity to deer K = random prop btwn X and 1
					alpha.start = runif(1, .5, 1)											# year 1 availability of deer to wolves = random prop btwn 0.5 and 1
					if(t==1){alpha.t = runif(1, .5, 1)}										# after yr 1, avail. deer is rand proportion from .5 to 1
					if(t==1){U.t= start.prop*K.t}											# for t = 1, deer 
					if(t>1){U.t = (as.numeric(p.prev$spring.deer) 							# after year 1, deer pop in pack i at start of interval t = ... 
							+ as.numeric(p.prev$recruit.deer) 								# spring deer in pack area i plus deer recruits in pack area i...
							- as.numeric(p.prev$hunt.deer)									# minus 	deer hunted in pack area i...
							- as.numeric(p.prev$pred.deer) 									# minus deer predated in pack area i by wolves...
							- as.numeric(p.prev$bear.deer))}								# minus adult deer predated in pack are i by bears, all in t-1
					if(U.t<0){U.t=0}														# no negative deer numbers allowed
					U.t = round(U.t, digits=0)												# no decimal deer numbers allowed		
		
				##### Wolf predation on deer
					C.t = rnorm(1, C, C.sd)													# deer eaten per wolf per year in pack i is a random #
					if(C.t <= 3.4){C.t = 3.4}												# 10% deer: wolves eat at least 3.4 deer per wolf per year
					if(C.t >= 34){C.t = 34}													# 100% deer: wolves no more than 34 deer per yr; Person & Bowyer 1996 app.3
					CP.t = round(P.t*C.t, 0)												# Total deer eaten by wolf pack i in time t; no decimal deer

				##### Calculate deer and wolf Density-dependent states	
				if(t<=wolf.lag){															# For wolf pop, cannot use wolf time-lagged response in 1st 2 yrs...
						rat = CP.t/(alpha.t*U.t + 1)}else{									# so instead, calculate "ratio" (deer needed:deer available) w/out lag...
						s.deer = names(packs.time[i,,(t-wolf.lag)]);						# where are spring deer stored
						s.deer = which(s.deer=="spring.deer");								# in which specific column
						rat = CP.t/(alpha.t*as.numeric(packs.time[i,s.deer,(t-wolf.lag)]) + 1)}	# but otherwise, ratio is based on deer available in (t-2), a 2-yr lag			
					if(rat > 1){															# if more deer needed than are available...
						ratio = 1}else{														# set ratio to max (which is 1)...otherwise,
						ratio = rat}														# ratio is just rat, which we calculated above
					DD.t <- (1-(U.t/K.t))													# For deer pop, DD via pop size (t-1) relative to deer K(t)

				##### Wolf reproduction (DD)
					pups <- litter*(1-ratio)												# DD wolf pup production, smaller litters as wolves approach K
					if(P.t >=2){L.t = rpois(1, pups)} else{L.t = 0}							# litter size is rand drawn from a poisson distribution
					if(L.t >= litter.max){L.t==litter.max}									# cap litter size, set by parameter litter.max
					if(P.t < 2){L.t = 0}														# cannot have pups if fewer than 2 wolves in a pack at start of year t
					if(L.t < 0){L.t = 0}														# cannot have fewer than 0 pups in a year
					Rp.t = round(L.t, 0)														# recruitment depends on DD pup production, and DD pup survival

				##### Deer population changes in year t
					Pa.t = (2*P.t + Rp.t)/2													# average wolf population size for t, used for predation on deer
					bear.fawn.t = betaval(bear.fawn, bear.fawn.sd)							# bear predation on fawns (DD)
					if(DD.t < bear.dd){														# if deer pop at > 0.5 of K...
						BF.t = bear.dd*bear.fawn.t}else{										# bear pred on fawns increasingly compensatory...
						BF.t = DD.t*bear.fawn.t}												# otherwise (deer pop < 0.5 of K), bear pred on fawns is additive
					BA.t = round(betaval(bear.adult, bear.adult.se)*U.t)						# predation rate of bears on adult deer
					R.t = round((1-BF.t)*U.t*rmax*(1-(U.t/K.t)^theta))						# DD fawn recruitment (theta logistic eq.), after bear predation
					if(winter==1){R.t=0}														# if it's a severe winter, deer recruitment goes to zero
					H.t = (U.t + R.t - BA.t)*h*(1 + 0.038*road.t)							# hunting ~ road length. No deer predation in here, hunting takes precedence		
					if(H.t> (U.t + R.t - BA.t)){												# can't hunt more deer than exist (spring deer + recruits - bear pred adults)...
						H.t= (U.t + R.t - CP.t*Pa.t - BA.t + 2)}								# but can hunt up to that number, but keep 2 deer  								
					H.t = round(H.t, 0)														# no decimal deer allowed to be hunted
					if(H.t<0){H.t=2}															# deer hunting cannot be < 2 
					spring.deer = U.t 														# Deer numbers to remember for the start of the (next) time step	
					if(spring.deer < 50){spring.deer = 50}			
					
					Deer.dt = U.t/pack.area													# Deer per km^2 per pack area
					Deer.K.prop = U.t/K.t													# How close are deer to K?
					if(Deer.K.prop>=deer.high.K){browse=1}else{browse=0}					# Record if deer density is high enough for conifefer browsing effects
					if(Deer.K.prop<=deer.low.K){recover=1}else{recover=0}					# Record if deer density is low enough for conifer recovery
																	
					
				##### Wolf population changes in year t
				
				#### Wolf harvest
					if(CP.t/(alpha.t*(U.t+1)) >= 50){deer.gone=T}else{deer.gone=F}			# the "threshold" below which wolf packs break up: 50 deer
					
					if(Tab.5.regression == FALSE){
						if(road.t/pack.area <= 0.9){											# predict wolf reported harvest per 100 km^2, is a function of road den
						HR.t = (0.073 + 1.126*(road.t/pack.area))^2}else{					# if rd density is <= 0.9 km/km^2, use this eq... 
						HR.t = (0.952 - 0.009*(pack.ocean[i]))^2}}else{						# but if rd density above 0.9, use this eq (Person & Russell 2008, p. 1547)	
							HR.t = (1.010 - 0.005*(pack.ocean[i]) + 0.207*(road.t/pack.area))^2}	 # Or use this single eq., (Person & Russell, Table 5) for all rd. densities
					
					wolf.reported = HR.t*pack.area/100										# wolf reported harvest predicted for pack i's area
					if(Scenario$wolf.harv ==FALSE)	{wolf.reported = 0}						# for sensitivity analysis, if we are not allowing any wolf harvest at all		
					if(wolf.reported>P.t){wolf.reported = P.t}								# if predicted harvest is greater than the pack size, harvest the whole pack
					Trap.t = wolf.reported*unreported.scalar									# multiply be unreported harvest scalar to get total wolves trapped
					T.t = round(Trap.t,0)													# no decidaml wolves trapped
					if(T.t>P.t){T.t= P.t}													# can't trap wolves that don't exist	 in pack i
					wolf.rep = round(wolf.reported, 0)										# round this number to "whole wolves" harvested
						
				#### Wolf dispersal
					disp.t = disp*ratio														# mean dispersal probability is DD for wolves
					disp.t = betaval(disp.t, disp.sd*ratio)									# generate a random betaval around this mean, with sd scaled by DD too
					D.t = (P.t + Rp.t - T.t)*(disp.t)										# wolves disperse only if they survive to recruitment and survive trapping
					D.t = round(D.t, 0)	
					if(deer.gone==T){D.t= P.t + Rp.t-T.t}									# if not enough deer to support 2 wolves, all disperse (pack breakup)
					if((P.t + Rp.t - T.t) > 20){D.t = P.t + Rp.t -T.t + D.t -20}				# any wolves over 20 individuals, add to dispersing pack members
					if(D.t<0){D.t = 0}
									
				#### Wolf chronic mortality
					cmort.t = cmort*ratio													# chronic mortality rate is DD for wolves
					cmort.t = betaval(cmort.t, cmort.sd*ratio)								# draw dd-adjusted mort rate from a beta dist, w/scaled SD
					M.t = (P.t + Rp.t - T.t - D.t)*(cmort.t)									# chronic mortality is applied to wolves not trapped and not dispersing
					M.t = round(M.t, 0)
					if(M.t < 0){M.t = 0}
					
				#### Wolf immigration
					if(P.t < 2 & CP.t/(alpha.t*(U.t+1)) <= 1/3){I.m = 2-P.t}else{I.m=0}		# Immigration can occur only if fewer than 2 wolves, and  deer for 3 wolves
					if(t > 1){I.avail = time.out$dispersers[t-1]}else{I.avail=0}				# How many potential immigrants in dispersers pool?
					if(I.avail>0 & I.m >=I.avail){I.t = I.avail}								# if space in pack greater than avail dispersers, imm =avail dispersers
					if(I.avail>0 & I.m < I.avail){I.t = I.m}									# if space in pack less than avail dispersersers, imm =space in pack
					if(I.m==0){I.t = 0}														# if no space in pack, imm= 0											
					if(I.avail == 0){I.t = 0}												# if no avail dispersers, imm=[]
					if(t==1){I.t = 0}														# no dispersal pool for t=1
					if(t>1){time.out$dispersers[t-1] = time.out$dispersers[t-1]-I.t}			# remove successful immigrants from disperser pool

					spring.wolf = P.t														# Spring wolf in year t	(includes repro, cmort, trap, imm, and emm from t-1)	
					fall.wolf = P.t + Rp.t													# Fall wolf, with young of the year, but w/out other demog. proceses yet				
					if(fall.wolf > 0){w.harv.perc.rep = wolf.rep/fall.wolf}else{				# record % of pack population harvested and reported
						w.harv.perc.rep = wolf.rep/0.00001}									# don't allow division by 0
					if(fall.wolf > 0){w.harv.perc.tot = T.t/fall.wolf}else{					# record % of pack population harvested total
						w.harv.perc.tot = T.t/0.0001}										# don't allow division by 0
					if(spring.wolf == 0){pack.extinct = 1}else{pack.extinct = 0}				# record when a pack goes extinct
					
					# Store each pack's results in rows, w/columns for each var.

					p.out = data.frame(	pack = 			i,								##### STORE RESULTS FOR EACH PACK (i) IN TIME t
										k =  			K.t,									# Deer K for pack i 
										road = 			road.t,									# road length for pack i
										spring.wolf = 	spring.wolf,							# spring wolf pack size for pack i
										fall.wolf = 		fall.wolf,							# fall wolf pack size for pack i
										wolf.harv.rep = wolf.rep,								# reported wolf harvest for pack i (NOT YET CAPPED!)
										wolf.harv.tot = T.t,									# total wolf harvest for pack i
										w.harv.perc.rep= w.harv.perc.rep	,					# percent of reported harvest out of fall pack size
										w.harv.perc.tot = w.harv.perc.tot,						# perc. of total harvest out of fall pack size
										cmort = 			M.t,								# chronic wolf mortality for pack i
										dispersers = 	D.t,									# wolves dispersing from pack i
										immigrants = 	I.t,									# wolves that immigrated to pack i
										spring.deer = 	spring.deer,							# deer spring population size in pack area i (at start of t)
										hunt.deer = 		H.t,								# deer that were hunted in pack area i
										recruit.deer = 	R.t,									# deer that recruited in pack area i
										pred.deer = 		CP.t,								# deer that were predated by wolves in pack area i
										bear.deer = 		BA.t,								# adult deer that were predated by bears in pack area i
										time = 			t,										# time (t)
										ratio = 			ratio,								# ratio (deer needed in t : deer available in (t-2)) for wolf DD
										pack.extinct =	pack.extinct,							# recond when a pack goes extinct in time t
										deer.dens = 	Deer.dt,								# Deer density per km^2
										browse =  		browse,									# record of if deer density exceeds browse impact threshold
										recover = 		recover)								# record if deer density is below browse recovery threshold
					if(i==1){pack.out=p.out}else{pack.out=rbind(pack.out, p.out)}
					

					
					}	############# End of loop through packs (i)


		disperse.new = sum(pack.out$dispersers)												# add up new dispersers from packs in year t (they've already been trapped)
		if(t==1){disperse.old = 0}else{disperse.old = time.out$dispersers[t-1]}				# dispersers left over from previous time step (will be trapped this yr)
		disp.surv = betaval(surv.disp, surv.disp.sd)											# generate survival rate for dispersers (chronic plus trapping)
		disperse.die = round(disperse.old*(1-disp.surv), 0)									# new dispersers have already survived trapping before dispersing in year t
		disperse.trap.tot = round(trap.prop.mort*disperse.die, 0)							# dispersers trapped, reported and unreported; Person & Russell 2008 p. 1545
		disperse.trap.rep = round(trap.prop.rep.mort*disperse.die, 0)						# dispersers trapped, only reported proportion; Person & Russell 2008 p. 1545
		disperse.pop = disperse.new + round(disperse.old,0)									# the pop of dispersers, for calculating proportion trapped
		disperse.pool = disperse.new + round(disperse.old, 0) - disperse.die					# the new disperser pool for this year						
		if(C.t*disperse.pool >= (sum(pack.out$spring.deer)-sum(pack.out$pred.deer))){disperse.pool = 0}				# if not enough deer to support dispersers, no dispersers
		if(disperse.pool < 0){disperse.pool = 0}

		harv.rep =  sum(pack.out$wolf.harv.rep, disperse.trap.rep)							# number of wolves harvested, reported (pack members plus dispersers) 
		if(sum(pack.out$fall.wolf, disperse.pop)==0){harv.rep=0}								# prevent NaN's when 0/0
		
																							# THIS IS WHERE THE HARVEST CAP IS IMPLEMENTED
		if(t==1){trap.trigger = harv.rep/sum(pack.out$fall.wolf, disperse.pop);				# calculate proportion of reported harvest out of fall wolf pop size, t = 1
				 untrap.val = round(trap.cap*sum(pack.out$fall.wolf))}else{					# number of wolves to be untrapped from packs if cap is exceeded, t = 1
				 	if(time.out$fall.wolf.tot[t-1]==0){										# if t >1, and wolf pop = 0...										
				 		trap.trigger = harv.rep/0.0001}else{									# trap.trigger is very large...otherwise
				 			trap.trigger = harv.rep/time.out$fall.wolf.tot[t-1]};			# calculate proportion of reported harvest out of fall wolf pop size, t > 1	
				 			untrap.val = round(trap.cap*time.out$fall.wolf.tot[t-1])}		# number of wolves that should be trapped, if cap is exceeded, t > 1
		if(trap.trigger)
		if(trap.trigger > trap.cap){															# IF reported harvest %(t) > harvest cap %...
			over.cap = 1;																	# record if cap is exceeded;
			untrap =  harv.rep -untrap.val;													# wolves to be un-trapped (pack members plus dispersers)	;
			if(untrap<0){untrap=0};															# make sure no negative numbers;
			if(disperse.trap.rep > disperse.old*trap.cap){									# IF more dispersers to be reported harvested than the trap cap %...;
				untrap.disp = disperse.trap.rep-disperse.old*trap.cap}else{					# undo disperser harvest up to the trap cap percent...;
				untrap.disp = 0};															# OTHERWISE undo none of disperser harvest;
			untrap.pack = untrap- untrap.disp;												# substract any "undone disperser harvest" from remaining reported harvest 			
			harvest.wts = pack.out$wolf.harv.rep/sum(pack.out$wolf.harv.rep);				# weights based on predicted reported harvest (uncapped) for each pack in time t
			w.untrap.pack = round(untrap.pack*harvest.wts)}else{								# undo harvest for each pack based on it's predicted reported harvest value wt.			
				over.cap = 0		;															# record if cap is not exceeded...
				untrap = 0																	# otherwise (if harest cap not exceeded)...
				untrap.disp = 0		;														# no wolf harvest will be undone for dispersers...
				w.untrap.pack = 0															# and no wolf harvest will be undone for packs
		}
																							# NOW WE WILL ACTUALLY UNDO THE HARVEST ABOVE THE HARVEST CAP
				pack.out$wolf.harv.tot = pack.out$wolf.harv.tot-	w.untrap.pack				# adjust total harvest by un-traped number of wolves																	
				pack.out$wolf.harv.rep = pack.out$wolf.harv.rep-w.untrap.pack				# un-trap the right number of pack members to meet cap
				pack.out$w.harv.perc.rep = pack.out$wolf.harv.rep/(pack.out$fall.wolf+0.0001	)# re-calc percent reported harvested
				pack.out$w.harv.perc.tot = pack.out$wolf.harv.tot/(pack.out$fall.wolf+0.0001	) # re-calc percent total harvested
				disperse.trap.rep = disperse.trap.rep-untrap.disp							# un-trap the right number of dispersers to meet cap
																							
		perc.harv.rep =  sum(pack.out$wolf.harv.rep, disperse.trap.rep)/sum(pack.out$fall.wolf, disperse.pop)		# calculate % pop harvested, reported (pack members plus dispersers) 
		if(sum(pack.out$fall.wolf, disperse.pop)==0){perc.harv.rep=0}												# prevent NaN's when 0/0
		perc.harv.tot = sum(pack.out$wolf.harv.tot, disperse.trap.tot)/sum(pack.out$fall.wolf, disperse.pop)			# % pop harvested, reported and unreported (pack plus disp.)
		if(sum(pack.out$fall.wolf, disperse.pop)==0){perc.harv.tot=0}												# prevent NaN's when 0/0
		perc.harv.pack.rep =  sum(pack.out$wolf.harv.rep)/sum(pack.out$fall.wolf)									# % pop harvested, reported (only pack members) 
		if(sum(pack.out$fall.wolf)==0){perc.harv.pack.rep=0}															# prevent NaN's when 0/0
		perc.harv.pack.tot = sum(pack.out$wolf.harv.tot)/sum(pack.out$fall.wolf)										# % pop harvested, reported and unreported (only pack members)
		if(sum(pack.out$fall.wolf)==0){perc.harv.pack.tot=0}															# prevent NaN's when 0/0
		
		if(sum(pack.out$spring.wolf, disperse.old)==0){extinct = 1}else{extinct = 0}									# total population size goes to zero?
		if(sum(pack.out$spring.wolf, disperse.old)<= qextinct){quas.extinct = 1}else{quas.extinct = 0}				# quasi-extinction threshold exceeded?				
		
		if(t==1){packs.time = pack.out}else{packs.time=abind(packs.time, pack.out, along=3)}		# store a 3-D stack of pack results through time


		pop.totals = data.frame(year = 					years[t],												##### STORE RESULTS OF ENTIRE POPULATION IN YEAR t
								k = 					sum(pack.out$k), 											# total study area deer K
								road = 					sum(pack.out$road), 											# total study area road length		
								spring.deer =			sum(pack.out$spring.deer),									# total study area deer pop size
								hunt.deer = 			sum(pack.out$hunt.deer), 									# total study area deer hunted
								spring.wolf=			sum(pack.out$spring.wolf),									# total study area spring wolves, packs only
								spring.wolf.tot = 		sum(pack.out$spring.wolf, disperse.old),						# total study area spring wolves, packs and dispersers
								fall.wolf = 			sum(pack.out$fall.wolf), 									# total study area fall wolves, packs only
								fall.wolf.tot = 		sum(pack.out$fall.wolf, disperse.pop),						# total study area fall wolves, packs and dispersers, 
								wolf.recruit = 			sum(pack.out$fall.wolf)-sum(pack.out$spring.wolf),			# total study area wolf recruits
								dispersers = 			disperse.pool, 												# total study area wolf disperser pool at end of t
								winter = 				winter,														# winter severity		
								wolf.pack.report.harv= 	sum(pack.out$wolf.harv.rep), 								# total study area reported harvest, packs only
								wolf.pack.tot.harv = 	sum(pack.out$wolf.harv.tot), 								# total study area all harvest, packs only
								wolf.pop.report.harv = 	sum(pack.out$wolf.harv.rep, disperse.trap.rep), 				# total study area reported harvest, packs and dispersers 
								wolf.pop.tot.harv  = 	sum(pack.out$wolf.harv.tot, disperse.trap.tot),				# total study area all harvest, packs and dispersers 
								perc.harv.pack.rep = 	perc.harv.pack.rep,											# percentage of pop harvested and reported (packs only) 
								perc.harv.pack.tot =  	perc.harv.pack.tot,											# percentage of pop harvested total (packs only)
								perc.harv.rep= 			perc.harv.rep,												# percentage of pop harvested and reported (packs and dispersers)
								perc.harv.total = 		perc.harv.tot,												# percentage of pop harvested total (packs and dispersers)
								cap.exceeded = 			over.cap,													# record if population-level harvest cap exceeded or not
								extinct.pop.t =			extinct,
								q.extinct.pop.t = 		quas.extinct,
								extinct.pack.perc = 	sum(pack.out$pack.extinct)/31,
								browse.pack.perc = 		sum(pack.out$browse)/31)
								

	if(t==1){time.out = pop.totals}else{time.out = rbind(time.out, pop.totals)}
	
	


	
	

	}########## End of loop through years (t)
	


if(m==1){sim.out = time.out}else{sim.out = abind(sim.out, time.out, along = 3)}

# Basic plot to track results
par(mar=c(6,4,5,5)+.1)
ymax = 60000
ymax2 = 600
col.wolf = "grey40"
col.deer = "tan3"


if(m ==1){
plot(time.out$year, time.out$spring.deer, type = "l", col=col.deer, xlab= "Years", ylab = "Deer pop size", main = "", ylim=c(0, ymax));
par(new=T);
plot(time.out$year, time.out$spring.wolf.tot, type="l", lwd=1, lty=1, col=col.wolf, xaxt="n",yaxt="n",xlab="",ylab="", ylim=c(0, ymax2));
axis(4, las=1);
abline(h=10, col="red", lty=2, lwd=2);
mtext("Wolf pop size", side=4, line=3);
legend("topright", col=c(col.deer, col.wolf), lty=c(1,1), lwd=1, legend=c("Deer pop","Wolf pop"), xjust=0.5, cex=1, bty="n");
}else{
	par(new=T);
	plot(time.out$year, time.out$spring.wolf.tot, type="l", lwd=1, lty=1, col=col.wolf, xaxt="n", yaxt="n", xlab="",ylab="", ylim=c(0, ymax2), frame.plot=F);
	abline(h=10, col="red", lty=2, lwd=2);
	par(new=T);
	plot(time.out$year, time.out$spring.deer, type = "l", col=col.deer, xaxt="n", yaxt="n", xlab = "", ylab = "", ylim=c(0, ymax), frame.plot=F)
	
}

print(m)

} ############# end of loop through Monte Carlo itterations (m)

dev.off()

##### now calculate metrics to save as outputs

m.CI = function(data, col.name, time.name, alph, dim){			# a function to empirically derive means, medians, and CIs. dim = which dim is time?, 

for(t in 1:dim(data)[dim]){
	yr = data[t,which(colnames(data[,,1])==time.name),1]
	s = sort(data[t,which(colnames(data[,,1])==col.name),])
	count = dim(data)[3]
	lcl.ind = s[floor(count*alph)]
	ucl.ind = s[ceiling(count-count*alph)]
	mn = mean(s)
	med = median(s)
	out.dat = data.frame(year = yr, mean = mn, median = med, lcl = lcl.ind, ucl = ucl.ind)
	if(t==1){out.data = out.dat}else{out.data= rbind(out.data, out.dat)}
}
return(out.data)
}

deer.mn 					<- m.CI(sim.out, "spring.deer", "year", 0.05, 1)
wolf.spring.mn 			<- m.CI(sim.out, "spring.wolf.tot", "year", 0.05, 1)
wolf.fall.mn 			<- m.CI(sim.out, "fall.wolf.tot", "year", 0.05, 1)
hunt.mn 					<- m.CI(sim.out, "hunt.deer", "year", 0.05, 1)
perc.harv.rep.mn 		<- m.CI(sim.out, "perc.harv.rep", "year", 0.05, 1)
perc.harv.tot.mn 		<- m.CI(sim.out, "perc.harv.total", "year", 0.05, 1)
cap.trap.triggered.mn	<- m.CI(sim.out, "cap.exceeded", "year", 0.05, 1)
extinct.pop.t.mn 		<- m.CI(sim.out, "extinct.pop.t", "year", 0.05, 1)
q.extinct.pop.t.mn 		<- m.CI(sim.out, "q.extinct.pop.t", "year", 0.05, 1)
extinct.pack.perc.mn 	<- m.CI(sim.out, "extinct.pack.perc", "year", 0.05, 1)
browse.pack.perc.mn 	<- m.CI(sim.out, "browse.pack.perc", "year", 0.05, 1)

wolf.spring.monte= sim.out[,which(names(time.out)=="spring.wolf.tot"),]
colnames(wolf.spring.monte) = paste("wolf", seq(1:monte.length), sep="_")

for(m in 1:monte.length){
	e.dat = data.frame(sim.out[,,m])
	if(1%in%e.dat$extinct.pop.t){e.monte=1}else{e.monte=0}
	if(1%in%e.dat$q.extinct.pop.t){q.e.monte=1}else{q.e.monte=0}
	e.monte.dat = data.frame(e.monte= e.monte, q.e.monte = q.e.monte)
	if(m==1){extinct.record = e.monte.dat}else{extinct.record = rbind(extinct.record, e.monte.dat)}
}
extinct.record=(colSums(extinct.record))/monte.length
extinct.record=data.frame(rate.extinct = extinct.record[1], rate.quasi.extinct = extinct.record[2])

write.table(deer.mn, paste(DirOut,Which.Scen,"deer.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(wolf.spring.mn, paste(DirOut,Which.Scen,"wolf.spring.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(wolf.fall.mn, paste(DirOut,Which.Scen,"wolf.fall.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(hunt.mn, paste(DirOut,Which.Scen,"hunt.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(perc.harv.rep.mn, paste(DirOut,Which.Scen,"perc.harv.rep.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(perc.harv.tot.mn, paste(DirOut,Which.Scen,"perc.harv.tot.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(cap.trap.triggered.mn, paste(DirOut,Which.Scen,"cap.trap.triggered.csv", sep="/"), sep=",", row.names=FALSE)
write.table(extinct.pop.t.mn, paste(DirOut,Which.Scen,"extinct.pop.t.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(q.extinct.pop.t.mn, paste(DirOut,Which.Scen,"q.extinct.pop.t.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(extinct.pack.perc.mn, paste(DirOut,Which.Scen,"extinct.pack.perc.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(browse.pack.perc.mn, paste(DirOut,Which.Scen,"browse.pack.perc.mn.csv", sep="/"), sep=",", row.names=FALSE)
write.table(extinct.record, paste(DirOut,Which.Scen,"extinct.rate.monte.csv", sep="/"), sep=",", row.names=FALSE)
write.table(wolf.spring.monte, paste(DirOut, Which.Scen, "wolf.spring.monte.csv", sep="/"), sep=",", row.names=FALSE)



} # end of F loop through Scenarios or Sensitivities

# end code

