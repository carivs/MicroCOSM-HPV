#Plotting Cancer Incidence:
suppressPackageStartupMessages({
	library(tidyverse)
	library(openxlsx)
})  
options(scipen = 100)
setwd("~/cplusplus/unchanged/MicroCOSMHPV")
#### Fixed inputs ####
start_year <- 2001

# Run script called Thembisa Populations.R
#source("~/cplusplus/Thembisa Populations.R")
THEM_POP <- FemalePop  %>% select("2001":"2120") 
THEM_ART <- FemaleART %>% select("2001":"2120") 
THEM_NOART <- FemaleNOART %>% select("2001":"2120")
THEM_NEG <- FemaleHIVneg %>% select("2001":"2120") 
WORLD3 <- c(8895,8508,8082,7850,7974,8191,7444,6756,6565,6198,5510,4701,4115,3092,2249,1763,1154,954)

#### Read in the model output ####
setwd("~/cplusplus/")
path <- paste0("/ClusterOutput01SEPT/output/")
path1 <- paste0("output/")
n <- 24
yrs <- 2120-start_year+1
ALL_SEGI=matrix(0, nrow=n, ncol=yrs)
Types_SEGI=matrix(0, nrow=n, ncol=yrs)
ALL_SEGI2=matrix(0, nrow=n, ncol=yrs)
Types_SEGI2=matrix(0, nrow=n, ncol=yrs)
new_cases=matrix(0, nrow=n, ncol=yrs)

ii=1
for(ii in c(1:n)){
	# import incident cancer
	neg=pos=art=0
	for(jj in 1:13){
		neg <- neg+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[1:18,(start_year-1985+2):137]
		pos <- pos+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[37:54,(start_year-1985+2):137]
		art <- art+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[19:36,(start_year-1985+2):137]
	}  

	# import HPV 16 and 18 incident cancer
	neg1618=pos1618=art1618=0
	for(jj in 1:2){
		neg1618 <- neg1618+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[1:18,(start_year-1985+2):137]
		pos1618 <- pos1618+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[37:54,(start_year-1985+2):137]
		art1618 <- art1618+read.table(paste0(path, ii, "_",jj-1,"NewCC.txt"))[19:36,(start_year-1985+2):137]
	}  
	
	#import MicroCOSM population
	MicroCOSM_POP <- read.table(paste0(path, ii, "_PopulationPyramid.txt"), nrows = 162)
	MicroCOSM_POP$num <- rep(1:9, each=18)
	MicroCOSM_POPL <- split(MicroCOSM_POP, MicroCOSM_POP$num)
	negPOP <- MicroCOSM_POPL[[1]][,(start_year-1985+1):136]
	ARTPOP <- MicroCOSM_POPL[[2]][,(start_year-1985+1):136]
	noARTPOP <- MicroCOSM_POPL[[3]][,(start_year-1985+1):136]
	
	# calculate age-standardised incidence, weighted by Thembisa pop
	x <- THEM_NEG*(neg/negPOP)
	y <- THEM_NOART*(pos/noARTPOP)
	z <- THEM_ART*(art/ARTPOP)
	x[is.na(x)] <- 0
	y[is.na(y)] <- 0
	z[is.na(z)] <- 0
	x[is.infinite(as.matrix(x))] <- 0
	y[is.infinite(as.matrix(y))] <- 0
	z[is.infinite(as.matrix(z))] <- 0
	ALL_SEGI[ii,] <-  apply(WORLD3*(x+y+z)/THEM_POP, 2, sum)
	
	# calculate age-standardised incidence for HPV 16 & 18, weighted by Thembisa pop
	x1 <- THEM_NEG*(neg1618/negPOP)
	y1<- THEM_NOART*(pos1618/noARTPOP)
	z1 <- THEM_ART*(art1618/ARTPOP)
	x1[is.na(x1)] <- 0
	y1[is.na(y1)] <- 0
	z1[is.na(z1)] <- 0
	x1[is.infinite(as.matrix(x1))] <- 0
	y1[is.infinite(as.matrix(y1))] <- 0
	z1[is.infinite(as.matrix(z1))] <- 0
	Types_SEGI[ii,] <-  apply(WORLD3*(x1+y1+z1)/THEM_POP, 2, sum)
	
	neg1=pos1=art1=0
	for(jj in 1:13){
		neg1 <- neg1+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[1:18,(start_year-1985+2):137]
		pos1 <- pos1+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[37:54,(start_year-1985+2):137]
		art1 <- art1+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[19:36,(start_year-1985+2):137]
	}  
	
	# import HPV 16 and 18 incident cancer
	neg16181=pos16181=art16181=0
	for(jj in 1:2){
		neg16181 <- neg16181+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[1:18,(start_year-1985+2):137]
		pos16181 <- pos16181+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[37:54,(start_year-1985+2):137]
		art16181 <- art16181+read.table(paste0(path1, ii, "_",jj-1,"NewCC.txt"))[19:36,(start_year-1985+2):137]
	}  
	
	#import MicroCOSM population
	MicroCOSM_POP1 <- read.table(paste0(path1, ii, "_PopulationPyramid.txt"), nrows = 162)
	MicroCOSM_POP1$num <- rep(1:9, each=18)
	MicroCOSM_POP1L <- split(MicroCOSM_POP1, MicroCOSM_POP1$num)
	negPOP1 <- MicroCOSM_POP1L[[1]][,(start_year-1985+1):136]
	ARTPOP1 <- MicroCOSM_POP1L[[2]][,(start_year-1985+1):136]
	noARTPOP1 <- MicroCOSM_POP1L[[3]][,(start_year-1985+1):136]
	
	
	# calculate age-standardised incidence, weighted by Thembisa pop
	x2 <- THEM_NEG*(neg1/negPOP1)
	y2 <- THEM_NOART*(pos1/noARTPOP1)
	z2 <- THEM_ART*(art1/ARTPOP1)
	x2[is.na(x2)] <- 0
	y2[is.na(y2)] <- 0
	z2[is.na(z2)] <- 0
	x2[is.infinite(as.matrix(x2))] <- 0
	y2[is.infinite(as.matrix(y2))] <- 0
	z2[is.infinite(as.matrix(z2))] <- 0
	ALL_SEGI2[ii,] <-  apply(WORLD3*(x2+y2+z2)/THEM_POP, 2, sum)
	
	# calculate age-standardised incidence for HPV 16 & 18, weighted by Thembisa pop
	x21 <- THEM_NEG*(neg16181/negPOP1)
	y21<- THEM_NOART*(pos16181/noARTPOP1)
	z21 <- THEM_ART*(art16181/ARTPOP1)
	x21[is.na(x21)] <- 0
	y21[is.na(y21)] <- 0
	z21[is.na(z21)] <- 0
	x21[is.infinite(as.matrix(x21))] <- 0
	y21[is.infinite(as.matrix(y21))] <- 0
	z21[is.infinite(as.matrix(z21))] <- 0
	Types_SEGI2[ii,] <-  apply(WORLD3*(x21+y21+z21)/THEM_POP, 2, sum)
	
	print(ii)
}  

plot.ASIR <- data.frame(
	year = rep(2001:2120, 4),
	ASIR = c(apply(ALL_SEGI, 2, median),
					 apply(Types_SEGI, 2, median),
					 apply(ALL_SEGI2, 2, median),
					 apply(Types_SEGI2, 2, median)), 
	Scenario = rep(c("TxV administered to those on ART (All)", "TxV administered to those on ART (Types16 and 18) ", "Status Quo (All)", "Status Quo (Types 16 and 18)"), each = 120)
)

plot.ASIR %>% 
	ggplot() +
	geom_line(aes(x = year, y = ASIR, color = Scenario), linewidth = 1) +
	geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
	geom_hline(yintercept = 4, linetype = "dashed", color = "red") +
	ylab("ASIR per 100,000 women") +
	scale_y_continuous(breaks = seq(0, 55, 5), limits = c(0, 55)) + 
	scale_x_continuous(breaks = seq(2020, 2120, 20)) +
	theme_bw() +
	theme(
		axis.title.x = element_blank(),
		legend.title = element_text(face = "bold", size = 12),
		legend.text = element_text(size = 10),
		legend.position = c(0.8, 0.8)
	)

