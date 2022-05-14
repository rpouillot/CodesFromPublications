####
# R Code for the article
#Risk Anal. 2012 Apr 10. doi: 10.1111/j.1539-6924.2012.01796.x.
#A Risk Assessment of Campylobacteriosis and Salmonellosis Linked to Chicken Meals Prepared in Households in Dakar, Senegal.
#Pouillot R, Garin B, Ravaonindrina N, Diop K, Ratsitorahina M, Ramanantsoa D, Rocourt J.
# By R.Pouillot

#Updated 3/12/2014 for new versions of Poisson and Binom
#Change Your Working Directory


library(mc2d)
ndvar(1001) # Choose the number of iterations in the variability dimension 
ndunc(101)  # Choose the number of iterations in the uncertainty dimension 


#Choose here if you want the "Campy" or "Salmo" model
Bacteria <- "Campy"
#Bacteria <- "Salmo"

#Choose the mitigation strategy you like (baseline = 0)
Mitigation <- 0

#rm(list=ls(all=TRUE))
#memory.limit(7900)
gc()
set.seed(666)
cat(Bacteria, "Mitigation:",Mitigation,"\n")

## FUNCTIONS


Cardinal <- function(T, muopt=0.7320, Tmin=5.699, Topt=40.01, Tmax=49.26){
	D <- (T-Tmax)*(T-Tmin)^2
	E <- (Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2*T))
	EGR <- ifelse(T < Tmin | T > Tmax, 0, 	muopt * D / E)
	return(EGR)
	}

# PROTECTED VERSIONS OF RPOIS AND RBINOM
rbinomProt <- function(n, size, prob){
  if(any(is.na(size)) | any(is.na(prob)) | any(prob > 1) | any(prob < 0) | any(size < 0)) stop("bad argument in rbinomProt")
  tmp <- rbinom(n,size,prob)
  quel <- is.na(tmp)
  if(any(quel)){ 
    size <- rep(size,length.out=n) #To deal with size or prob of size < element of quel
    prob <- rep(prob,length.out=n)
    tmp[quel] <- pmin(pmax(
      round(rnorm(sum(quel),size[quel]*prob[quel],
                  sqrt(size[quel]*prob[quel]*(1-prob[quel])))),
      0,size[quel]))
  }
  return(tmp)
}

rpoisProt <- function(n, lambda){
  if(any(is.na(lambda)) | any(lambda < 0)) stop("bad argument in rpoisProt")
  tmp <- rpois(n,lambda)
  quel <- is.na(tmp)
  if(any(quel)) {
    lambda <- rep(lambda, length.out=n)
    tmp[quel] <- pmin(round(rnorm(sum(quel),lambda[quel],sqrt(lambda[quel]))),0)
  }
  return(tmp)
}

# Temperature Profiles
tempProf <- as.matrix(read.csv("AllSN.csv"))
tempProf <- tempProf[-1, tempProf[1,]==1]
StartPre <- tempProf[1,]
StartPost <- tempProf[2,]
summary(StartPost-StartPre)
tempProf <- tempProf[-c(1,2),]

nProf <- ncol(tempProf)
ProfilesPost <- ProfilesPre <- NULL
for(i in 1:nProf){
  l <- max(which(!is.na(tempProf[,i])))
  ProfilesPost <- c(ProfilesPost, list(tempProf[StartPost[i]:l,i]))
  ProfilesPre <- c(ProfilesPre,   list(tempProf[1:StartPre[i],i]))
}

#####
#####
# SPECIFICITIES CAMPY
#####
#####
if(Bacteria == "Campy"){
	T0boot <- dget("f3BootCampy")
	T0boot <- T0boot[rep(1:nrow(T0boot),length=ndunc()),]
	Prev <- mcdata(1,"0")
	GrowthPre <- GrowthPost <- mcdata(0,"0")
	
	DRpar1 <- read.delim("campylobacter-updated.txt",sep="\t", header=FALSE)
	DRpar2 <- read.delim("eta-r-black.txt",sep="\t", header=FALSE)
	DRpar3 <- read.delim("eta-r-outbreak.txt",sep="\t", header=FALSE)
	DRpar <- cbind(DRpar1,DRpar2,DRpar3)
	colnames(DRpar) <- c("alpha","beta","etaAdults","rAdults","etaKids","rKids")
	DRpar <- DRpar[sample(nrow(DRpar)),]											# Shuffle
	
	DRpar <- DRpar[rep(1:nrow(DRpar),length=ndunc()),]								# Take the first ndunc()
	alpha <- mcdata(DRpar$alpha,"U")
	beta <- mcdata(DRpar$beta,"U")

	rAdults <-  mcdata(DRpar$rAdults,"U")
	etaAdults <-  mcdata(DRpar$etaAdults,"U")
	powerAdults <- -1 # decrease with dose

	rKids <-  mcdata(DRpar$rKids,"U")
	etaKids <- mcdata(DRpar$etaKids,"U")
	powerKids <- 1 # increase with dose

} 

#####
#####
# SPECIFICITIES SALMO
#####
#####


if(Bacteria == "Salmo"){
	T0boot <- dget("f3BootSalm")
	T0boot[,"prev"] <- exp(T0boot[,"prev"]) / (1+exp(T0boot[,"prev"]))
	T0boot <- T0boot[rep(1:nrow(T0boot),length=ndunc()),]
	
	Prev <- mcdata(T0boot[,"prev"],"U")
	
	#Growth model
	#curve(Cardinal(x), from=5, to=55)

	muopt <- mcstoc(rnorm,"U",mean=0.7320,sd=(0.7498-.7143)/3.92) # From CI Oscar 2002
	Tmin  <- mcstoc(rnorm,"U",mean=5.699,sd=(7.308-4.090)/3.92)   # From CI Oscar 2002
	Tmax  <- mcstoc(rnorm,"U",mean=49.26,sd=(49.64-48.89)/3.92)   # From CI Oscar 2002
	Topt  <- mcstoc(rnorm,"U",mean=40.01,sd=(40.51-39.51)/3.92)   # From CI Oscar 2002

	GrowthPre <- GrowthPost <- array(NA, dim=c(ndvar(),ndunc(),1))
	rProf <-  sample(1:nProf, ndvar(), replace=TRUE)	# Pick one profiles by Simu

	for(i in 1:ndunc()){
		basePre <- sapply(1:nProf, function(x) sum(Cardinal(ProfilesPre[[x]],muopt[i],Tmin[i],Topt[i],Tmax[i]) * (3/60)))
		GrowthPre[,i,1] <- basePre[rProf]
		basePost <- sapply(1:nProf, function(x) sum(Cardinal(ProfilesPost[[x]],muopt[i],Tmin[i],Topt[i],Tmax[i]) * (3/60)))
		GrowthPost[,i,1] <- basePost[rProf]
	}	
	GrowthPost <- mcdata(GrowthPost,"VU")
	GrowthPre <- mcdata(GrowthPre,"VU")
	
	if(Mitigation == 7) GrowthPre <- GrowthPost <- mcdata(0,"0")
	
	DRpar <- read.delim("parameters salmo.txt",sep="\t", header=TRUE)
	DRpar <- DRpar[sample(nrow(DRpar)),]											# Shuffle
	DRpar <- DRpar[rep(1:nrow(DRpar),length=ndunc()),]								# Take the first ndunc()
	alpha <- mcdata(DRpar$alpha,"U")
	beta <- mcdata(DRpar$beta,"U")

	rAdults <- rKids <- mcdata(DRpar$r,"U")
	etaAdults <- etaKids <- mcdata(DRpar$eta,"U")
	powerAdults <- powerKids <- 1 # Increase with dose

} 

llMax <- 8
lMax <- 10^llMax

if(Mitigation == 1) Prev <- Prev / 2 	
prevT0 <- mcstoc(rbern,"VU",Prev)
meanlevelT0 <- mcdata(T0boot[,"mean"], type="U")

if(Mitigation == 2) meanlevelT0 <- meanlevelT0-2 
if(Mitigation == 3) meanlevelT0 <- meanlevelT0-5 

sdlevelT0   <- mcdata(T0boot[,"sd"], type="U")
levelT0Pos     <- 10^mcstoc(rnorm, type="VU",mean=meanlevelT0, sd=sdlevelT0)
levelT0 <- levelT0Pos * prevT0



# Transfer Coefficients

TChicRTE <- pmin(10^mcstoc(rnorm,"V",mean=-1.72,sd=1.07),1)	# Hand Vegetable!
TRTEChic <- pmin(10^mcstoc(rnorm,"V",mean=-1.72,sd=1.07),1)	# Hand Vegetable!

TChicHand <- pmin(10^mcstoc(rnorm,"V",mean=-1.69,sd=0.81),1)	# From Hoeltzer 
THandChic <- pmin(10^mcstoc(rnorm,"V",mean=-4.96,sd=0.37),1) 	# From Hoeltzer

TRTEHand <- pmin(10^mcstoc(rnorm,"V",mean=-1.72,sd=1.07),1)
THandRTE <- pmin(10^mcstoc(rnorm,"V",mean=-1.72,sd=1.07),1)

TChicBoard <- pmin(10^mcstoc(rnorm,"V",mean=-1.45,sd=1.39),1)	
TBoardChic <- pmin(10^mcstoc(rnorm,"V",mean=-0.15,sd=0.07),1)	

TRTEBoard <- pmin(10^mcstoc(rnorm,"V",mean=-1.42,sd=0.52),1)	
TBoardRTE <- pmin(10^mcstoc(rnorm,"V",mean=-1.42,sd=0.52),1)	

TRTEKnife <- pmin(10^mcstoc(rnorm,"V",mean=-2.43,sd=0.69),1)	
TKnifeRTE <-  pmin(10^mcstoc(rnorm,"V",mean=-2.43,sd=0.69),1)	
TChicKnife <- pmin(10^mcstoc(rnorm,"V",mean=-2.43,sd=0.69),1)	
TKnifeChic <-   pmin(10^mcstoc(rnorm,"V",mean=-2.43,sd=0.69),1)	

THandKnife <- TKnifeHand <- 0	# 0 because touch the HANDLE

TDishChic <- pmin(10^mcstoc(rnorm,"V",mean=-2.70,sd=0.45),1)	
TChicDish <- pmin(10^mcstoc(rnorm,"V",mean=-2.70,sd=0.45),1)	

WashBoard <-  10^mcstoc(rpert, type="V", min=-7, mode=-4.5, max=-1)
RinseBoard <- 10^mcstoc(rpert, type="V", min=-1.5, mode=-0.5, max=0)

if(Mitigation == 4) RinseBoard <- WashBoard

WashDish <-  10^mcstoc(rpert, type="V", min=-7, mode=-4.5, max=-1)

WashHand <-  mcstoc(rbeta, type="V", shape1=0.24, shape2=6.67)
RinseHand <- 10^mcstoc(rpert, type="V", min=-1.5, mode=-0.5, max=0)

if(Mitigation == 4) RinseHand <- WashHand


Hand <- RTE <- Board <- Knife <- Dish <- 0

############################################################################# Stage 0
# Chicken neck size: 10 cm long * 1.5 diameter. Half for campy. 
surfNeck <- round(1.5 * pi * 10 / 2, 1)
surfChicken <- mcstoc(rnorm, type="V", mean=1232, sd=165)

conc <- levelT0 / surfNeck * surfChicken
Chicken <-  mcstoc(rpoisProt, type="VU", conc)
Chicken <- pmin(Chicken, lMax * surfChicken)  

ChickT0 <- Chicken

############################################################################# Stage 1
# Growth Pre cook
Chicken <- round(10^(log10(Chicken/surfChicken) + GrowthPre) * surfChicken)  
Chicken <- ChickT1 <- round(pmin(Chicken,  lMax * surfChicken))

ChickT1 <- Chicken 
HandT1 <- RTET1 <- BoardT1 <- KnifeT1 <- DishT1 <- mcdata(0,"0")

############################################################################# Stage 2
# Chicken RTE

pChickRTE <- mcstoc(rbeta, type="U", shape1= 2+1, shape2= 72-2+1)
PChickRTE <- mcstoc(rbern, type="VU", p = pChickRTE)

NChicken <- mcstoc(rbinomProt, type="VU", size=Chicken, prob=1-TChicRTE) # ASSUME NO BACTERIA ON RTE, Transfer 2 objects
NRTE <- Chicken - NChicken

Chicken <- mcprobtree(PChickRTE, list("0" = Chicken, "1" = NChicken),type="VU") 
RTE    <- NRTE * PChickRTE # 0 if no contact


# Chicken Cooks Hand
pChickHand <- mcstoc(rbeta, type="U", shape1= 72+1, shape2= 1)
PChickHand <- mcstoc(rbern, type="VU", p = pChickHand)

NChicken <- mcstoc(rbinomProt, type="VU", size=Chicken, prob=1-TChicHand) # ASSUME NO BACTERIA ON HAND
NHand <- Chicken - NChicken

Chicken <- mcprobtree(PChickHand, list("0" = Chicken, "1" = NChicken),type="VU") 
Hand    <- NHand * PChickHand # 0 if no contact

# Chicken, Hand, Board and Knife
pCutBefore <- mcstoc(rbeta, type="U", shape1=55+1, shape2=72-55+1)
PCutBefore <- mcstoc(rbern, type="VU", p = pCutBefore)

#Chicken, Hands, Board, Knife 

temp <- Chicken	# To Draw the Multinomial
fromChickenC <- mcstoc(rbinomProt, type="VU", size=temp, prob=(1-TChicHand) * (1-TChicBoard) * (1-TChicKnife)) # Stay on the Chicken
temp <- temp - fromChickenC
fromChickenH <- mcstoc(rbinomProt, type="VU", size=temp, prob= TChicHand / (TChicHand + TChicBoard + TChicKnife))
temp <- temp - fromChickenH
fromChickenB <- mcstoc(rbinomProt, type="VU", size=temp, prob= TChicBoard / (TChicBoard+TChicKnife))
fromChickenK <-  temp - fromChickenB

temp <- Hand
fromHandH <- mcstoc(rbinomProt, type="VU", size=temp, prob=1-THandChic)
fromHandC <- temp - fromHandH

NChicken <- fromChickenC + fromHandC
NHand <- fromChickenH + fromHandH
NBoard <- fromChickenB
NKnife <- fromChickenK

Chicken <- mcprobtree(PCutBefore, list("0" = Chicken, "1" = NChicken),type="VU") 
Hand    <- mcprobtree(PCutBefore, list("0" = Hand, "1" = NHand),type="VU") 
Board    <- NBoard * PCutBefore 	# ASSUME NO BACTERIA ON BOARD
Knife <- NKnife * PCutBefore		# ASSUME NO BACTERIA ON KNIFE
 
# Chicken in Dish

Dish <-  mcstoc(rbinomProt, type="VU", size=Chicken, prob=TChicDish)
Chicken <-  Chicken - Dish 	# ASSUME NO BACTERIA ON DISH

HandT2 <- Hand
RTET2 <- RTE
BoardT2 <- Board
KnifeT2 <- Knife
DishT2 <- Dish

############################################################################# Stage 3


# Wash Board and Knife
pWashBoard <- mcstoc(rdirichlet, type="U", alpha= c(9+1,41+1,5+1),nvariates=3)
if(Mitigation == 5) pWashBoard <- c(0,0,1) 
PWashBoard <- mcstoc(rempiricalD, type="VU", values=0:2, prob=pWashBoard)

Persistence <- mcprobtree(PWashBoard, type="VU", list("0"= mcdata(1,"VU"),"1"= RinseBoard, "2"=WashBoard))

Hand  <- mcstoc(rbinomProt, type="VU", size=Hand, prob=Persistence)
Board <- mcstoc(rbinomProt, type="VU", size=Board, prob=Persistence)
Knife <- mcstoc(rbinomProt, type="VU", size=Knife, prob=Persistence)

# Wash Dish
pWashDish <- mcstoc(rdirichlet, type="U", alpha= c(13+1,9+1,48+1),nvariates=3)	# 1: Do not Wash, 2:Wash, 3:Another Dish
if(Mitigation == 5) pWashDish[,,1] <- 0
PWashDish <- mcstoc(rempiricalD, type="VU", values=0:2, prob=pWashDish)

Persistence <- mcprobtree(PWashDish, type="VU", list("0"= mcdata(1,"VU"),"1"= WashDish, "2"=mcdata(0,"VU")))
Dish <- mcstoc(rbinomProt, type="VU", size=Dish, prob=Persistence)
Hand  <- mcstoc(rbinomProt, type="VU", size=Hand, prob=Persistence)

# Wash Hand
pWashHand <- mcstoc(rdirichlet, type="U", alpha= c(44+1,18+1,10+1),nvariates=3)
if(Mitigation == 5 | Mitigation == 6) pWashHand <- c(0,0,1)
PWashHand <- mcstoc(rempiricalD, type="VU", values=0:2, prob=pWashHand)

Persistence <- mcprobtree(PWashHand, type="VU", list("0"= mcdata(1,"VU"),"1"= RinseHand , "2"=WashHand))

Hand <- mcstoc(rbinomProt, type="VU", size=Hand, prob=Persistence)

if(Mitigation == 6) Board <- Knife  <- Dish <- mcdata(0,"0")


# Cooking
ChickT2 <- Chicken
Chicken <- ChickT3 <- 0 * Chicken

HandT3 <- Hand
RTET3 <- RTE
BoardT3 <- Board
KnifeT3 <- Knife
DishT3 <- Dish

############################################################################# Stage 4

# Touch RTE

pHandRTE <- mcstoc(rbeta, type="U", shape1= 36+1, shape2= 44-36+1)
PHandRTE <- mcstoc(rbern, type="VU", p = pHandRTE)

NHand <- mcstoc(rbinomProt, type="VU", size=Hand, prob=1-THandRTE) + mcstoc(rbinomProt, type="VU", size=RTE, prob=TRTEHand)
NRTE <-  Hand + RTE - NHand

pIsRTE <- mcstoc(rbeta,"U",shape1=44+1,shape2=72-44+1)
isRTE <- mcstoc(rbern,"VU",p=pIsRTE)

Hand <- mcprobtree(PHandRTE*isRTE, list("0" = Hand, "1" = NHand),type="VU") 
RTE <- mcprobtree(PHandRTE*isRTE, list("0" = RTE, "1" = NRTE),type="VU") 


# Back on the Board
pBackBoard <- mcstoc(rbeta, type="U", shape1=12+1, shape2=72-12+1)
PBackBoard <- mcstoc(rbern, type="VU", p = pBackBoard)

NChicken <- mcstoc(rbinomProt, type="VU", size=Board, prob=TBoardChic) 
NBoard <-  Board - NChicken

Chicken <- mcprobtree(PBackBoard, list("0" = Chicken, "1" = NChicken), type="VU") 
Board <- mcprobtree(PBackBoard, list("0" = Board, "1" = NBoard), type="VU") 

# Chicken in Dish

All <- Chicken + Dish
Chicken <- mcstoc(rbinomProt, type="VU", size=Dish, prob=TDishChic) +   mcstoc(rbinomProt, type="VU", size=Chicken, prob=TChicDish)
Dish <-  All - Chicken


#RTE, Hands, and Board 
pCutRTE <- mcstoc(rbeta, type="U", shape1= 14+1, shape2= 72-14+1)
PCutRTE <- mcstoc(rbern, type="VU", p = pCutRTE) * isRTE

temp <- RTE	
fromRTER <- mcstoc(rbinomProt, type="VU", size=temp, prob=(1-TRTEHand)*(1-TRTEBoard))	# Stay on the RTE
temp <- temp - fromRTER
fromRTEH <- mcstoc(rbinomProt, type="VU", size=temp, prob= TRTEHand / (TRTEHand+TRTEBoard))
fromRTEB <- temp - fromRTEH

temp <- Hand
fromHandH <- mcstoc(rbinomProt, type="VU", size=temp, prob=1-THandRTE)
fromHandR <- temp - fromHandH

temp <- Board
fromBoardB <- mcstoc(rbinomProt, type="VU", size=temp, prob=1-TBoardRTE)
fromBoardR <- temp - fromBoardB

temp <- Knife
fromKnifeK <- mcstoc(rbinomProt, type="VU", size=temp, prob=1-TKnifeRTE)
fromKnifeR <- temp - fromKnifeK


NRTE <- fromRTER + fromHandR + fromBoardR
NHand <- fromRTEH + fromHandH
NBoard <- fromRTEB + fromBoardB 

RTE <- mcprobtree(PCutRTE, list("0" = RTE, "1" = NRTE),type="VU") 
Hand    <- mcprobtree(PCutRTE, list("0" = Hand, "1" = NHand),type="VU") 
Board    <-  mcprobtree(PCutRTE, list("0" = Board, "1" = NBoard),type="VU") 


ChickT4 <- Chicken
HandT4 <- Hand
RTET4 <- RTE
BoardT4 <- Board
KnifeT4 <- Knife
DishT4 <- Dish

############################################################################# Stage 5
# Growth

Chicken <- round(10^(log10(Chicken/surfChicken) + GrowthPost) * surfChicken)  
Chicken <- ChickT5 <-round(pmin(Chicken,  lMax * surfChicken))

RTET5 <- RTE

####
####
# CONSO ET DR
###
###

pNbEaters <- mcstoc(rdirichlet,"U",alpha=1+c(0,0,1,10,12,6,9,13,10,1,4,4,1,1),nvariates=14)
nbEaters <- mcstoc(rempiricalD,"VU",values=c(1:14),prob=pNbEaters)

Dose <- mcstoc(rbinomProt, type= "VU", size = isRTE*RTE + Chicken, prob= 1/nbEaters)
pinf <- mcstoc(rbeta,"VU", shape1=alpha, shape2=beta)
PInf <- 1-(1-pinf)^Dose


tauKids <- mcstoc(rgamma, "VU", shape = rKids, scale = ifelse(Dose == 0, 0, etaKids*(Dose)^powerKids))	# if Dose = 0, Pinf = 0
pIllgInfKids <- 1-exp(-tauKids)


tauAdults <- mcstoc(rgamma, "VU", shape = rAdults, scale = ifelse(Dose == 0, 0, etaAdults*(Dose)^powerAdults))	# if Dose = 0, Pinf = 0
pIllgInfAdults <- 1-exp(-tauAdults)

# Alt beta Poisson

PIllKids <- PInf * pIllgInfKids					# For Salmo
PIllAdults <- PInf * pIllgInfAdults 

Risk <- mc(ChickT0, Chicken0Pos=ChickT0>0,Chicken,ChickenPos=Chicken>0,RTE,RTEPos=RTE>0,Dose,DosePos=Dose>0,PInf,PIllKids)
if(Bacteria == "Campy") Risk <- mc(Risk,PIllAdults) 


print(sR <- summary(Risk))
#dput(PIllKids,paste("Risk",Bacteria,"Mitigation",Mitigation,sep="-"))
plot(Risk)

