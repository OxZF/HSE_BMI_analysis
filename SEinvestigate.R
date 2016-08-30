# --------------------------------------
# DPhil: working on standard errors
# SEinvestigate.R
# 30/08/16
# Zoe Fannon, University of Oxford
# --------------------------------------

# SET UP
library(apc)
setwd("C:/Users/me/Documents/DPhil/HSE_BMI_analysis/Functions/")
source("write.functions.v2.R")                                     # gets functions required for apc on indiv data
setwd("C:/Users/me/Documents/R/APCmphil.thesis")
test.data <- read.csv("artificial_individ_AP.csv")                 # has three ages, three cohorts, five periods

# Refresh knowledge of indiv apc commands
test.d.collin <- apc.get.design.indiv.collinear(test.data, datatype = "AP")
test.d.mod    <- apc.get.design.indiv.model(test.d.collin, dep.var= "continuous")
test.fit      <- apc.fit.indiv(test.d.mod, model.family="gaussian")
var.apc.plot.fit(test.fit, main.outer="test.ssdd", cex=1.5, cex.axis=.9, type="sum.sum",
                 mgp=c(1.8, .4, 0), theight = 1.3, cex.main = 2, mar=c(3.4, 3, 3.4, 1), oma=c(2, 2, 4, 2))   # Specs stolen from thesis_example_results
var.apc.plot.fit(test.fit, main.outer="test.detrend", cex=1.5, cex.axis=.9, type="detrend",
                 mgp=c(1.8, .4, 0), theight = 1.3, cex.main = 2, mar=c(3.4, 3, 3.4, 1), oma=c(2, 2, 4, 2))


# Look more closely at some of the elements
View(test.d.collin$structure.design.collinear)
View(test.fit$coefficients.canonical)


#### Steal parts from apc.identify to understand them ####
apc.fit.model <- test.fit
##############################
#	get values
coefficients	<- apc.fit.model$coefficients.canonical	#	4 columns
covariance		<- apc.fit.model$covariance.canonical
slopes			<- apc.fit.model$slopes
difdif			<- apc.fit.model$difdif
dates			<- apc.fit.model$dates
index.age		<- apc.fit.model$index.age
index.per		<- apc.fit.model$index.per
index.coh		<- apc.fit.model$index.coh
age1			<- apc.fit.model$age1
per1			<- apc.fit.model$per1
coh1			<- apc.fit.model$coh1
unit			<- apc.fit.model$unit
per.zero		<- apc.fit.model$per.zero
per.odd			<- apc.fit.model$per.odd
U				<- apc.fit.model$U
age.max			<- apc.fit.model$age.max
per.max			<- apc.fit.model$per.max
coh.max			<- apc.fit.model$coh.max
model.design	<- apc.fit.model$model.design
##############################
#	derived values
det.max	<- 1+sum(slopes)
det.sub	<- 1+sum(slopes)-sum(difdif)
xi.max	<- det.max+difdif[1]*age.max+difdif[2]*per.max+difdif[3]*coh.max
xi.sub	<- det.sub+difdif[1]*age.max+difdif[2]*per.max+difdif[3]*coh.max
xi.dif	<- det.sub+difdif[1]*(age.max-1)+difdif[2]*(per.max-1)+difdif[3]*(coh.max-1)
xi		<- det.max+difdif[1]*(age.max-2)+difdif[2]*(per.max-2)+difdif[3]*(coh.max-2)
##############################
#	construct for indices for double sums of double difference parameters
#	first for use with (double difference) canonical parameters
index.age.max	<- NULL
index.per.max	<- NULL
index.coh.max	<- NULL
start		<- det.max
if(difdif[1])	{	index.age.max	<- start+seq(1,age.max);	start	<- start+age.max	}
if(difdif[2])	{	index.per.max	<- start+seq(1,per.max);	start	<- start+per.max	}
if(difdif[3])	{	index.coh.max	<- start+seq(1,coh.max);	start	<- start+coh.max	}
#	then for use with in reparametrised submodels in terms of sums of differences 
index.age.sub	<- NULL
index.per.sub	<- NULL
index.coh.sub	<- NULL
if(model.design != "APC")
{
  start		<- det.sub
  if(difdif[1])	{	index.age.sub	<- start+seq(1,age.max);	start	<- start+age.max	}
  if(difdif[2])	{	index.per.sub	<- start+seq(1,per.max);	start	<- start+per.max	}
  if(difdif[3])	{	index.coh.sub	<- start+seq(1,coh.max);	start	<- start+coh.max	}
}
#	then for use with in reparametrised submodels in terms of differences 
index.age.dif	<- NULL
index.per.dif	<- NULL
index.coh.dif	<- NULL
if(model.design != "APC")
{
  start		<- det.sub
  if(difdif[1])	{	index.age.dif	<- start+seq(1,age.max-1);	start	<- start+age.max-1	}
  if(difdif[2])	{	index.per.dif	<- start+seq(1,per.max-1);	start	<- start+per.max-1	}
  if(difdif[3])	{	index.coh.dif	<- start+seq(1,coh.max-1);	start	<- start+coh.max-1	}
}
##############################
#	construct dates for for double difference parameters 
#	first for use with (double difference) canonical parameters
dates.max		<- matrix(data=NA,nrow=xi.max,ncol=1)
if(difdif[1])	dates.max[index.age.max,1]	<- age1+seq(0,age.max-1)*unit	
if(difdif[2])	dates.max[index.per.max,1]	<- per1+seq(0,per.max-1)*unit
if(difdif[3])	dates.max[index.coh.max,1]	<- coh1+seq(0,coh.max-1)*unit
#	then for use with in reparametrised submodels of sums of differences 
dates.sub		<- NULL
if(model.design != "APC")
{
  dates.sub		<- matrix(data=NA,nrow=xi.sub,ncol=1)
  if(difdif[1])	dates.sub[index.age.sub,1]	<- age1+seq(0,age.max-1)*unit	
  if(difdif[2])	dates.sub[index.per.sub,1]	<- per1+seq(0,per.max-1)*unit
  if(difdif[3])	dates.sub[index.coh.sub,1]	<- coh1+seq(0,coh.max-1)*unit
}	
#	then for use with in reparametrised submodels of sums of differences 
dates.dif		<- NULL
if(model.design != "APC")
{
  dates.dif		<- matrix(data=NA,nrow=xi.dif,ncol=1)
  if(difdif[1])	dates.dif[index.age.dif,1]	<- age1+seq(1,age.max-1)*unit	
  if(difdif[2])	dates.dif[index.per.dif,1]	<- per1+seq(1,per.max-1)*unit
  if(difdif[3])	dates.dif[index.coh.dif,1]	<- coh1+seq(1,coh.max-1)*unit
}	
##############################
#	get linear transformation matrix
#		from canonical parameter
#		to standard representation
#	level + slope_age (i-1) + slope_coh (k-1)
#		+ sum sum DD age 	[padded with zeros]
#		+ sum sum DD period	[padded with zeros]
#		+ sum sum DD cohort	[padded with zeros]
#	a summation function is needed
#	for sum sum DD age:		use with U=U
#	for sum sum DD cohort:	use with U=U
#	for sum sum DD period, L=per.odd=TRUE:	use with U=2
#	for sum sum DD period, L=per.odd=FALSE:	use with U=1	
function.ssdd	<- function(n,U)
  #	BN, 4 mar 2015
  #	U is the anchoring point in the summation
{	#	function.ssdd
  m	<- matrix(data=0,nrow=n+2,ncol=n)
  if(U>1)
    for(row in 1:(U-1))
      m[row,row:(U-1)]	<- 1:(U-row)	
    if(U<n+1)		
      for(row in (U+2):(n+2))
        m[row,U:(row-2)]	<- (row-U-1):1	
      return(m)	
}	#	function.ssdd
##############################
#	declare linear transformation matrix
m.ssdd	<-	matrix(data=0,nrow=xi.max,ncol=xi)									
m.ssdd[1:det.max,1:det.max]	<- diag(det.max)							#	level/trend terms
if(difdif[1])	m.ssdd[index.age.max,index.age]	<- function.ssdd(age.max-2,U)			#	alpha
if(difdif[2])	m.ssdd[index.per.max,index.per]	<- function.ssdd(per.max-2,per.odd+1)	#	beta 
if(difdif[3])	m.ssdd[index.coh.max,index.coh]	<- function.ssdd(coh.max-2,U)			# 	gamma
##############################
#	get linear transformation matrix
#		from standard representation
#		to detrended representation
#	this matrix more complicated because
#	linear trends are moved from the
#	time effects to the slopes.
#
#	a detrending function is needed
function.detrend	<- function(n)
  #	BN, 3 apr 2015
  #	in:		n			is the dimension
  #	Out		m			matrix of dimension n x n
  #						for detrending an n vector,
  #						takes an identity matrix
  #						replaces first column with (col-n)/(n-1)
  #						replaces  last column with (1-col)/(n-1)
{	#	function.detrend
  #	m defines the detrending
  m			<- diag(n);
  m[1:n,1]	<- (seq(1:n)-n)/(n-1);
  m[1:n,n]	<- (1-seq(1:n))/(n-1);
  m[1,1]		<- 0;
  m[n,n]		<- 0;
  return(m)
}	#	function.detrend


#### not part of the code for apc_identify: trying to get a handle on function.detrend ####
dog     <- matrix(seq(1,11), nrow=1, ncol=11)
time    <- as.matrix(seq(1,11)) 
det.dog <- dog %*% function.detrend(11)

reg1 <- glm.fit(time, t(dog))
summary.glm(reg1)                                # perfect fit
reg2 <- glm.fit(time, t(det.dog))               
summary.glm(reg2)                                # no fit, exactly

dogsqd <- dog^2
det.dogsqd <- dogsqd %*% function.detrend(11)

reg3 <- glm.fit(time, t(dogsqd))
summary.glm(reg3)                                # good fit
reg4 <- glm.fit(time, t(det.dogsqd))
summary.glm(reg4)                                # no fit, near exactly

noisydog <- dog + rnorm(11)
det.noisydog <- noisydog %*% function.detrend(11)
noisydogsqd <- dogsqd + rnorm(11)
det.noisydogsqd <- noisydogsqd %*% function.detrend(11)

reg5 <- glm.fit(time, t(noisydog))
summary.glm(reg5)                                # good fit
reg6 <- glm.fit(time, t(det.noisydog))
summary.glm(reg6)                                # no fit, near exactly
reg7 <- glm.fit(time, t(noisydogsqd))
summary.glm(reg7)                                # good fit
reg8 <- glm.fit(time, t(det.noisydogsqd))
summary.glm(reg8)                                # no fit, near exactly

#### back to apc.identify ####

#	declare linear transformation matrix
m.detrend	<-	diag(xi.max)
m.detrendA <- m.detrend ########################added to help myself
#	move anchoring of linear trend from U to 1.
if(sum(slopes)==2)
  m.detrend[1,2:3]	<- 1-U
if(sum(slopes)==1 && !slopes[2])
  m.detrend[1,2]		<- 1-U
m.detrendB <- m.detrend ########################added to help myself
#	detrend age effects, move linear trend to deterministics
if(difdif[1])
{	m.detrend[1,index.age.max[1]]	<- 1
m.detrend[2,index.age.max[c(1,age.max)]]	<- c(-1,1)/(age.max-1)
m.detrend[index.age.max,index.age.max]	<- function.detrend(age.max)		
}
m.detrendC <- m.detrend ########################added to help myself
if(difdif[3])
{	# there are 2 slopes if slopes=c(1,0,1)
  # there is  1 slope  if slopes=c(0,0,1)
  # recall det.max <- 1+sum(slopes)
  m.detrend[1,index.coh.max[1]]	<- 1
  m.detrend[det.max,index.coh.max[c(1,coh.max)]]	<- c(-1,1)/(coh.max-1)
  m.detrend[index.coh.max,index.coh.max]	<- function.detrend(coh.max)		
}
m.detrendD <- m.detrend ########################added to help myself
if(difdif[2])
{	# if slopes=c(1,0,1) the period slope gives age & cohort slopes with equal weight
  # if slopes=c(0,1,0) the period slope gives a period     slope
  if(!slopes[2])  m.detrend[1,index.per.max[c(1,per.max)]] <- c(1,0)+c(1,-1)*per.zero/(per.max-1)
  if(slopes[2]) 	m.detrend[1,index.per.max[c(1,per.max)]] <- c(1,0)
  if(slopes[1])	m.detrend[2,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
  if(slopes[2])	m.detrend[2,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
  if(slopes[3])	m.detrend[3,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
  m.detrend[index.per.max,index.per.max]	<- function.detrend(per.max)		
}
m.detrendE <- m.detrend ########################added to help myself
##############################	


#### not part of the code for apc_identify: trying to get a handle on m.detrend ####

View(m.detrendA) # literally just the diagonal matrix
m.detrendA == m.detrendB # all that has changed is two elements in the top row
m.detrendB == m.detrendC # changes in the rows 1+2 + 4-6, cols 4+6 area
m.detrendC == m.detrendD # changes in the rows 1+3 + 12-18, cols 12+18 area
m.detrendE == m.detrendD # changes in the rows 1-3 + 7-11, cols 7+11 area

#### back to apc_identif ####

##############################
#	Now, manipulate estimates using m.ssdd and m.detrend
##############################
#	get estimates
coefficients.ssdd		<- m.ssdd 		%*% coefficients
coefficients.detrend	<- m.detrend 	%*% coefficients.ssdd
##############################
#	construct row names 
names.ssdd		<- c("level")
names.detrend	<- c("level")
if(slopes[1])	{	names.ssdd		<- c(names.ssdd	  ,"age slope")
names.detrend	<- c(names.detrend,"age slope")		}		
if(slopes[2])	{	names.ssdd		<- c(names.ssdd	  ,"period slope")
names.detrend	<- c(names.detrend,"period slope")	}
if(slopes[3])	{	names.ssdd		<- c(names.ssdd	  ,"cohort slope")
names.detrend	<- c(names.detrend,"cohort slope")	}		
if(difdif[1])
  for(i in 1:age.max)
  {	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_age_"   ,as.character((dates.max[index.age.max,1])[i]),sep=""))
  names.detrend	<- c(names.detrend,paste("SS_DD_age_"   ,as.character((dates.max[index.age.max,1])[i]),sep=""))	}
if(difdif[2])															      
  for(i in 1:per.max)															      
  {	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_period_",as.character((dates.max[index.per.max,1])[i]),sep=""))
  names.detrend	<- c(names.detrend,paste("SS_DD_period_",as.character((dates.max[index.per.max,1])[i]),sep=""))	}		
if(difdif[3])																	      
  for(i in 1:coh.max)															      
  {	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_cohort_",as.character((dates.max[index.coh.max,1])[i]),sep=""))
  names.detrend	<- c(names.detrend,paste("SS_DD_cohort_",as.character((dates.max[index.coh.max,1])[i]),sep="")) }
rownames(coefficients.ssdd	 )	<- names.ssdd		
rownames(coefficients.detrend)	<- names.detrend