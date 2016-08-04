# Writes four functions required for APC individual data analysis
# is compatible with AP trapezoid structure
# combines write.functions.R and changes.for.AP.trap.R (both now redundant)


library(apc)
library(plyr)

##### write apc.get.design.indiv.collinear #####
apc.get.design.indiv.collinear <- function(data, unit=1, datatype, 
                                           n.coh.excl.start=0, n.coh.excl.end=0){
  ##### check data has appropriate varnames #####
  if(isTRUE(("age" %in% colnames(data)) 
            && ("period" %in% colnames(data))
            && ("cohort" %in% colnames(data))) == FALSE)
    stop("apc.error: data missing one of age, period, cohort")
  
  datatype.list <- c("AP", "AC rectangular", "CP", "AC triangular", "AP trapezoid")
  if (!isTRUE(datatype %in% datatype.list))
    stop("datatype not recognised")
  
  ##### get contextual information about data #####  
  age1 <- min(data$age)
  per1 <- min(data$period)  
  coh1 <- min(data$cohort)
  
  ##### set dimensions, create indexation matrix #####
  #set dimensionality of data
  I <- as.numeric(length(unique(data$age)))
  J <- as.numeric(length(unique(data$period)))
  K <- as.numeric(length(unique(data$cohort)))
  
  #construct L: {L+1 < j < L+J}
  #create index matrix for data in sample
  #both depend on structure of data (AC rectangular, AP, CP, AC triangular)
  if (datatype == "AC triangular")
  { 
    L=0       
    i.index <- 1:I
    k.index <- 1:K
    i.value.mu <- rep(i.index, each=K)
    k.value.mu <- rep(k.index, I)
    j.value.mu <- i.value.mu + k.value.mu - 1
    mu.index <- as.data.frame(cbind(i.value.mu, j.value.mu, k.value.mu))
  } 
  # AC triangular
  if (datatype == "AC rectangular")
  {
    L=0       
    k.value.mu <- vector(length=(K*(K+1)/2))
    startpoint <- 1
    for(k in 0:(K-1)){
      length <- K-k
      k.value.mu[startpoint:(startpoint+length-1)] <- c(1:(K-k))
      startpoint <- startpoint+length
    }
    i.value.mu <- rep(1:I, seq(K, 1))
    j.value.mu <- i.value.mu + k.value.mu - 1
    mu.index <- as.data.frame(cbind(i.value.mu, j.value.mu, k.value.mu))
  }
  # AC rectangular
  if (datatype == "AP")
  {
    L <- I-1  
    i.index <- 1:I
    j.index <- (L+1):(L+J)
    i.value.mu <- rep(i.index, each=J)
    j.value.mu <- rep(j.index, I)
    k.value.mu <- j.value.mu - i.value.mu + 1
    mu.index <- as.data.frame(cbind(i.value.mu, j.value.mu, k.value.mu))
  } 
  # AP
  if (datatype == "CP")
  {
    L=K-1 
    k.index <- 1:K
    j.index <- (L+1):(L+J)
    k.value.mu <- rep(k.index, each=J)
    j.value.mu <- rep(j.index, K)
    i.value.mu <- j.value.mu - k.value.mu + 1
    mu.index <- as.data.frame(cbind(i.value.mu, j.value.mu, k.value.mu))
  } 
  # CP 
  if (datatype == "AP trapezoid")
  {
    L <- I-1  
    i.index <- 1:I
    j.index <- (L+1):(L+J)
    i.value.mu <- rep(i.index, each=J)
    j.value.mu <- rep(j.index, I)
    k.value.mu <- j.value.mu - i.value.mu + 1
    mu.index <- as.data.frame(cbind(i.value.mu, j.value.mu, k.value.mu))
    
    max.k <- max(k.value.mu)  
    mu.index1 <- mu.index[mu.index$k.value.mu > n.coh.excl.start, ]
    mu.index2 <- mu.index1[mu.index1$k.value.mu <= (max.k - n.coh.excl.end), ]
    mu.index <- mu.index2
  }
  # AP with corners cut out
  
  stopifnot(exists("L")==TRUE)
  #construct U, the reference point of the linear plane
  U <- as.integer((L+3)/2)
  
  ##### create and fill design matrix #####
  #create space of design matrix
  n.rows <- nrow(mu.index)
  n.param.design <- I+J+K-2 # assuming we want all of them
  design <- matrix(data=0, nrow=n.rows, ncol=n.param.design+3) #add additional 3 columns for index
  
  #begin by defining L.odd, max.i, max.j, max.k
  
  L.odd <- !L %% 2 == 0 # iintc that L modulo 2 is 0 (true when odd)
  age.max <- I
  coh.max <- K + n.coh.excl.start + n.coh.excl.end 
  per.max <- J
  
  #fill row-by-row
  for (row in 1:n.rows) {
    i <- mu.index[row, 1]
    j <- mu.index[row, 2]
    k <- mu.index[row, 3]
    design[row, 1] <- 1
    design[row, 2] <- i - U
    design[row, 4] <- k - U
    design[row, 3] <- design[row, 2] + design[row, 4]
    
    #age DDs
    if (i < U) 
      design[row, (4 + i):(4 + U - 1)] <- seq(1, U - i)
    if (i > U + 1) 
      design[row, (4 + U):(4 + U + i - U - 2)] <- seq(i - U - 1, 1)
    
    #period DDs
    if (L.odd && j == 2 * (U - 1)) 
      design[row, (2 + age.max + 1)] <- 1
    if (j > 2 * U) 
      design[row, (2 + age.max + L.odd + 
                     1):(2 + age.max + L.odd + j - 2 * U)] <- seq(j - 2 * U, 1)
    
    #cohort DDs
    if (k < U) 
      design[row, (age.max + per.max + k - n.coh.excl.start):(age.max +per.max + U - 1 - n.coh.excl.start)] <- seq(1, U - k)
    if (k > U + 1) 
      design[row, (age.max + per.max + U - n.coh.excl.start):(age.max + 
                                                                per.max + U + k - U - 2 - n.coh.excl.start)] <- seq(k - U - 1, 1)
    
    #index
    design[row, n.param.design+1] <- i
    design[row, n.param.design+2] <- j
    design[row, n.param.design+3] <- k
  }
  
  designdf <- as.data.frame(design)
  #assign names to variables in design dataframe
  firstfour <- c("level", "age slope", "period slope", "cohort slope")
  age.names <- paste("DD_age", (age1+2):(age1+2+I-2-1), sep="_")
  per.names <- paste("DD_period", (per1+2):(per1+2+J-2-1), sep="_")
  coh.names <- paste("DD_cohort", (coh1+2):(coh1+2+K-2-1), sep="_")
  index <- c("i.value", "j.value", "k.value")
  colnames(designdf) <- c(firstfour, age.names, per.names, coh.names, index)
  
  ##### add ijk values and thus design matrix to dataset#####
  
  #define adjustments between apc and ijk
  age.adjustment <- max(data$age) - max(mu.index$i.value.mu)
  period.adjustment <- max(data$period) - max(mu.index$j.value.mu)
  cohort.adjustment <- max(data$cohort) - max(mu.index$k.value.mu)
  
  #redefine these vectors to be length of data rather than mu.index
  i.value <- data$age - age.adjustment
  j.value <- data$period - period.adjustment
  k.value <- data$cohort - cohort.adjustment
  
  ijk.data <- as.data.frame(cbind(data, i.value, j.value, k.value))
  
  #associates APC design matrix with each individual datapoint as needed
  data.with.design <- join(ijk.data, designdf, by=c("i.value", "j.value", "k.value"))
  
  ##### create final design matrix based on model choice#####
  
  exclude.all <- c("X", "indiv.ID", "age", "period", "cohort", "i.value", "j.value", "k.value")
  
  first.exclusion <- names(data.with.design) %in% c(exclude.all)
  full.design.collinear <- data.with.design[!first.exclusion]
  
  ##### return #####
  valuables <- list(structure.design.collinear = designdf, 
                    full.design.collinear = full.design.collinear, 
                    unit = unit,
                    age1 = age1,
                    per1 = per1,
                    coh1 = coh1,
                    age.max = age.max,
                    per.max = per.max,
                    coh.max = coh.max,
                    per.zero = L,
                    per.odd = L.odd,
                    U = U)
  
  return (valuables)
}


##### write apc.get.design.indiv.model #####
apc.get.design.indiv.model <- function(apc.get.design.indiv.collinear, model.design = "APC",
                                       dep.var = NULL, covariates = NULL){
  
  ##### label important components #####
  structure.design.collinear <- apc.get.design.indiv.collinear$structure.design.collinear 
  full.design.collinear <- apc.get.design.indiv.collinear$full.design.collinear
  
  set.coh.DDs <- colnames(structure.design.collinear)[grep("^DD_cohort_", colnames(structure.design.collinear))]
  set.age.DDs <- colnames(structure.design.collinear)[grep("^DD_age_", colnames(structure.design.collinear))]
  set.per.DDs <- colnames(structure.design.collinear)[grep("^DD_period_", colnames(structure.design.collinear))]
  
  ##### model inclusions #####
  
  if(model.design=="APC")	{	slopes <- c(1,0,1); difdif <- c(1,1,1);	}
  if(model.design=="AP" )	{	slopes <- c(1,0,1); difdif <- c(1,1,0);	}
  if(model.design=="AC" )	{	slopes <- c(1,0,1); difdif <- c(1,0,1);	}
  if(model.design=="PC" )	{	slopes <- c(1,0,1); difdif <- c(0,1,1);	}
  if(model.design=="Ad" )	{	slopes <- c(1,0,1); difdif <- c(1,0,0);	}
  if(model.design=="Pd" )	{	slopes <- c(1,0,1); difdif <- c(0,1,0);	}
  if(model.design=="Cd" )	{	slopes <- c(1,0,1); difdif <- c(0,0,1);	}
  if(model.design=="A"  )	{	slopes <- c(1,0,0); difdif <- c(1,0,0);	}
  if(model.design=="P"  )	{	slopes <- c(0,1,0); difdif <- c(0,1,0);	}
  if(model.design=="C"  )	{	slopes <- c(0,0,1); difdif <- c(0,0,1);	}
  if(model.design=="t"  )	{	slopes <- c(1,0,1); difdif <- c(0,0,0);	}
  if(model.design=="tA" )	{	slopes <- c(1,0,0); difdif <- c(0,0,0);	}
  if(model.design=="tP" )	{	slopes <- c(0,1,0); difdif <- c(0,0,0);	}
  if(model.design=="tC" )	{	slopes <- c(0,0,1); difdif <- c(0,0,0);	}
  if(model.design=="1"  )	{	slopes <- c(0,0,0); difdif <- c(0,0,0);	}
  
  incl <- vector(mode="character")
  incl <- c(incl, "level")
  if (slopes[1]) {incl <- c(incl, "age slope")}
  if (slopes[2]) {incl <- c(incl, "period slope")}
  if (slopes[3]) {incl <- c(incl, "cohort slope")}
  if (difdif[1]) {incl <- c(incl, set.age.DDs)}
  if (difdif[2]) {incl <- c(incl, set.per.DDs)}
  if (difdif[3]) {incl <- c(incl, set.coh.DDs)}
  
  full.inclusion <- names(full.design.collinear) %in% c(incl, covariates)
  final.design <- full.design.collinear[full.inclusion]
  
  DV <- full.design.collinear[dep.var]
  
  
  ##### return #####
  valuables <- list(structure.design.collinear = structure.design.collinear, 
                    full.design.collinear = full.design.collinear,
                    full.design = final.design,
                    dep.var = DV,
                    slopes = slopes,
                    difdif = difdif,
                    xi.dim = length(incl),
                    model.design = model.design,
                    
                    unit = apc.get.design.indiv.collinear$unit,
                    age1 = apc.get.design.indiv.collinear$age1,
                    per1 = apc.get.design.indiv.collinear$per1,
                    coh1 = apc.get.design.indiv.collinear$coh1,
                    age.max = apc.get.design.indiv.collinear$age.max,
                    per.max = apc.get.design.indiv.collinear$per.max,
                    coh.max = apc.get.design.indiv.collinear$coh.max,
                    per.zero = apc.get.design.indiv.collinear$per.zero,
                    per.odd = apc.get.design.indiv.collinear$per.odd,
                    U = apc.get.design.indiv.collinear$U)
  return(valuables)
}


##### write apc.fit.indiv #####
apc.fit.indiv <- function (apc.get.design.indiv.model, model.family=NULL, n.coh.excl.start = 0, n.coh.excl.end = 0){
  
  ##### take out what is needed #####
  design <- apc.get.design.indiv.model$full.design
  DV <- apc.get.design.indiv.model$dep.var
  xi.dim1 <- apc.get.design.indiv.model$xi.dim
  model.design <- apc.get.design.indiv.model$model.design
  
  difdif <- apc.get.design.indiv.model$difdif
  slopes <- apc.get.design.indiv.model$slopes
  age.max <- apc.get.design.indiv.model$age.max
  per.max <- apc.get.design.indiv.model$per.max
  coh.max <- apc.get.design.indiv.model$coh.max
  age1 <- apc.get.design.indiv.model$age1
  per1 <- apc.get.design.indiv.model$per1
  coh1 <- apc.get.design.indiv.model$coh1
  unit <- apc.get.design.indiv.model$unit
  
  ##### run regression #####
  # need to extend this once other models are understood
  if (model.family == "binomial")
    fit <- glm.fit(design, DV[, 1], family=binomial())
  if (model.family == "gaussian")
    fit <- glm.fit(design, DV[, 1], family=gaussian())
  
  ##### get coeffs, covariance #####
  coefficients	<- summary.glm(fit)$coefficients
  covariance	<- summary.glm(fit)$cov.scaled
  
  n.coeff <- length(coefficients[, 1])
  n.coeff.canonical <- xi.dim1
  
  if (model.design == "1"){
    coefficients.canonical	<- t(coefficients[(n.coeff-n.coeff.canonical+1):n.coeff, ])
    covariance.canonical	<- summary.glm(fit)$cov.scaled[(n.coeff-n.coeff.canonical+1):n.coeff, 
                                                        (n.coeff-n.coeff.canonical+1):n.coeff]  
    
    coefficients.covariates <- t(coefficients[1:(n.coeff-n.coeff.canonical), ])
    
    #### WHY are these happening here when the programme already does it? ####
    #	get standard errors 
    coefficients.canonical[,2]	<- sqrt(covariance.canonical)
    #	get t-statistics
    coefficients.canonical[,3]	<- coefficients.canonical[,1] 	/ coefficients.canonical[,2]
    #	get p-values
    coefficients.canonical[,4]	<- 2*pnorm(abs(coefficients.canonical[  ,3]),lower.tail=FALSE)
  }
  
  if (!model.design == "1"){
    coefficients.canonical	<- summary.glm(fit)$coefficients[(n.coeff-n.coeff.canonical+1):n.coeff, ]
    covariance.canonical	<- summary.glm(fit)$cov.scaled[(n.coeff-n.coeff.canonical+1):n.coeff, 
                                                        (n.coeff-n.coeff.canonical+1):n.coeff]  
    
    coefficients.covariates <- coefficients[1:(n.coeff-n.coeff.canonical), ]
    
    #### WHY are these happening here when the programme already does it? ####
    #	get standard errors 
    coefficients.canonical[,2]	<- sqrt(diag(covariance.canonical))
    #	get t-statistics
    coefficients.canonical[,3]	<- coefficients.canonical[,1] 	/ coefficients.canonical[,2]
    #	get p-values
    coefficients.canonical[,4]	<- 2*pnorm(abs(coefficients.canonical[  ,3]),lower.tail=FALSE)
  }
  
  ##### Apparently we also need these for plotting #####
  index.age	<- NULL
  index.per	<- NULL
  index.coh	<- NULL
  
  # redefine coh.max for this section (as per options A, c in possibility table)
  coh.max.from.earlier <- coh.max
  coh.max <- coh.max.from.earlier - n.coh.excl.start - n.coh.excl.end
  
  start		<- 1+sum(slopes)
  if(difdif[1])	{	index.age	<- start+seq(1,age.max-2);	start	<- start+age.max-2	}
  if(difdif[2])	{	index.per	<- start+seq(1,per.max-2);	start	<- start+per.max-2	}
  if(difdif[3])	{	index.coh	<- start+seq(1,coh.max-2);	start	<- start+coh.max-2	}
  xi.dim2		<- start
  
  dates		<- matrix(data=NA,nrow=xi.dim2,ncol=1)			
  if(difdif[1])	dates[index.age,1]	<- age1+seq(2,age.max-1)*unit	
  if(difdif[2])	dates[index.per,1]	<- per1+seq(2,per.max-1)*unit
  if(difdif[3])	dates[index.coh,1]	<- coh1+seq(2,coh.max-1)*unit
  
  ##### return #####
  valuables <- c(fit, list(coefficients.canonical = coefficients.canonical,
                           covariance.canonical = covariance.canonical,
                           dates = dates,
                           index.age = index.age,
                           index.coh = index.coh,
                           index.per = index.per,
                           
                           coefficients.covariates = coefficients.covariates, 
                           
                           model.family = model.family,
                           
                           difdif = difdif,
                           slopes = slopes,
                           age1 = age1,
                           per1 = per1,
                           coh1 = coh1,
                           unit = unit,
                           age.max = age.max,
                           per.max = per.max,
                           coh.max = coh.max,
                           
                           model.design = apc.get.design.indiv.model$model.design,
                           per.zero = apc.get.design.indiv.model$per.zero,
                           per.odd = apc.get.design.indiv.model$per.odd,
                           U = apc.get.design.indiv.model$U))
  return(valuables)
}


##### write	apc.fit.indiv.table #####
apc.fit.indiv.table	<- function(apc.get.design.indiv.collinear, model.family, DV, covariates = NULL)
{	#	apc.fit.table
  ######################
  #	model families
  model.family.list		<- c("binomial.dose.response","poisson.response","od.poisson.response",
                          "poisson.dose.response","gaussian.rates","gaussian.response",
                          "log.normal.response", "binomial", "gaussian")                                  # added binomial
  model.family.gaussian	<- c("gaussian.rates","gaussian.response","log.normal.response", "gaussian")
  model.family.od			<- c("od.poisson.response")
  ######################
  #	check input
  if(isTRUE(model.family %in% model.family.list)==FALSE)
    return(cat("apc.fit.table error: model.family has wrong argument \n"))
  ######################
  #	get index
  #  if(is.null(apc.index)==TRUE)
  #    apc.index	<- apc.get.index(apc.data.list)
  ######################
  #	Function to get one line of table from two fits
  fit.tab.line.glm	<- function(fit.U,fit.R,gaussian)
    #	BN 20 Sep 2013
  {
    dev.U	<- fit.U$deviance
    dev.R	<- fit.R$deviance
    df.U	<- fit.U$df.residual
    df.R	<- fit.R$df.residual
    LR	<- dev.R-dev.U
    df	<- df.R-df.U
    aic	<- fit.R$aic
    if(gaussian)		
      return(round(c(dev.R,df.R,LR,df,pchisq(LR,df,lower.tail=FALSE),aic),digits=3))
    else	
      return(round(c(dev.R,df.R,pchisq(dev.R,df.R,lower.tail=FALSE),LR,df,pchisq(LR,df,lower.tail=FALSE),aic),digits=3))
  }
  ######################	
  model.design.list	<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC", "1")
  ######################	
  #	number of columns
  ncol <- 7
  if(isTRUE(model.family %in% model.family.gaussian))	ncol <- 6
  if(isTRUE(model.family %in% model.family.od))		ncol <- 9
  ######################
  #	declare table
  fit.tab		<- matrix(nrow=length(model.design.list),ncol=ncol,data=NA)
  
  
  ###########THIS IS THE SECTION THAT NEEDS EDITING
  
  
  #	unrestricted apc model
  design.apc <- apc.get.design.indiv.model(apc.get.design.indiv.collinear, model.design = "APC", dep.var = DV, covariates = covariates) 
  fit.apc	<- apc.fit.indiv(design.apc, model.family = model.family)             #check call to DV and model.family
  #	model list
  for(i in 1:length(model.design.list))
  {
    design.sub <- apc.get.design.indiv.model(apc.get.design.indiv.collinear, model.design = model.design.list[i], dep.var = DV, covariates = covariates)
    fit.sub			<- apc.fit.indiv(design.sub, model.family = model.family)       # again check call to model.family
    if(isTRUE(model.family %in% model.family.gaussian))
      fit.tab[i,1:6]		<- fit.tab.line.glm(fit.apc, fit.sub,1)
    else	
      fit.tab[i,1:7]		<- fit.tab.line.glm(fit.apc, fit.sub,0)         # check whether we are getting the expected outputs 
    # from my fit code to enter this - see above e.g. df, aic
  }
  
  
  ##################
  
  
  
  #	OD F-test
  if(isTRUE(model.family %in% model.family.od))
    for(i in 2:length(model.design.list))
    {	fit.tab[i,8]	= (fit.tab[i,4]/fit.tab[i,5])/(fit.tab[1,1]/fit.tab[1,2])
    fit.tab[i,9]	= round(pf(fit.tab[i,8],fit.tab[i,5],fit.tab[1,2],lower.tail=FALSE),digits=3)
    }
  #	insert NAs
  insert			<- c(4,5,6)
  if(isTRUE(model.family %in% model.family.gaussian))
    insert			<- c(3,4,5)
  fit.tab[1,insert] <- NA
  #	row names
  rownames(fit.tab)	<- model.design.list
  #	column names
  colnames		<- c("-2logL","df.residual","prob(>chi_sq)","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic")
  if(isTRUE(model.family %in% model.family.gaussian))
    colnames		<- c("-2logL","df.residual","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic")
  if(isTRUE(model.family %in% model.family.od))
    colnames		<- c("-2logL","df.residual","prob(>chi_sq)","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic","F","prob(>F)")	
  colnames(fit.tab)	<- colnames		
  ######################
  return(fit.tab)
}	#	apc.fit.table

##### write model.DVcont.manyEVcont.TS #####
model.DVcont.manyEVcont.TS <- function(data, depvar, evars=NULL){
  
  # evars needs to be a c(, ,)
  
  #get n.group
  get.n.group <-  ddply(data, .variables = c("cell.name"), 
                        function(dfr, colnm){length(dfr[,colnm])}, depvar)
  colnames(get.n.group)[2] <- "n.group"  
  
  if(!is.null(evars)){
    # 1 reg continuous and otherEV on TS. Get EV excluding TS (residual)
    #first stage
    
    part1 <- get.n.group
    # get all groupsum evars
    for (i in evars){
      get.sum.c2 <- ddply(data, .variables = c("cell.name"),
                          function(dfr, colnm){sum(dfr[, colnm])}, i)
      colnames(get.sum.c2)[2] <- paste("groupsum", i, sep=".")
      part1 <- join(part1, get.sum.c2, by="cell.name")
    }
    
    # get all the means
    groupmeans.matrix <- as.data.frame(matrix(nrow=nrow(part1), ncol=length(evars)))
    for (i in 1:length(evars)){
      cbar <- as.matrix(part1[, 2+i] / part1[,2]) # this divides each groupsum by n.group
      groupmeans.matrix[, i] <- cbar
      colnames(groupmeans.matrix)[i] <- paste("groupmean", evars[i], sep = ".")
    }
    
    ### groupmean is the estimated value of that regressor based on TS dummies
    # now need to get the residual for each evar once that estimated value is subtracted out
    
    part2 <- cbind(part1, groupmeans.matrix)
    data.and.first.coeff <- join(data, part2, by="cell.name")
    
    # residuals of c from x
    residuals.matrix <- as.data.frame(matrix(nrow=nrow(data.and.first.coeff), ncol=length(evars)))
    for (i in 1:length(evars)){
      true <- data.and.first.coeff[evars[i]]
      predicted.name <- paste("groupmean", evars[i], sep=".")
      predicted <- data.and.first.coeff[predicted.name]
      resid <- true - predicted
      residuals.matrix[, i] <- resid
      colnames(residuals.matrix)[i] <- paste("resid", evars[i], sep=".")
    }
    
    data.resid.c.on.x <- cbind(data.and.first.coeff, residuals.matrix)  
    
    
    # second stage: regress y on c residuals
    stage2.depvar <- data.resid.c.on.x[depvar]
    stage2.evar <- residuals.matrix
    stage2.regress <- glm.fit(as.matrix(stage2.evar), as.matrix(stage2.depvar), family = gaussian())
    truebetahat <- stage2.regress$coefficients
    
    
    #third stage: get y less the evars times their truebetahats (ie predicted y)
    predicted.y.from.evars <- as.matrix(data.resid.c.on.x[evars]) %*% as.matrix(truebetahat)
    y.less.truec <- as.matrix(data.resid.c.on.x[depvar] - predicted.y.from.evars)
    colnames(y.less.truec) <- "y.less.truec"
    
    #fourth stage: add these y.less.truec to the dataset, then get their means (which become the predictors from TS of them)
    data.y.less.truec <- cbind(data.resid.c.on.x, y.less.truec)
    get.sum.y <-  ddply(data.y.less.truec, .variables = c("cell.name"), 
                        function(dfr, colnm){sum(dfr[,colnm])}, "y.less.truec")
    colnames(get.sum.y)[2] <- paste("groupsum.resids", depvar, sep=".")
    
    data.and.second.coeff <- join(data.y.less.truec, get.sum.y, by="cell.name")
    sumycall <- paste("groupsum.resids", depvar, sep=".")
    truekappahat <- data.and.second.coeff[sumycall]/data.and.second.coeff["n.group"]
    colnames(truekappahat) <- "truekappahat"
    data.kappa <- cbind(data.and.second.coeff, truekappahat)
  }  else{
    # we skip straight to the fourth stage! 
    data.ngroup <- join(data, get.n.group, by="cell.name")
    get.sum.y <- ddply(data.ngroup, .variables=c("cell.name"),
                       function(dfr, colnm){sum(dfr[,colnm])}, depvar)
    colnames(get.sum.y)[2] <- paste("groupsum", depvar, sep=".")
    data.and.coeff <- join(data.ngroup, get.sum.y, by="cell.name")
    sumycall <- paste("groupsum", depvar, sep=".")
    truekappahat <- data.and.coeff[sumycall]/data.and.coeff["n.group"]
    colnames(truekappahat) <- "truekappahat"
    data.kappa <- cbind(data.and.coeff, truekappahat)
    # and just for the stage 5 and the output
    predicted.y.from.evars <- 0
    truebetahat <- NA
  }
  
  #fifth stage : calculate final residuals
  final.residuals <- data.kappa[depvar] - predicted.y.from.evars - data.kappa["truekappahat"]
  colnames(final.residuals) <- "final.residuals"
  final.data <- cbind(data.kappa, final.residuals)
  
  #return
  
  valuables <- list(final.data = final.data,
                    TSS = sum(final.data[depvar]^2),
                    trueRSS = sum(final.data["final.residuals"]^2),
                    truebetahat = truebetahat,
                    truekappahat= truekappahat)
  return(valuables)
}


##### write ftest.apc.vs.ts #####
ftest.apc.vs.ts <- function(data, datatype, model.design="APC", depvar, covariates=NULL,
                            model.family="gaussian", unit=1, n.coh.excl.start=0, n.coh.excl.end = 0,
                            existing.apc.fit.indiv=NULL, existing.TS=NULL){
  
  #fit the two models - continuous as DV, dichot as EV
  #first the APC one
  if(is.null(existing.apc.fit.indiv)){
    collinear <- apc.get.design.indiv.collinear(data, datatype=datatype, unit = unit,
                                                n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end)
    design <- apc.get.design.indiv.model(collinear, model.design = model.design, dep.var=depvar, covariates=covariates)
    fitted <- apc.fit.indiv(design, model.family=model.family, n.coh.excl.start = n.coh.excl.start, 
                            n.coh.excl.end = n.coh.excl.end)
  } else 
    fitted <- existing.apc.fit.indiv
  
  
  TSS.APC <- sum(fitted$y^2)
  RSS.APC <- sum(fitted$residuals^2)
  ESS.APC <- sum(fitted$fitted.values^2)
  
  nparam.APC <- length(fitted$coefficients)
  sample.N.APC <- length(fitted$y)
  
  #second the dummy one
  if(is.null(existing.TS)){
    dummymodel <- model.DVcont.manyEVcont.TS(data, depvar = depvar,
                                             evars = covariates)
  } else
    dummymodel <- existing.TS
  
  RSS.dummy <- dummymodel$trueRSS
  if (!is.na(dummymodel$truebetahat[1])){
    nparam.dummy <- length(unique(dummymodel$final.data$cell.name))+length(dummymodel$truebetahat)
  } else{
    nparam.dummy <- length(unique(dummymodel$final.data$cell.name))
  }
  
  #elements of F-test
  big.model.RSS <- RSS.dummy
  small.model.RSS <- RSS.APC
  big.model.nparams <- nparam.dummy
  small.model.nparams <- nparam.APC
  sample.N <- sample.N.APC
  
  df.num <- big.model.nparams - small.model.nparams
  df.denom <- sample.N - big.model.nparams
  df <- paste("(", paste(as.character(df.num), as.character(df.denom), sep=", "), ")", sep="")
  
  fstat.num <- (small.model.RSS - big.model.RSS)/df.num
  fstat.denom <- (big.model.RSS - 0)/df.denom
  fstat <- fstat.num/fstat.denom
  
  #get pvalue
  p.value <- pf(fstat, df1 = df.num, df2 = df.denom, lower.tail = FALSE)
  
  #valuables
  valuables <- list(fstat = fstat,
                    df = df,
                    df.num = df.num,
                    df.denom = df.denom,
                    p.value = p.value,
                    aic.APC = fitted$aic)
  return(valuables)
}

##### write ftest.apcreduced.vs.apc #####
ftest.apcreduced.vs.apc <- function(data, datatype, model.design="APC", reduced.model.design,
                                    depvar, covariates=NULL, model.family="gaussian", unit=1, 
                                    n.coh.excl.start=0, n.coh.excl.end=0,
                                    existing.apc.fit.indiv=NULL, existing.apcreduced.fit.indiv=NULL){
  
  #fit the two models - continuous as DV, dichot as EV
  #first the reduced APC one
  if(is.null(existing.apcreduced.fit.indiv)){
    collinear.R <- apc.get.design.indiv.collinear(data, datatype=datatype, unit = unit,
                                                  n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end)
    design.R <- apc.get.design.indiv.model(apc.get.design.indiv.collinear = collinear.R, 
                                           model.design = reduced.model.design, dep.var=depvar,
                                           covariates=covariates)
    fitted.R <- apc.fit.indiv(design.R, model.family=model.family, n.coh.excl.start = n.coh.excl.start,
                              n.coh.excl.end = n.coh.excl.end)
  } else 
    fitted.R <- existing.apcreduced.fit.indiv
  
  TSS.APC.R <- sum(fitted.R$y^2)
  RSS.APC.R <- sum(fitted.R$residuals^2)
  ESS.APC.R <- sum(fitted.R$fitted.values^2)
  
  nparam.APC.R <- length(fitted.R$coefficients)
  sample.N.APC.R <- length(fitted.R$y)
  
  #second the full APC one
  if(is.null(existing.apc.fit.indiv)){
    collinear.U <- apc.get.design.indiv.collinear(data, datatype=datatype, unit = unit,
                                                  n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end)
    design.U <- apc.get.design.indiv.model(apc.get.design.indiv.collinear = collinear.U, 
                                           model.design = model.design, dep.var=depvar,
                                           covariates=covariates)
    fitted.U <- apc.fit.indiv(design.U, model.family=model.family, n.coh.excl.start = n.coh.excl.start,
                              n.coh.excl.end = n.coh.excl.end)
  } else 
    fitted.U <- existing.apc.fit.indiv
  
  
  TSS.APC.U <- sum(fitted.U$y^2)
  RSS.APC.U <- sum(fitted.U$residuals^2)
  ESS.APC.U <- sum(fitted.U$fitted.values^2)
  
  nparam.APC.U <- length(fitted.U$coefficients)
  sample.N.APC.U <- length(fitted.U$y)
  
  #elements of F-test
  big.model.RSS <- RSS.APC.U
  small.model.RSS <- RSS.APC.R
  big.model.nparams <- nparam.APC.U
  small.model.nparams <- nparam.APC.R
  sample.N <- sample.N.APC.U
  
  df.num <- big.model.nparams - small.model.nparams
  df.denom <- sample.N - big.model.nparams
  df <- paste("(", paste(as.character(df.num), as.character(df.denom), sep=", "), ")", sep="")
  
  fstat.num <- (small.model.RSS - big.model.RSS)/df.num
  fstat.denom <- (big.model.RSS - 0)/df.denom
  fstat <- fstat.num/fstat.denom
  
  #get pvalue
  p.value <- pf(fstat, df1 = df.num, df2 = df.denom, lower.tail = FALSE)
  
  #valuables
  valuables <- list(fstat = fstat,
                    df = df,
                    df.num = df.num,
                    df.denom = df.denom,
                    p.value = p.value,
                    aic.U = fitted.U$aic,
                    aic.R = fitted.R$aic)
  return(valuables)
}

##### write Ftest.comparison.table #####
Ftest.comparison.table <- function(data, datatype, depvar, covariates=NULL, unit=1, 
                                   n.coh.excl.start=0, n.coh.excl.end=0){
  get.table.line <- function(ftest1, ftest2, first=FALSE){
    if(isTRUE(first))
      line <- c(round(ftest1$fstat, 3), ftest1$df.num, round(ftest1$p.value,3), NA,
                NA, NA, round(ftest1$aic.APC,3))
    else
      line <- c(round(ftest1$fstat, 3), ftest1$df.num, round(ftest1$p.value,3), round(ftest2$fstat, 3),
                ftest2$df.num, round(ftest2$p.value, 3), round(ftest2$aic.R,3))
    return(line)
  }
  
  model.design.list	<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC", "1")
  
  fit.tab <- matrix(nrow=length(model.design.list),ncol=7,data=NA)
  
  fullAPCcollinear <- apc.get.design.indiv.collinear(data, datatype=datatype, unit=unit,
                                                     n.coh.excl.start=n.coh.excl.start, n.coh.excl.end=n.coh.excl.end)
  fullAPCdesign <- apc.get.design.indiv.model(fullAPCcollinear, dep.var = depvar, covariates = covariates)
  fullAPCfit <- apc.fit.indiv(fullAPCdesign, model.family="gaussian", n.coh.excl.start = n.coh.excl.start,
                              n.coh.excl.end = n.coh.excl.end)
  TSfit <- model.DVcont.manyEVcont.TS(data, depvar=depvar, evars=covariates)
  
  ftest.full.vs.TS <- ftest.apc.vs.ts(existing.apc.fit.indiv = fullAPCfit, existing.TS = TSfit)
  first.line <- get.table.line(ftest1 = ftest.full.vs.TS, first=TRUE)
  fit.tab[1, 1:7] <- first.line
  
  for (i in 2:length(model.design.list)){
    againstTS <- ftest.apc.vs.ts(data, datatype = datatype, model.design=model.design.list[i],
                                 depvar = depvar,covariates = covariates, unit = unit,
                                 n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end,
                                 existing.TS = TSfit)
    againstAPC <- ftest.apcreduced.vs.apc(data, datatype=datatype, reduced.model.design = model.design.list[i],
                                          depvar=depvar, covariates=covariates, unit = unit, 
                                          n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end,
                                          existing.apc.fit.indiv = fullAPCfit)
    fit.tab[i, 1:7] <- get.table.line(againstTS, againstAPC)
  }
  
  df.againstTS <- paste("DF ( * , ", as.character(againstTS$df.denom), ")", sep="")
  df.againstAPC <- paste("DF ( * , ", as.character(againstAPC$df.denom), ")", sep="")
  
  table.colnames <- c("F-test vs TS", df.againstTS, "p-value", "F-test vs APC", 
                      df.againstAPC, "p-value", "AIC")
  
  colnames(fit.tab) <- table.colnames
  rownames(fit.tab) <- model.design.list
  
  return(fit.tab)
}


##### write var.apc.plot.fit
var.apc.plot.fit	<- function(apc.fit.model,scale=FALSE,sdv.at.zero=TRUE,type="detrend",sub.plot=NULL,main.outer=NULL,main.sub=NULL,cex=NULL,cex.axis=NULL, cex.main=2, mgp=c(2, 1, 0), theight=1, mar=c(4, 3, 2, 0), oma = c(0, 0, 5, 1))
# incorporates alterations to some of the visuals  
{	#	apc.plot.fit
  #################
  #	change type
  if(type=="ss.dd")	type<-"sum.sum"
  #################
  #	get values from fit
  coefficients.canonical	<- apc.fit.model$coefficients.canonical	
  slopes					<- apc.fit.model$slopes					
  difdif					<- apc.fit.model$difdif					
  index.age				<- apc.fit.model$index.age 				
  index.per				<- apc.fit.model$index.per 				
  index.coh				<- apc.fit.model$index.coh 				
  dates					<- apc.fit.model$dates
  model.design			<- apc.fit.model$model.design
  model.family			<- apc.fit.model$model.family
  age1					<- apc.fit.model$age1
  per1					<- apc.fit.model$per1
  coh1					<- apc.fit.model$coh1
  unit					<- apc.fit.model$unit
  age.max					<- apc.fit.model$age.max
  per.max					<- apc.fit.model$per.max
  coh.max					<- apc.fit.model$coh.max	
  per.odd					<- apc.fit.model$per.odd
  U						<- apc.fit.model$U
  #################
  # 	identify fit
  apc.id	<- apc.identify(apc.fit.model)
  index.age.max			<- apc.id$index.age.max 				
  index.per.max			<- apc.id$index.per.max 				
  index.coh.max			<- apc.id$index.coh.max 				
  dates.max				<- apc.id$dates.max
  index.age.sub	 		<- apc.id$index.age.sub	  		 
  index.per.sub			<- apc.id$index.per.sub	  		 
  index.coh.sub			<- apc.id$index.coh.sub	  		 
  dates.sub				<- apc.id$dates.sub		  		 
  index.age.dif	 		<- apc.id$index.age.dif	 		 
  index.per.dif			<- apc.id$index.per.dif	 		 
  index.coh.dif			<- apc.id$index.coh.dif	 		 
  dates.dif				<- apc.id$dates.dif		 		 
  coefficients.ssdd		<- apc.id$coefficients.ssdd	
  coefficients.detrend	<- apc.id$coefficients.detrend	
  coefficients.demean		<- apc.id$coefficients.demean	
  coefficients.dif		<- apc.id$coefficients.dif		
  ##############################
  #	check model design
  if(isTRUE(type %in% c("dif")))
  {	if(model.design=="APC")
    return(cat("apc.error: differences not identified when model.design is APC.  Type cannot be demean or dif \n"))
    if(model.design %in% c("Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1"))	
      return(cat("apc.error: types demean and dif not implemented for model designs At, Pt, Ct, A, P, C, t, tA, tP, tC, 1 \n"))
  }
  ###########################################
  #	construct ingredients to plot depending on type
  mixed	<- FALSE
  if(model.family=="poisson.response")	mixed	<- TRUE
  ###########################################
  #	declare variables
  v.do.plot		<- vector(length=9)
  l.dates			<- list(a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1,i=1)
  l.coefficients	<- list(a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1,i=1)
  v.main.sub		<- vector(length=9) 
  v.xlab			<- vector(length=9)								  
  v.intercept		<- vector(length=9)
  v.tau			<- vector(length=9)
  ###########################################
  #	type is "detrend" or "sum.sum"	
  if(type %in% c("detrend","sum.sum"))
  {	if(type=="detrend")
    main	<- paste("APC canonical parameters & detrended representation","\n","model.design=",model.design, "(1/2 std blue/red)")
  if(type=="sum.sum")
    main	<- paste("APC canonical parameters & standard representation","\n","model.design=",model.design, "(1/2 std blue/red)")	
  ##################
  #	do plot?
  v.do.plot[1:3]	<- difdif	
  v.do.plot[4]	<- isTRUE(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","t","tA","tP"))
  v.do.plot[5]	<- TRUE
  v.do.plot[6]	<- isTRUE(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","C","t","tC"))
  v.do.plot[7:9]	<- difdif	
  ##################
  #	sub.main
  v.main.sub[1]	<- expression(paste("(a) ",Delta^2,alpha))
  v.main.sub[2]	<- expression(paste("(b) ",Delta^2,beta))
  v.main.sub[3]	<- expression(paste("(c) ",Delta^2,gamma))
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t","A","tA"))
    v.main.sub[4]	<- "(d)  first linear trend"
  if(model.design %in% c("A","tA"))
    v.main.sub[4]	<-"(d)  age linear trend"
  if(model.design %in% c("P","tP"))									
    v.main.sub[4]	<- "(d)  period linear trend"	
  if(!mixed)
    v.main.sub[5]	<- "(e)  level"
  if(mixed)
    v.main.sub[5]	<- "(e)  aggregate mean"
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t"))
    v.main.sub[6]	<- "(f)  second linear trend"
  if(model.design %in% c("C","tC"))									
    v.main.sub[6]	<- "(f)  cohort linear trend"
  if(type=="detrend")
  {	v.main.sub[7]	<- expression(paste("(g) detrended ",Sigma^2,Delta^2,alpha))
  v.main.sub[8]	<- expression(paste("(h) detrended ",Sigma^2,Delta^2,beta))
  v.main.sub[9]	<- expression(paste("(i) detrended ",Sigma^2,Delta^2,gamma))
  }	
  if(type=="sum.sum")
  {	v.main.sub[7]	<- expression(paste("(g) ",Sigma^2,Delta^2,alpha))
  v.main.sub[8]	<- expression(paste("(h) ",Sigma^2,Delta^2,beta))
  v.main.sub[9]	<- expression(paste("(i) ",Sigma^2,Delta^2,gamma))
  }
  ##################
  #	dates
  l.dates[[1]]	<- dates[index.age,1]
  l.dates[[2]]	<- dates[index.per,1]
  l.dates[[3]]	<- dates[index.coh,1]
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t","A","tA"))
    l.dates[[4]]	<- age1+seq(0,age.max-1)*unit
  if(model.design %in% c("P","tP"))									
    l.dates[[4]]	<- per1+seq(0,per.max-1)*unit
  l.dates[[5]]	<- c(0,1)  # matrix(data=c(0,1)		     ,nrow=2	  ,ncol=1)
  l.dates[[6]]	<- coh1+seq(0,coh.max-1)*unit
  l.dates[[7]]	<- dates.max[index.age.max,1]
  l.dates[[8]]	<- dates.max[index.per.max,1]
  l.dates[[9]]	<- dates.max[index.coh.max,1]		
  ##################
  #	coefficients
  l.coefficients[[1]]	<- coefficients.canonical[index.age,]
  l.coefficients[[2]]	<- coefficients.canonical[index.per,]
  l.coefficients[[3]]	<- coefficients.canonical[index.coh,]
  if(type=="detrend")
  {	coefficients.sum.sum	<- coefficients.detrend
  UU	<- 1
  }
  if(type=="sum.sum")
  {	coefficients.sum.sum	<- coefficients.ssdd
  UU <- U
  }
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t","A","tA"))
    l.coefficients[[4]]	<- matrix(data=seq(1,age.max)-UU  		,nrow=age.max,ncol=1) %*% coefficients.sum.sum[2,]
  if(model.design %in% c("P","tP"))									
    l.coefficients[[4]]	<- matrix(data=seq(1,per.max)-per.odd-1	,nrow=per.max,ncol=1) %*% coefficients.sum.sum[2,]
  l.coefficients[[5]]		<- matrix(data=c(1,1)		       		,nrow=2	     ,ncol=1) %*% coefficients.sum.sum[1,]
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t"))
    l.coefficients[[6]]	<- matrix(data=seq(1,coh.max)-UU 		,nrow=coh.max,ncol=1) %*% coefficients.sum.sum[3,]
  if(model.design %in% c("C","tC"))									
    l.coefficients[[6]]	<- matrix(data=seq(1,coh.max)-UU 	 	,nrow=coh.max,ncol=1) %*% coefficients.sum.sum[2,]		
  l.coefficients[[7]]	<- coefficients.sum.sum[index.age.max,]
  l.coefficients[[8]]	<- coefficients.sum.sum[index.per.max,]
  l.coefficients[[9]]	<- coefficients.sum.sum[index.coh.max,]
  ####################
  #	xlab
  v.xlab[1:3]	<- c("age","period","cohort")
  if(model.design %in% c("APC","AP","AC","PC","Ad","Pd","Cd","t","A","tA"))
    v.xlab[4]	<- "age"
  if(model.design %in% c("P","tP"))									
    v.xlab[4]	<- "period"
  if(!mixed)
    v.xlab[5]	<- "age, period, cohort"
  if(mixed)
    v.xlab[5]		<- ""
  v.xlab[6]	<- "cohort"
  v.xlab[7:9]	<- c("age","period","cohort")
  ####################
  #	intercept
  v.intercept[5]	<- TRUE
  ####################
  #	tau
  if(mixed)
    v.tau[5]	<- TRUE
  } 
  ###########################################
  #	type is "dif"	
  if(type %in% c("dif"))
  {	main	<- paste("Difference parameters & demeaned representation","\n","model.design=",model.design, "(1/2 std blue/red)")
  ##################
  #	do plot?
  v.do.plot[1:3]	<- difdif
  v.do.plot[5]	<- TRUE
  v.do.plot[7:9]	<- difdif
  ##################
  #	sub.main
  v.main.sub[1]	<- expression(paste("(a) ",Delta,alpha))
  v.main.sub[2]	<- expression(paste("(b) ",Delta,beta))
  v.main.sub[3]	<- expression(paste("(c) ",Delta,gamma))
  if(!mixed)	v.main.sub[5]	<- "(e)  level"
  if(mixed)	v.main.sub[5]	<- "(e)  aggregate mean"
  v.main.sub[7]	<- expression(paste("(g) demeaned ",Sigma,Delta,alpha))
  v.main.sub[8]	<- expression(paste("(h) demeaned ",Sigma,Delta,beta))
  v.main.sub[9]	<- expression(paste("(i) demeaned ",Sigma,Delta,gamma))
  ##################
  #	dates
  l.dates[[1]]	<- dates.dif[index.age.dif,1]
  l.dates[[2]]	<- dates.dif[index.per.dif,1]
  l.dates[[3]]	<- dates.dif[index.coh.dif,1]
  l.dates[[5]]	<- c(0,1)  # matrix(data=c(0,1)		     ,nrow=2	  ,ncol=1)
  l.dates[[7]]	<- dates.sub[index.age.sub,1]
  l.dates[[8]]	<- dates.sub[index.per.sub,1]
  l.dates[[9]]	<- dates.sub[index.coh.sub,1]		
  ##################
  #	coefficients
  l.coefficients[[1]]	<- coefficients.dif[index.age.dif,]
  l.coefficients[[2]]	<- coefficients.dif[index.per.dif,]
  l.coefficients[[3]]	<- coefficients.dif[index.coh.dif,]
  l.coefficients[[5]]		<- matrix(data=c(1,1)		     ,nrow=2	  ,ncol=1) %*% coefficients.demean[1,]
  l.coefficients[[7]]	<- coefficients.demean[index.age.sub,]
  l.coefficients[[8]]	<- coefficients.demean[index.per.sub,]
  l.coefficients[[9]]	<- coefficients.demean[index.coh.sub,]
  ####################
  #	xlab
  v.xlab[1:3]	<- c("age","period","cohort")
  if(!mixed)	v.xlab[5]	<- "age, period, cohort"
  if(mixed) 	v.xlab[5]		<- ""
  v.xlab[7:9]	<- c("age","period","cohort")
  ####################
  #	intercept
  v.intercept[5]	<- TRUE
  ####################
  #	tau
  if(mixed)
    v.tau[5]	<- TRUE
  } 
  #############################################
  #	use arguments of function
  if(is.null(main.outer)==FALSE)
    main	<- main.outer
  if(is.null(main.sub)==FALSE)
    v.main.sub	<- main.sub
  if(scale==1 && model.family=="binomial.dose.response")	scale <- 2	
  #######################################################
  #	Internal function to plot estimates with sdv
  #######################################################   
  function.plot.est.sdv	<- function(dates,coefficients,xlab="",main="",scale=0,sdv.at.zero=FALSE,intercept=FALSE,tau=FALSE,cex=NULL,cex.axis=NULL)
    #	BN 2 Dec 2013
  {	#	function.plot.est.sdv
    #	define function that can move to exponential scale
    function.scale	<- function(x,scale=0)
    {	if(scale==0)	x.scale <- x
    if(scale==1)	x.scale <- exp(x)
    if(scale==2)	x.scale <- exp(x)/(1+exp(x))
    return(x.scale)
    }
    ################
    #	IF MORE THAN ONE OBSERVATION
    if(length(dates)>1)
    {
      ################
      #	get estimates and sdv
      dat	<- as.vector(dates)
      est	<- as.vector(coefficients[,1])
      sdv <- as.vector(coefficients[,2])
      #	get center for sdv	
      sdv0	<- (1-sdv.at.zero)*est		
      ################
      #	set ylim
      y.lower	<- min(0,min(function.scale(est,scale)),max(2*min(function.scale(est,scale)),max(function.scale(sdv0-sdv,scale),na.rm=TRUE)))
      y.upper	<- max(0,max(function.scale(est,scale)),min(2*max(function.scale(est,scale)),min(function.scale(sdv0+sdv,scale),na.rm=TRUE)))
      if(max(est)-min(est)<min(sdv,na.rm=TRUE)/2)
        cat("apc.plot.fit warning: sdv large in for plot",i,"- possibly not plotted\n")
      ################
      
      #	plot
      plot(dat,function.scale(est,scale),type="l",ann=FALSE, axes=FALSE,ylim=c(y.lower,y.upper),cex.axis=cex.axis)
      if(intercept == FALSE) mtext(side = 1, text = xlab, line = mgp[1]  ,cex=cex.axis) 
      title(main=list(main, cex=cex), line=theight) 
      if(intercept == FALSE) axis(1, cex.axis=cex.axis, mgp = mgp)
      axis(2, las=2, cex.axis=cex.axis, at=round(c(y.lower, mean(c(y.lower, y.upper)), y.upper), 1))
      box()
      if(tau==FALSE)
      {	lines(dat,function.scale(sdv0+  sdv,scale),lty=2,col="blue" )
        lines(dat,function.scale(sdv0-1*sdv,scale),lty=2,col="blue" )
        lines(dat,function.scale(sdv0+2*sdv,scale),lty=3,col="red",lwd=2)
        lines(dat,function.scale(sdv0-2*sdv,scale),lty=3,col="red",lwd=2)
      }	
      abline(0,0)
    }
    ################
    #	IF ONLY ONE OBSERVATION
    if(length(dates)==1)
    {
      ################
      #	get estimates and sdv
      dat	<- dates
      est	<- coefficients[1]
      sdv <- coefficients[2]
      #	get center for sdv	
      sdv0	<- (1-sdv.at.zero)*est		
      ################
      #	set ylim
      y.lower	<- min(0,min(function.scale(est,scale)),max( 2*min(function.scale(est,scale)),max(function.scale(sdv0-sdv,scale),na.rm=TRUE) ))
      y.upper	<- max(0,max(function.scale(est,scale)),min( 2*max(function.scale(est,scale)),min(function.scale(sdv0-sdv,scale),na.rm=TRUE) ))
      ################
      #	plot
      plot(dat,function.scale(est,scale),type="p",ann=FALSE, axes=FALSE,xaxt=xaxt,ylim=c(y.lower,y.upper),xlim=c(dat-1,dat+1),pch=19,cex.axis=cex.axis)
      if (intercept == FALSE) mtext(side = 1, text = xlab, line = mgp[1]  ,cex=cex.axis) 
      title(main=list(main, cex=cex), line=theight) 
      if (intercept == FALSE) axis(1, cex.axis=cex.axis, mgp = mgp)
      axis(2, las=2, cex.axis=cex.axis, at=round(c(y.lower, mean(c(y.lower, y.upper)), y.upper), 1))
      if(tau==FALSE)
      {	points(dat,function.scale(sdv0+  sdv,scale),col="blue" )
        points(dat,function.scale(sdv0-1*sdv,scale),col="blue" )
        points(dat,function.scale(sdv0+2*sdv,scale),col="red")
        points(dat,function.scale(sdv0-2*sdv,scale),col="red")
      }
      abline(0,0)
    }	
  }	#	function.plot.est.sdv
  #############################################
  #	plots
  if(is.null(sub.plot)==TRUE)
  {
    par(mfrow=c(3,3)); par(mar=mar,oma=oma);
    if(is.null(cex)==TRUE) 		cex 		<- 1
    if(is.null(cex.axis)==TRUE) cex.axis 	<- 1		
    for(i in 1:9)
    {		
      if(v.do.plot[i]==TRUE)
        function.plot.est.sdv(l.dates[[i]],l.coefficients[[i]],xlab=v.xlab[i],main=v.main.sub[i],intercept=v.intercept[i],tau=v.tau[i],scale=scale,sdv.at.zero=sdv.at.zero,cex=cex,cex.axis=cex.axis)
      else
        frame()
    }		
    title(main=list(main, cex=cex.main), outer=TRUE)
  }
  else
  {	if(sub.plot=="a")	i <- 1
  if(sub.plot=="b")	i <- 2
  if(sub.plot=="c")	i <- 3
  if(sub.plot=="d")	i <- 4
  if(sub.plot=="e")	i <- 5
  if(sub.plot=="f")	i <- 6
  if(sub.plot=="g")	i <- 7
  if(sub.plot=="h")	i <- 8
  if(sub.plot=="i")	i <- 9
  par(mfrow=c(1,1));	par(mar=c(4,5,3,1),oma=c(0,0,0,0)); cex <- NULL
  main	<- main.sub
  if(is.null(main)==TRUE)	main <- v.main.sub[i]
  if(v.do.plot[i]==TRUE)
    function.plot.est.sdv(l.dates[[i]],l.coefficients[[i]],xlab=v.xlab[i],main=main,scale=scale,sdv.at.zero=sdv.at.zero,cex=cex)
  else
    return(cat("apc.plot.fit error: cannot draw this sub.plot. Check sub.plot is correct \n"))
  }
}	#	apc.plot.fit


##### alternative approach to estimating TS model (BN) #####
estimate.TS <- function(data, depvar, evars=NULL){
  # evars needs to be a c(, ,)
  
  # get n.group, ie the number of people in each ik cell
  # this will be needed when regressing on X_TS 
  get.n.group <-  ddply(data, .variables = c("cell.name"), 
                        function(dfr, colnm){length(dfr[,colnm])}, depvar)
  colnames(get.n.group)[2] <- "n.group"  
  
  # need to distinguish between cases with and without evars. First the case with.
  if(!is.null(evars)){
    
    # Step 1: Regress covariates on X_TS to get residuals and coefficient
    
    part1 <- get.n.group
    
    # begin by summing the covariate values for observations in each ik cell
    for (i in evars){
      get.sum.c2 <- ddply(data, .variables = c("cell.name"),
                          function(dfr, colnm){sum(dfr[, colnm])}, i)
      colnames(get.sum.c2)[2] <- paste("groupsum", i, sep=".")
      part1 <- join(part1, get.sum.c2, by="cell.name")
    }
    
    # now use the sums to get the mean covariate values in each ik cell
    groupmeans.matrix <- as.data.frame(matrix(nrow=nrow(part1), ncol=length(evars)))
    for (i in 1:length(evars)){
      cbar <- as.matrix(part1[, 2+i] / part1[,2]) # this divides each groupsum by n.group
      groupmeans.matrix[, i] <- cbar
      colnames(groupmeans.matrix)[i] <- paste("groupmean", evars[i], sep = ".")
    }
    
    # this matrix of group means is in fact equal to phi hat, see the matrix expression
    # for kappa hat in section 4.3 of the thesis
    
    # because X_TS is a matrix of unit vectors, the predicted value of C is just the
    # appropriate element of the estimated coefficient vector: C_hat is a set of phi_hats
    
    part2 <- cbind(part1, groupmeans.matrix)
    data.and.first.coeff <- join(data, part2, by="cell.name")
    # this now contains all the data, the vector of group counts, and the matrix of 
    # phi-hats
    
    # now construct the residuals of C based on X_TS ie C-tilde
    residuals.matrix <- as.data.frame(matrix(nrow=nrow(data.and.first.coeff), ncol=length(evars)))
    for (i in 1:length(evars)){
      true <- data.and.first.coeff[evars[i]]
      predicted.name <- paste("groupmean", evars[i], sep=".")
      predicted <- data.and.first.coeff[predicted.name]
      resid <- true - predicted
      residuals.matrix[, i] <- resid
      colnames(residuals.matrix)[i] <- paste("resid", evars[i], sep=".")
    }
    
    # this contains all the data, the group counts, the matrix of phi-hats, the matrix
    # of C-tildes
    data.resid.c.on.x <- cbind(data.and.first.coeff, residuals.matrix)  
    
    # Step 2: regress Y on C-tildes
    stage2.depvar <- data.resid.c.on.x[depvar]
    stage2.evar <- residuals.matrix
    stage2.regress <- glm.fit(as.matrix(stage2.evar), as.matrix(stage2.depvar), family = gaussian())
    truelambdahat <- stage2.regress$coefficients
    # here truelambdahat is the lambda estimator as described in thesis
    
    # Step 3: regress Y on X_TS
    # Begin by summing Y for observations in each ik cell
    get.sum.y <-  ddply(data.resid.c.on.x, .variables = c("cell.name"), 
                        function(dfr, colnm){sum(dfr[,colnm])}, depvar)
    colnames(get.sum.y)[2] <- paste("groupsum", depvar, sep=".")    
    
    # Construct rhohat
    ysum.and.count <- join(get.n.group, get.sum.y, by="cell.name")
    truerhohat <- ysum.and.count[, 3]/ysum.and.count[, 2]
    yinfo.and.rho <- cbind(ysum.and.count, truerhohat)
    
    # Step 4: use all existing estimated coefficient vectors to construct kappa-hat
    conformable.lambda <- (as.matrix(as.numeric(truelambdahat)))
    conformable.phi <- as.matrix(groupmeans.matrix)
    phi.times.lambda <- conformable.phi %*% conformable.lambda
    
    truekappahat <- truerhohat - phi.times.lambda
    
    # associate truekappahat with cell.name
    yinfo.and.kappa <- cbind(yinfo.and.rho, truekappahat)
    # join data matrix and truekappahat, so that we have predicted y from kappahat
    data.all.coeffs <- join(data.resid.c.on.x, yinfo.and.kappa, by="cell.name")
    
    # Step 5: calculate residuals
    predicted.y.from.evars <- as.matrix(data.all.coeffs[evars]) %*% as.matrix(truelambdahat)
    y.residuals <- data.all.coeffs[depvar] - predicted.y.from.evars - data.all.coeffs["truekappahat"]
    colnames(y.residuals) <- "y.residuals"
    
    # Alternative: y.residuals without kappahat
    predicted.y.from.ctilde <- as.matrix(residuals.matrix) %*% as.matrix(truelambdahat)
    y.residuals.nokappa <- data.all.coeffs[depvar] - predicted.y.from.ctilde - data.all.coeffs["truerhohat"]
    colnames(y.residuals.nokappa) <- "y.residuals.nokappa"
    
    final.data <- cbind(data.all.coeffs, y.residuals, y.residuals.nokappa)
    
  } else{
    # need to deal with the case where there are no covariates. Basically just regress
    # y on TS - ie use only rhohat
    
    # construct sums of y for each ik cell
    get.sum.y <-  ddply(data, .variables = c("cell.name"), 
                        function(dfr, colnm){sum(dfr[,colnm])}, depvar)
    colnames(get.sum.y)[2] <- paste("groupsum", depvar, sep=".")    
    
    # Construct rhohat
    ysum.and.count <- join(get.n.group, get.sum.y, by="cell.name")
    truerhohat <- ysum.and.count[, 3]/ysum.and.count[, 2]
    yinfo.and.rho <- cbind(ysum.and.count, truerhohat)
    
    # Join data and rho by cell.name
    data.all.coeffs <- join(data, yinfo.and.rho, by="cell.name")
    y.residuals <- data.all.coeffs[depvar] - data.all.coeffs["truerhohat"]
    colnames(y.residuals) <- "y.residuals"
    
    final.data <- cbind(data.all.coeffs, y.residuals)
    
    truelambdahat <- NA
    truekappahat <- NA
  }
  
  # Return valuables
  
  valuables <- list(final.data = final.data,
                    TSS = sum(final.data[depvar]^2),
                    trueRSS = sum(final.data["y.residuals"]^2),
                    truelambdahat = truelambdahat,
                    truerhohat = truerhohat,
                    truekappahat= truekappahat)
}

##### alternative write ftest.apc.vs.ts #####
alt.ftest.apc.vs.ts <- function(data, datatype, model.design="APC", depvar, covariates=NULL,
                            model.family="gaussian", unit=1, n.coh.excl.start=0, n.coh.excl.end = 0,
                            existing.apc.fit.indiv=NULL, existing.TS=NULL){
  
  #fit the two models - continuous as DV, dichot as EV
  #first the APC one
  if(is.null(existing.apc.fit.indiv)){
    collinear <- apc.get.design.indiv.collinear(data, datatype=datatype, unit = unit,
                                                n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end)
    design <- apc.get.design.indiv.model(collinear, model.design = model.design, dep.var=depvar, covariates=covariates)
    fitted <- apc.fit.indiv(design, model.family=model.family, n.coh.excl.start = n.coh.excl.start, 
                            n.coh.excl.end = n.coh.excl.end)
  } else 
    fitted <- existing.apc.fit.indiv
  
  
  TSS.APC <- sum(fitted$y^2)
  RSS.APC <- sum(fitted$residuals^2)
  ESS.APC <- sum(fitted$fitted.values^2)
  
  nparam.APC <- length(fitted$coefficients)
  sample.N.APC <- length(fitted$y)
  
  #second the dummy one
  if(is.null(existing.TS)){
    dummymodel <- estimate.TS(data, depvar = depvar, evars = covariates)
  } else
    dummymodel <- existing.TS
  
  RSS.dummy <- dummymodel$trueRSS
  if (!is.na(dummymodel$truelambdahat[1])){
    nparam.dummy <- length(unique(dummymodel$final.data$cell.name))+length(dummymodel$truelambdahat)
  } else{
    nparam.dummy <- length(unique(dummymodel$final.data$cell.name))
  }
  
  #elements of F-test
  big.model.RSS <- RSS.dummy
  small.model.RSS <- RSS.APC
  big.model.nparams <- nparam.dummy
  small.model.nparams <- nparam.APC
  sample.N <- sample.N.APC
  
  df.num <- big.model.nparams - small.model.nparams
  df.denom <- sample.N - big.model.nparams
  df <- paste("(", paste(as.character(df.num), as.character(df.denom), sep=", "), ")", sep="")
  
  fstat.num <- (small.model.RSS - big.model.RSS)/df.num
  fstat.denom <- (big.model.RSS - 0)/df.denom
  fstat <- fstat.num/fstat.denom
  
  #get pvalue
  p.value <- pf(fstat, df1 = df.num, df2 = df.denom, lower.tail = FALSE)
  
  #valuables
  valuables <- list(fstat = fstat,
                    df = df,
                    df.num = df.num,
                    df.denom = df.denom,
                    p.value = p.value,
                    aic.APC = fitted$aic)
  return(valuables)
}

##### alternative write Ftest.comparison.table #####
alt.Ftest.comparison.table <- function(data, datatype, depvar, covariates=NULL, unit=1, 
                                   n.coh.excl.start=0, n.coh.excl.end=0){
  get.table.line <- function(ftest1, ftest2, first=FALSE){
    if(isTRUE(first))
      line <- c(round(ftest1$fstat, 3), ftest1$df.num, round(ftest1$p.value,3), NA,
                NA, NA, round(ftest1$aic.APC,3))
    else
      line <- c(round(ftest1$fstat, 3), ftest1$df.num, round(ftest1$p.value,3), round(ftest2$fstat, 3),
                ftest2$df.num, round(ftest2$p.value, 3), round(ftest2$aic.R,3))
    return(line)
  }
  
  model.design.list	<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC", "1")
  
  fit.tab <- matrix(nrow=length(model.design.list),ncol=7,data=NA)
  
  fullAPCcollinear <- apc.get.design.indiv.collinear(data, datatype=datatype, unit=unit,
                                                     n.coh.excl.start=n.coh.excl.start, n.coh.excl.end=n.coh.excl.end)
  fullAPCdesign <- apc.get.design.indiv.model(fullAPCcollinear, dep.var = depvar, covariates = covariates)
  fullAPCfit <- apc.fit.indiv(fullAPCdesign, model.family="gaussian", n.coh.excl.start = n.coh.excl.start,
                              n.coh.excl.end = n.coh.excl.end)
  TSfit <- estimate.TS(data, depvar=depvar, evars=covariates)
  
  ftest.full.vs.TS <- alt.ftest.apc.vs.ts(existing.apc.fit.indiv = fullAPCfit, existing.TS = TSfit)
  first.line <- get.table.line(ftest1 = ftest.full.vs.TS, first=TRUE)
  fit.tab[1, 1:7] <- first.line
  
  for (i in 2:length(model.design.list)){
    againstTS <- alt.ftest.apc.vs.ts(data, datatype = datatype, model.design=model.design.list[i],
                                 depvar = depvar,covariates = covariates, unit = unit,
                                 n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end,
                                 existing.TS = TSfit)
    againstAPC <- ftest.apcreduced.vs.apc(data, datatype=datatype, reduced.model.design = model.design.list[i],
                                          depvar=depvar, covariates=covariates, unit = unit, 
                                          n.coh.excl.start = n.coh.excl.start, n.coh.excl.end = n.coh.excl.end,
                                          existing.apc.fit.indiv = fullAPCfit)
    fit.tab[i, 1:7] <- get.table.line(againstTS, againstAPC)
  }
  
  df.againstTS <- paste("DF ( * , ", as.character(againstTS$df.denom), ")", sep="")
  df.againstAPC <- paste("DF ( * , ", as.character(againstAPC$df.denom), ")", sep="")
  
  table.colnames <- c("F-test vs TS", df.againstTS, "p-value", "F-test vs APC", 
                      df.againstAPC, "p-value", "AIC")
  
  colnames(fit.tab) <- table.colnames
  rownames(fit.tab) <- model.design.list
  
  return(fit.tab)
}
