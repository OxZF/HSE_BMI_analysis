# apc.fit.indiv taken from write.functions.v2 08.08.16 in attempt to trace SE


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

