# var.apc.plot.fit as taken from write.functions.v2 08.08.16 to establish how SE are plotted
# will also want to work out how detrending is achieved



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