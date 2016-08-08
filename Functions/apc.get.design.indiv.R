# apc.get.design.indiv.model taken from write.functions.v2 08.08.16


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