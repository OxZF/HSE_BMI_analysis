# --------------------------
# DPhil
# revising.functions.R
# 13/09/16
# --------------------------

# script to investigate existing functions and look
# into revising them.

# --------------------------

setwd("C:/Users/me/Documents")

#### Access the basic material I need ####

# write functions as they currently are
source(file = "DPhil/HSE_BMI_analysis/Functions/write.functions.v2.R")

# call in all artificial data
d.AP <- read.csv("DPhil/HSE_BMI_analysis/Inputs/Artificial.data/artificial_individ_AP.csv")
d.AP.cc <- read.csv("DPhil/HSE_BMI_analysis/Inputs/Artificial.data/artificial_individ_AP_cohscut.csv")
d.CP <- read.csv("DPhil/HSE_BMI_analysis/Inputs/Artificial.data/artificial_individ_CP.csv")
d.AC.r <- read.csv("DPhil/HSE_BMI_analysis/Inputs/Artificial.data/artificial_individ_AC.rect.csv")
d.AC.t <- read.csv("DPhil/HSE_BMI_analysis/Inputs/Artificial.data/artificial_individ_AC.tri.csv")

#### 1: cell.name introduction ####
cell.name <- paste("ij", as.character(educdata$age), as.character(educdata$period), sep="_")
educdata <- cbind(cell.name, educdata)

design.AP.collin <- apc.get.design.indiv.collinear(d.AP, datatype="AP")
View(design.AP.collin$full.design.collinear)

f.AP <- Ftest.comparison.table(d.AP, datatype="AP")   # throws cell.name error

cell.name <- paste("ij", as.character(d.AP$age), as.character(d.AP$period), sep="_")
d.AP.wcn <- cbind(cell.name, d.AP)

f.AP <- Ftest.comparison.table(d.AP.wcn, datatype = "AP")   # throws glm.fit error, not clear that problem is missing depvar, revise
f.AP <- Ftest.comparison.table(d.AP.wcn, datatype = "AP", depvar = "continuous")

# try ik cell name (should work fine)
rm(cell.name, f.AP, d.AP.wcn)
cell.name <- paste("ik", as.character(d.AP$age), as.character(d.AP$cohort), sep="_")
d.AP.wcn2 <- cbind(cell.name, d.AP)

f.AP2 <- Ftest.comparison.table(d.AP.wcn, datatype = "AP", depvar = "continuous")

# show equal outcome using ik/ij 
all(f.AP == f.AP2, na.rm=TRUE)

# add cell.name internally to apc.get.design.indiv.collinear and investigate
design.AP.collin2 <- apc.get.design.indiv.collinear(d.AP, datatype="AP")
View(design.AP.collin2$full.design.collinear)
# check all elements of both matrices (aside from cell.name line) are equal
all(design.AP.collin2$full.design.collinear[, -1]==design.AP.collin$full.design.collinear)  # yep, they are :)

f.AP3 <- Ftest.comparison.table(d.AP, datatype = "AP", depvar = "continuous") # doesn't work. Problem inside model.DVcont.manyEVcont.TS
# probs will also need to add to estimate.TS but one step at a time...

# add cell.name internally to model.DVcont.manyEVcont using this structure:
# if (!"cell.name" %in% colnames(data)) {
#   cell.name <- paste("ik", as.character(data$age), as.character(data$cohort), sep="_")
#   data <- cbind(cell.name, data)
# }
# also changed it inside apc.get.design.indiv.collinear to this

f.AP3 <- Ftest.comparison.table(d.AP, datatype = "AP", depvar = "continuous") # works
f.AP4 <- Ftest.comparison.table(d.AP.wcn, datatype = "AP", depvar = "continuous") # works
f.AP5 <- Ftest.comparison.table(d.AP.wcn2, datatype = "AP", depvar = "continuous") # works

all(f.AP3==f.AP5, na.rm=TRUE) # checked w/ variety of numbers and all good

design.AP.collin3 <- apc.get.design.indiv.collinear(d.AP, datatype = "AP")
design.AP.collin4 <- apc.get.design.indiv.collinear(d.AP.wcn, datatype = "AP")

all(design.AP.collin2$structure.design.collinear == design.AP.collin3$structure.design.collinear) # 2, 3, 4 all the same

f.AP6 <- alt.Ftest.comparison.table(d.AP, datatype = "AP", depvar="continuous") # so this is where estimate.TS comes in. 
# Add same cell.name construction to estimate.TS

f.AP6 <- alt.Ftest.comparison.table(d.AP, datatype = "AP", depvar="continuous") # so this is where estimate.TS comes in. 

all(f.AP6==f.AP5, na.rm=TRUE) # all equal!

# Now want to get it to throw a warning if cell.name already defined using structure below
data <- d.AP
data <- d.AP.wcn

# if (!"cell.name" %in% colnames(data)) {
#   cell.name <- paste("ik", as.character(data$age), as.character(data$cohort), sep="_")
#   data <- cbind(cell.name, data)
# } else {
#   warning("Variable cell.name must be unique identifier of age-cohort cells; if preexisting variable cell.name has alternate meaning please rename")
# }
 
# and now check how it runs in code
design.AP.collin5 <- apc.get.design.indiv.collinear(d.AP, datatype="AP")
design.AP.collin6 <- apc.get.design.indiv.collinear(d.AP.wcn, datatype="AP")

f.AP7 <- Ftest.comparison.table(d.AP, datatype = "AP", depvar = "continuous") # works
f.AP8 <- Ftest.comparison.table(d.AP.wcn, datatype = "AP", depvar = "continuous") # works
f.AP9 <- alt.Ftest.comparison.table(d.AP, datatype = "AP", depvar = "continuous") # works
f.AP10 <- alt.Ftest.comparison.table(d.AP.wcn, datatype = "AP", depvar = "continuous") # works
# all works - but we are generating 30 warnings! Which means we're running 30 internal estimations - that can't be efficient.

#### 2: addressed case of using wrong datatype or getting n.coh.excl wrong ####

#### 3: not specifying dep.var ####
# dropped some bloating in apc.get.design.indiv.model output
# warning in apc.get.design.indiv.model if dep.var not specified, option to specify dep.var in apc.fit.indiv
# don't specify DV anywhere
design.CP.collin <- apc.get.design.indiv.collinear(d.CP, datatype = "CP")
design.CP.model <- apc.get.design.indiv.model(design.CP.collin, model.design="APC")
fit.CP.model <- apc.fit.indiv(design.CP.model, model.family="gaussian")
# specify DV in fit
design.CP.model2 <- apc.get.design.indiv.model(design.CP.collin, model.design="APC")
regressand <- as.matrix(d.CP$continuous)
fit.CP.model2 <- apc.fit.indiv(design.CP.model2, model.family="gaussian", DV = regressand)
# specify DV in design
design.CP.model3 <- apc.get.design.indiv.model(design.CP.collin, model.design = "APC", dep.var="continuous")
fit.CP.model3 <- apc.fit.indiv(design.CP.model3, model.family="gaussian")
# specify same DV both fit and design
fit.CP.model4 <- apc.fit.indiv(design.CP.model3, model.family="gaussian", DV=regressand)

regressand2 <- as.matrix(rnorm(96, 0, 1))
# specify different DV fit, old one design
fit.CP.model5 <- apc.fit.indiv(design.CP.model3, model.family="gaussian", DV=regressand2)
# specify different DV fit only
fit.CP.model6 <- apc.fit.indiv(design.CP.model2, model.family="gaussian", DV=regressand2)

all(fit.CP.model2$fitted.values==fit.CP.model3$fitted.values) # true, good
all(fit.CP.model2$fitted.values==fit.CP.model4$fitted.values) # true, good
all(fit.CP.model2$fitted.values==fit.CP.model5$fitted.values) # false, good
all(fit.CP.model5$fitted.values==fit.CP.model6$fitted.values) # true, good

rm(list=ls(pattern = "fit.CP.model*"))
rm(list=ls(pattern = "design.CP*"))
rm(list=ls(pattern="regressand*"))

#### 4: model.design/model.family specification ####
design.ACt.collin <- apc.get.design.indiv.collinear(d.AC.t, datatype = "AC triangular")
design.ACt.model <- apc.get.design.indiv.model(design.ACt.collin, model.design="APC", dep.var="continuous")
fit.ACt.model <- apc.fit.indiv(design.ACt.model, model.family="gaussian")

design.ACt.model2 <- apc.get.design.indiv.model(design.ACt.collin, model.design="cat", dep.var="continuous")
# OK so it should through an "invalid model design" error 
# and now it does!

fit.ACt.model2 <- apc.fit.indiv(design.ACt.model, model.family="turkey")
# and now this throws a family error. Good

# next step is to check this happens in all the Ftest functions as well.