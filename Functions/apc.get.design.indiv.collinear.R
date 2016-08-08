# apc.get.design.indiv.collinear taken from write.functions.v2 08.08.16


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
