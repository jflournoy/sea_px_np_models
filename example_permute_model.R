# #Uncomment to test
# load('../chavg.fear/debug.Rdata')
# attach(designmat)
# permutationNumber <- 1
# nperm <- 1000
# v=10000

packages <- list('nlme', 'clubSandwich', 'permute', 'neuropointillist')
loaded <- lapply(packages, library, character.only = TRUE)

processVoxel <- function(v) {
  #This processVoxel function will produce the group map for 
  #a single permutation.
  
  permutationRDS = 'permutation_set-free.RDS'
  targetDV = 'WCEN_CHRONICAVG'
  formula = BRAIN ~ 1 + TIMECENTER + GCEN_CHRONICAVG + WCEN_CHRONICAVG
  
  BRAIN <- voxeldat[,v]
  NOVAL <- 999
  retvals <- numeric()
  model_data <- data.frame(
    BRAIN = BRAIN, 
    TIMECENTER = TIMECENTER, 
    GCEN_CHRONICAVG = GCEN_CHRONICAVG, 
    WCEN_CHRONICAVG = WCEN_CHRONICAVG, 
    idnum = idnum)
  
  if(v%%5e3 == 0){
    apercent <- sprintf("%3.0f%%", v/dim(voxeldat)[2]*100)
    status_message <- paste0(apercent, ": Voxel ", v, " of ", dim(voxeldat)[2], "...")
    message(status_message)
  }

  # `method = REML` for unbiased estimate of variance parameters.
  # See: 
  # Luke, S. G. (2017). Evaluating significance in linear mixed-effects
  # models in R. Behavior Research Methods, 49(4), 1494–1502. 
  # https://doi.org/10.3758/s13428-016-0809-y
  
  permuteModel <- neuropointillist::npointLmePermutation(permutationNumber = permutationNumber, 
                                                         permutationRDS = permutationRDS,
                                                         targetDV = targetDV,
                                                         z_sw = TRUE, vcov = 'CR2',
                                                         formula = formula,
                                                         random = ~1 | idnum,
                                                         data = model_data,
                                                         lmeOpts = list(method = "REML", na.action=na.omit))
  if(is.null(permuteModel$Z)){
    Z <- NOVAL
  } else {
    Z <- permuteModel$Z
  }
  names(Z) <- paste0(targetDV,'-z_sw')
  return(Z)
}

# (timetime <- system.time({
# perms <- lapply(1:1000, function(p){
#   permutationNumber <<- p
#   aperm <- processVoxel(1e5+1)
# })
# }))
