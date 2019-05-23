packages <- list('nlme', 'clubSandwich', 'reghelper')
loaded <- lapply(packages, library, character.only = TRUE)

# #
# load('~/data/SEA/fMRI_sensitivity/debug.Rdata')
# attach(designmat)
# v <- 1e4
# #

processVoxel <- function(v) {
  BRAIN <- voxeldat[,v]
  NOVAL <- 999
  retvals <- numeric()
  retnames <- c()
  
  # Some effects may be allowed to vary within-person (anything
  # for which we have multiple observations per person):
  # - Linear effect of time, TIMECENTER
  # - Linear effect of chronic stress severity deviations, WGEN_CHRONICSEV
  # - Their interaction, WCEN_CHRONICSEV:TIMECENTER
  model_data <- data.frame(BRAIN = BRAIN,
                           TIMECENTER = TIMECENTER,
                           GCEN_CHRONICSEV = GCEN_CHRONICSEV,
                           WCEN_CHRONICSEV = WCEN_CHRONICSEV,
                           idnum = idnum)
  model_formula <- BRAIN ~ TIMECENTER + GCEN_CHRONICSEV + WCEN_CHRONICSEV
  e <- try( CHRONIC_BETWEEN_WITHIN2 <-
              lme(model_formula,
                  random=~1 | idnum, data = model_data, 
                  method = "ML", na.action=na.omit) )
  
  e_sandwich <- try(CHRONIC_BETWEEN_WITHIN2_sw <- coef_test(e, vcov = 'CR2'))
  
  if (inherits(e, "try-error")){ 
    message("error thrown at voxel ", v)
    message(e);
    CHRONIC_BETWEEN_WITHIN2  <- NULL 
  }
  if (inherits(e_sandwich, "try-error")){
    message("can't compute corrected standard errors", v)
    message(e_sandwich);
    CHRONIC_BETWEEN_WITHIN2_sw <- NULL
  }
  
  name_col_correspondence <- c('-est' = 'Value', '-t' = 't-value', '-p' = 'p-value')
  name_col_correspondence_sw <- c('-t_sw' = 'tstat', '-df_sw' = 'df', '-p_satt_sw' = 'p_Satt')
  
  if (!(is.null(CHRONIC_BETWEEN_WITHIN2 ))) {
    coefs <- coef(summary(CHRONIC_BETWEEN_WITHIN2))
    reg_retvals <- unlist(lapply(name_col_correspondence, function(colname){
      if(colname == name_col_correspondence[[3]]){
        return(1-coefs[, colname])
      } else {
        return(coefs[, colname])
      }
    }))
    names(reg_retvals) <- gsub('[\\(\\)]', '', gsub('(.*?)\\.(.*)', '\\2\\1', names(reg_retvals)))
    
    std_coefs <- coef(beta(CHRONIC_BETWEEN_WITHIN2))
    std_reg_retvals <- std_coefs[, 'Value']
    names(std_reg_retvals) <- paste0(gsub('\\(*(.*?)\\)*\\.*\\z*$', '\\1', names(std_reg_retvals)), '-beta')
    
  } else {
    model_terms <- c('Intercept', attr(terms(model_formula), 'term.labels'))
    reg_retvals <- rep(NOVAL, length(model_terms) * length(name_col_correspondence))
    names(reg_retvals) <- paste0(rep(model_terms, length(name_col_correspondence)), 
                                 rep(names(name_col_correspondence), each = length(model_terms)))
    
    model_terms <- c('Intercept', attr(terms(model_formula), 'term.labels'))
    std_reg_retvals <- rep(NOVAL, length(model_terms))
    names(std_reg_retvals) <- paste0(model_terms, '-beta')
  }
  if(!(is.null(CHRONIC_BETWEEN_WITHIN2_sw))){
    coefs_sw <- CHRONIC_BETWEEN_WITHIN2_sw
    reg_retvals_sw <- unlist(lapply(name_col_correspondence_sw, function(colname){
      if(colname == name_col_correspondence_sw[[3]]){
        rval <- 1-coefs_sw[, colname] 
      } else {
        rval <- coefs_sw[, colname]
      }
      names(rval) <- rownames(coefs_sw)
      return(rval)
    }))
    names(reg_retvals_sw) <- gsub('[\\(\\)]', '', gsub('(.*?)\\.(.*)', '\\2\\1', names(reg_retvals_sw)))
  } else {
    
  }
  retvals <- c(reg_retvals, std_reg_retvals)
  return(retvals)
}
