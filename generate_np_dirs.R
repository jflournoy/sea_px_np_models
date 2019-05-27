
library(argparse)

create_np_model_dir <- function(basedir, model_name, overwrite = FALSE){
  if(dir.exists(basedir)){
    model_dir <- file.path(basedir, model_name)
    if(dir.exists(model_dir) & !overwrite){
      stop('Model directory already exists. See help.')
    } else if(dir.exists(model_dir) & overwrite) {
      message('Overwriting model files in ', model_dir, '...')
      message('Note, this does not _delete_ any other files!')
    } else {
      dir.create(model_dir)  
    }
  } else {
    stop('Directory does not exist: ', basedir)
  }
  return(model_dir)
}

make_model_header <- function(y, x, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  
  x_b_name <- paste0(between_prefix, x)
  x_w_name <- paste0(within_prefix, x)
  
  header <- '
packages <- list(\'nlme\', \'clubSandwich\', \'reghelper\')
loaded <- lapply(packages, library, character.only = TRUE)

processVoxel <- function(v) {
  BRAIN <- voxeldat[,v]
  NOVAL <- 999
  retvals <- numeric()'
  
  if(x == 'BRAIN'){
    header <- paste0(header,'
  brain_df <- data.frame(BRAIN, idnum) #idnum is defined by neuropoint
  
  # Compute Within-person Mean of brain data
  win_mean <- aggregate(brain_df, by = list(idnum), FUN = mean)[,1:2]
  colnames(win_mean) <- c("idnum", "WIN_BRAIN_MEAN")
  brain_df <- merge(brain_df, win_mean, by="idnum")

  grandmean <- mean(ag$win_mean)

  # This is a bit of an extra step... Tara had this here so I\'m leaving it
  brain_df[, "', x_b_name,'"] <- brain_df$WIN_BRAIN_MEAN - grandmean
  brain_df[, "', x_w_name,'"] <- BRAIN - brain_df$WIN_BRAIN_MEAN
  
  # attach to these new variables
  ', x_b_name,' <- brain_df[, "', x_b_name,'"]
  ', x_w_name,' <- brain_df[, "', x_w_name,'"]
') 
  }
  allvars <- c(y, covariates, x_b_name, x_w_name, id_var)
  data_frame_args <- paste(
    paste0(c(y, covariates, x_b_name, x_w_name), ' = ', c(y, covariates, x_b_name, x_w_name)),
    collapse = ', \n    ')
  header <- paste0(header,'
  model_data <- data.frame(\n    ', data_frame_args,')')
  return(header)
}

make_lme_model_syntax <- function(y, x, model_name, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  x_b_name <- paste0(between_prefix, x)
  x_w_name <- paste0(within_prefix, x)
  model_name_sw <- paste0(model_name, '_sw')
  
  model_formula <- paste0(y, ' ~ ', paste(covariates, collapse = ' + '), ' + ', x_b_name, ' + ', x_w_name)
  model_syntax <- paste0('
  model_formula <- ', model_formula, '

  # `method = REML` for unbiased estimate of variance parameters.
  # See: 
  # Luke, S. G. (2017). Evaluating significance in linear mixed-effects
  # models in R. Behavior Research Methods, 49(4), 1494–1502. 
  # https://doi.org/10.3758/s13428-016-0809-y
  e <- try(', model_name, ' <-
              nlme::lme(model_formula,
                        random = ~1 | ', id_var, ', data = model_data, 
                        method = "REML", na.action=na.omit) )

  # Compute cluster corrected standard errors to account for, e.g., residual autocorrelation.
  # See coef_test help for description of "CR2"
  e_sandwich <- try(', model_name_sw, ' <- clubSandwich::coef_test(e, vcov = "CR2"))
')
  return(model_syntax)
}

make_model_footer <- function(y, x, model_name, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  x_b_name <- paste0(between_prefix, x)
  x_w_name <- paste0(within_prefix, x)
  model_name_sw <- paste0(model_name, '_sw')
  
  footer <- paste0('
  if (inherits(e, "try-error")){ 
    message("error thrown at voxel ", v)
    message(e);
    ', model_name, ' <- NULL 
  }
  if (inherits(e_sandwich, "try-error")){
    message("can\'t compute corrected standard errors at voxel ", v)
    message(e_sandwich);
    ', model_name_sw, ' <- NULL
  }
  name_col_correspondence <- c("-est" = "Value", "-t" = "t-value", "-p" = "p-value")
  name_col_correspondence_sw <- c("-t_sw" = "tstat", "-df_sw" = "df", "-p_satt_sw" = "p_Satt")
  if (!(is.null(', model_name, '))) {
    coefs <- coef(summary(', model_name, '))
    reg_retvals <- unlist(lapply(name_col_correspondence, function(colname){
      if(colname == name_col_correspondence[[3]]){
        return(1-coefs[, colname])
      } else {
        return(coefs[, colname])
      }
    }))
    names(reg_retvals) <- gsub("[\\(\\)]", "", gsub("(.*?)\\.(.*)", "\\2\\1", names(reg_retvals)))
    
    std_coefs <- coef(beta(', model_name, '))
    std_reg_retvals <- std_coefs[, "Value"]
    names(std_reg_retvals) <- paste0(gsub("\\(*(.*?)\\)*\\.*\\z*$", "\\1", names(std_reg_retvals)), "-beta")
  } else {
    model_terms <- c("Intercept", attr(terms(model_formula), "term.labels"))
    reg_retvals <- rep(NOVAL, length(model_terms) * length(name_col_correspondence))
    names(reg_retvals) <- paste0(rep(model_terms, length(name_col_correspondence)), 
                                 rep(names(name_col_correspondence), each = length(model_terms)))
    std_reg_retvals <- rep(NOVAL, length(model_terms))
    names(std_reg_retvals) <- paste0(model_terms, "-beta")
  }
  if(!(is.null(', model_name_sw, '))){
    coefs_sw <- ', model_name_sw, '
    reg_retvals_sw <- unlist(lapply(name_col_correspondence_sw, function(colname){
      if(colname == name_col_correspondence_sw[[3]]){
        rval <- 1-coefs_sw[, colname] 
      } else {
        rval <- coefs_sw[, colname]
      }
      names(rval) <- rownames(coefs_sw)
      return(rval)
    }))
    names(reg_retvals_sw) <- gsub("[\\(\\)]", "", gsub("(.*?)\\.(.*)", "\\2\\1", names(reg_retvals_sw)))
  } else {
    model_terms_sw <- c("Intercept", attr(terms(model_formula), "term.labels"))
    reg_retvals_sw <- rep(NOVAL, length(model_terms_sw) * length(name_col_correspondence_sw))
    names(reg_retvals_sw) <- paste0(rep(model_terms_sw, length(name_col_correspondence_sw)), 
                                    rep(names(name_col_correspondence_sw), each = length(model_terms_sw)))
  }
  retvals <- c(reg_retvals, std_reg_retvals, reg_retvals_sw)
  return(retvals)
}')
  return(footer)
}

make_readargs <- function(mask_fname, set1, setlabels, model_file, testvoxel = '1000', output,
                          debugfile = 'debug.Rdata', hpc = 'slurm', jobs = '40'){
  output <- file.path(output, paste0(output, '.'))
  setlabel_args <- paste0('"--setlabels', 1:length(setlabels), '", "', setlabels, '", ')
  readargs <- paste0('
cmdargs <- c("-m","', mask_fname, '", "--set1", "', set1, '",
             ', setlabel_args, '
             "--model", "', model_file, '",
             "--testvoxel", "', testvoxel, '",
             "--output", "', output, '",
             "--debugfile", "', debugfile, '",
             "--', hpc, 'N", "', jobs, '")
')
  return(readargs)
}

write_model_script <- function(model_dir, y, x, model_name, 
                               within_prefix = 'WCEN_', between_prefix = 'GCEN_', 
                               covariates = NULL, id_var = 'idnum',  overwrite = FALSE){
  header <- make_model_header(y = y, x = x, 
                              within_prefix = within_prefix, between_prefix = between_prefix, 
                              covariates = covariates, id_var = id_var)
  model_syntax <- make_lme_model_syntax(y = y, x = x, model_name = model_name, 
                                        within_prefix = within_prefix, between_prefix = between_prefix, 
                                        covariates = covariates, id_var = id_var)
  footer <- make_model_footer(y = y, x = x, model_name = model_name, 
                              within_prefix = within_prefix, between_prefix = between_prefix, 
                              covariates = covariates, id_var = id_var)
  model_text <- paste(header, model_syntax, footer, sep = '\n')
  model_file <- file.path(model_dir, 'model.R')
  message('Model file: ', model_file)
  if(file.exists(model_file) & !overwrite){
    stop('Model file exists. See help.')
  } else {
    f <- file(model_file, open = 'w')
    writeLines(model_text, f)
    close(f)
  }
  return(model_file)
}

move_mask <- function(mask, model_dir){
  message('Copying ', mask, ' to ', model_dir)
  mask_file <- file.path(model_dir, basename(mask))
  if(!file.exists(mask)){
    stop('Mask file does not exist')
  } else {
    
    file.copy(mask, mask_file)
  }
  return(mask_file)
}

parser <- ArgumentParser(description='Create a neuropoint model directory to be run')
parser$add_argument('base_dir', type="character",
                    help='This is where the model directory will be created')
parser$add_argument('model_name', type="character", help='Name of the model')
parser$add_argument('--IV', type="character",
                    help = 'IV variable name. Use "BRAIN" to specify BOLD contrast as the IV',
                    required = TRUE)
parser$add_argument('--DV', type="character",
                    help = 'DV variable name. Use "BRAIN" to specify BOLD contrast as the DV',
                    required = TRUE)
parser$add_argument('--mask', type="character", help='Path to mask (will be copied to model directory)', required = TRUE)
parser$add_argument('--set1', type="character", help='Path to input set of imaging files', required = TRUE)
parser$add_argument('--setlabels', type="character", help='Paths to imaging file labels (and covariates)', 
                    nargs = '+', required = TRUE)
parser$add_argument('--output', type="character", help='String will become "string/string." prefix to NP output', required = TRUE)
parser$add_argument('--overwrite', action='store_true', help = 'Overwrite existing model files?')
parser$add_argument('--win_pre', type="character", help='Within-person IV name prefix', default = 'WCEN_')
parser$add_argument('--bw_pre', type="character", help='Between-person IV name prefix', default = 'GCEN_')
parser$add_argument('--id_var', type="character", help='Column name of id variable', default = 'idnum')
parser$add_argument('--covariates', type="character", nargs = '+', help='Column names of covariates', default = NULL)
parser$add_argument('--testvoxel', type="character", help='Test voxel number', default = '10000')
parser$add_argument('--job', type="character", help='Number of jobs', default = '40')
parser$print_help()

args <- parser$parse_args(c('~/', 'episodic_auc', '--IV', 'EPISODICTOT', 
                            '--mask', 'mask.nii.gz',
                            '--DV', 'BRAIN', 
                            '--covariates', 'TIMECENTER', '--overwrite'))

if(! 'BRAIN' %in% c(args$DV, args$IV)){
  stop('Either --DV or --IV must be "BRAIN"... else why are you using neuropoint?')
}

model_name <- make.names(args$model_name)

model_dir <- create_np_model_dir(args$base_dir, model_name, overwrite = args$overwrite)
mask_file <- move_mask(args$mask, model_dir)
model_script <- write_model_script(model_dir = model_dir,
                                   y = args$DV, 
                                   x = args$IV, 
                                   model_name = model_name, 
                                   within_prefix = args$win_pre, 
                                   between_prefix = args$bw_pre,
                                   covariates = args$covariates, 
                                   id_var = args$covariates,
                                   overwrite = args$overwrite)

make_readargs(mask_fname = basename(mask_file), 
              set1 = args$set1, 
              setlabels = args$setlabels, 
              model_file = basename(model_file), 
              testvoxel = args$testvoxel, 
              output = args$output,
              debugfile = 'debug.Rdata', hpc = 'slurm', jobs = args$jobs)
