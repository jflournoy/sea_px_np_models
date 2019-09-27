
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

make_brain_centering_code <- function(x_w_name, x_b_name, id_var){
  brain_centering <- paste0('
  brain_df <- data.frame(BRAIN, ', id_var, ') #the id variable is defined in the setlabels file
  
  # Compute Within-person Mean of brain data
  win_mean <- aggregate(brain_df, by = list(', id_var, '), FUN = mean)[,1:2]
  colnames(win_mean) <- c("', id_var, '", "WIN_BRAIN_MEAN")
  brain_df <- merge(brain_df, win_mean, by="', id_var, '")
  
  grandmean <- mean(win_mean$WIN_BRAIN_MEAN)
  
  brain_df[, "', x_b_name,'"] <- brain_df$WIN_BRAIN_MEAN - grandmean
  brain_df[, "', x_w_name,'"] <- BRAIN - brain_df$WIN_BRAIN_MEAN
  
  # attach to these new variables
  ', x_b_name,' <- brain_df[, "', x_b_name,'"]
  ', x_w_name,' <- brain_df[, "', x_w_name,'"]
  ')
  return(brain_centering)
}

make_model_formula <- function(y, x_w_name, x_b_name, covariates){
  model_formula <- paste0(y, ' ~ 1 + ', paste(covariates, collapse = ' + '), ' + ', x_b_name, ' + ', x_w_name)
  return(model_formula)
}

disagg_x_names <- function(x, between_prefix, within_prefix){
  x_b_name <- paste0(between_prefix, x)
  x_w_name <- paste0(within_prefix, x)
  return(c(x_b_name = x_b_name, x_w_name = x_w_name))
}

make_model_data <- function(y, x, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  disagg_names <- disagg_x_names(x, between_prefix, within_prefix)
  x_b_name <- disagg_names['x_b_name']
  x_w_name <- disagg_names['x_w_name']
  
  mode_data_code <- '\n
  BRAIN <- voxeldat[,v]'
  
  if(x == 'BRAIN'){
    mode_data_code <- paste0(mode_data_code, 
                     make_brain_centering_code(x_w_name = x_w_name, 
                                               x_b_name = x_b_name, 
                                               id_var = id_var)) 
  }
  allvars <- c(y, covariates, x_b_name, x_w_name, id_var)
  data_frame_args <- paste(
    paste0(c(y, covariates, x_b_name, x_w_name, id_var), ' = ', c(y, covariates, x_b_name, x_w_name, id_var)),
    collapse = ', \n    ')
  mode_data_code <- paste0(mode_data_code,'
  model_data <- na.omit(data.frame(\n    ', data_frame_args,'))')
  return(mode_data_code)
}

make_model_header <- function(y, x, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum', is_permute = FALSE, permutationRDS = NULL){
  
  disagg_names <- disagg_x_names(x, between_prefix, within_prefix)
  x_b_name <- disagg_names['x_b_name']
  x_w_name <- disagg_names['x_w_name']
  
  header <- '
packages <- list(\'nlme\', \'clubSandwich\', \'reghelper\')
loaded <- lapply(packages, library, character.only = TRUE)

processVoxel <- function(v) {

  NOVAL <- 999
  retvals <- numeric()'
  if(is_permute){
    model_formula <- make_model_formula(y = y, x_w_name = x_w_name, x_b_name = x_b_name, covariates = covariates)
    header <- paste0(header, '\n
  permutationRDS = \'', permutationRDS,'\'
  targetDV = \'', x_w_name, '\'
  formula = ', model_formula)
  }
  header <- paste0(header, 
                   make_model_data(y = y, x = x, 
                                   within_prefix = within_prefix, 
                                   between_prefix = between_prefix, 
                                   covariates = covariates, 
                                   id_var = id_var))
  if(is_permute){
    header <- paste0(header, '\n
  if(v%%5e3 == 0){
    apercent <- sprintf("%3.0f%%", v/dim(voxeldat)[2]*100)
                     status_message <- paste0(apercent, ": Voxel ", v, " of ", dim(voxeldat)[2], "...")
                     message(status_message)
  }')
  }
  return(header)
}

make_lme_model_syntax <- function(y, x, model_name, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  disagg_names <- disagg_x_names(x, between_prefix, within_prefix)
  x_b_name <- disagg_names['x_b_name']
  x_w_name <- disagg_names['x_w_name']
  model_name_sw <- paste0(model_name, '_sw')
  
  model_formula <- make_model_formula(y = y, covariates = covariates, x_w_name = x_w_name, x_b_name = x_b_name)
  model_syntax <- paste0('
  model_formula <- ', model_formula, '

  # `method = REML` for unbiased estimate of variance parameters.
  # See: 
  # Luke, S. G. (2017). Evaluating significance in linear mixed-effects
  # models in R. Behavior Research Methods, 49(4), 1494–1502. 
  # https://doi.org/10.3758/s13428-016-0809-y
  e <- try(', model_name, ' <-
              nlme::lme(', model_formula, ',
                        random = ~1 | ', id_var, ', data = model_data, 
                        method = "REML", na.action=na.omit) )

  # Compute cluster corrected standard errors to account for, e.g., residual autocorrelation.
  # See coef_test help for description of "CR2"
  e_sandwich <- try(', model_name_sw, ' <- clubSandwich::coef_test(e, vcov = "CR2"))
')
  return(model_syntax)
}

make_permute_tail <- function(){
  permute_tail <- paste0("
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
                                                         lmeOpts = list(method = 'REML', na.action=na.omit))
  if(is.null(permuteModel$Z)){
    Z <- NOVAL
  } else {
    Z <- permuteModel$Z
  }
  names(Z) <- paste0(targetDV,'-z_sw')
  return(Z)
}")
  return(permute_tail)
}

make_model_footer <- function(y, x, model_name, within_prefix = 'WCEN_', between_prefix = 'GCEN_', covariates = NULL, id_var = 'idnum'){
  disagg_names <- disagg_x_names(x, between_prefix, within_prefix)
  x_b_name <- disagg_names['x_b_name']
  x_w_name <- disagg_names['x_w_name']
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
    names(reg_retvals) <- gsub("[\\\\(\\\\)]", "", gsub("(.*?)\\\\.(.*)", "\\\\2\\\\1", names(reg_retvals)))
    
    std_coefs <- coef(beta(', model_name, '))
    std_reg_retvals <- std_coefs[, "Value"]
    names(std_reg_retvals) <- paste0(gsub("\\\\(*(.*?)\\\\)*\\\\.*z*$", "\\\\1", names(std_reg_retvals)), "-beta")
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
    names(reg_retvals_sw) <- gsub("[\\\\(\\\\)]", "", gsub("(.*?)\\\\.(.*)", "\\\\2\\\\1", names(reg_retvals_sw)))
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
                          debugfile = 'debug.Rdata', hpc = 'slurm', jobs = '40',
                          is_permute = FALSE, permuteN){
  output <- file.path(output, paste0(output, '.'))
  setlabel_args <- paste0('"--setlabels', 1:length(setlabels), '", "', setlabels, '", ')
  readargs <- paste0('
cmdargs <- c("-m","', mask_fname, '", "--set1", "', set1, '",
             ', setlabel_args, '
             "--model", "', model_file, '",
             "--testvoxel", "', testvoxel, '",
             "--output", "', output, '",
             "--debugfile", "', debugfile, '",
             "--', hpc, 'N", "', jobs, '"')
  if(is_permute){
    readargs <- paste0(readargs, ',
             "--permute", "', permuteN,'"')
  }
  readargs <- paste0(readargs, ')
')
  return(readargs)
}

make_writeperms <- function(y, x, within_prefix, between_prefix, covariates, id_var, model_dir, debug_data, permuteN, permutationRDS, seed = NULL){
  model_data <- make_model_data(y = y, x = x, 
                                within_prefix = within_prefix, 
                                between_prefix = between_prefix, 
                                covariates = covariates, 
                                id_var = id_var)
  writeperms <- paste0("
library(permute)

message('Loading data...')
#This presumes you've arleady set up the target model
load('", debug_data, "')
attach(designmat)
nperm <- ", permuteN, " #Make sure this number matches input to NeuroPointillist or is bigger
v <- 1 #to get example brain data

", model_data, "

set.seed(", seed, ")
message('Generating ', nperm, ' permutations...')
ctrl.free <- how(within = Within(type = 'free'), nperm = nperm, blocks = ", id_var, ")
perm_set.free <- shuffleSet(n = ", id_var, ", control = ctrl.free)
permpath <- file.path('", model_dir,"', '", permutationRDS, "')
message('Saving permutations to ', permpath)
saveRDS(perm_set.free, permpath)")
  return(writeperms)
}

write_model_script <- function(model_dir, y, x, model_name, 
                               within_prefix = 'WCEN_', between_prefix = 'GCEN_', 
                               covariates = NULL, id_var = 'idnum',  overwrite = FALSE, 
                               is_permute = FALSE, permutationRDS = NULL){
  
  if(is_permute){
    header <- make_model_header(y = y, x = x, 
                                within_prefix = within_prefix, between_prefix = between_prefix, 
                                covariates = covariates, id_var = id_var, 
                                is_permute = is_permute, permutationRDS = permutationRDS)
    permute_tail <- make_permute_tail()
    model_text <- paste(header, permute_tail, sep = '\n')
    model_file <- file.path(model_dir, 'permute_free_model.R')
  } else {
    header <- make_model_header(y = y, x = x, 
                                within_prefix = within_prefix, between_prefix = between_prefix, 
                                covariates = covariates, id_var = id_var, is_permute = FALSE)
    model_syntax <- make_lme_model_syntax(y = y, x = x, model_name = model_name, 
                                          within_prefix = within_prefix, between_prefix = between_prefix, 
                                          covariates = covariates, id_var = id_var)
    footer <- make_model_footer(y = y, x = x, model_name = model_name, 
                                within_prefix = within_prefix, between_prefix = between_prefix, 
                                covariates = covariates, id_var = id_var)
    model_text <- paste(header, model_syntax, footer, sep = '\n')
    model_file <- file.path(model_dir, 'model.R')
  }
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

write_read_args <- function(model_dir, mask_fname, set1, setlabels, model_file, testvoxel = '1000', output,
                            debugfile = 'debug.Rdata', hpc = 'slurm', jobs = '40', overwrite = FALSE,
                            is_permute = FALSE, permuteN){
  readargs_text <- make_readargs(mask_fname = mask_fname, 
                                 set1 = set1, 
                                 setlabels = setlabels, 
                                 model_file = model_file, 
                                 testvoxel = testvoxel, 
                                 output = output,
                                 debugfile = debugfile, 
                                 hpc = hpc, 
                                 jobs = jobs,
                                 is_permute = is_permute, 
                                 permuteN = permuteN)
  readargs_file <- file.path(model_dir, 'readargs.R')
  message('readargs file: ', readargs_file)
  if(file.exists(readargs_file) & !overwrite){
    stop('readargs file exists. See help.')
  } else {
    f <- file(readargs_file, open = 'w')
    writeLines(readargs_text, f)
    close(f)
  }
  return(readargs_file)
}
 
write_writeperms <- function(y, x, within_prefix, between_prefix, covariates, id_var, model_dir, debug_data, permuteN, permutationRDS, overwrite = FALSE, seed = NULL){
  writeperms_text <- make_writeperms(y = y, x = x, 
                                     within_prefix = within_prefix, 
                                     between_prefix = between_prefix, 
                                     covariates = covariates, 
                                     id_var = id_var,
                                     model_dir = model_dir, 
                                     debug_data = debug_data, 
                                     permuteN = permuteN, 
                                     permutationRDS = permutationRDS,
                                     seed = seed)
  writeperms_file <- file.path(model_dir, 'write_permutations.R')
  message('permutation generation file: ', writeperms_file)
  if(file.exists(writeperms_file) & !overwrite){
    stop('permutation generation file exists. See help.')
  } else {
    f <- file(writeperms_file, open = 'w')
    writeLines(writeperms_text, f)
    close(f)
    message('Remember to move it to the npoint-generated directory after running npoint.')
  }
  return(writeperms_file)
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

run_writeperms <- function(writeperms_file){
  source(writeperms_file)
  return(NULL)
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
parser$add_argument('--jobs', type="character", help='Number of jobs', default = '40')
parser$add_argument('--permute', action='store_true', help = 'Are we building permutation models?')
parser$add_argument('--permuteN', type="character", help = 'Number of permutations', default = '1000')
parser$add_argument('--debugdata', type="character", help = 'Path to debug data to use in generating permutation matrix', default = '')
parser$add_argument('--seed', type="integer", help="Seed for permutation randomization.", default=NULL)
parser$add_argument('--permutationRDS', type="character", help='Name of RDS file containing permutation matrix', default = 'permutation_set-free.RDS')

args <- parser$parse_args()

if(! 'BRAIN' %in% c(args$DV, args$IV)){
  stop('Either --DV or --IV must be "BRAIN"... else why are you using neuropoint?')
}
if(args$permute & args$permutationRDS == ''){
  stop('Permutation file RDS not given -- required if `--permute` is set.')
}
if(args$permute & is.null(args$seed)){
  stop('Randome seed required if `--permute` is set. Use `--seed INT`')
}

model_name <- make.names(args$model_name)
message('Creating model ', model_name)

model_dir <- create_np_model_dir(args$base_dir, model_name, overwrite = args$overwrite)
mask_file <- move_mask(args$mask, model_dir)
model_script <- write_model_script(model_dir = model_dir,
                                   y = args$DV, 
                                   x = args$IV, 
                                   model_name = model_name, 
                                   within_prefix = args$win_pre, 
                                   between_prefix = args$bw_pre,
                                   covariates = args$covariates, 
                                   id_var = args$id_var,
                                   overwrite = args$overwrite,
                                   is_permute = args$permute,
                                   permutationRDS = args$permutationRDS)

readargs_file <- write_read_args(model_dir = model_dir,
                                 mask_fname = basename(mask_file), 
                                 set1 = args$set1, 
                                 setlabels = args$setlabels, 
                                 model_file = basename(model_script), 
                                 testvoxel = args$testvoxel, 
                                 output = args$output,
                                 debugfile = 'debug.Rdata', 
                                 hpc = 'slurm', 
                                 jobs = args$jobs,
                                 overwrite = args$overwrite,
                                 is_permute = args$permute,
                                 permuteN = args$permuteN)

if(args$permute){
  writeperms_file <- write_writeperms(y = args$DV, 
                                      x = args$IV, 
                                      within_prefix = args$win_pre, 
                                      between_prefix = args$bw_pre,
                                      covariates = args$covariates, 
                                      id_var = args$id_var,
                                      model_dir = model_dir, 
                                      debug_data = args$debugdata, 
                                      permuteN = args$permuteN, 
                                      permutationRDS = args$permutationRDS, 
                                      overwrite = args$overwrite,
                                      seed = args$seed)
  message('Generating permutations...')
  nada <- run_writeperms(writeperms_file)
}

### TESTING
# 
# args <- parser$parse_args(c(
#   "--IV", "BRAIN",
#   "--DV", "GAD7_TOT",
#   "--mask", "~/NewNeuropoint/dep.fear/make_perms/mask.nii.gz",
#   "--set1", "/mnt/stressdevlab/stress_pipeline/Group/FaceReactivity/NewNeuropoint/datafiles/setfilenames_happyGTcalm.txt",
#   "--setlabels", "/mnt/stressdevlab/stress_pipeline/Group/FaceReactivity/NewNeuropoint/datafiles/depanxcov-midpoint5.csv",
#   "--output", "test.test.gad.happy",
#   "--win_pre", "WCEN_",
#   "--bw_pre", "GCEN_",
#   "--id_var", "idnum",
#   "--covariates", "TIMECENTER",
#   "--jobs", "60",
#   "--overwrite",
#   "--permute",
#   "--permuteN", "1000",
#   "--debugdata", "~/NewNeuropoint/gad.fear/gad.fear/debug.Rdata",
#   "--seed", "8982",
#   "--permutationRDS", "permute_free.RDS",
#   getwd(),
#   "test.test.gad.happy"))
