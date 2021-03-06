#' Run Latin hypercube sampling
#'@description
#'Run lake models using Latin hypercube sampling for model parameters.
#'
#' @name run_LHC
#' @param config_file filepath; to LakeEnsemblr yaml master config file
#' @param num integer; the number of random parameter sets to generate. If param file is provided num = number of parameters in that file.
#' @param param_file filepath; to previously created parameter file set. If NULL creates a new parameter set. Defaults to NULL
#' #' @param method character; Method for calibration. Can be 'met', 'model' or 'both'. Needs to be specified by the user.
#' @param obs_file filepath; to LakeEnsemblR standardised observed water temperature profile data. If included adds observed data to netCDF and list if they are set to TRUE. Defaults to NULL.
#' @param config_file filepath; to LakeEnsemblr yaml master config file
#' @param model vector; model to export driving data. Options include c('GOTM', 'GLM', 'Simstrat', 'FLake')
#' @param folder filepath; to folder which contains the model folders generated by export_config()
#' @param meteo_file filepath; to met file which is in the standardised LakeEnsemblR format.
#' @param folder filepath; to folder which contains the model folders generated by export_config()
#' @param spin_up numeric; Number of days to disregard as spin-up for analysis.
#'
#' @examples
#' \dontrun{
#'pars <- c('wind_factor', 'swr_factor', 'lw_factor')
#'mat <- matrix(data = c(0.5,2,0.5,1.5,0.5,1.5), nrow = 3, byrow = T)
#'df <- as.data.frame(mat)
#'rownames(df) <- pars
#'run_LHC(parRange = parRange, num = 100, obs_file = 'LakeEnsemblR_wtemp_profile_standard.csv', config_file = 'Feeagh_master_config.yaml', model = 'FLake', meteo_file = 'LakeEnsemblR_meteo_standard.csv')
#' }
#'pars <- c('wind_factor', 'swr_factor', 'lw_factor')
#'mat <- matrix(data = c(0.5,2,0.5,1.5,0.5,1.5), nrow = 3, byrow = T)
#'df <- as.data.frame(mat)
#'rownames(df) <- pars
#'run_LHC(parRange = parRange, num = 100, obs_file = 'LakeEnsemblR_wtemp_profile_standard.csv', config_file = 'Feeagh_master_config.yaml', model = 'FLake', meteo_file = 'LakeEnsemblR_meteo_standard.csv')
#'@importFrom FME Latinhyper
#'@importFrom gotmtools get_yaml_value calc_cc input_nml sum_stat input_yaml get_vari
#'@importFrom glmtools get_nml_value
#'@importFrom reshape2 dcast
#'@importFrom FLakeR run_flake
#'@importFrom GLM3r run_glm
#'@importFrom GOTMr run_gotm
#'@importFrom SimstratR run_simstrat
#'@importFrom lubridate round_date seconds_to_period
#'
#' @export


run_LHC <- function(config_file, num = NULL, param_file = NULL, method, model = c('FLake', 'GLM', 'GOTM', 'Simstrat'), folder = '.', spin_up = NULL){

  # It's advisable to set timezone to GMT in order to avoid errors when reading time
  original_tz = Sys.getenv("TZ")
  Sys.setenv(TZ="GMT")
  tz = "UTC"
  # Set working directory
  oldwd <- getwd()

  # this way if the function exits for any reason, success or failure, these are reset:
  on.exit({
    setwd(oldwd)
    Sys.setenv(TZ=original_tz)
  })



  # get lat and lon - currently hack getting from GOTM but maybe could be in global config file?
  yaml = file.path(folder,config_file)

  # Function to be added to gotmtools
  lat <- get_yaml_value(file = yaml, label = 'location', key = 'latitude')
  lon <- get_yaml_value(file = yaml, label = 'location', key = 'longitude')
  depth <- get_yaml_value(file = yaml, label = 'location', key = 'depth')
  start <- get_yaml_value(file = yaml, label = 'time', key = 'start')
  stop <- get_yaml_value(file = yaml, label = 'location', key = 'stop')
  meteo_file <- get_yaml_value(file = yaml, label = 'meteo', key = 'file')
  obs_file <- get_yaml_value(file = yaml, label = 'observations', key = 'file')
  time_unit <- get_yaml_value(config_file, "output", "time_unit")
  time_step <- get_yaml_value(config_file, "output", "time_step")
  met_timestep <- get_yaml_value(config_file, "meteo", "time_step")

  # Create output time vector
  if(is.null(spin_up)){
    out_time <- seq.POSIXt(as.POSIXct(start, tz = tz), as.POSIXct(stop, tz = tz), by = paste(time_step, time_unit))
  }else{
    start <- as.POSIXct(start, tz = tz) + spin_up * 24 * 60 * 60
    stop <- as.POSIXct(stop, tz = tz)
    out_time <- seq.POSIXt(as.POSIXct(start, tz = tz), as.POSIXct(stop, tz = tz), by = paste(time_step, time_unit))
  }
  out_time <- data.frame(datetime = out_time)

  if(met_timestep == 86400){
    out_hour <- hour(start)
  }else{
    out_hour = 0
  }

  # read in Observed data
  message('Loading observed wtemp data...')
  obs <- read.csv(file.path(folder, obs_file), stringsAsFactors = FALSE)
  message('Finished!')
  obs$datetime <- as.POSIXct(obs$datetime, tz = tz)

  # Susbet to out_time
  obs <- obs[obs$datetime %in% out_time$datetime,]

  obs_deps <- unique(obs$Depth_meter)

  # change data format from long to wide
  obs_out <- dcast(obs, datetime ~ Depth_meter, value.var = 'Water_Temperature_celsius')
  str_depths <- colnames(obs_out)[2:ncol(obs_out)]
  colnames(obs_out) <- c('datetime',paste('wtr_',str_depths, sep=""))
  obs_out$datetime <- as.POSIXct(obs_out$datetime)



  # Which hemisphere?
  if(lat > 0){
    NH = TRUE
  }else{
    NH = FALSE
  }


  ### Import data
  # I'd prefer to use a function that can read both comma and tab delimited. data.table::fread does this, but then it's data.table
  message('Loading met data...')
  met = read.csv(file.path(folder, meteo_file), stringsAsFactors = F)
  message('Finished!')
  met[,1] <- as.POSIXct(met[,1])
  # Check time step
  tstep <- diff(as.numeric(met[,1]))

  if((mean(tstep) - 86400)/86400 < -0.05){
    daily = FALSE
    subdaily = TRUE
  } else {
    daily = TRUE
    subdaily = FALSE
  }

  if(is.null(param_file)){
    param_file <- sample_LHC(config_file = config_file, num = num, method = method, folder = folder)
  }
  params <- read.csv(param_file, stringsAsFactors = FALSE)
  num = nrow(params)

  all_pars <- NULL
  all_strat <- NULL

  # FLake
  #####
  if('FLake' %in% model){

    # Format met file
    fla_met <- format_met(met = met, model = 'FLake', daily = daily, config_file = config_file)

    # Select nml file for running FLake
    nml_file <- get_yaml_value(config_file, "config_files", "flake")
    nml_file <- file.path(folder, nml_file)
    # Select nml file again
    nml_file_run <- basename(get_yaml_value(config_file, "config_files", "flake"))

    mean_depth <- suppressWarnings(get_nml_value(arg_name = 'depth_w_lk', nml_file = nml_file))
    depths <- seq(0,mean_depth,by = get_yaml_value(config_file,"output", "depths"))

    # Input values to nml
    nml_file <- list.files(file.path(folder, 'FLake'))[grep('nml', list.files(file.path(folder, 'FLake')))]
    nml_file <- file.path(folder, 'FLake', nml_file)

    input_nml(nml_file, 'SIMULATION_PARAMS', 'time_step_number', nrow(fla_met))
    input_nml(nml_file, 'METEO', 'meteofile', paste0("'",'LHS_meteo_file.dat',"'"))

    tmp_met_file <- file.path(folder, 'FLake', 'LHS_meteo_file.dat')

    # Add in obs depths which are not in depths and less than mean depth
    add_deps <- obs_deps[!(obs_deps %in% depths)]
    add_deps <- add_deps[which(add_deps < mean_depth)]
    depths <- c(add_deps, depths)
    depths <- depths[order(depths)]

    for(i in 1:nrow(params)){

      if(method %in% c('met', 'both')){

        scale_met(met = fla_met, pars = params[i,], model = 'FLake', out_file = tmp_met_file)
      }

      run_flake(sim_folder = file.path(folder, 'FLake'), nml_file = nml_file_run)

      # Extract output

      out <- calc_stats(obs, model = 'FLake', depths = depths, NH = NH, flake_nml = nml_file, out_time = out_time, out_hour = out_hour)

      fit <- out$fit
      strat <- out$strat

      fit$par_id <- params$par_id[i]
      strat$par_id <- params$par_id[i]


      if(i == 1){
        fit_stats <- fit
        strat_stats <- strat
      }else{
        fit_stats <- rbind.data.frame(fit_stats, fit)
        strat_stats <- rbind.data.frame(strat_stats, strat)
      }
      print(paste0('[',i,'/', nrow(params),']'))
    }

    fit_file <- gsub('params', 'fitness', param_file)
    strat_file <- gsub('params', 'strat', param_file)


    write.csv(fit_stats, file.path(folder, 'FLake', 'output', fit_file), quote = FALSE, row.names = FALSE)
    write.csv(strat_stats, file.path(folder, 'FLake', 'output', strat_file), quote = FALSE, row.names = FALSE)


    fit_stats$model <- 'FLake'
    strat_stats$model <- 'FLake'

    if(is.null(all_pars)){
      all_pars <- fit_stats
    }else{
      all_pars <- rbind.data.frame(all_pars, fit_stats)
    }

    if(is.null(all_strat)){
      all_strat <- strat_stats
    }else{
      all_strat <- rbind.data.frame(all_strat, strat_stats)
    }

    message('FLake: Finished Latin Hypercube Sampling calibration')

  }

  # GLM
  #####
  if('GLM' %in% model){
    glm_met <- format_met(met = met, model = 'GLM', config_file = config_file, daily = daily)

    if("LongWave" %in% colnames(glm_met)){
      lw_type = 'LW_IN'
    }else{
      lw_type = 'LW_IN' ### Needs to be developed catch if no LW
    }

    # Input to nml file
    nml_path <- file.path(folder, get_yaml_value(config_file, "config_files", "glm"))
    nml <- glmtools::read_nml(nml_path)

    nml_list <- list('subdaily' = subdaily, 'lw_type' = lw_type, 'meteo_fl' = 'temp_meteo_file.csv')
    nml <- glmtools::set_nml(nml, arg_list = nml_list)

    glmtools::write_nml(nml, nml_path)

    # Input values to nml
    nml_file <- file.path(folder, 'GLM', 'glm3.nml')

    input_nml(nml_file, 'meteorology', 'meteo_fl', paste0("'",'LHS_meteo_file.csv',"'"))

    # Get depths for comparison
    depths <- obs_deps

    tmp_met_file <- file.path(folder, 'GLM', 'LHS_meteo_file.csv')


    for(i in 1:nrow(params)){


      scale_met(met = glm_met, pars = params[i,], model = 'GLM', out_file = tmp_met_file)

      run_glm(sim_folder = file.path(folder, 'GLM'))

      # Extract output
      # Add in obs depths which are not in depths and less than mean depth

      # Extract output
      out <- calc_stats(obs, model = 'GLM', depths = depths, NH = NH)

      fit <- out$fit
      strat <- out$strat

      fit$par_id <- params$par_id[i]
      strat$par_id <- params$par_id[i]


      if(i == 1){
        fit_stats <- fit
        strat_stats <- strat
      }else{
        fit_stats <- rbind.data.frame(fit_stats, fit)
        strat_stats <- rbind.data.frame(strat_stats, strat)
      }
      print(paste0('[',i,'/', nrow(params),']'))
    }

    fit_file <- gsub('params', 'fitness', param_file)
    strat_file <- gsub('params', 'strat', param_file)


    write.csv(fit_stats, file.path(folder, 'GLM', 'output', fit_file), quote = FALSE, row.names = FALSE)
    write.csv(strat_stats, file.path(folder, 'GLM', 'output', strat_file), quote = FALSE, row.names = FALSE)

    fit_stats$model <- 'GLM'
    strat_stats$model <- 'GLM'

    if(is.null(all_pars)){
      all_pars <- fit_stats
    }else{
      all_pars <- rbind.data.frame(all_pars, fit_stats)
    }

    if(is.null(all_strat)){
      all_strat <- strat_stats
    }else{
      all_strat <- rbind.data.frame(all_strat, strat_stats)
    }

    message('GLM: Finished Latin Hypercube Sampling calibration')

  }

  ## GOTM
  if('GOTM' %in% model){

    met_got <- format_met(met = met, model = 'GOTM', daily = daily, config_file = config_file)

    got_yaml <- file.path(folder,get_yaml_value(config_file, "config_files", "gotm"))

    met_outfile <- 'LHS_meteo_file.dat'

    # Get depths for comparison
    depths = -obs_deps
    obs_got <- obs
    obs_got[,2] <- -obs_got[,2]

    out_file <- file.path(folder, 'GOTM', met_outfile)

    yaml_file <- file.path(folder, get_yaml_value(config_file, "config_files", "gotm"))

    for(i in 1:nrow(params)){

      scale_met(met = met_got, pars = params[i,], model = 'GOTM', out_file = out_file)

      if(i == 1){
        # Helper function
        set_met_config_yaml(met = out_file, yaml_file = got_yaml)
      }

      run_gotm(sim_folder = file.path(folder, 'GOTM'), yaml_file = basename(yaml_file))


      # Extract output
      out <- calc_stats(obs, model = 'GOTM', depths = depths, NH = NH)

      fit <- out$fit
      strat <- out$strat

      fit$par_id <- params$par_id[i]
      strat$par_id <- params$par_id[i]


      if(i == 1){
        fit_stats <- fit
        strat_stats <- strat
      }else{
        fit_stats <- rbind.data.frame(fit_stats, fit)
        strat_stats <- rbind.data.frame(strat_stats, strat)
      }

      print(paste0('[',i,'/', nrow(params),']'))
    }

    fit_file <- gsub('params', 'fitness', param_file)
    strat_file <- gsub('params', 'strat', param_file)

    write.csv(fit_stats, file.path(folder, 'GOTM', 'output', fit_file), quote = FALSE, row.names = FALSE)
    write.csv(strat_stats, file.path(folder, 'GOTM', 'output', strat_file), quote = FALSE, row.names = FALSE)

    fit_stats$model <- 'GOTM'
    strat_stats$model <- 'GOTM'

    if(is.null(all_pars)){
      all_pars <- fit_stats
    }else{
      all_pars <- rbind.data.frame(all_pars, fit_stats)
    }

    if(is.null(all_strat)){
      all_strat <- strat_stats
    }else{
      all_strat <- rbind.data.frame(all_strat, strat_stats)
    }

    message('GOTM: Finished Latin Hypercube Sampling calibration')

  }

  ## Simstrat
  if('Simstrat' %in% model){

    # par file for running Simstrat
    par_file <- basename(get_yaml_value(config_file, "config_files", "simstrat"))
    par_fpath <- file.path(folder, 'Simstrat', par_file)

    met_simst <- format_met(met = met, model = 'Simstrat', config_file = config_file, daily = daily)

    met_outfile <- 'LHS_meteo_file.dat'


    input_json(file = par_fpath, label = 'Input', key = 'Forcing', paste0('"', met_outfile, '"'))

    # Need to input start and stop into json par file
    timestep <- get_json_value(par_fpath, "Simulation", "Timestep s")
    reference_year <- get_json_value(par_fpath, "Simulation", "Start year")

    # Set times
    reference_year <- year(as.POSIXct(start))
    input_json(par_fpath, "Simulation", "Start year", reference_year)
    start_date_simulation <- as.POSIXct(start, format = '%Y-%m-%d %H:%M:%S', tz = tz)
    end_date_simulation <- as.POSIXct(stop, format = '%Y-%m-%d %H:%M:%S', tz = tz)
    input_json(par_fpath, "Simulation", "Start d", as.numeric(difftime(start_date_simulation, as.POSIXct(paste0(reference_year,"-01-01")), units = "days")))
    input_json(par_fpath, "Simulation", "End d", as.numeric(difftime(end_date_simulation, as.POSIXct(paste0(reference_year,"-01-01")), units = "days")))
    input_json(par_fpath, "Output", "Times", time_step)

    #depths
    depths = -obs_deps


    out_file <- file.path(folder, 'Simstrat', met_outfile)

    for(i in 1:nrow(params)){

      scale_met(met = met_simst, pars = params[i,], model = 'Simstrat', out_file = out_file)

      run_simstrat(sim_folder = file.path(folder, 'Simstrat'), par_file = par_file, verbose = FALSE)

      # Extract output
      out <- calc_stats(obs, model = 'Simstrat', depths = depths, NH = NH, par_file = par_fpath, start = start)


      fit <- out$fit
      strat <- out$strat

      fit$par_id <- params$par_id[i]
      strat$par_id <- params$par_id[i]


      if(i == 1){
        fit_stats <- fit
        strat_stats <- strat
      }else{
        fit_stats <- rbind.data.frame(fit_stats, fit)
        strat_stats <- rbind.data.frame(strat_stats, strat)
      }


      print(paste0('[',i,'/', nrow(params),']'))
    }

    fit_file <- gsub('params', 'fitness', param_file)
    strat_file <- gsub('params', 'strat', param_file)

    write.csv(fit_stats, file.path(folder, 'Simstrat', 'output', fit_file), quote = FALSE, row.names = FALSE)
    write.csv(strat_stats, file.path(folder, 'Simstrat', 'output', strat_file), quote = FALSE, row.names = FALSE)

    fit_stats$model <- 'Simstrat'
    strat_stats$model <- 'Simstrat'

    if(is.null(all_pars)){
      all_pars <- fit_stats
    }else{
      all_pars <- rbind.data.frame(all_pars, fit_stats)
    }

    if(is.null(all_strat)){
      all_strat <- strat_stats
    }else{
      all_strat <- rbind.data.frame(all_strat, strat_stats)
    }

    message('Simstrat: Finished Latin Hypercube Sampling calibration')
  }

  dir.create(file.path(folder,'output'), showWarnings = FALSE)

  fit_file <- gsub('LHS_params', 'LHS_fit', param_file)
  strat_file <- gsub('LHS_params', 'LHS_strat', param_file)

  write.csv(all_pars, file.path(folder,'output', fit_file), quote = FALSE, row.names = FALSE)

  write.csv(all_strat, file.path(folder,'output', strat_file), quote = FALSE, row.names = FALSE)

}


