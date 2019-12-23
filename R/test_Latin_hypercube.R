#initial clean up
rm(list = ls())
graphics.off()
cat("\f")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('../data/feeagh')

# Load libraries
library(GOTMr);library(SimstratR);library(GLM3r);library(FLakeR);library(gotmtools);library(glmtools)
library(lubridate);library(plyr);library(ncdf4); library(ggplot2);library(reshape2)

# Load functions
source('../../R/export_config.R')
source('../../R/export_meteo.R')
source('../../R/export_init_cond.R')
source('../../R/run_ensemble.R')
source('../../R/run_Latin_hypercube.R')


# Load helper functions
source('../../R/helper_functions/input_json.R') # Potential function for 'simstrattools'
source('../../R/helper_functions/get_json_value.R') # Potential function for 'simstrattools'
source('../../R/helper_functions/input_nml.R') # This versions preserves comments in the nml
source('../../R/helper_functions/get_wtemp_df.R') # Potential function for flaketools
source('../../R/helper_functions/analyse_strat.R') # Potential function for flaketools

masterConfigFile <- 'Feeagh_master_config.yaml'


pars <- c('wind_factor', 'swr_factor', 'lw_factor')
mat <- matrix(data = c(0.5,2,0.5,1.5,0.5,1.5), nrow = 3, byrow = T)
df <- as.data.frame(mat)
rownames(df) <- pars

# Run Latin_hypercube sample
run_Latin_hypercube(parRange = df, num = 300, obs_file = 'LakeEnsemblR_wtemp_profile_standard.csv', param_file = 'latin_hypercube_params_FLake_GLM_GOTM_Simstrat_201912231225.csv', config_file = 'Feeagh_master_config.yaml', model = c('GLM', 'GOTM', 'Simstrat'), meteo_file = 'LakeEnsemblR_meteo_standard.csv')

# 2. Create meteo driver files
export_meteo(masterConfigFile, model = c('FLake', 'GLM', 'GOTM', 'Simstrat'),
             meteo_file = 'LakeEnsemblR_meteo_standard.csv',
             lhc_file = 'output/LHC_calibration_results_p300_201912231355.csv',
             metric = 'RMSE')

# 4. Run ensemble lake models
wtemp_list <- run_ensemble(config_file = masterConfigFile, model = c('FLake', 'GLM', 'GOTM', 'Simstrat'), return_list = TRUE,
                           create_netcdf = TRUE, obs_file = 'LakeEnsemblR_wtemp_profile_standard.csv')


####
# Plot model output using gotmtools/ggplot2
####
#Extract names of all the variables in netCDF
ens_out <- 'output/ensemble_output.nc4'
vars <- gotmtools::list_vars(ens_out)
vars # Print variables

plist <- list() # Initialize empty list for storing plots of each variable
for(i in 1:(length(vars)-1)){
  p1 <- gotmtools::plot_vari(ncdf = ens_out,
                             var = vars[i],
                             incl_time = FALSE,
                             limits = c(0,22),
                             zlab = 'degC')
  p1 <- p1 + scale_y_reverse() + #Reverse y-axis
    ggtitle(vars[i]) + # Add title using variable name
    xlab('')+ # Remove x-label
    theme_bw(base_size = 18) # Increase font size of plots
  plist[[i]] <- p1
}

# Plot all model simulations
# install.packages('ggpubr')
g1 <- ggpubr::ggarrange(plotlist = plist, ncol = 1, common.legend = TRUE, legend = 'right')
g1
ggsave('output/model_ensemble_watertemp_LHC.png', g1,  dpi = 300,width = 384,height = 300, units = 'mm')
####


###
# Model diagnostics plot
obs <- na.exclude(wtemp_list[[length(wtemp_list)]])
obs <- reshape2::melt(obs, id.vars = 1)
obs[,2] <- as.character(obs[,2])
obs[,2] <- as.numeric(gsub('wtr_','',obs[,2]))
colnames(obs) <- c('datetime','Depth_meter','Water_Temperature_celsius')
obs <- obs[order(obs[,1], obs[,2]),]

# Loop through each model and calculate diagnostics and generate diagnostic plot
for(i in 1:(length(wtemp_list)-1)){
  
  # Convert from wide format to long format - could be a function i.e. incorporated into gotmtools:wide2long()
  mod <- reshape2::melt(wtemp_list[[i]], id.vars = 1)
  mod[,2] <- as.character(mod[,2])
  mod[,2] <- as.numeric(gsub('wtr_','',mod[,2]))
  colnames(mod) <- c('datetime','Depth_meter','Water_Temperature_celsius')
  mod <- mod[order(mod[,1], mod[,2]),] # Reorder so datetime and depth increasing 
  # Check if same dimensions and merge by date and depth
  if(nrow(mod) != nrow(obs)){
    mod <- merge(obs, mod, by = c(1,2), all.x = T)
    mod <- mod[order(mod[,1], mod[,2]),]
    mod <- mod[,c(1,2,4)]
    colnames(mod) <- c('datetime','Depth_meter','Water_Temperature_celsius')
  }
  g2 <- diag_plots(mod, obs, colourblind = F, na.rm = T) # Needs to be sorted
  ggsave(paste0('output/diag_plot_', names(wtemp_list[i]), '_LHC.png'), g2,  dpi = 220,width = 384,height = 216, units = 'mm') 
  
  # Print name and output statistics
  print(names(wtemp_list)[i])
  print(sum_stat(mod, obs, depth = T))
}
###



pars <- read.csv('latin_hypercube_params_FLake_GLM_GOTM_Simstrat_201912231225.csv')



## FLake
res <- read.csv('FLake/output/latin_hypercube_calibration_results_p300_201912231237.csv')

dat <- merge(res, pars, by = 'par_id')
dat$model <- 'FLake'
all_par <- dat
fla_par <- dat[which.min(dat$RMSE), c(1,2,9:14)]

my.cols = RColorBrewer::brewer.pal(11, "Spectral")
sub <- which(mlt$variable %in% c('NSE', 'RMSE', 'Pearson_r'))
p1 <- ggplot(dat, aes(wind_factor, swr_factor, colour = RMSE))+
  geom_point(size =2)+
  geom_point(data = dat[which.min(dat$RMSE),], size =4, shape = 21)+
  scale_color_gradientn(colours = (my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1

p2 <- ggplot(dat, aes(wind_factor, swr_factor, colour = NSE))+
  geom_point(size =2)+
  geom_point(data = dat[which.max(dat$NSE),], size =4, shape = 21)+
  scale_color_gradientn(colours = rev(my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(dat, aes(wind_factor, swr_factor, colour = Pearson_r))+
  geom_point(size =2)+
  geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  scale_color_gradientn(colours = rev(my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g1 <- ggpubr::ggarrange(p1,p2,p3,nrow=3, align = 'v')
g1
ggsave('output/FLake_LHC_plot.png', plot = g1, dpi = 200,width = 324,height = 312, units = 'mm')


## GLM
library(plotly)
res <- read.csv('GLM/output/latin_hypercube_calibration_results_p300_201912231250.csv')

dat <- merge(res, pars, by = 'par_id')
dat$model <- 'GLM'
all_par <- rbind.data.frame(all_par, dat)
glm_par <- dat[which.min(dat$RMSE), c(1,2,9:14)]

p <- plot_ly(dat, x = ~wind_factor, y = ~swr_factor, z = ~lw_factor, color = ~RMSE, colors = my.cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'wind_factor'),
                      yaxis = list(title = 'swr_factor'),
                      zaxis = list(title = 'lw_factor')))

p

p <- plot_ly(dat, x = ~wind_factor, y = ~swr_factor, z = ~lw_factor, color = ~NSE, colors = rev(my.cols)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'wind_factor'),
                      yaxis = list(title = 'swr_factor'),
                      zaxis = list(title = 'lw_factor')))

p

## GOTM
res <- read.csv('GOTM/output/latin_hypercube_calibration_results_p300_201912231304.csv')
# pars <- read.csv('latin_hypercube_params_GOTM.csv')


dat <- merge(res, pars, by = 'par_id')
dat$model <- 'GOTM'
all_par <- rbind.data.frame(all_par, dat)
got_par <- dat[which.min(dat$RMSE), c(1,2,9:14)]

my.cols = RColorBrewer::brewer.pal(11, "Spectral")
sub <- which(mlt$variable %in% c('NSE', 'RMSE', 'Pearson_r'))
p1 <- ggplot(dat, aes(wind_factor, swr_factor, colour = RMSE))+
  geom_point(size =2)+
  geom_point(data = dat[which.min(dat$RMSE),], size =4, shape = 21)+
  scale_color_gradientn(colours = (my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1

p2 <- ggplot(dat, aes(wind_factor, swr_factor, colour = NSE))+
  geom_point(size =2)+
  geom_point(data = dat[which.max(dat$NSE),], size =4, shape = 21)+
  scale_color_gradientn(colours = rev(my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(dat, aes(wind_factor, swr_factor, colour = Pearson_r))+
  geom_point(size =2)+
  geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  scale_color_gradientn(colours = rev(my.cols))+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g1 <- ggpubr::ggarrange(p1,p2,p3,nrow=3, align = 'v')
g1
ggsave('output/GOTM_LHC_plot.png', plot = g1, dpi = 200,width = 324,height = 312, units = 'mm')


## Simstrat
library(plotly)
res <- read.csv('Simstrat/output/latin_hypercube_calibration_results_p300_201912231316.csv')

dat <- merge(res, pars, by = 'par_id')
dat$model <- 'Simstrat'
all_par <- rbind.data.frame(all_par, dat)
sim_par <- dat[which.min(dat$RMSE), c(1,2,9:14)]

my.cols = RColorBrewer::brewer.pal(11, "Spectral")

p <- plot_ly(dat, x = ~wind_factor, y = ~swr_factor, z = ~lw_factor, color = ~RMSE, colors = my.cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'wind_factor'),
                      yaxis = list(title = 'swr_factor'),
                      zaxis = list(title = 'lw_factor')))

p

fla_par
glm_par
got_par
sim_par


# Compare parameters
p1 <- ggplot(all_par, aes(wind_factor, RMSE, colour = model))+
  geom_point(size =2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1

p2 <- ggplot(all_par, aes(swr_factor, RMSE, colour = model))+
  geom_point(size =2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2

p3 <- ggplot()+
  geom_point(data = all_par[!(all_par$model %in% c('GLM', 'Simstrat')),], aes(lw_factor, RMSE, colour = model),size =0.001)+
  geom_point(data = all_par[(all_par$model %in% c('GLM', 'Simstrat')),], aes(lw_factor, RMSE, colour = model), size = 2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
g1 <- ggpubr::ggarrange(p1,p2,p3,nrow=3, align = 'v')
g1

ggsave('output/all_models_RMSE_LHC_plot.png', plot = g1, dpi = 200,width = 324,height = 312, units = 'mm')

# Compare parameters - Zoom
p1 <- ggplot(all_par, aes(wind_factor, RMSE, colour = model))+
  geom_point(size =2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(0,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1

p2 <- ggplot(all_par, aes(swr_factor, RMSE, colour = model))+
  geom_point(size =2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(0,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2

p3 <- ggplot()+
  geom_point(data = all_par[!(all_par$model %in% c('GLM', 'Simstrat')),], aes(lw_factor, RMSE, colour = model),size =0.001)+
  geom_point(data = all_par[(all_par$model %in% c('GLM', 'Simstrat')),], aes(lw_factor, RMSE, colour = model), size = 2)+
  # geom_point(data = dat[which.max(dat$Pearson_r),], size =4, shape = 21)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(0,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
g1 <- ggpubr::ggarrange(p1,p2,p3,nrow=3, align = 'v')
g1

ggsave('output/all_models_RMSE_zoom_LHC_plot.png', plot = g1, dpi = 200,width = 324,height = 312, units = 'mm')


mlt <- melt(all_par, id.vars = 'model', measure.vars = c('RMSE'))
ggplot(mlt, aes(x = value, colour = model))+
  geom_density()+
  theme_bw(base_size = 24)+
  xlab('RMSE')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1a <- ggplot(all_par, aes(x = wind_factor, colour = model))+
  geom_density()+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(0,3), xlim = c(0.5,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1b <- ggplot(all_par[all_par$RMSE <2,], aes(x = wind_factor, colour = model))+
  geom_density()+
  coord_cartesian(ylim = c(0,3), xlim = c(0.5,2))+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggpubr::ggarrange(p1a, p1b, nrow = 1, labels = c('Uninformed', 'Informed'))

ggplot(all_par[all_par$RMSE <2,], aes(x = swr_factor, colour = model))+
  geom_density()+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(all_par[all_par$RMSE <2,], aes(x = lw_factor, colour = model))+
  geom_density()+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
