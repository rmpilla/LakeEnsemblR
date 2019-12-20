
source('../../R/run_Latin_hypercube.R')
pars <- c('wind_factor', 'swr_factor', 'lw_factor')
mat <- matrix(data = c(0.5,2,0.5,1.5,0.5,1.5), nrow = 3, byrow = T)
df <- as.data.frame(mat)
rownames(df) <- pars

# Run Latin_hypercube sample
run_Latin_hypercube(parRange = parRange, num = 100, obs_file = 'LakeEnsemblR_wtemp_profile_standard.csv', param_file = 'latin_hypercube_params_FLake_GLM.csv', config_file = 'Feeagh_master_config.yaml', model = c('GOTM'), meteo_file = 'LakeEnsemblR_meteo_standard.csv')


pars <- read.csv('latin_hypercube_params_FLake_GLM.csv')



## FLake
res <- read.csv('FLake/output/latin_hypercube_calibration_results.csv')

dat <- merge(res, pars, by = 'par_id')
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
res <- read.csv('GLM/output/latin_hypercube_calibration_results.csv')

dat <- merge(res, pars, by = 'par_id')
glm_par <- dat[which.min(dat$RMSE), c(1,2,9:14)]

my.cols = RColorBrewer::brewer.pal(11, "Spectral")

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
res <- read.csv('GOTM/output/latin_hypercube_calibration_results.csv')
# pars <- read.csv('latin_hypercube_params_GOTM.csv')


dat <- merge(res, pars, by = 'par_id')
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

fla_par
glm_par
got_par
