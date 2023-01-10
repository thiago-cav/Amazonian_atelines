## biomod2: Multi species modelling ----
## Most of the R code below for multi-species distribution modelling is
## available in the vignette of the biomod package
## (https://rdrr.io/rforge/biomod2/src/inst/doc/Multi_species_computation.R),
## except for the parameter implementation settings

rm(list = ls())

## setup environment ----
setwd('workdir')

## load the required packages
library(biomod2)
library(raster)
library(rasterVis)
library(gridExtra)
library(reshape2)

## read data ----
## species occurences data
data <- read.csv('../data/atelins_occ.csv', stringsAsFactors = FALSE)
head(data)
table(data$species)
spp_to_model <- unique(data$species)

## curent climatic variables
stk_current <-
  stack(paste0("../data/stk_current.gri"))

## 2050 climatic variables
stk_2050_BC_45 <-
  stack(paste0("../data/stk_2070_v2.gri"))


## 2070 climatic variables
stk_2070_BC_45 <-
  stack(paste0("../data/stk_2100_v2.gri"))


stk_current
stk_2050_BC_45
stk_2070_BC_45



## build species modelling wrapper ----
biomod2_wrapper <- function(sp){
  cat("\n> species : ", sp)

  ## get occurrences points
  sp_dat <- data[data$species == sp, ]

  ## formating the data
  sp_format <-
    BIOMOD_FormatingData(
      resp.var = rep(1, nrow(sp_dat)),
      expl.var = stk_current,
      resp.xy = sp_dat[, c("long", "lat")],
      resp.name = sp,
      PA.strategy = "random",
      PA.nb.rep = 2,
      PA.nb.absences = 1000
    )

  ## print formatting summary
  sp_format

  ## save image of input data summary
  if(!exists(sp)) dir.create(sp)
  pdf(paste(sp, "/", sp ,"_data_formated.pdf", sep="" ))
  try(plot(sp_format))
  dev.off()

  ## define models options
  ## define models options
  sp_opt <- BIOMOD_ModelingOptions(RF = list(do.classif = TRUE,
                                             ntree = 500,
                                             mtry = 3,
                                             nodesize = 5,
                                             maxnodes = NULL),
                                   GAM = list(method = "REML"),
                                   GBM = list(distribution = 'bernoulli',
                                              n.trees = 1000,
                                              interaction.depth = 7,
                                              n.minobsinnode = 5,
                                              shrinkage = 0.001,
                                              bag.fraction = 0.75,
                                              train.fraction = 1,
                                              cv.folds = 5,# 5-fold cross-validation
                                              keep.data = FALSE,
                                              verbose = FALSE,
                                              perf.method = 'cv'))


  ## model species
  sp_model <- BIOMOD_Modeling(
    sp_format,
    models = c('GBM', 'GAM', 'RF'),
    models.options = sp_opt,
    NbRunEval = 2,
    DataSplit = 75,
    Yweights = NULL,
    VarImport = 3,
    models.eval.meth = c('TSS', 'ROC'),
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = FALSE,
    modeling.id = "demo2"
  )

  ## save some graphical outputs
  #### models scores
  pdf(paste0(sp, "/", sp , "_models_scores.pdf"))
  try(gg1 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'models', plot = FALSE))
  try(gg2 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'data_set', plot = FALSE))
  try(gg3 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'cv_run', plot = FALSE))
  try(grid.arrange(gg1, gg2, gg3))
  dev.off()

  ## build ensemble models
  sp_ens_model <-
    BIOMOD_EnsembleModeling(
      modeling.output = sp_model,
      em.by = 'all',
      eval.metric = 'TSS',
      eval.metric.quality.threshold = 0.6,
      models.eval.meth = c('TSS','ROC'),
      prob.mean = FALSE,
      prob.mean.weight = TRUE,
      VarImport = 0
    )

  ## do projections
  proj_scen <- c("current", "2050_BC_45", "2070_BC_45")

  for(scen in proj_scen){
    cat("\n> projections of ", scen)

    ## single model projections
    sp_proj <-
      BIOMOD_Projection(
        modeling.output = sp_model,
        new.env = get(paste0("stk_", scen)),
        proj.name = scen,
        selected.models = 'all',
        binary.meth = "TSS",
        filtered.meth = NULL,
        compress = TRUE,
        build.clamping.mask = FALSE,
        do.stack = FALSE,
        output.format = ".img"
      )

    ## ensemble model projections
    sp_ens_proj <-
      BIOMOD_EnsembleForecasting(
        EM.output = sp_ens_model,
        projection.output = sp_proj,
        binary.meth = "TSS",
        compress = TRUE,
        do.stack = FALSE,
        output.format = ".img"
      )
  }

  return(paste0(sp," modelling completed !"))
}















#####

require(snowfall)
sfInit(parallel = TRUE, cpus = 8) ## here we only require 4 cpus
sfExportAll()
sfLibrary(biomod2)

ptm = proc.time()
sf_out <- sfLapply(spp_to_model, biomod2_wrapper)
processing.time = proc.time() - ptm
processing_time_in_minutes <- processing.time[3]/60
processing_time_in_minutes















## launch the spiecies modelling wrapper over species list ----
if(require(snowfall)){ ## parallel computation
  ## start the cluster
  sfInit(parallel = TRUE, cpus = 5) ## here we only require 4 cpus
  sfExportAll()
  sfLibrary(biomod2)
  ## launch our wrapper in parallel
  sf_out <- sfLapply(spp_to_model, biomod2_wrapper)
  ## stop the cluster
  sfStop()
} else { ## sequencial computation
  for (sp in spp_to_model){
    biomod2_wrapper(sp)
  }
  ## or with a lapply function in sequential model
  ## all_species_bm <- lapply(spp_to_model, biomod2_wrapper)
}












## Post-processing results exploration ------

## current conditons
### load binary projections
f_em_wmean_bin_current <-
  paste0(
    spp_to_model,
    "/proj_current/individual_projections/",
    spp_to_model, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
  )

### sum all projections
if(length(f_em_wmean_bin_current) >= 2){
  ## initialisation
  taxo_alpha_div_current <- raster(f_em_wmean_bin_current[1])
  for(f in f_em_wmean_bin_current[-1]){
    taxo_alpha_div_current <- taxo_alpha_div_current + raster(f)
  }
}

## 2050 conditons
### load binaries projections
f_em_wmean_bin_2050 <-
  paste0(
    spp_to_model,
    "/proj_2050_BC_45/individual_projections/",
    spp_to_model, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
  )

### sum all projections
if(length(f_em_wmean_bin_2050) >= 2){
  ## initialisation
  taxo_alpha_div_2050 <- raster(f_em_wmean_bin_2050[1])
  for(f in f_em_wmean_bin_2050[-1]){
    taxo_alpha_div_2050 <- taxo_alpha_div_2050 + raster(f)
  }
}

## 2070 conditons
### load binaries projections
f_em_wmean_bin_2070 <-
  paste0(
    spp_to_model,
    "/proj_2070_BC_45//individual_projections/",
    spp_to_model, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
  )

### sum all projections
if(length(f_em_wmean_bin_2070) >= 2){
  ## initialisation
  taxo_alpha_div_2070 <- raster(f_em_wmean_bin_2070[1])
  for(f in f_em_wmean_bin_2070){
    taxo_alpha_div_2070 <- taxo_alpha_div_2070 + raster(f)
  }
}

## plot the alpha-div maps
levelplot(
  stack(
    c(
      current = taxo_alpha_div_current,
      in_2050 = taxo_alpha_div_2050,
      in_2070 = taxo_alpha_div_2070
    )
  ),
  main = expression(paste("Larus ", alpha, "-diversity")),
  par.settings = BuRdTheme
)



dev.off()


par(mfrow=c(1,3),oma = c(0, 0, 5, 0))
plot(taxo_alpha_div_current, xlim=c(-85,-40), ylim=c(-20, 10), main = "Current")
plot(taxo_alpha_div_2050, xlim=c(-85,-40), ylim=c(-20, 10), main = "in 2070")
plot(taxo_alpha_div_2070, xlim=c(-85,-40), ylim=c(-20, 10), main = "in 2100")
mtext(expression(paste("Atelins ", alpha, "-diversity")), outer = TRUE,cex = 2)


##### reload the BIOMOD_Modelling() output #####
## keep a track
(bm_out_file <- load("Lagothrix.cana/Lagothrix.cana.demo2.models.out"))

ls()

## lets store this object with another name and free workspace
my_bm_model <- get(bm_out_file)
rm(list = c(bm_out_file, 'bm_out_file'))
ls()


## check our model is well loaded
my_bm_model

#plotting model metrics scores
models_scores_graph(my_bm_model, metrics = c("TSS", "ROC"), by = 'models', plot = TRUE)

#getting evaluation metrics as a dataframe
evatualtion_results<-as.data.frame(get_evaluations(my_bm_model))


#checking if some models fail
my_bm_model@models.failed

#checking variables importance in each model
my_bm_model@variables.importances





##### reloading the ENSEMBLE output #####
## keep a track
(bm_out_file <- load("Lagothrix.cana/Lagothrix.cana.demo2ensemble.models.out"))

ls()

## lets store this object with another name and free workspace
my_ensemble_model <- get(bm_out_file)
rm(list = c(bm_out_file, 'bm_out_file'))
ls()


## check our model is well loaded
my_ensemble_model



x<-as.data.frame(get_evaluations(my_ensemble_model))




##### loading the individual projections for L. cana #####
raster_data <- list.files(path="Ateles.paniscus/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_paniscus<-current[[1]]

dev.off()
par(mfrow=c(1,2),oma = c(0, 0, 5, 0))
plot(ensemble_current_paniscus, xlim=c(-85,-40), ylim=c(-20, 10))

ensemble_binary_paniscus<-current[[2]]
plot(ensemble_binary_paniscus, xlim=c(-85,-40), ylim=c(-20, 10))
mtext(paste("Gray woollies_Current habitat suitability and binary predictions"), outer = TRUE,cex = 2)


##### loading the individual projections for Lagothrix.poeppigii #####
raster_data <- list.files(path="C:/Users/thiag/OneDrive/Documentos/Projetos/biomod2_video_multi_species_modelling/workdir/Lagothrix.poeppigii/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_poeppigii<-current[[1]]

dev.off()
par(mfrow=c(1,2),oma = c(0, 0, 5, 0))
plot(ensemble_current_poeppigii, xlim=c(-85,-40), ylim=c(-20, 10))

ensemble_binary_poeppigii<-current[[2]]
plot(ensemble_binary_poeppigii, xlim=c(-85,-40), ylim=c(-20, 10))
mtext(paste("L_poeppigii_Current habitat suitability and binary predictions"), outer = TRUE,cex = 2)





##### loading the individual projections for A. belzebuth #####
raster_data <- list.files(path="Ateles.belzebuth/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_belzebuth<-current[[1]]

dev.off()
par(mfrow=c(1,2),oma = c(0, 0, 5, 0))
plot(ensemble_current_belzebuth, xlim=c(-85,-40), ylim=c(-20, 10))

ensemble_belzebuth<-current[[2]]
plot(ensemble_belzebuth, xlim=c(-85,-40), ylim=c(-20, 10))
mtext(paste("Gray woollies_Current habitat suitability and binary predictions"), outer = TRUE,cex = 2)


raster_data <- list.files(path="Lagothrix.flavicauda/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_flav<-current[[1]]
plot(ensemble_current_flav, xlim=c(-85,-40), ylim=c(-20, 10))
points(data$long, data$lat, pch=19, col=factor(data$species), cex=0.3)

head(data)


#saving the richness rasters
writeRaster(taxo_alpha_div_current, filename="Alpha diversity current", format="GTiff", overwrite=TRUE)
writeRaster(taxo_alpha_div_2050, filename="Alpha diversity in 2070", format="GTiff", overwrite=TRUE)
writeRaster(taxo_alpha_div_2070, filename="Alpha diversity in 2100", format="GTiff", overwrite=TRUE)



raster_data <- list.files(path="Ateles.chamek/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_chamek<-current[[1]]
plot(ensemble_current_chamek, xlim=c(-85,-40), ylim=c(-20, 10))
points(data$long, data$lat, pch=19, col=factor(data$species), cex=0.3)



raster_data <- list.files(path="Ateles.marginatus/proj_current/individual_projections",
                          pattern="*img$", full.name=TRUE)

current<-stack(raster_data)
rm(raster_data)

ensemble_current_marg<-current[[1]]
plot(ensemble_current_marg, xlim=c(-85,-40), ylim=c(-20, 10))
points(data$long, data$lat, pch=19, col=factor(data$species), cex=0.5)


plot(taxo_alpha_div_current, xlim=c(-85,-40), ylim=c(-20, 10))


ensemble_current_marg_2<-current[[2]]
plot(ensemble_current_marg_2, xlim=c(-85,-40), ylim=c(-20, 10))
plot(taxo_alpha_div_2070, xlim=c(-85,-40), ylim=c(-20, 10))
