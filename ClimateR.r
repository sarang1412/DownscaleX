# Load necessary libraries
options(java.parameters = "-Xmx20g")
library(loadeR)     
library(transformeR)
library(downscaleR) 
library(visualizeR) 
library(climate4R.value) 
library(magrittr)
library(gridExtra)
library(RColorBrewer)
library(sp)            
library(downscaleR.keras)

# Load climate data from NetCDF files
ncum_r <- loadGridData(dataset = "/home/trg1/DATA/ncumr_day1rf_jjas2021-24.nc",
                  var = "APCP_surface", lonLim = c(62, 106),
                  latLim = c(6, 41), 
                  years = 2021:2024)

obs <- loadGridData(dataset = "/home/trg1/candy/IMD_MSG-2020-24-jjas.nc",
                  var = "rf", lonLim = c(62, 106),
                  latLim = c(6, 41), 
                  years = 2020:2024)

# Subset the data for training and testing periods
xT <- subsetGrid(ncum_r, years = 2021:2022)
xt <- subsetGrid(ncum_r, years = 2023:2024)
yT <- subsetGrid(obs, years = 2021:2022)
yt <- subsetGrid(obs, years = 2023:2024)

# Convert observation data to binary grids
yT_bin <- binaryGrid(yT, threshold = 1, condition = "GT")
yt_bin <- binaryGrid(yt, threshold = 1, condition = "GT")

# Define color palette for plots
colsindex <- brewer.pal(n = 9, "YlGnBu")
cb <- colorRampPalette(colsindex)
cb2 <- colorRampPalette(colsindex)
# Create coordinate grids for plotting
coords_x <- expand.grid(xt$xyCoords$x,xt$xyCoords$y) ; names(coords_x) <- c("x","y")
coords_y <- expand.grid(yt$xyCoords$x,yt$xyCoords$y) ; names(coords_y) <- c("x","y")


# Create spatial plots for the data
pplot <- list()
pplot[[1]] <- spatialPlot(climatology(subsetGrid(xt)), backdrop.theme = "coastline",
                          main = "NCUM-R",
                          col.regions = cb2,
                          at = seq(0, 100, 10),
                          set.min = -3, set.max = 300, colorkey = TRUE, 
                          #sp.layout = list(list(SpatialPoints(coords_x), first = FALSE, col = rgb(0, 0, 0, alpha = 0.3), pch = 20, cex = 0.01))
                          )
pplot[[2]] <- spatialPlot(climatology(yt), backdrop.theme = "coastline", 
                          main = "Observation_IMD_MSG",
                          col.regions = cb,
                          at = seq(0, 100, 10),
                          set.min = -3, set.max = 300, colorkey = TRUE, 
                          #sp.layout = list(list(SpatialPoints(coords_y), first = FALSE, col = rgb(0, 0, 0, alpha = 0.3), pch = 20, cex = 0.01))
                          )

# Arrange the plots side by side
lay = rbind(c(1, 2))
grid.arrange(grobs = pplot, layout_matrix = lay)

# Define deep learning model architectures
deepName <- c("CNN-LM", "CNN1", "CNN10", "CNN-PR", "CNNdense")
architectures <- function(architecture, input_shape, output_shape) {
  if (architecture == "CNN-LM") {
    inputs <- layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x, filters = 50, kernel_size = c(3, 3), activation = 'linear', padding = "same")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3, 3), activation = 'linear', padding = "same")
    l3 = layer_conv_2d(l2, filters = 1, kernel_size = c(3, 3), activation = 'linear', padding = "same")
    l4 = layer_flatten(l3)
    parameter1 = layer_dense(l4, units = output_shape, activation = "sigmoid")
    parameter2 = layer_dense(l4, units = output_shape)
    parameter3 = layer_dense(l4, units = output_shape)
    outputs = layer_concatenate(list(parameter1, parameter2, parameter3))
    model <- keras_model(inputs = inputs, outputs = outputs)
  }
  
  if (architecture == "CNN1") {
    inputs <- layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x, filters = 50, kernel_size = c(3, 3), activation = 'relu', padding = "same")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3, 3), activation = 'relu', padding = "same")
    l3 = layer_conv_2d(l2, filters = 1, kernel_size = c(3, 3), activation = 'relu', padding = "same")
    l4 = layer_flatten(l3)
    parameter1 = layer_dense(l4, units = output_shape, activation = "sigmoid")
    parameter2 = layer_dense(l4, units = output_shape)
    parameter3 = layer_dense(l4, units = output_shape)
    outputs = layer_concatenate(list(parameter1, parameter2, parameter3))
    model <- keras_model(inputs = inputs, outputs = outputs)
  }
  
  if (architecture == "CNN10") {
    inputs <- layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x, filters = 50, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l3 = layer_conv_2d(l2, filters = 10, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l4 = layer_flatten(l3)
    parameter1 = layer_dense(l4, units = output_shape, activation = "sigmoid")
    parameter2 = layer_dense(l4, units = output_shape)
    parameter3 = layer_dense(l4, units = output_shape)
    outputs = layer_concatenate(list(parameter1, parameter2, parameter3))
    model <- keras_model(inputs = inputs, outputs = outputs)
  }
  
  if (architecture == "CNN-PR") {
    inputs <- layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x, filters = 10, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l3 = layer_conv_2d(l2, filters = 50, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l4 = layer_flatten(l3)
    parameter1 = layer_dense(l4, units = output_shape, activation = "sigmoid")
    parameter2 = layer_dense(l4, units = output_shape)
    parameter3 = layer_dense(l4, units = output_shape)
    outputs = layer_concatenate(list(parameter1, parameter2, parameter3))
    model <- keras_model(inputs = inputs, outputs = outputs)
  }
  
  if (architecture == "CNNdense") {
    inputs <- layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x, filters = 50, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l3 = layer_conv_2d(l2, filters = 10, kernel_size = c(3, 3), activation = 'relu', padding = "valid")
    l4 = layer_flatten(l3)
    l5 = layer_dense(l4, units = 50, activation = "relu")
    l6 = layer_dense(l5, units = 50, activation = "relu")
    parameter1 = layer_dense(l6, units = output_shape, activation = "sigmoid")
    parameter2 = layer_dense(l6, units = output_shape)
    parameter3 = layer_dense(l6, units = output_shape)
    outputs = layer_concatenate(list(parameter1, parameter2, parameter3))
    model <- keras_model(inputs = inputs, outputs = outputs)
  }
  
  return(model)
}

# Prepare data for training and testing
xy.T <- prepareData.keras(xT, binaryGrid(gridArithmetics(yT, 1, operator = "-"),
                                        condition = "GE",
                                        threshold = 0,
                                        partial = TRUE),
                          first.connection = "conv",
                          last.connection = "dense",
                          channels = "last")
xy.tT <- prepareNewData.keras(xT, xy.T)
xy.t <- prepareNewData.keras(xt, xy.T)

# Define simulation parameters
simulateName <- c("deterministic", "stochastic")
simulateDeep <- c(FALSE, TRUE)

# Train and evaluate models
lapply(1:length(deepName), FUN = function(z) {
  model <- architectures(architecture = deepName[z],
                         input_shape = dim(xy.T$x.global)[-1],
                         output_shape = dim(xy.T$y$Data)[2])
  downscaleTrain.keras(obj = xy.T,
             model = model,
             clear.session = TRUE,
             compile.args = list("loss" = bernouilliGamma.loss_function(last.connection = "dense"),
                                           "optimizer" = optimizer_adam(lr = 0.0001)),
                       fit.args = list("batch_size" = 100,
                            "epochs" = 1000,
                            "validation_split" = 0.1,
                            "verbose" = 1,
                            "callbacks" = list(callback_early_stopping(patience = 30),
                                callback_model_checkpoint(filepath=paste0('./models/precip/',deepName[z],'.h5'),
                                monitor='val_loss', save_best_only=TRUE))))
  lapply(1:length(simulateDeep), FUN = function(zz) {
    pred_ocu_train <- downscalePredict.keras(newdata = xy.tT,
                                             model = list("filepath" = 
                                                      paste0("./models/precip/",deepName[z],".h5"), 
                                                     "custom_objects" = 
                                                     c("custom_loss" = 
                                                      bernouilliGamma.loss_function(
                                                        last.connection = "dense"))),
                                             C4R.template = yT,
                                             clear.session = TRUE) %>% 
      subsetGrid(var = "pr1")
    pred <- downscalePredict.keras(newdata = xy.t,
                                   model = list("filepath" = 
                                              paste0("./models/precip/",deepName[z],".h5"), 
                                              "custom_objects" = 
                                              c("custom_loss" = 
                                              bernouilliGamma.loss_function(last.connection = "dense"))),
                                   C4R.template = yT,
                                   clear.session = TRUE) 
    pred <- bernouilliGamma.statistics(p = subsetGrid(pred, var = "pr1"),
                                       alpha = subsetGrid(pred, var = "pr2"),
                                       beta = subsetGrid(pred, var = "pr3"),
                                       simulate = simulateDeep[zz],
                                       bias = 1)
    pred_ocu <- subsetGrid(pred, var = "probOfRain") %>% redim(drop = TRUE)
    pred_amo <- subsetGrid(pred, var = "amountOfRain") %>% redim(drop = TRUE)
    pred_bin <- binaryGrid(pred_ocu, ref.obs = yT_bin, ref.pred = pred_ocu_train); rm(pred_ocu_train)
    save(pred_bin, pred_ocu, pred_amo, file = 
            paste0("./Data/precip/predictions_", simulateName[zz], "_", deepName[z], ".rda"))
  })
})

# Define validation parameters
simulateName <- c(rep("deterministic", 5), "stochastic", rep("deterministic", 3))
models <- c("glm1", "glm4",
               "CNN-LM", "CNN1", "CNN10",
               "CNN-PR", "CNNdense")
measures <- c("ts.rocss", "ts.RMSE", "ts.rs", rep("biasRel", 3), rep("bias", 3))
index <- c(rep(NA, 3), "Mean", rep("P98", 2), "AnnualCycleRelAmp",
           "WetAnnualMaxSpell", "DryAnnualMaxSpell")

# Evaluate models
validation.list <- lapply(1:length(measures), FUN = function(z) {
  lapply(1:length(models), FUN = function(zz) {
    args <- list()
    load(paste0("./Data/precip/predictions_", simulateName[z], "_", models[zz], ".rda"))
    if (simulateName[z] == "deterministic") {
      pred <- gridArithmetics(pred_bin, pred_amo, operator = "*")
      if (measures[z] == "ts.rocss") {
        args[["y"]] <- yt_bin; args[["x"]] <- pred_ocu
      } else if (measures[z] == "ts.RMSE") {
        args[["y"]] <- yt; args[["x"]] <- pred_amo
        args[["condition"]] = "GT"; args[["threshold"]] = 1; args[["which.wetdays"]] = "Observation"  
      } else {
        args[["y"]] <- yt; args[["x"]] <- pred
      }
    } else {
      pred <- gridArithmetics(pred_ocu, pred_amo, operator = "*")
      args[["y"]] <- yt; args[["x"]] <- pred
    }
    args[["measure.code"]] <- measures[z]
    if (!is.na(index[z])) args[["index.code"]] <- index[z]
    do.call("valueMeasure", args)$Measure
  }) %>% makeMultiGrid()
})
save(validation.list, file = "./Data/precip/validation.rda")

# Plot validation results
par(mfrow = c(3, 3)) 
ylabs <- c("ROCSS", "RMSE (wet days, mm)",
           "Spearman Corr.", "biasRel(%)",
           "biasRel P98 (DET, %)", "biasRel P98 (STO, %)",
           "Annual Cycle Rel. Amplitude", "biasRel WetAMS (days)",
           "biasRel DryAMS (days)")

lapply(1:length(validation.list), FUN = function(z) {
  # Set y-axis limits based on the measure
  if (z == 1) {ylim <- c(0.65, 0.9)}
  if (z == 2) {ylim <- c(3, 6.5)}
  if (z == 3) {ylim <- c(0.5, 0.8)}
  if (z == 4) {ylim <- c(-0.2, 0.2)}
  if (z == 5) {ylim <- c(-0.4, 0.0)}
  if (z == 6) {ylim <- c(-0.2, 0.2)}
  if (z == 7) {ylim <- c(-1, 1)}
  if (any(z == c(8, 9))) {ylim <- c(-1, 1)}
  
  # Extract and reshape the validation data
  index <- (validation.list[[z]] %>% redim(drop = TRUE))$Data
  dim(index) <- c(nrow(index), prod(dim(index)[2:3]))
  indLand <- (!apply(index, MARGIN = 2, anyNA)) %>% which()
  index <- index[, indLand] %>% t()
  
  # Calculate median and quantiles for the boxplot
  mglm4 <- median(index[, 2], na.rm = TRUE)
  perc <- apply(index, MARGIN = 2, FUN = function(z) quantile(z, probs = c(0.1, 0.9)))
  
  # Create the boxplot
  boxplot(index, outline = FALSE, ylim = ylim, range = 0.0001, ylab = ylabs[z], asp = 1)
  lines(c(0, 8), c(mglm4, mglm4), col = "red")
  for (i in 1:ncol(index)) lines(c(i, i), perc[, i], lty = 2)
})
