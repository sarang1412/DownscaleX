library(zoo)
library(hydroTSM)
library(hydroGOF)
library(sf)
library(raster)
library(RFmerge)
library(terra)


#data
precip_data <- rast("E:/Dissertation/data/RainGG/RF25_ind2023_rfp25_slice.nc")
IMD_data    <- rast("E:/Dissertation/obs_JJAS2023.nc")
DEM_data    <- raster("E:/Dissertation/data/SHP&DEM/end.tif")
#SHP        <- st_read("E:/Dissertation/data/SHP&DEM/.shp")

#coords
coords <- terra::xyFromCell(precip_data, 1:ncell(precip_data))
lon <- coords[, 1]
lat <- coords[, 2]
cell_ids <- 1:ncell(precip_data)

# Extract rainfall values for all layers and grid cells and make dataset
ID <- names(precip_data)
station <- data.frame(ID = cell_ids,lon = lon,lat = lat)
rainfall_values <- data.frame(cell_ids, lon, lat, terra::values(precip_data))

#utmz43n.p4s <- CRS("+init=epsg:32643")
IMD.utm<- projectRaster(from=IMD_data, crs=)
DEM.utm<- projectRaster(from=DEM_data, crs=)
stations.utm <- sf::st_transform(stations, crs=)
SHP.utm <- sf::st_transform(SHP, crs=)

covariates.utm <- list(IMD=IMD.utm, DEM=DEM.utm)

drty.out <- "C:/Users/candy/Desktop/vi"

rfmep <- RFmerge(x=precip_data, metadata=station.utm, cov=covariates.utm,
                 id="ID", lat="lat", lon="lon", mask=SHP.utm,
                 training=0.8, write2disk=TRUE, drty.out= drty.out, 
                 parallel="parallel")
