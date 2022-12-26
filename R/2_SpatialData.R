#### Script to read and manipulate the spatial files needed to perform the JSDM ####
# edit: 24.jan.22
# edit: 10/05/22 - define new cell sizes (1km2 of cell size gives a lot of information to estimate), trying with 1km2
# edit: 08/09/22 - final objects

# Required packages
require(raster)
require(rgdal)
require(ggplot2)

# put all spatial files in the same resolution (1km2)----------------------------
wd <- getwd()

# read the downloaded bioclimatic variables in 30s (~1km2)
# create a list of .tif files that exist in the wd
bio.dir <- file.path(here::here("data/raw/WGS84/Bioclim_wc2.1"))
setwd(bio.dir)

# import all files at once
list.files(bio.dir)

bioclim_30s <- raster::stack(list.files(pattern = ".tif"))
bioclim_30s@layers 

# renaming
names(bioclim_30s) <- gsub("wc2.1_30s_", "", names(bioclim_30s))
names(bioclim_30s)
res(bioclim_30s)

# redefine the working directory (root)
setwd(wd)

# read the raster of elevation (SRTM elevation data)
elev_30s <- raster::raster(here::here("data/raw/WGS84/wc2.1_30s_elev.tif"))
elev_30s

# read the raster of land use and land cover by MapBiomas project
AFtrinac.r <- raster::raster(here::here("data/raw/WGS84/MapBiomas_c1/mapbiomas-atlantic-forest-collection-10-2019.tif"))
plot(AFtrinac.r)

# read the shapefile of AF trinational limits
AFtrinac.shp <- rgdal::readOGR(here::here("data/raw/WGS84/MapBiomas_c1/Limits_AFtrinacional.shp"))
AFtrinac.shp
# raster::plot(AFtrinac.shp)

# prepare to resampling

# Resampling the pixel to 1 km2 
# 1 sec = ~30.87m
# 30 sec = 926.1m
# 926.1m = 0.008333333 degrees
# 1000m = 0.008998308 degrees
# 2000m = 0.01799662 degrees

AF.1km <- raster::raster(AFtrinac.shp, res = 0.008998308)

# resampling bioclim variables
bio.AF.1km <- raster::resample(bioclim_30s, AF.1km, method = "bilinear")
res(bio.AF.1km)

saveRDS(bio.AF.1km, here::here("data/processed/spatial_objects/bio_AF_1km.rds"))

# resampling elevation data
elev.AF.1km <- raster::resample(elev_30s, AF.1km, method = "bilinear")
res(elev.AF.1km)

saveRDS(elev.AF.1km, here::here("data/processed/spatial_objects/elev_AF_1km.rds"))

# Croping raster files to study area
# crop bioclimatic by the AF shapefile
bio.AF.1km.mask <- mask(bio.AF.1km, AFtrinac.shp)
bio.AF.1km.crop <- crop(bio.AF.1km.mask, extent(AFtrinac.shp))
plot(bio.AF.1km.crop$bio_2)

df.bio.AF.1km <- as.data.frame(bio.AF.1km.crop, xy = T)

# writeRaster(bio.AF.1km.crop, here::here("data/processed/spatial_objects/bio_AF_1km_crop.tif"),
#                                         overwrite = TRUE)

saveRDS(bio.AF.1km.crop, here::here("data/processed/spatial_objects/bio_AF_1km_crop.rds"))
saveRDS(df.bio.AF.1km, here::here("data/processed/spatial_objects/df_bio_AF_1km.rds"))


# crop the elevation by the AF shapefile
elev.AF.1km.mask <- mask(elev.AF.1km, AFtrinac.shp)
elev.AF.1km.crop <- crop(elev.AF.1km.mask, extent(AFtrinac.shp))
plot(elev.AF.1km.crop)

df.elev.AF.1km <- as.data.frame(elev.AF.1km.crop, xy = T)

saveRDS(elev.AF.1km.crop, here::here("data/processed/spatial_objects/elev_AF_1km_crop.rds"))
saveRDS(df.elev.AF.1km, here::here("data/processed/spatial_objects/df_elev_AF_1km.rds"))

# resampling the mapbiomas LULC
lulc.AF.1km <- raster::resample(AFtrinac.r, AF.1km, method = "ngb")
res(lulc.AF.1km)

# Reclassify the categories of MapBiomas 

mtx.reclass <- matrix(c(0, 0, NA, # NA
                        1, 3, 1, # natural forest formation
                        3, 4, 2, # Natural open habitats
                        4, 9, 3, # Agricultural use
                        9, 11, 5, # Non Forest Natural Formation
                        11, 12, 2,
                        12, 13, 5, 
                        13, 21, 3,
                        21, 22, 4, # Urban area
                        22, 29, 2,
                        29, 33, 6), # Water
                      ncol = 3, byrow = TRUE)

lulc.AF.1km.reclass <- reclassify(lulc.AF.1km, mtx.reclass)

is.factor(lulc.AF.1km.reclass)

lulc.AF.1km.reclass <- as.factor(lulc.AF.1km.reclass)

x <- levels(lulc.AF.1km.reclass)[[1]]
x$code <- c("NA", "Forest_Formation", "Open_Formations",
            "Agriculture", "Urban_area", "NonForest_Formation", "Water")

levels(lulc.AF.1km.reclass) <- x
levels(lulc.AF.1km.reclass)

class(lulc.AF.1km.reclass)

#writeRaster(lulc.AF.1km.reclass, here::here("data/processed/lulc_AF_1km_reclass.tif"), 
#            format = "GTiff", overwrite = T) # we used this object to perform reprojection

saveRDS(lulc.AF.1km.reclass, here::here("data/processed/spatial_objects/lulc_AF_1km_reclass.rds"))

df.lulc.AF.1km.reclass <- as.data.frame(lulc.AF.1km.reclass, xy = T)
head(df.lulc.AF.1km.reclass)
colnames(df.lulc.AF.1km.reclass)[3] <- "code"

saveRDS(df.lulc.AF.1km.reclass, here::here("data/processed/spatial_objects/df_lulc_AF_1km_reclass.rds"))

# We will include the information about ecoregions, for all grid cells
# we use the TNC's terrestrial ecoregions of the world 
# (https://geospatial.tnc.org/datasets/b1636d640ede4d6ca8f5e369f2dc368b/about)

TNC <- readOGR(here::here("data/raw/WGS84/Ecoregions_AF/tnc_terr_ecoregions.shp"))
AFtrinac.shp <- readOGR(here::here("data/raw/WGS84/MapBiomas_c1/Limits_AFtrinacional.shp"))
elev.AF.1km.crop <- readRDS(here::here("data/processed/spatial_objects/elev_AF_1km_crop.rds"))

## TNC data is a shapefile and we need a raster format
tnc.rast <- raster::rasterize(TNC, elev.AF.1km.crop)
plot(tnc.rast)

# masking and cropping TNC from AF extension 
tnc.mask <- raster::mask(tnc.rast, AFtrinac.shp)
tnc.crop <- raster::crop(tnc.mask, extent(AFtrinac.shp))
plot(tnc.crop)

df.tnc <- as.data.frame(tnc.crop, xy = T)
head(df.tnc)

saveRDS(tnc.crop, here::here("data/processed/spatial_objects/TNC_1km.rds"))
saveRDS(df.tnc, here::here("data/processed/spatial_objects/df_TNC_1km.rds"))


################################################################################


# Extract environmental variables for butterflies communities -------------
library(raster)

# read spatial files in the same resolution
bio.AF.1km.crop <- readRDS(here::here("data/processed/spatial_objects/bio_AF_1km_crop.rds"))
df.bio.AF.1km <- readRDS(here::here("data/processed/spatial_objects/df_bio_AF_1km.rds"))

elev.AF.1km.crop <- readRDS(here::here("data/processed/spatial_objects/elev_AF_1km_crop.rds"))
df.elev.AF.1km <- readRDS(here::here("data/processed/spatial_objects/df_elev_AF_1km.rds"))
plot(elev.AF.1km.mask)

lulc.AF.1km.reclass <- readRDS(here::here("data/processed/spatial_objects/lulc_AF_1km_reclass.rds"))
df.lulc.AF.1km.reclass  <- readRDS(here::here("data/processed/spatial_objects/df_lulc_AF_1km_reclass.rds"))

# community information
env.bfly <- readRDS(here::here("data/processed/comm_objects/env.bfly.rds"))
head(env.bfly)

# separate the coordinates of each community in long lat format (important)
tmp <- env.bfly[, c("Longitude", "Latitude")] 

# put the coordinates in the same projection of data
points <- SpatialPoints(tmp, proj4string = lulc.AF.1km.reclass@crs) 
points

# plot(lulc.AF.1km.reclass)
# plot(points, add = T, pch = 19)

# 1 - extract the proportion of each LULC class for each community

e <- raster::extract(lulc.AF.1km.reclass, points, buffer = 1000)

class.counts <- lapply(e, table)
class.prop <- lapply(e, FUN = function(x) { prop.table(table(x)) })

source(here::here("R/functions/rbindFill_function.R"))
p.prop <- rbind.fill(class.prop)
prop.replace <- replace(x = p.prop, list = is.na(p.prop), values = 0)

colnames(prop.replace) 
# 1 = Natural Forest Formation
# 3 = Agriculture
# 2 = Open natural formations
# 0 = non informed
# 4 = Urbanization
# 6 = Water 
# 5 = Other non forest formations

colnames(prop.replace) <- c("Forest_Formation", "Agriculture", "Open_Formation",
                            "Non_informed", "Urbanization", "Water", "NonForest_Formation")
head(prop.replace)

# 2 - extracting the bioclim variables
require(raster)
points <- SpatialPoints(tmp, proj4string = bio.AF.1km.mask@crs) # put the coordinates in the same projection of data
points

values <- extract(bio.AF.1km.crop, points) # Extracting the bioclim of interest and for the coords sites
values

elev <- raster::extract(elev.AF.1km.crop, points) # Extracting the bioclim of interest and for the coords sites
elev


bio.sel <- values[,c("bio_1", "bio_4", "bio_12", "bio_15", "bio_18")]
colnames(bio.sel) <- c("Annual_temp", "Temp_seasonality", "Annual_prec",
                       "Prec_seasonality", "Prec_hot_qtr")
head(bio.sel)

# put all variable in one object
env.bfly <- cbind(env.bfly, prop.replace[,-4], bio.sel, elev)
colnames(env.bfly)

head(env.bfly)

rownames(env.bfly) <- env.bfly$Sites_ID

# 4 - TNC ecoregions

tnc.crop <- readRDS(here::here("data/processed/spatial_objects/TNC_1km.rds"))
df.tnc <- readRDS(here::here("data/processed/spatial_objects/df_TNC_1km.rds"))

ecor <- raster::extract(tnc.crop, points, buffer = 10) # Extracting the bioclim of interest and for the coords sites
(ecor)
# not all ecoregions have sampling communities:(

saveRDS(env.bfly, here::here("data/processed/comm_objects/Envir_covariates.rds"))

################################################################################


# Clean data --------------------------------------------------------------
# here we will remove NA and non-meaningful values for raster layers
# we will use this data to predict butterflies distribution in TAF

df.bio.1km <- readRDS(here::here("data/processed/spatial_objects/df_bio_AF_1km.rds"))
df.elev.1km <- readRDS(here::here("data/processed/spatial_objects/df_elev_AF_1km.rds"))
df.lulc.1km <- readRDS(here::here("data/processed/spatial_objects/df_lulc_AF_1km_reclass.rds"))

dim(df.bio.1km)
df.bio.1km <- na.omit(df.bio.1km) # remove NA
summary(df.bio.1km.s$bio_1)

# remove cells that have lower non-meaningful values
df.bio.1km.s <- subset(df.bio.1km, bio_1 > 0) 
head(df.bio.1km.s)

# elev
df.elev.1km.s <- df.elev.1km[match(rownames(df.bio.1km.s), rownames(df.elev.1km)), ]
colnames(df.elev.1km.s)[3] <- "Elevation"
summary(df.elev.1km.s$Elevation)

# LULC - by classes
df.lulc.1km.s <- df.lulc.1km[match(rownames(df.bio.1km.s), rownames(df.lulc.1km)),]
colnames(df.lulc.1km.s)

df.envir <- cbind(df.bio.1km.s[,c(1,2,3,12)], elev = df.elev.1km.s$Elevation, LULC = df.lulc.1km.s$code)
head(df.envir)

# LULC - buffer for all grid cells
lulc.AF.1km <- readRDS(here::here("data/processed/spatial_objects/lulc_AF_1km_reclass.rds"))

# separate the coordinates in long lat (important)
tmp <- df.envir[, c("x", "y")] 

# put the coordinates in the same projection of data
points <- SpatialPoints(tmp, proj4string = lulc.AF.1km@crs) 
points

e <- extract(lulc.AF.1km, points, buffer = 1000)

class.counts <- lapply(e, table)
class.prop <- lapply(e, FUN = function(x) { prop.table(table(x)) })


source(here::here("R/functions/rbindFill_function.R"))
p.prop <- rbind.fill(class.prop)
prop.replace <- replace(x = p.prop, list = is.na(p.prop), values = 0)

colnames(prop.replace)
prop.replace <- prop.replace[, c("1", "3", "2", "0", "4", "6", "5")]
colnames(prop.replace) <- c("Forest_Formation", "Agriculture", "Open_Formations",
                           "Non_informed", "Urbanization", "Water", "NonForest_Formation")

saveRDS(prop.replace, here::here("data/processed/spatial_objects/prop_LULC.rds"))

df.envir2 <- cbind(df.envir, prop.replace)
head(df.envir2)


# and put the ecoregion information in the environmental file
colnames(df.tnc)
unique(df.tnc$layer_ECO_NAME)

df.envir3 <- cbind(df.envir2, Ecor = df.tnc[match(rownames(df.envir2), rownames(df.tnc)), 5])
head(df.envir3)

df.envir3 <- na.omit(df.envir3)

saveRDS(df.envir3, here::here("data/processed/spatial_objects/Clean_data_envir.rds"))