# Script for read and manipulate the fruit-feeding butterfly data --------
# 1.community data (sites x species)
# 2.sites data (lat/long, and metadata for each community)
# 3. phylogenetic data for all species sampled
# 4. traits data (FWL)

# required packages
library(reshape)
library(ggplot2)

# spatial
library(adespatial)
library(rgdal)
library(raster)
library(sp)

#phylogenetic
library(ape)
library(geiger)
library(phytools)

## Fruit-feeding butterflies dataset ---------------------------------------
# information of species names for each community
bfly.species <- read.csv(here::here("data/raw/bfly/ATLANTIC_BFLY_species.csv"),
                         header = TRUE, sep = ";")
str(bfly.species)
unique(bfly.species$Study)

# information of metadata for each community (methodology employed, trap hours, coords)
bfly.sites <- read.csv(here::here("data/raw/bfly/ATLANTIC_BFLY_sites.csv"),
                       header = TRUE, sep = ";")
str(bfly.sites)

unique(bfly.sites$Study)

points.tmp <- data.frame(long = unique(bfly.sites$Longitude), lat = unique(bfly.sites$Latitude), 
                         Sites_ID = unique(bfly.sites$Sites_ID))
head(points.tmp)


# cleaning community data ----------------------------------------------------

# removing communities out of range of the AF trinational
af.shape <- readOGR(here::here("data/raw/WGS84/MapBiomas_c1/Limits_AFtrinacional.shp"))

# ggplot() + geom_polygon(data = af.shape, aes(x = long, y = lat, group = group), 
#                         colour = "black", fill = NA)

comm.points <- SpatialPoints(as.matrix(points.tmp[,1:2]), proj4string = af.shape@proj4string) # put the coordinates in the same projection of data

source(here::here("R/functions/find_outliers_function.R"))
rm.points <- find.outliers(points = comm.points, poly = af.shape)
rm.points

points.tmp$ID <- cbind(seq(1:nrow(points.tmp)))
points.tmp$status <- ifelse(points.tmp$ID %in% rm.points$row_idx, "out", "in")
head(points.tmp)
# plot(af.shape)

p.in <- subset(points.tmp, status == "in")
p.in <- SpatialPoints(p.in[,1:2], proj4string = af.shape@proj4string)
# plot(p.in, add = T, col = "blue", pch = 19)

p.out <- subset(points.tmp, status == "out")
p.out <- SpatialPoints(p.out[,1:2], proj4string = af.shape@proj4string)
# plot(p.out, add = T, col = "red", pch = 1)

env.bfly <- bfly.sites[subset(points.tmp, status == "in")$ID, -2]

saveRDS(env.bfly, here::here("data/processed/data/Final_objects/env.bfly.rds"))


# creating the community data ---------------------------------------------
# merging both data
match(unique(bfly.species$Study), unique(bfly.sites$Study))

bfly.data <- merge(x = bfly.species, y = env.bfly, by = c("Sites_ID", "Study"), all.y = T)
head(bfly.data)

# create a community matrix
bfly.data$Occ <- rep(1, dim(bfly.data)[1])

junk.melt <- melt(bfly.data, id.var = c("Species", "Sites_ID"), 
                  measure.var = "Occ")
Y <- cast(junk.melt, Sites_ID ~ Species)
tail(Y)
colnames(Y)
rownames(Y) <- Y$Sites_ID

# remove the first column
Y <- Y[,-1]
# check if are only 0/1 values
which(Y > 1, arr.ind = T)
rowSums(Y)
table(colSums(Y))

comm.bfly <- ifelse(Y >= 1, 1, 0)
max(comm.bfly)
head(comm.bfly)


#####################################################################################
# read traits data ----

bfly.traits <- read.csv(here::here("data/raw/bfly/Bfly_TraitsBase.csv"), sep = ";")
head(bfly.traits)

# removing spp without traits information
# FWL is the forewing length and is a proxy for dispersion and body length
bfly.traits <- subset(bfly.traits, !is.na(FWL))
str(bfly.traits)
unique(bfly.traits$Species)
sort(bfly.traits$WTR)

library(dplyr)
bfly_grupo <- group_by(bfly.traits, Species)
bfly_grupo %>% str()
bfly.fwl <- data.frame(FWL = summarise(bfly_grupo, mean(FWL)),
                       FWW = summarise(bfly_grupo, mean(FWW, na.rm = T))[,2],
                       TW = summarise(bfly_grupo, mean(TW, na.rm = T))[,2],
                       AW = summarise(bfly_grupo, mean(AW, na.rm = T))[,2],
                       AR = summarise(bfly_grupo, mean(as.numeric(AR), na.rm = T))[,2],
                       WTR = summarise(bfly_grupo, mean(as.numeric(WTR), na.rm = T))[,2])
bfly.fwl
colnames(bfly.fwl) <- c("Species", "FWL", "FWW", "TW", "AW", "AR", "WTR")

rm.tr <- setdiff(bfly.fwl$Species, colnames(comm.bfly)) # present in traits but not in comm
rm.cm <- setdiff(colnames(comm.bfly), bfly.fwl$Species) # present in comm but not in traits

comm.bfly <- comm.bfly[,-match(rm.cm, colnames(comm.bfly))]
dim(comm.bfly)

bfly.fwl <- bfly.fwl[-match(rm.tr, bfly.fwl$Species), ]
rownames(bfly.fwl) <- bfly.fwl$Species
bfly.fwl <- bfly.fwl[,-1]
dim(bfly.fwl)

match(rownames(bfly.fwl), colnames(comm.bfly))

saveRDS(comm.bfly, here::here("data/processed/comm_objects/comm.bfly.rds"))
saveRDS(bfly.fwl, here::here("data/processed/comm_objects/bfly_trait.rds"))


# Phylogenetic tree for butterflies (Chazot et al. 2021) ------------------

bfly.tree <- read.tree(here::here("data/raw/bfly/Nymphalidae_tree_Chazot2021.new"))
(bfly.tree)

# Manipulate the phylogenetic tree to add genus Amiga -----------------------
# based on Espeland et al 2019

# find the position to add the genus as sister-group
find.node <- fastMRCA(bfly.tree, "Megeuptychia_monopunctata", "Taydebis_peculiaris")
plot(extract.clade(bfly.tree, find.node))
posit <- (bfly.tree$edge.length[sapply(find.node, function(x, y) which(y == x),
                                       y = bfly.tree$edge[,2])])/2 # add a new genus to tree cutting the branch in two equal parts

## adding the genus "Amiga"
bfly.tree <- bind.tip(bfly.tree, "Amiga_arnaca", where = find.node, position = posit)

# add the genus Stegosatyrus as sister group "Cissia" and "Megisto"
# find the position to add the genus as sister-group
bfly.tree$tip.label

find.node <- fastMRCA(bfly.tree, "Megisto_rubricata", "Moneuptychia_soter")
plot(extract.clade(bfly.tree, find.node))

## adding the genus "Stegosatyrus"
bfly.tree <- bind.tip(bfly.tree, "Stegosatyrus_periphas", where = find.node, position = 1)


# Updating some genus names ----

tips <- bfly.tree$tip.label

sep.names <- strsplit(tips, split = "_")

df.tips <- matrix(unlist(sep.names), ncol = 2, nrow = length(tips), byrow = T)

df.tips[which(df.tips[,1] == "Moneuptychia"), ]

# Paulogramma to Catagramma
df.tips[which(df.tips[,1] == "Paulogramma"), 1] <- "Catagramma"

# Paryphthimoides phronius to Cissia phronius
df.tips[which(df.tips[,1] == "Paryphthimoides"), 1][5] <- "Cissia"

# Guaianaza pronophila to Frosterinaria pronophila
df.tips[which(df.tips[,1] == "Guaianaza"), 1] <- "Forsterinaria"

# Moneuptychia paeon and M. griseldis to Carminda paeon and C. griseldis
df.tips[which(df.tips[,1] == "Moneuptychia"), 1][8:9] <- "Carminda"

tips <- paste(df.tips[,1], "_", df.tips[,2], sep = "")

# Save this tree as the current bfly.tree
bfly.tree$tip.label <- tips

bfly.tree$tip.label[which(bfly.tree$tip.label == "Carminda_paeon")]

# Add missing species to backbone tree
df.bfly <- data.frame(spp = colnames(comm.bfly), row.names = colnames(comm.bfly))

source(here::here("R/functions/function_treedata_modif.R"))

# Insert the another species to the tree, as polytomies
spp.insert <- treedata_modif(bfly.tree, df.bfly)$nc$data_not_tree

for (i in 1:length(spp.insert)) {
  bfly.tree <- phytools::add.species.to.genus(bfly.tree, spp.insert[i])
}
bfly.tree

# saving in the folder
write.tree(bfly.tree, here::here("data/processed/Nymphalidae_2952spp.new"))

# Cutting the bfly.tree to my specific sampling pool to Atlantic trinational forest
phy.bfly <- treedata_modif(bfly.tree, df.bfly)$phy
phy.bfly$tip.label

saveRDS(phy.bfly, here::here("data/processed/comm_objects/tree_bfly_AF.rds"))
