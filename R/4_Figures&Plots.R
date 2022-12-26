# Figure plots
library(ggtreeExtra)
library(ggtree)
library(Hmsc)
library(tidyverse)
library(phytools)
library(corrplot)
library(ggcorrplot)
library(ggstance)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)

# Fig 1 - TAF limits and comminity location ----

grid <- readRDS(here::here("data/processed/spatial_objects/Clean_data_envir.rds"))
head(grid)
comm.points <- readRDS(here::here("data/processed/comm_objects/env.bfly.rds"))

# Drawing maps programmatically with R, sf and ggplot2 
theme_set(theme_bw())

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

nt <- subset(world, region_wb == "Latin America & Caribbean")

gg.nt <- ggplot(data = nt) +
  geom_sf(fill = "gray90") +
  geom_rect(xmin = -58, xmax = -34, ymin = -34, ymax = -2.5, 
            fill = NA, colour = "darkred") +
  annotate(geom = "text", x = -50, y = 20, label = "Atlantic Ocean", 
           fontface = "italic", color = "grey22", angle = -20) +
  annotate(geom = "text", x = -90, y = -20, label = "Pacific Ocean", 
           fontface = "italic", color = "grey22") + 
  coord_sf(datum = NA) + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA))
gg.nt


# map of the AF location
p1 <- ggplot(data = nt) +
  geom_sf(fill = "gray90", col = "black") +
  coord_sf(xlim = c(-60, -34), ylim = c(-35, 5), expand = T, 
           crs = st_crs(4326), )
p1

p2 <- p1 + annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -45, y = -28, label = "Atlantic Ocean", 
           fontface = "italic", color = "grey22", angle = 45) +
  labs(x = "Longitude", y = "Latitude")#, title = "Trinational Atlantic Forest")
p2

p3 <- p2 + geom_tile(data = grid, aes(x = x, y = y, fill = Annual_temp), fill = "forestgreen", 
                     alpha = 0.7) +
  geom_point(data = comm.points, aes(x = Longitude, y = Latitude), col = "darkred") 
p3

p.AF <- p3 + annotation_custom(
  grob = ggplotGrob(gg.nt),
  xmin = -45,
  xmax = -61,
  ymin = -10,
  ymax = 6)
p.AF
cowplot::save_plot(filename = here::here("output/figures/Fig1_TAF_limits.png"), plot = p.AF,
                   base_height = 8, base_width = 8)


# Preparing the figure 2 - Environmental variables
colnames(grid)

p.temp <-  ggplot(data = nt) +
  geom_sf(fill = "gray90", col = "black") +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = grid, aes(x = x, y = y, fill = bio_1)) + #, alpha = .7) +
  scale_fill_viridis_c(name = "Annual Mean\nTemperature") +
  labs(x = "Longitude", y = "Latitude", tag = "a)") +
  theme(legend.position = c(.8, .2), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.background = element_blank(), legend.key.size = unit(1, "lines"))

p.prec <-  ggplot(data = nt) +
  geom_sf(fill = "gray90", col = "black") +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = grid, aes(x = x, y = y, fill = bio_18)) + #, alpha = .7) +
  scale_fill_viridis_c(name = "Precipitation of \nWarmest Quarter") +
  labs(x = "Longitude", y = "Latitude", tag = "b)") +
  theme(legend.position = c(.8, .2), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.background = element_blank(), legend.key.size = unit(1, "lines"))

p.elev <-  ggplot(data = nt) +
  geom_sf(fill = "gray90", col = "black") +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = grid, aes(x = x, y = y, fill = elev)) + #, alpha = .7) +
  scale_fill_viridis_c(name = "Elevation") +
  labs(x = "Longitude", y = "Latitude", tag = "c)") +
  theme(legend.position = c(.8, .2), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.background = element_blank(), legend.key.size = unit(1, "lines"))

# we will plot the landscape variables as LULC classes
df.lulc.1km <- readRDS(here::here("data/processed/spatial_objects/df_lulc_AF_1km_reclass.rds"))
df.lulc.1km <- df.lulc.1km[-which(df.lulc.1km$code == "NA"),]
df.lulc.1km <- df.lulc.1km[-which(df.lulc.1km$code == "Water"),]
df.lulc.1km <- df.lulc.1km[-which(df.lulc.1km$code == "NonForest_Formation"),]

colnames(df.lulc.1km)
unique(df.lulc.1km$code)

p.lulc <-  ggplot(data = nt) +
  geom_sf(fill = "gray90", col = "black") +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = df.lulc.1km, aes(x = x, y = y, fill = code)) + #, alpha = .7) +
  scale_fill_viridis_d(direction = -1, name = "Land Use and \nLand Cover", labels = c("Agriculture", 
                                                                                      "Natural Forest Formation",
                                                                                      "Natural NonForest\nFormations",
                                                                                      "Natural Open Formations",
                                                                                      "Urban Area", "Water")) +
  labs(x = "Longitude", y = "Latitude", tag = "d)") +
  theme(legend.position = c(.8, .2), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.background = element_blank(), legend.key.size = unit(1, "lines"))
windows()
cowplot::plot_grid(p.temp, p.prec, p.elev, p.lulc, ncol = 2)

cowplot::save_plot(here::here("output/figures/Fig2_EnvirCovariates.png"),
                   cowplot::plot_grid(p.temp, p.prec, p.elev, p.lulc, ncol = 2),
                  base_height = 12, base_width = 12)

# read the HMSC object
m.FULL <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))[[1]]

# read the species data
bfly.species <- read.csv(here::here("data/raw/bfly/ATLANTIC_BFLY_species.csv"),
                         header = TRUE, sep = ";")
colnames(bfly.species)

bfly.species$Colors <- ifelse(bfly.species$Subfamily == "Biblidinae", "#0072B2",
                              ifelse(bfly.species$Subfamily == "Charaxinae", "#44AA99",
                                     ifelse(bfly.species$Subfamily == "Nymphalinae", "#CC79A7", "#E69F00")))
levels(as.factor(bfly.species$Colors))

tree <- m.FULL$phyloTree

# FIG S1: prevalence and species richness -----------------------------------
P <- readRDS(here::here("output/Final_objects/Summary_obs.rds"))[[1]]
S <- readRDS(here::here("output/Final_objects/Summary_obs.rds"))[[2]]

p.prev <- ggplot(as.data.frame(P), aes(x = P)) + 
  geom_density(alpha = .5,  color = "#009E73", fill = "#009E73") +
  labs(x = "Species prevalence", tag = "a)")

p.rich <- ggplot(as.data.frame(S), aes(x = S)) + 
  geom_density(alpha = .5,  color = "#CC79A7", fill = "#CC79A7") +
  labs(x = "Species Richness", tag = "b)")

cowplot::save_plot(here::here("output/figures/FigS1_Summary_res.png"), 
                   cowplot::plot_grid(p.prev, p.rich),
                   base_height = 4, base_width = 8)

# FIG S2 - MCMC convergence ------------------------------------------------
mcmc.all <- readRDS(here::here("output/Final_objects/psrf_all.rds"))
reshape2::melt(mcmc.all$Beta[[1]]) # model full

mypal <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", 
           "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")

scales::show_col(mypal)

mcmc.beta <- ggplot(reshape2::melt(mcmc.all$Beta[[1]]), aes(x = value, colour = Var2, fill = Var2)) + 
  geom_density(alpha = .5) +
  scale_color_manual(values = mypal[c(3,8)], name = NULL) +
  scale_fill_manual(values = mypal[c(3,8)], name = NULL) +
  geom_vline(xintercept = 1.1, color = "red", linetype = 2) +
  labs(x = expression(paste("Potential scale reduction factor ", beta, sep = "")), 
       tag = "a)") + theme(legend.position = "none")
mcmc.beta

mcmc.gamma <- ggplot(reshape2::melt(mcmc.all$Gamma[[1]]), aes(x = value, colour = Var2, fill = Var2)) + 
  geom_density(alpha = .5) +
  scale_color_manual(values = mypal[c(3,8)], name = NULL) +
  scale_fill_manual(values = mypal[c(3,8)], name = NULL) +
  geom_vline(xintercept = 1.1, color = "red", linetype = 2) +
  labs(x = expression(paste("Potential scale reduction factor ", gamma, sep = "")), 
       tag = "b)") + theme(legend.position = "none")
mcmc.gamma

mcmc.omega <- ggplot(reshape2::melt(mcmc.all$Omega[[1]]), aes(x = value, colour = Var2, fill = Var2)) + 
  geom_density(alpha = .5) +
  scale_color_manual(values = mypal[c(3,8)], name = NULL) +
  scale_fill_manual(values = mypal[c(3,8)], name = NULL) +
  geom_vline(xintercept = 1.1, color = "red", linetype = 2) +
  labs(x = expression(paste("Potential scale reduction factor ", Omega, sep = "")), 
       tag = "c)") + theme(legend.position = "none")
mcmc.omega

leg <- cowplot::get_legend(mcmc.omega + 
                      theme(legend.position = "left"))

cowplot::save_plot(here::here("output/figures/FigS2_MCMC_convergence.png"),
                   cowplot::plot_grid(mcmc.beta, mcmc.gamma, mcmc.omega, leg,
                                      rel_widths = c(1,1,1, .3), ncol = 4),
                   base_height = 4, base_width = 10)

mcmc.all$Rho[[1]]
#         Point est. Upper C.I.
# [1,]    1.06391   1.214543
# good

mcmc.all$Alpha[[1]]

#                   Point est. Upper C.I.
# Alpha1[factor1]   1.298465   2.196262
# Alpha1[factor2]   1.021358   1.070310
# Alpha1[factor3]   1.068746   1.229371
# Alpha1[factor4]   1.222058   1.643737
# Alpha1[factor5]   1.264517   1.913944
# Alpha1[factor6]   1.440652   4.795983

# you should change the mcmc.all$Beta[[1]] for [[2]] or [[3]] to evaluate MCMC convergence
# for spatial and environmental models


# Environmental filter (beta parameters) ----------------------------------
# FIG 3: Beta parameters -------------------------------------------------

postBeta.full <- getPostEstimate(m.FULL, parName = "Beta")
plotBeta(m.FULL, postBeta.full, param = "Support", SpeciesOrder = "Tree")

betaP <- postBeta.full$support
supportLevel <- .7
toPlot <- 2 * betaP - 1
toPlot <- toPlot * ((betaP > supportLevel) + (betaP < 
                                                (1 - supportLevel)) > 0)

df.betas <- as.data.frame(t(toPlot))
colnames(df.betas) <- m.FULL$covNames
head(df.betas)

df.betas <- df.betas[match(tree$tip.label, rownames(df.betas)),]

df.betas$label <- rownames(df.betas)
df.betas$Subfamily <- bfly.species[match(rownames(df.betas), bfly.species$Species), "Subfamily"]
df.betas$Tribe <- bfly.species[match(rownames(df.betas), bfly.species$Species), "Tribe"]
head(df.betas)
unique(df.betas$Tribe)

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Satyrini = subset(df.betas, Tribe == "Satyrini")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Morphini = subset(df.betas, Tribe == "Morphini")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Brassolini = subset(df.betas, Tribe == "Brassolini")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Haeterini = subset(df.betas, Tribe == "Haeterini")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Charaxinae = subset(df.betas, Subfamily == "Charaxinae")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Biblidinae = subset(df.betas, Subfamily == "Biblidinae")[,"label"]))

tree <- makeNodeLabel(phy = tree, method = "u",
                      nodeList = list(Nymphalinae = subset(df.betas, Subfamily == "Nymphalinae")[,"label"]))

tree$node.label

phy <- dplyr::full_join(tree, df.betas, by = "label")

mypal <- c("#DC0000B2", "#52854C", "#4E84C4",  "#00A087B2", "#E64B35B2", 
           "#F4EDCA", "#D16103", "#999933")

scales::show_col(mypal)

# The circular layout tree.
# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
labdf <- data.frame(node = nodeid(tree, tree$node.label[nchar(tree$node.label)>1]), 
                    labs = tree$node.label[nchar(tree$node.label)>1])
labdf

p <- ggtree(tree, layout = "fan", size = 0.15, open.angle = 10) +
  geom_hilight(data = labdf, mapping = aes(node = node, fill = labs, colour = labs),
               alpha = 0.3,
               size = 0.05) + labs (fill = "Clade", colour = "Clade") +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal)
p 

colnames(df.betas)[1:8] <- c("Int", "For", "Open", "Agri", "Urb", "Temp", "Prec", "Elev")
dat2 <- reshape2::melt(df.betas)
head(dat2)

# 
head(df.betas)
df.betas$new <- ifelse(df.betas$Subfamily == "Satyrinae", paste(df.betas$Tribe), 
                       paste(df.betas$Subfamily))
df.betas$new <- as.factor(df.betas$new)
levels(df.betas$new) <- c("Nymphalinae", "Biblidinae", "Charaxinae", "Brassolini",
                          "Morphini", "Haeterini", "Melanitini", "Satyrini")
levels(df.betas$new)

ggplot(df.betas, aes(x = Int, y = For, colour = new)) + geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & For < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & For < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2)

int.for <- ggplot(df.betas, aes(x = Int, y = For, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & For < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & For < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Forest Formation")
int.for

int.open <- ggplot(df.betas, aes(x = Int, y = Open, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Open < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Open < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Natural Open Formation")
int.open

int.agri <- ggplot(df.betas, aes(x = Int, y = Agri, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Agri < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Agri < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Agriculture")
int.agri

int.urb <- ggplot(df.betas, aes(x = Int, y = Urb, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Urb < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Urb < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Urbanization")
int.urb

int.temp <- ggplot(df.betas, aes(x = Int, y = Temp, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Temp < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Temp < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Annual mean temperature")
int.temp

int.prec <- ggplot(df.betas, aes(x = Int, y = Prec, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Prec < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Prec < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Precipitation of the warmest quarter")
int.prec

int.elev <- ggplot(df.betas, aes(x = Int, y = Elev, fill = new, shape = new)) +
  geom_point() +
  geom_text(data = filter(df.betas, Int < 0 & Elev < 0), 
            aes(label = gsub( "_", " ", filter(df.betas, Int < 0 & Elev < 0)[,"label"])), 
            vjust = "inward", hjust = "inward", color = "black", size = 2) +
  scale_fill_manual(values = mypal, name = NULL) +
  scale_shape_manual(values = c(21,21,21,24,24,24,24,24), name = NULL) +
  theme(legend.position = "none") + labs(x = "Mean occurrence probability", 
                                         y = "Elevation")
int.elev
leg <- cowplot::get_legend(int.elev + theme(legend.position = "right"))

# Fig S7 - variation in responses of rares and common species
cowplot::plot_grid(int.for, int.open, int.agri, int.urb,
                   int.temp, int.prec, int.elev, leg, ncol = 4)

cowplot::save_plot(here::here("output/figures/FigS7_Int_betas.png"),
                   cowplot::plot_grid(int.for, int.open, int.agri, int.urb,
                                      int.temp, int.prec, int.elev, leg, ncol = 4),
                   base_height = 8, base_width = 12)


p.beta <- p + new_scale_fill() +
  geom_fruit(data = dat2, geom = geom_tile,
             mapping = aes(y = label, x = variable, fill = value),
             color = "white", offset = 0.03, size = 0.1, pwidth = .85,
             axis.params = list(axis = "x", text.size = 4, vjust = 1)) +
  scale_fill_gradient2() + labs (fill = "Support")
p.beta + layout_rectangular()

dat3 <- as.data.frame(m.FULL$Tr[,2:3])
dat3$labs <- rownames(dat3)
dat3 <- dat3[match(tree$tip.label, rownames(dat3)),]

p.beta.tr <- p.beta + new_scale_fill() +
  geom_fruit(data = dat3, geom = geom_bar,
             mapping = aes(y = labs, x = AR),
             pwidth = 0.4, offset = .1, color = "gray50",
             orientation = "y", fill = "white",
             stat="identity", axis.params = list(axis = "x", text.size = 4, vjust = 1)) +
  annotate(geom = "text", x = 155, y = 263, label = "AR",
           hjust = "center", size = 5)
p.beta.tr + layout_rectangular() 

pfim <- p.beta.tr + new_scale_fill() +
  geom_fruit(data = dat3, geom = geom_bar,
             mapping = aes(y = labs, x = FWL),
             pwidth = 0.5, offset = .1, color = "gray50",
             orientation = "y", fill = "white",
             stat="identity", axis.params = list(axis = "x", text.size = 4, vjust = 1)) +
  geom_treescale(fontsize = 2, linesize =.3, x = 4.9, y = 0.1, width = 20) + 
  annotate(geom = "text", x = 202, y = 263, label = "FWL",
           hjust = "center", size = 5)
p.niche <- pfim + layout_rectangular()
p.niche + theme(legend.position = "none")

############################################################################
# Traits and phylogenetic signal ------------------------------------------
############################################################################

# FIG 4: Gamma parameter  -------------------------------------------------

supportLevel = .9

post <- getPostEstimate(m.FULL, parName = "Gamma")
plotGamma(m.FULL, post = post, supportLevel = 0.9)

gammaP = post$support
toPlot = 2 * gammaP - 1
toPlot = toPlot * ((gammaP > supportLevel) + (gammaP < 
                                                (1 - supportLevel)) > 0)
colnames(toPlot) <- m.FULL$trNames
rownames(toPlot) <- c("Int", "For", "Open", "Agri", "Urb", "Temp", "Prec", "Elev")

post.gammas <- reshape2::melt(t(toPlot))

p.gamma <- ggplot(data = post.gammas, aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = value), size = 1, color = "black") +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradient2(limit = c(-1, 1)) +  coord_fixed(.5) +
  theme_classic()
p.gamma

leg <- cowplot::get_legend(p)
leg2 <- cowplot::get_legend(p.gamma)

f.leg <- cowplot::plot_grid(leg, leg2)

p2 <- cowplot::plot_grid(p.gamma + theme(legend.position = "none") + labs(tag = "b)"), f.leg, ncol = 1,
                         rel_heights = c(2,1))

p.BetaTr <- cowplot::plot_grid(p.niche + theme(legend.position = "none") + labs(tag = "a)"), p2,
                   rel_widths = c(4, 1))
p.BetaTr

cowplot::save_plot(here::here("output/figures/Figure4_nicheTR.png"), p.BetaTr,
                   base_height = 10, base_width = 14)


# RHO - phylogenetic signal and 
# Alpha - spatial signal
mpost <- convertToCodaObject(m.FULL)
rho <- round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)
alpha <- round(summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]],2)


rho
alpha


# Residual species association --------------------------------------------
# FIG S3: Omega parameters --------------------------------------------

OmegaCor <- computeAssociations(m.FULL)
supportLevel <- 0.9
toPlot <- ((OmegaCor[[1]]$support > supportLevel)
           + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
toPlot <- sign(toPlot)
plotOrder <- match(m.FULL$phyloTree$tip.label, rownames(toPlot))
toPlot <- toPlot[plotOrder,plotOrder]
colnames(toPlot)

tmp <- as.data.frame(matrix(unlist(strsplit(colnames(toPlot), "_")), ncol = 2, 
       nrow = ncol(toPlot), byrow = T))
tmp$new <- paste(substr(tmp$V1, 1, 1), ". ", tmp$V2, sep = "")
head(tmp)

colnames(toPlot) <- tmp$new
rownames(toPlot) <- tmp$new

rm.names <- match(rownames(toPlot)[which(rowSums(toPlot) == 1)], rownames(toPlot))

p.teste <- ggcorrplot(toPlot[-rm.names, -rm.names], hc.order = TRUE, type = "lower",
         outline.col = "white", ggtheme = ggplot2::theme_classic, tl.cex = 6)

p.omega <- p.teste + scale_fill_gradient2(limit = c(-1, 1)) +# coord_fixed() +
  theme(axis.text.y = element_text(size = 4, hjust = 1), panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90, size = 4, hjust = 1, vjust = 0.5), legend.position = "bottom") 
p.omega

cowplot::save_plot(here::here("output/figures/FigS4_Omega.png"), 
          p.omega,
          base_width = 10, base_height = 10)


# FIG 5: Variance partitioning ---------------------------------------------------
# To be able to group the environmental variables, we look at the design matrix X that 
# Hmsc has constructed by applying the XFormula to the XData.

colnames(m.FULL$X)

# We observe that the columns 2-5 relate to land use and land cover (LULC)
# and columns 6-7 to macroclimatic variation, 8 is a topographic variable.
# Arbitrarily, we include the intercept in the LULC variables, and thus perform the 
# variance partitioning with the following grouping of the columns of the X matrix.

groupnames <- c("LULC","Bioclim", "elev")
group <- c(1,1,1,1,1,2,2,3)

VP.full <- computeVariancePartitioning(m.FULL, group = group, groupnames = groupnames)
plotVariancePartitioning(m.FULL, VP.full)

VP.full$vals <- VP.full$vals[,match(m.FULL$phyloTree$tip.label, 
                                    colnames(VP.full$vals))]

head(VP.full$vals)

df.VP <- reshape2::melt(VP.full$vals)
df.VP <- df.VP[,c(2,3,1)]
head(df.VP)

p <- ggtree(tree, size = 0.15) +
  geom_treescale(fontsize = 2, linesize =.3, x = 0, y = 0.1, width = 20)

p.facet <- facet_plot(p,
                 panel = "Variance Partitioning",
                 data = df.VP,
                 geom = geom_barh,
                 mapping = aes(x = value, fill = Var1),
                 stat = "identity")

prop.var <- round(rowMeans(VP.full$vals), 4)*100
new.labs <- c(paste("LULC (", prop.var[1], "%)", sep = ""),
              paste("Climatic (", prop.var[2], "%)", sep = ""),
              paste("Elevation (", prop.var[3], "%)", sep = ""),
              paste("Random (", prop.var[4], "%)", sep = ""))

mypal <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", 
  "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
scales::show_col(mypal)

p3 <- p.facet + scale_fill_manual(name = NULL, values = mypal[c(1,3,6,9)], labels = new.labs) + 
  labs(x = NULL, y = NULL) 

p.VP <- p3 + theme(strip.background = element_blank(), strip.text = element_blank())


cowplot::save_plot(here::here("output/figures/FigS3_VarPart.png"), p.VP, 
          base_height = 8, base_width = 10)

# Accessing the proportion of explanation of traits on covariates and occurrence/abundance of species
VP.full$R2T$Beta
round(sum(VP.full$R2T$Beta)*100, 2)
VP.full$R2T$Y*100


# Community-level responses -----------------------------------------------

# We start by making gradient plots that visualize how the communities vary among the 
# environmental variables.
source(here::here("R/functions/plotGradient_modify_function.R"))

m.FULL$XFormula

# Proportion of Forest Formation
Gradient <- constructGradient(m.FULL, focalVariable = "Forest_Formation",
                              non.focalVariables = list("Open_Formation" = list(1), 
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

# Species richness
s.for <- plotGradient_modify(hM = m.FULL, Gradient = Gradient, predY = predY, 
                    measure = "S")

fwl.for <- plotGradient_modify(hM = m.FULL, Gradient = Gradient, predY = predY, 
                    measure = "T", index = 2)

ar.for <- plotGradient_modify(hM = m.FULL, Gradient = Gradient, predY = predY, 
                    measure = "T", index = 3)

p.for <- cowplot::plot_grid(s.for + xlab("Forest"), fwl.for + xlab("Forest"),
                            ar.for + xlab("Forest"), ncol = 1)

# Open formations
Gradient <- constructGradient(m.FULL, focalVariable = "Open_Formation",
                              non.focalVariables = list("Forest_Formation" = list(1), 
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

# Species richness
s.open <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", 
                              showData = TRUE)

# trait
fwl.open <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.open <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.open <- cowplot::plot_grid(s.open + ylab(NULL) + xlab("Open"), 
                             fwl.open + ylab(NULL) + xlab("Open"),
                             ar.open + ylab(NULL) + xlab("Open"), ncol = 1)

# urbanization
Gradient <- constructGradient(m.FULL, focalVariable = "Urbanization",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

s.urb <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

fwl.urb <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.urb <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.urb <- cowplot::plot_grid(s.urb + ylab(NULL), fwl.urb + ylab(NULL),
                            ar.urb + ylab(NULL), ncol = 1)

# Agriculture
Gradient <- constructGradient(m.FULL, focalVariable = "Agriculture",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

s.agri <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

fwl.agri <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.agri <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.agri <- cowplot::plot_grid(s.agri + ylab(NULL), fwl.agri + ylab(NULL),
                             ar.agri + ylab(NULL), ncol = 1)

# Annual_temp
Gradient <- constructGradient(m.FULL, focalVariable = "Annual_temp",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

s.temp <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

fwl.temp <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.temp <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.temp <- cowplot::plot_grid(s.temp + ylab(NULL) + xlab("Temperature"),
                             fwl.temp + ylab(NULL) + xlab("Temperature"), 
                             ar.temp + ylab(NULL) + xlab("Temperature"), ncol = 1)

# Prec_hot_qtr
Gradient <- constructGradient(m.FULL, focalVariable = "Prec_hot_qtr",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

s.prec <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

fwl.prec <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.prec <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.prec <- cowplot::plot_grid(s.prec + ylab(NULL) + xlab("Precipitation"), 
                             fwl.prec + ylab(NULL) + xlab("Precipitation"),
                             ar.prec + ylab(NULL) + xlab("Precipitation"), ncol = 1)

# elevation
Gradient <- constructGradient(m.FULL, focalVariable = "elev",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

s.elev <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

fwl.elev <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 2, showData = TRUE)

ar.elev <- plotGradient_modify(m.FULL, Gradient, pred = predY, measure = "T",
             index = 3, showData = TRUE)

p.elev <- cowplot::plot_grid(s.elev + ylab(NULL) + xlab("Elevation"),
                             fwl.elev + ylab(NULL) + xlab("Elevation"),
                             ar.elev + ylab(NULL) + xlab("Elevation"), ncol = 1)

cowplot::plot_grid(p.for, p.open, p.urb, p.agri, p.temp, p.prec, p.elev, nrow = 1)

cowplot::save_plot(here::here("output/figures/Fig4_Comm_response.png"), 
                   cowplot::plot_grid(p.for, p.open, p.urb, p.agri, p.temp,
                                      p.prec, p.elev, nrow = 1),
                   base_height = 6, base_width = 12)


# Spatial predictions -----------------------------------------------------

# by subfamily/tribe

mapData <- readRDS(here::here("data/processed/comm_objects/spatial_pred.rds"))
colnames(mapData)

# Mean occ for subfamilies (we will sum the mean occ for all species of the subfamily/tribe)
mapData$Biblidinae <- rowSums(mapData[, c(3:53)])/length(mapData[, c(3:53)])
mapData$Nymphalinae <- rowSums(mapData[, c(54:59)])/length(mapData[, c(54:59)])
mapData$Charaxinae <- rowSums(mapData[, c(60:93)])/length(mapData[, c(60:93)])
mapData$Morphini <- rowSums(mapData[, c(94:105)])/length(mapData[, c(94:105)])
mapData$Brassolini <- rowSums(mapData[, c(106:138)])/length(mapData[, c(106:138)])
mapData$Haeterini <- rowSums(mapData[, c(139:142)])/length(mapData[, c(139:142)])
mapData$Satyrini <- rowSums(mapData[, c(144:260)])/length(mapData[, c(144:260)])

# Drawing maps programmatically with R, sf and ggplot2 
theme_set(theme_bw())

world <- ne_countries(scale = "medium", returnclass = "sf")

nt <- subset(world, region_wb == "Latin America & Caribbean")

p.bib <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Biblidinae)) + 
  labs(title = "Biblidinae", x = " ", y = "Latitude") +
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.bib

p.nym <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Nymphalinae)) + 
  labs(title = "Nymphalinae", x = " ", y = " ") + 
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.nym

p.cha <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Charaxinae)) + 
  labs(title = "Charaxinae", x = " ", y = " ") +
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.cha

p.mor <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Morphini)) + 
  labs(title = "Morphini", x = " ", y = " ") +
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.mor

p.bra <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Brassolini)) + 
  labs(title = "Brassolini", x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.bra

p.hae <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Haeterini)) + 
  labs(title = "Haeterini", x = "Longitude", y = " ") +
 scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")
p.hae

p.sat <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = Satyrini)) + 
  labs(title = "Satyrini", x = "Longitude", y = " ") + 
  scale_fill_viridis_c(limits = c(0,1)) +
  theme(legend.position = "none")

leg <- cowplot::get_legend(p.sat + theme(legend.position = "right", 
                                         legend.title = element_blank()))

# Fig S4 - Subfamily/tribe mean occurrence probability
cowplot::plot_grid(p.mor, p.bra, p.hae, p.sat,
                   p.bib, p.nym, p.cha, leg, ncol = 4)

cowplot::save_plot(here::here("output/figures/FigS5_Spat_Sfam_Tribe.png"), 
                   cowplot::plot_grid(p.bib, p.nym, p.cha, leg,
                                      p.bra, p.mor, p.hae, p.sat,
                                      ncol = 4), base_height = 10, 
                   base_width = 12)

# We next plot predicted species richness

p.S <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = S)) + 
  labs(title = "Species Richness", x = "Longitude", y = "Latitude") + 
  scale_fill_viridis_c(name = NULL)
p.S

# We next plot the community-weighted mean forewing length
p.FWL <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = FWL)) + 
  labs(title = "Forewing length", x = "Longitude", y = " ", tag = "b)") + 
  scale_fill_viridis_c(name = NULL) 
p.FWL

# Aspect ratio
p.AR <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = mapData, aes(x = x, y = y, fill = AR)) + 
  labs(title = "Aspect ratio", x = "Longitude", y = " ", tag = "c)") + 
  scale_fill_viridis_c(name = NULL) 
p.AR

# taking into account the uncertainty
df.ric <- readRDS(here::here("data/processed/comm_objects/Mean_sd_richness.rds"))
head(df.ric)

# Posterior mean of the expected richness
p.ES <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = df.ric, aes(x = Longitude, y = Latitude, fill = ES)) + 
  labs(title = "Expected Richness", tag = "a)") + 
  scale_fill_viridis_c(name = NULL)
p.ES

# FIG 7 - predicted Richness and CWM for TAF
cowplot::save_plot(here::here("output/figures/Fig5_Pred_S_CWM.png"),
                   cowplot::plot_grid(p.ES, p.FWL, p.AR, ncol = 2),
                   base_height = 12, base_width = 10)


# Fig S5 - Standard deviation os posterior mean species richness
p.sdS <- ggplot(data = nt) +
  geom_sf(fill = "white", size = .3) +
  coord_sf(xlim = c(-58, -34), ylim = c(-34, -2.5), expand = T, 
           crs = st_crs(4326)) +
  geom_tile(data = df.ric, aes(x = Longitude, y = Latitude, fill = sdS)) + 
  labs(title = "Standand Deviation") + scale_fill_viridis_c(option = "B", name = NULL) 
p.sdS

range(df.ric$ES)
range(df.ric$sdS)

range(mapData$FWL)
range(mapData$AR)

cowplot::save_plot(here::here("output/figures/FigS6_sdS.png"),p.sdS,
                   base_height = 6, base_width = 8)
