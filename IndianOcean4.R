#' ---
#' title: "Create a Spatial Plan"
#' author: "Jason D. Everett, modified by Lea Fourchault"
#' affiliation: "UQ/CSIRO/UNSW/Tropimundo"
#' date: "Last compiled on April 8th, 2022"
#' output: 
#'   html_document
#' ---
#'

#' ## Overview
#' This code has been written to simplify the process for running a _prioritizr_ analysis on a given region. It is still a work in progress so feel free to submit pull requests with new features and code improvements.
#' 
#' The code depends on `sf`, `terra`, `tidyverse`, `rnaturalearth`, `prioritizr`, `stars`, `patchwork`.    
#' 
#' To use this code, you will need to download and expand `MME1DATA-Q1215/SpatialPlanning/Data.zip` to the directory `GitHub/SpatialPlanning/Data/`. Note that the download is only 2GB, but the expanded data is 35 GB in size. If you need help subsetting the data to your region due to memory or HD space constraints, contact Jason.
#' 
#' Script below modified by Lea Fourchault for the Cross-sectoral conservation of the Indian Ocean project
#' 
#' ## Preliminaries 
source("SpatPlan_Extras.R") # Load the extras, including functions and libraries
# source("SpatPlan_Process_AquaMaps.R") # This script reprocesses AquaMaps. WARNING: Lots of time (10s hrs) and memory (10s GB)
# source("SpatPlan_Process_MPAs.R") # Only run if you need to reprocess the MPA data. You will need the WDPA files

if (!file.exists("Data")) {
  stop("The Data folder does not exist at SpatialPlanning/Data. Please download from the RDM and then try again.")
}

reprocess<- TRUE # Do we want to reprocess the PUs/AquaMaps data or use saved versions

fSpatPlan_Get_Polyg <- function(filen, PUs){
  
  polyg_crs <- st_crs(PUs)
  Polyg1 <- st_read(filen) %>% 
    st_transform(polyg_crs) 
  
  Polyg1PUs <- PUs %>% 
    st_centroid() %>% # 50% overlap
    st_within(Polyg1, sparse = FALSE) %>% 
    rowSums() %>%
    as.logical()
  
  # Polyg1PUs <- PUs %>% 
  #   mutate(locked_in = Polyg1PUs) ##assign new column to PUs (call it name of shp)
  # 
  return(Polyg1PUs)
}

##' Define region and boundary

cCRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # Mollweide projection ESRI:54009
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(cCRS) #

Region <- c(xmin = 18, xmax = 120, ymin = -45, ymax = 34) 
Bndry <- fSpatPlan_Get_Boundary(Region, cCRS) # boundary

#' Intersect bndry and high seas

ABNJ <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/World_High_Seas_v1/High_Seas_v1.shp") 
ABNJRobin <- st_transform(ABNJ, cCRS) # high seas in mollweide proj, does not work when using pipe
ABNJRobBuff  <- st_buffer(ABNJRobin, 0) # need to create buffer of 0 to avoid Error in CPL_geos_op2(op, x, y) : Evaluation error: TopologyException: Input geom 0 is invalid: Self-intersection at or near point -> overwrite offending geometries
BndryABNJ <- st_intersection(ABNJRobBuff, Bndry)

#' Create PUs

PU_size <- 1000 # in km2, but check resolution for this project?
Shape <- "square" # shape of PUs
PUs <- fSpatPlan_Get_PlanningUnits(BndryABNJ, world, PU_size, Shape) # modified PUs with ABNJ only

#' Plot PUs

(ggPU <- fSpatPlan_PlotPUs(PUs, world)) 

#' ## Get the conservation features 
#' IBAs

#loading file 
IBA <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp" # no need to st_transform here because done in function when crossing with PUs

#crossing shp & PUs
PUxIBA <- fSpatPlan_Get_Polyg(IBA, PUs)  
PUs <- PUs %>% 
  mutate(IBA = PUxIBA) #add column to PU table with

#' IMMAs

#loading file
IMMA <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IMMAs/iucn-imma.shp"

#crossing shp & PUs
PUxIMMA <- fSpatPlan_Get_Polyg(IMMA, PUs)  
PUs <- PUs %>% 
  mutate(IMMA = PUxIMMA)

#' EBSAs

NEIO_09<- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp"
PUxNEIO_09 <- fSpatPlan_Get_Polyg(NEIO_09, PUs)  
PUs <- PUs %>% 
  mutate(NEIO_09 = PUxNEIO_09) # sumatra upwelling

NWIO_14 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp"
PUxNWIO_14 <- fSpatPlan_Get_Polyg(NWIO_14, PUs)  
PUs <- PUs %>% 
  mutate(NWIO_14 = PUxNWIO_14) # arabian oxygen min zone

SIO_11 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp"
PUxSIO_11 <- fSpatPlan_Get_Polyg(SIO_11, PUs)  
PUs <- PUs %>% 
  mutate(SIO_11 = PUxSIO_11) # agulhas upwelling 

SIO_19 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp"
PUxSIO_19 <- fSpatPlan_Get_Polyg(SIO_19, PUs)  
PUs <- PUs %>% 
  mutate(SIO_19 = PUxSIO_19)

SIO_22 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp"
PUxSIO_22 <- fSpatPlan_Get_Polyg(SIO_22, PUs)  
PUs <- PUs %>% 
  mutate(SIO_22 = PUxSIO_22)

SIO_23 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp"
PUxSIO_23 <- fSpatPlan_Get_Polyg(SIO_23, PUs)  
PUs <- PUs %>% 
  mutate(SIO_23 = PUxSIO_23)

SIO_30 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp"
PUx_SIO30 <- fSpatPlan_Get_Polyg(SIO_30, PUs)  
PUs <- PUs %>% 
  mutate(SIO_30 = PUx_SIO30)

SIO_32 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp"
PUxSIO_32 <- fSpatPlan_Get_Polyg(SIO_32, PUs)  
PUs <- PUs %>% 
  mutate(SIO_32 = PUxSIO_32)

SIO_35<- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp"
PUxSIO_35 <- fSpatPlan_Get_Polyg(SIO_35, PUs)  
PUs <- PUs %>% 
  mutate(SIO_35 = PUxSIO_35)

SIO_36 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp"
PUxSIO_36 <- fSpatPlan_Get_Polyg(SIO_36, PUs)  
PUs <- PUs %>% 
  mutate(SIO_36 = PUxSIO_36)

SIO_37  <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp"
PUx_SIO37 <- fSpatPlan_Get_Polyg(SIO_37, PUs)  
PUs <- PUs %>% 
  mutate(SIO_37 = PUx_SIO37)

#' Deep sea features

Seamounts <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Seamounts.shp"
PUxSeamounts <- fSpatPlan_Get_Polyg(Seamounts, PUs)  
PUs <- PUs %>% 
  mutate(Seamounts = PUxSeamounts)

Plateaus <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Plateaus.shp"
PUxPlateaus <- fSpatPlan_Get_Polyg(Plateaus, PUs)  
PUs <- PUs %>% 
  mutate(Plateaus = PUxPlateaus)

Vents <- read.csv("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/vents.csv")
VentsActive <- filter(Vents, Activity == "active, confirmed" | Activity == "active, inferred") 
VentsActiveSF <- st_as_sf(VentsActive, coords = c(x = "Longitude", y = "Latitude"), crs = 4326) %>% 
  st_transform(cCRS)

PUxVentsActive <- st_contains(PUs,VentsActiveSF, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(Active_Vents = PUxVentsActive)

VentsInactive <- filter(Vents, Activity == "inactive")
VentsInactiveSF <- st_as_sf(VentsInactive, coords = c(x = "Longitude", y = "Latitude"), crs = 4326) %>% 
  st_transform(cCRS)

PUxVentsInactive <- st_contains(PUs,VentsInactiveSF, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(Inactive_Vents = PUxVentsInactive)

#' ## Lock out shipping routes

# create shipping intensity column in PUs

Shipping <- terra::rast("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Shipping/shipping.tif")
ShippingRobin <- terra::project(Shipping, "ESRI:54009") #transform to mollweide proj
PUxShipping <- exactextractr::exact_extract(ShippingRobin, PUs, "sum") # attribute sum of values of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(ShippingIntensity = PUxShipping) # add column to PU data frame

# explore data

max(PUs$ShippingIntensity)
#373235.6
mean(PUs$ShippingIntensity)
#6323.317
median(PUs$ShippingIntensity)
#0
sd(PUs$ShippingIntensity)
#20811.08

#' ## Get the fishing plan with shipping routes locked out

# get fishing cost using extact extractr, as for mining - might be better to have consistent method to compare? also, max and mean value make more sense like this, I think, and we are sure to have the sum of values per PU

FishingRaster <- terra::rast("Data/Cost/Cost_Raster_Sum.grd")
FishingRobin <- terra::project(FishingRaster, "ESRI:54009") #transform to mollweide proj
PUxFishing <- exactextractr::exact_extract(FishingRobin, PUs, "sum") # attribute mean value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(Fishing2 = PUxFishing) # add column to PU data frame

max(PUs$Fishing2)
#[1] 9646.386 USD/PU/yr
mean(PUs$Fishing2)
#[1] 57.71586 USD/PU/yr

#(ggCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

#' ## Get the mining plan with mining exploration and reserved areas & shipping routes locked out 

## lock out exploration areas and reserved areas

ExplPMN <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/01_pmn_exploration_areas.shp"
PUxExplPMN <- fSpatPlan_Get_Polyg(ExplPMN, PUs)  
PUs <- PUs %>% 
  mutate(ExplPMN = PUxExplPMN)

ExplPMS <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/02_pms_exploration_areas.shp"
PUxExplPMS <- fSpatPlan_Get_Polyg(ExplPMS, PUs)  
PUs <- PUs %>% 
  mutate(ExplPMS = PUxExplPMS)

ResPMN <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/04_pmn_reserved_areas.shp"
PUxResPMN <- fSpatPlan_Get_Polyg(ResPMN, PUs)  
PUs <- PUs %>% 
  mutate(ResPMN = PUxResPMN)

# create locked_out column of exploration and reserved mining areas (logical) # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

PUs["locked_out_mining"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN

## get the mineral resources value

# load sf objects for mineral types

sfPMN <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")%>%
  st_transform(cCRS) #mollweide proj

sfCFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")%>%
  st_transform(cCRS) #mollweide proj

# attribute average cost to mineral type (based on area covered, interaction calculation explained in Master_Lea_v2)

sfCFC$value <- 14307812.057 #sqkm value

sfPMN$value <- 94569447.19  #sqkm value

# create template raster based on PU to matching resolution, and intersect with sf objects (also because high resolution eg res=100 generates error)

rPUs <- raster::raster(PUs, res=1000) #Formal class RasterLayer, resolution 1 sqkm
rCFC <- fasterize::fasterize(sfCFC, rPUs, field = "value") 
rPMN <- fasterize::fasterize(sfPMN, rPUs, field = "value") 

# extract sum of all pixel values contained in each PU to get sum of USD value of all minerals

PUxrCFC <- exactextractr::exact_extract(rCFC, PUs, "sum") # attribute mean value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(CFCValue = PUxrCFC)

PUxrPMN <- exactextractr::exact_extract(rPMN, PUs, "sum") # attribute mean value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(PMNValue = PUxrPMN)

# get total mineral resource value per PU

PUs$MiningCost <- PUs$CFCValue + PUs$PMNValue

#' ## Prioritizr problems

#' Fishing, shipping as linear constraint

p_FishPlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_Sol <- solve(p_FishPlan)

sum(PUs$ShippingIntensity) #[1] 210465273

p_FishShipPlan0.1 <- p_FishPlan%>%
  add_linear_constraints(threshold = 210465273*0.05, # plan should have less than 10% of total ships in selected PUs
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_FishShipPlan0.2 <- p_FishPlan%>%
  add_linear_constraints(threshold = 210465273*0.2, #20%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_FishShipPlan0.3 <- p_FishPlan%>%
  add_linear_constraints(threshold = 210465273*0.3, #30%
                         sense = "<=",
                         data = PUs$ShippingIntensity)
p_FishShipPlan0.4 <- p_FishPlan%>%
  add_linear_constraints(threshold = 210465273*0.3, #40%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_FishShipPlan0.1_Sol <- solve(p_FishShipPlan0.1, force = TRUE)
p_FishShipPlan0.2_Sol <- solve(p_FishShipPlan0.2, force = TRUE)
p_FishShipPlan0.3_Sol <- solve(p_FishShipPlan0.3, force = TRUE)
p_FishShipPlan0.4_Sol <- solve(p_FishShipPlan0.4, force = TRUE)

p_FishShipPlan0.1_Sol_Plot <- ggplot() +
  ggtitle("Fishing, Shipping Const = 10%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_FishShipPlan0.1_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FishShipPlan0.1_Sol_Plot

p_FishShipPlan0.2_Sol_Plot <- ggplot() +
  ggtitle("Fishing, Shipping = 20%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_FishShipPlan0.2_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FishShipPlan0.2_Sol_Plot

p_FishShipPlan0.3_Sol_Plot <- ggplot() +
  ggtitle("Fishing, Shipping = 30%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_FishShipPlan0.3_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FishShipPlan0.3_Sol_Plot

p_FishShipPlan0.4_Sol_Plot <- ggplot() +
  ggtitle("Fishing, Shipping = 40%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_FishShipPlan0.4_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FishShipPlan0.4_Sol_Plot


try <- cowplot::plot_grid(p_FishShipPlan0.1_Sol_Plot, p_FishShipPlan0.2_Sol_Plot, p_FishShipPlan0.3_Sol_Plot, p_FishShipPlan0.4_Sol_Plot, 
                          labels = c("A", "B", "C", "D"),
                          ncol = 2, nrow = 2)

#' MinePlan problem, shipping  as constraint

p_MinePlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_mining") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_Sol <- solve(p_MinePlan, force = TRUE) # need force = TRUE because planning units with very high (> 1e+6) cost values

p_MineShipPlan0.1 <- p_MinePlan%>%
  add_linear_constraints(threshold = 210465273*0.05, # plan should have less than 10% of total ships in selected PUs
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_MineShipPlan0.2 <- p_MinePlan%>%
  add_linear_constraints(threshold = 210465273*0.05, #15%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_MineShipPlan0.3 <- p_MinePlan%>%
  add_linear_constraints(threshold = 210465273*0.10, #20%
                         sense = "<=",
                         data = PUs$ShippingIntensity)
p_MineShipPlan0.4 <- p_MinePlan%>%
  add_linear_constraints(threshold = 210465273*0.20, #25%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_MineShipPlan0.1_Sol <- solve(p_MineShipPlan0.1, force = TRUE)
p_MineShipPlan0.2_Sol <- solve(p_MineShipPlan0.2, force = TRUE)
p_MineShipPlan0.3_Sol <- solve(p_MineShipPlan0.3, force = TRUE)
p_MineShipPlan0.4_Sol <- solve(p_MineShipPlan0.4, force = TRUE)

p_MineShipPlan0.1_Sol_Plot <- ggplot() +
  ggtitle("Mining, Shipping = 5%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_MineShipPlan0.1_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_MineShipPlan0.1_Sol_Plot

p_MineShipPlan0.2_Sol_Plot <- ggplot() +
  ggtitle("Mining, Shipping = 10%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_MineShipPlan0.2_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_MineShipPlan0.2_Sol_Plot

p_MineShipPlan0.3_Sol_Plot <- ggplot() +
  ggtitle("Mining, Shipping = 15%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_MineShipPlan0.3_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_MineShipPlan0.3_Sol_Plot

p_MineShipPlan0.4_Sol_Plot <- ggplot() +
  ggtitle("Mining, Shipping = 20%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_MineShipPlan0.4_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_MineShipPlan0.4_Sol_Plot

p_MineShipPlan0.1_Sol_Plot + p_MineShipPlan0.2_Sol_Plot



#' ## Ship problem, zero cost and ship as constraint

# create NoCost column  in PUs

PUs$NoCost <- 0

# create problems

p_ShipPlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "NoCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_Sol <- solve(p_ShipPlan, force = TRUE) # need force = TRUE because planning units with very high (> 1e+6) cost values

p_ShipPlan0.1 <- p_ShipPlan%>% #divide shipping by 1000 to scale it down
  add_linear_constraints(threshold = 210465.273*0.1, # plan should have less than 10% of total ships in selected PUs
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_ShipPlan0.2 <- p_ShipPlan%>%
  add_linear_constraints(threshold = 210465.273*0.2, #20%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_ShipPlan0.3 <- p_ShipPlan%>%
  add_linear_constraints(threshold = 210465.273*0.3, #30%
                         sense = "<=",
                         data = PUs$ShippingIntensity)
p_ShipPlan0.4 <- p_ShipPlan%>%
  add_linear_constraints(threshold = 210465.273*0.3, #40%
                         sense = "<=",
                         data = PUs$ShippingIntensity)

p_ShipPlan0.1_Sol <- solve(p_ShipPlan0.1, force = TRUE)
p_ShipPlan0.2_Sol <- solve(p_ShipPlan0.2, force = TRUE)
p_ShipPlan0.3_Sol <- solve(p_ShipPlan0.3, force = TRUE)
p_ShipPlan0.4_Sol <- solve(p_ShipPlan0.4, force = TRUE)

p_ShipPlan0.1_Sol_Plot <- ggplot() +
  ggtitle("NoCost, Shipping = 10%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_ShipPlan0.1_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ShipPlan0.1_Sol_Plot

p_ShipPlan0.2_Sol_Plot <- ggplot() +
  ggtitle("NoCost,Shipping = 20%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_ShipPlan0.2_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ShipPlan0.2_Sol_Plot

p_ShipPlan0.3_Sol_Plot <- ggplot() +
  ggtitle("NoCost, Shipping = 30%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_ShipPlan0.3_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ShipPlan0.3_Sol_Plot

p_ShipPlan0.4_Sol_Plot <- ggplot() +
  ggtitle("NoCost, Shipping = 40%") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_ShipPlan0.4_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ShipPlan0.4_Sol_Plot

#' Comparing 3 plans

ctrl <- p_FishShipPlan_Sol
rmv <- p_MineShipPlan_Sol
aaa <- p_ShipPlan_Sol

compare <- fCompareSolutions(ctrl, rmv, aaa)

fCompareSolutions <- function(ctrl, rmv, aaa){
  
  #Select only the ID column and the solution column; change the name of the columns in 'cellID' and 'ctrl_sol'  
  
  ctrl2 <- ctrl %>% #fish
    dplyr::select(geometry, solution_1) %>% #ID does not exist, select geometry column instead
    rename(geometry = geometry, ctrl_sol = solution_1)
  
  rmv2 <- rmv %>% #mine
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, rmv_sol = solution_1)
  
  rmv3 <- aaa %>% #ship
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, rmv_sol3 = solution_1)
  
  #Create an object with the three solutions
  
  soln <- ctrl2 # fish
  soln$rmv_sol <- rmv2$rmv_sol #mine
  soln$rmv_sol3 <- rmv3$rmv_sol3 #ship
  
  #Combine the result and call them differently if a PUs in present in both the solution, or only in one of them
  
  soln <- mutate(soln, Combined = as.numeric(ctrl_sol + rmv_sol + rmv_sol3)) %>%
    mutate(Compare = case_when(Combined == 3 ~ "Optimal for All Industries",
                               ctrl_sol == 1 & rmv_sol == 1 & rmv_sol3 == 0 ~ "Optimal for Fishing and Mining",                            
                               ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 1 ~ "Optimal for Mining and Shipping",
                               ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 1 ~ "Optimal for Fishing and Shipping",
                               ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 0 ~ "Optimal for Fishing",
                               ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 0 ~ "Optimal for Mining",
                               ctrl_sol == 0 & rmv_sol == 0 & rmv_sol3 == 1 ~ "Optimal for Shipping",
                               Combined == 0 ~ "Not Selected"),
           Compare = factor(Compare, levels = c("Optimal for All Industries",
                                                "Optimal for Fishing and Mining",                                    
                                                "Optimal for Mining and Shipping",
                                                "Optimal for Fishing and Shipping",
                                                "Optimal for Fishing",
                                                "Optimal for Mining",
                                                "Optimal for Shipping", 
                                                "Not Selected"))) %>%
    filter(!is.na(Compare))
  
  label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%\nCost: ",round((sum(rmv$solution_1*rmv$cost) - sum(ctrl$solution_1*ctrl$cost))/ sum(ctrl$solution_1*ctrl$cost)*100, 1), "%\n") #Error in label_loc$label <- paste0("Area: ", round((sum(rmv$solution_1) -  : 
  # object 'label_loc' not found
  label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%")
  
  #Plot the results
  
  gg <- ggplot() +
    geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE) +
    geom_sf(data = PUs, fill = "lightsteelblue2", color = "grey64", size = 0.05, show.legend = FALSE) +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
    geom_sf(data = soln, aes(fill = Compare), colour = "black", size = 0.0001) +
    theme_bw() +
    theme(axis.title = element_blank(), text = element_text(size = 20), plot.title = element_text(size = 12)) +
    scale_fill_manual(values = c("Optimal for All Industries" = "red",
                                 "Optimal for Fishing and Mining" = "magenta",
                                 "Optimal for Mining and Shipping" = "blue",
                                 "Optimal for Fishing and Shipping" = "darkorchid",
                                 "Optimal for Fishing" = "limegreen",
                                 "Optimal for Mining" = "orange",
                                 "Optimal for Shipping" = "lightsalmon1",
                                 "Not Selected" = "white"))
  
  return(gg)
}

gg

#' Comparing 3 plans

ctrl <- p_ShipPlan0.1_Sol
rmv <- p_ShipPlan0.2_Sol
aaa <- p_ShipPlan0.3_Sol

compare <- fCompareSolutions(ctrl, rmv, aaa)

fCompareSolutions <- function(ctrl, rmv, aaa){
  
  #Select only the ID column and the solution column; change the name of the columns in 'cellID' and 'ctrl_sol'  
  
  ctrl2 <- ctrl %>% #fish
    dplyr::select(geometry, solution_1) %>% #ID does not exist, select geometry column instead
    rename(geometry = geometry, ctrl_sol = solution_1)
  
  rmv2 <- rmv %>% #mine
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, rmv_sol = solution_1)
  
  rmv3 <- aaa %>% #ship
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, rmv_sol3 = solution_1)
  
  #Create an object with the three solutions
  
  soln <- ctrl2 # 0.1
  soln$rmv_sol <- rmv2$rmv_sol # 0.2
  soln$rmv_sol3 <- rmv3$rmv_sol3 # 0.3
  
  #Combine the result and call them differently if a PUs in present in both the solution, or only in one of them
  
  soln <- mutate(soln, Combined = as.numeric(ctrl_sol + rmv_sol + rmv_sol3)) %>%
    mutate(Compare = case_when(Combined == 3 ~ "Optimal for All Constraints",
                               ctrl_sol == 1 & rmv_sol == 1 & rmv_sol3 == 0 ~ "0.1 and 0.2",                            
                               ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 1 ~ "0.2 and 0.3",
                               ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 1 ~ "0.1 and 0.3",
                               ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 0 ~ "0.1",
                               ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 0 ~ "0.2",
                               ctrl_sol == 0 & rmv_sol == 0 & rmv_sol3 == 1 ~ "0.3",
                               Combined == 0 ~ "Not Selected"),
           Compare = factor(Compare, levels = c("Optimal for All Constraints",
                                                "0.1 and 0.2",                                    
                                                "0.2 and 0.3",
                                                "0.1 and 0.3",
                                                "0.1",
                                                "0.2",
                                                "0.3", 
                                                "Not Selected"))) %>%
    filter(!is.na(Compare))
  
  label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%\nCost: ",round((sum(rmv$solution_1*rmv$cost) - sum(ctrl$solution_1*ctrl$cost))/ sum(ctrl$solution_1*ctrl$cost)*100, 1), "%\n") #Error in label_loc$label <- paste0("Area: ", round((sum(rmv$solution_1) -  : 
  # object 'label_loc' not found
  label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%")
  
  #Plot the results
  
  gg <- ggplot() +
    geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE) +
    geom_sf(data = PUs, fill = "lightsteelblue2", color = "grey64", size = 0.05, show.legend = FALSE) +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
    geom_sf(data = soln, aes(fill = Compare), colour = "black", size = 0.0001) +
    theme_bw() +
    theme(axis.title = element_blank(), text = element_text(size = 20), plot.title = element_text(size = 12)) +
    scale_fill_manual(values = c("Optimal for All Constraints" = "red",
                                 "0.1 and 0.2" = "magenta",
                                 "0.2 and 0.3" = "blue",
                                 "0.1 and 0.3" = "darkorchid",
                                 "0.1" = "limegreen",
                                 "0.2" = "orange",
                                 "0.3" = "lightsalmon1",
                                 "Not Selected" = "white"))
  
  return(gg)
}

gg


#' ## Ferrier scores

df7 <- p_FishShipPlan0.1_Sol%>% #with 5% threshold for shipping
  dplyr::select(geometry, solution_1)
FerrierFishD <- eval_ferrier_importance(p_FishShipPlan0.1, df7)

df8 <- p_MineShipPlan0.1_Sol%>%
  dplyr::select(geometry, solution_1)
FerrierMineD <- eval_ferrier_importance(p_MineShipPlan0.1, df8)

df9 <- p_ShipPlan_Sol%>%
  dplyr::select(geometry, solution_1)
FerrierShipD <- eval_ferrier_importance(p_ShipPlan, df9)

## plotting Ferrier scores

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierFish4, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for FishPlan, Approach D")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierMine, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for MinePlan, Approach D")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierShip, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for ShipPlan, Approach D")

#' ## Replacement cost scores
# https://www.sciencedirect.com/science/article/pii/S0006320706001820?casa_token=nbc0AGdxJ8IAAAAA:4XoXBy88sLnIzJMfRAq6gz5FkI5UMCNA7oRcxgXAcew0_dpDCJkbdzJI7eIq8WjL-sKaRsFuN-Qi
# A replacement-cost value of zero tells us that there exists an alternative solution with the same properties as the current (best) solution has, i.e. same cost, and same biodiversity value (although probably obtained with representation levels that are different from those in the original optimal selection). A replacement cost larger than zero means that any alternative solution including/excluding the focal site will have either a lower biodiversity value or a larger cost than the optimal one.

ReplacementFish <- eval_replacement_importance(p_FishShipPlan0.1, df7) #redo
ReplacementMine <- eval_replacement_importance(p_MineShipPlan0.1, df8, force = TRUE) #with shipping threshold 5%
ReplacementShip <- eval_replacement_importance(p_ShipPlan, df3, force = TRUE)


## plotting Replacement scores 

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementFish, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for FishShipPlan, Approach D")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementMine, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for MinePlan, Approach A")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementShip, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for ShipPlan, Approach A")


