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


#' ## Get the mining plan

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

PUs["locked_out_mining"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN # create locked_out column of exploration and reserved mining areas (logical). If =<1 is TRUE, then locked_out is TRUE # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

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

PUxrCFC <- exactextractr::exact_extract(rCFC, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% #seems like NaN become 0s on their own?? after you click to check sfCFC and sfPMN?
  mutate(CFCValue = PUxrCFC)

PUxrPMN <- exactextractr::exact_extract(rPMN, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(PMNValue = PUxrPMN) #seems like NaN become 0s on their own?? after you click to check sfCFC and sfPMN?

# get total mineral resource value per PU

PUs$MiningCost <- PUs$CFCValue + PUs$PMNValue

#' ## Get the fishing plan

#FishingCost <- terra::rast("Data/Cost/Cost_Raster_Sum.grd") %>% 
 # terra::as.polygons(trunc = FALSE, dissolve = FALSE, na.rm=FALSE) %>% # Convert to polygon data
  #st_as_sf() %>% # Convert to sf
  #st_transform(cCRS) %>% # transform to robinson
  #st_interpolate_aw(PUs, extensive = TRUE) %>% ## intersect with PUs # set to = TRUE to have sum per PU rather than mean?
  #rename(FishingCost = layer)

#PUs <- PUs %>% 
 # mutate(Fishing = FishingCost$FishingCost) # be clear about units

#> min(PUs$Fishing) #with extensive = TRUE, and PU as square (but still 1000km2 area), which is odd because it's supposed to be the sum of all values rather than the average, but the values are lower for the sum than for the average 
#[1] 0.09251705
#> max(PUs$Fishing)
#[1] 1888.157
#> mean(PUs$Fishing)
#[1] 10.15606

# the method below is the same as for shipping and mining - might be better to have consitent method to compare? also, max and mean value make more sense like this, I think, and we are sure to have the sum of values per PU

FishingRaster <- terra::rast("Data/Cost/Cost_Raster_Sum.grd")
FishingRobin <- terra::project(FishingRaster, "ESRI:54009") #transform to mollweide proj
PUxFishing <- exactextractr::exact_extract(FishingRobin, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(Fishing = PUxFishing) # add column to PU data frame

max(PUs$Fishing)
#[1] 9646.386 USD/PU/yr
mean(PUs$Fishing)
#[1] 57.71586 USD/PU/yr

#(ggCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

#' ## Get the shipping plan

Shipping <- terra::rast("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Shipping/shipping.tif")
ShippingRobin <- terra::project(Shipping, "ESRI:54009") #transform to mollweide proj
PUxShipping <- exactextractr::exact_extract(ShippingRobin, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(ShippingIntensity = PUxShipping) # add column to PU data frame

max(PUs$ShippingIntensity)
#[1] 373235.6
mean(PUs$ShippingIntensity)
#[1] 6323.317


#' ## Prioritizr problems

#' FishPlan problem (no longhurst, no shipping)

p_FishPlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
  add_binary_decisions() %>%
  #add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_Sol <- solve(p_FishPlan)

#' MinePlan problem (no longhurst, no shipping)

PUs["locked_out_mining"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN # create locked_out column of exploration and reserved mining areas (logical). If =<1 is TRUE, then locked_out is TRUE # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

p_MinePlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_mining") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_Sol <- solve(p_MinePlan, force = TRUE) # need force = TRUE because planning units with very high (> 1e+6) cost values


#' ShipPlan problems (ShipPlan = shipping intensity as cost; ShipPlan2 = shipping intensity locked out, cost = 0)


p_ShipPlan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
  add_binary_decisions() %>%
  #add_locked_out_constraints(locked_out = "locked_out_fishing") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_Sol <- solve(p_ShipPlan, force = TRUE)

#' ## Ferrier scores

df <- p_FishPlan_Sol%>%
  dplyr::select(geometry, solution_1)
FerrierFishA <- eval_ferrier_importance(p_FishPlan, df)

df2 <- p_MinePlan_Sol%>%
  dplyr::select(geometry, solution_1)
FerrierMineA <- eval_ferrier_importance(p_MinePlan, df2)

df3 <- p_ShipPlan_Sol%>%
  dplyr::select(geometry, solution_1)
FerrierShipA <- eval_ferrier_importance(p_ShipPlan, df3)

#' ## Replacement cost scores
# https://www.sciencedirect.com/science/article/pii/S0006320706001820?casa_token=nbc0AGdxJ8IAAAAA:4XoXBy88sLnIzJMfRAq6gz5FkI5UMCNA7oRcxgXAcew0_dpDCJkbdzJI7eIq8WjL-sKaRsFuN-Qi
# A replacement-cost value of zero tells us that there exists an alternative solution with the same properties as the current (best) solution has, i.e. same cost, and same biodiversity value (although probably obtained with representation levels that are different from those in the original optimal selection). A replacement cost larger than zero means that any alternative solution including/excluding the focal site will have either a lower biodiversity value or a larger cost than the optimal one.

ReplacementFishA <- eval_replacement_importance(p_FishPlan, df)
ReplacementMineA <- eval_replacement_importance(p_MinePlan, df2, force = TRUE)
ReplacementShipA <- eval_replacement_importance(p_ShipPlan, df3, force = TRUE)


#' ## Plotting


#' Plotting conservation features


IMMA <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IMMAs/iucn-imma.shp")%>% 
  st_transform(cCRS)
IBA <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp")%>% 
  st_transform(cCRS)
NEIO_09<- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp")%>% 
  st_transform(cCRS)
NWIO_14 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_11 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_19 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_22 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_23 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_30 <- st_read  ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_32 <- st_read  ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_35<- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_36 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp")%>% 
  st_transform(cCRS)
SIO_37  <-st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp")%>% 
  st_transform(cCRS)
Seamounts <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Seamounts.shp")%>% 
  st_transform(cCRS)
Plateaus <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Plateaus.shp")%>% 
  st_transform(cCRS)

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = NEIO_09,  aes(color="A"), show.legend = "line")+ 
  geom_sf(data = NWIO_14,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_11,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_19,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_22,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_23,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_30,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_32,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_35,  aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_36, aes(color="A"), show.legend = "line")+
  geom_sf(data = SIO_37,  aes(color="A"), show.legend = "line")+
  geom_sf(data = Seamounts,  aes(color="B"), show.legend = "line")+
  geom_sf(data = Plateaus,  aes(color="C"), show.legend = "line")+
  geom_sf(data = IBA,  aes(color="D"), show.legend = "line") + 
  geom_sf(data = IMMA,  aes(color="E"), show.legend = "line")+ 
  geom_sf(data = VentsActiveSF,  aes(color="F"), show.legend = "line")+
  geom_sf(data = VentsInactiveSF,  aes(color="G"), show.legend = "line")+
  ggtitle("IBA, IMMA, EBSAs, Seamounts, Plateaus, Vents")+ 
  #coord_sf(xlim = c(18, 130), ylim = c(34, -60), expand = FALSE)+
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)+
  scale_color_manual(values = c("A" = "grey", "B" = "black", "C" = "blue", "D" = "green", "E" = "purple", "F" = "red", "G" = "brown"), 
                     labels = c("EBSAs", "Seamounts", "Plateaus", "IBAs", "IMMAs", "Active vents", "Inactive Vents"),
                     name = "Legend") 

#' plotting locked out areas

world <- ne_countries(scale = "medium", returnclass = "sf")%>% 
  st_transform(cCRS)
#IndianO <- c(xmin = 20, xmax = 120, ymin = -60, ymax = 34)
Expl_PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/01_pmn_exploration_areas.shp")%>% 
  st_transform(cCRS)
Expl_PMS <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/02_pms_exploration_areas.shp")%>% 
  st_transform(cCRS)
Res_PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/04_pmn_reserved_areas.shp")%>% 
  st_transform(cCRS)

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = Expl_PMN, aes(color="A"), show.legend = "line") + 
  geom_sf(data = Expl_PMS, aes(color="B"), show.legend = "line") +
  geom_sf(data = Res_PMN, aes(color="C"), show.legend = "line") +
  scale_color_manual(values = c("A" = "red", "B" = "yellow", "C" = "black"),
                     labels = c("Expl_PMN", "Expl_PMS", "Res_PMN"),
                     name = "Legend") +
  ggtitle("Exploration & Reserved areas") + 
  #coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)


#' Plotting fishing cost + log10 fishing cost

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = PUs, aes(color=Log10_Fishing)) + 
  ggtitle("Fishing Cost (Log10)") + 
  #coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$Fishing2)) +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  ggtitle("Fishing Opportunity Cost using 'sum' (USD/PU/yr)") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)
#coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$Fishing)) +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  ggtitle("Fishing Opportunity Cost using 'mean' (USD/PU/yr)") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)
#coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

#' Plotting mining layer

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = CFC, aes(color="A"), show.legend = "line") + 
  geom_sf(data = PMN, aes(color="B"), show.legend = "line") + 
  scale_color_manual(values = c("A" = "grey", "B" = "brown"),
                     labels = c("CFC", "PMN"),
                     name = "Legend") +
  ggtitle("Polymetallic Nodule and Cobalt-rich Ferromanganese Crust areas") + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$MiningCost)) +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  ggtitle("Mining Opportunity Cost using 'sum' (USD/PU)") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

#' Plotting shipping layer

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = PUs, aes(color=ShippingIntensity)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Shipping Intensity (ships/PU/yr)")



#' Plotting FishPlan

p_FishPlan_Sol_Plot <- ggplot() +
  ggtitle("Optimal Conservation Plan for the Fishing Industry") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_FishPlan_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FishPlan_Sol_Plot

#' PLotting MinePlan


p_MinePlan_Sol_Plot <- ggplot() +
  ggtitle("Optimal Conservation Plan for the Deep-sea Mining Industry") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_MinePlan_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_MinePlan_Sol_Plot

#' Plotting ShipPlan

p_ShipPlan_Sol_Plot <- ggplot() +
  ggtitle("Optimal Conservation Plan for the Shipping Industry") +
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = p_ShipPlan_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ShipPlan_Sol_Plot

#' Comparing 3 plans

ctrl <- p_FishPlan_Sol
rmv <- p_MinePlan_Sol
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

## plotting Ferrier scores

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierFishA, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for FishPlan, Approach A")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierMineA, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for MinePlan, Approach A")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = FerrierShipA, aes(color=total)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Ferrier Scores for ShipPlan, Approach A")

## plotting Replacement scores 

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementFishA, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for FishPlan, Approach A")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementMineA, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for MinePlan, Approach A")

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE)+
  geom_sf(data = ReplacementShipA, aes(color=rc)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Replacement Scores for ShipPlan, Approach A")

