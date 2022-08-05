#' ---
#' title: "Cross-sectoral spatial planning"
#' author: "Lea Fourchault"
#' affiliation: "UQ/CSIRO/UNSW/Tropimundo"
#' date: "Last compiled on July 12th, 2022"
#' output: 
#'   html_document
#' ---
#'

#' ## Overview
#' 
#' This code has been written to simplify the process for running a _prioritizr_ analysis on a given region. It is still a work in progress so feel free to submit pull requests with new features and code improvements.
#' 
#' The code depends on `sf`, `terra`, `tidyverse`, `rnaturalearth`, `prioritizr`, `stars`, `patchwork`, `exactextractr`.    
#' 
#' To use this code, you will need to download and expand `MME1DATA-Q1215/SpatialPlanning/Data.zip` to the directory `GitHub/SpatialPlanning/Data/`. Note that this version does not use Aquamaps.
#' 
#' Contact: Lea Fourchault or Jason D. Everett or Camille K. V. Buenafe from the MME at UQ.
#' 
#' ## Preliminaries 
#' 
source("SpatPlan_Extras.R") # Load the extras, including functions and libraries

if (!file.exists("Data")) {
  stop("The Data folder does not exist at SpatialPlanning/Data. Please download from the RDM and then try again.")
}


##' Define region and boundary

cCRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # Mollweide projection ESRI:54009 to ensure equal area of PUs across latitude

world <- ne_countries(scale = "medium", returnclass = "sf")%>%
  st_transform(cCRS) # apply Mollweide projection to background worldmap

saveRDS(world, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "world.rds"))

Region <- c(xmin = 18, xmax = 120, ymin = -45, ymax = 34) # Indian Ocean, adjusted to include ecological data of interest and exclude low-quality fishing cost data
Bndry <- fSpatPlan_Get_Boundary(Region, cCRS) # boundary of plannin region

#' Intersect bndry and high seas (ABNJ)

ABNJ <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/World_High_Seas_v1/High_Seas_v1.shp") 
ABNJMoll <- st_transform(ABNJ, cCRS) # Apply Mollweide projection to ABNJ. In subsequent lines, does not work when using %>%
ABNJMollBuff  <- st_buffer(ABNJMoll, 0) # need to create buffer of 0 to avoid Error in CPL_geos_op2(op, x, y) : Evaluation error: TopologyException: Input geom 0 is invalid: Self-intersection at or near point -> overwrite offending geometries
BndryABNJ <- st_intersection(ABNJMollBuff, Bndry)

#' Create PUs

reprocess<- TRUE # We want to reprocess the PUs, not use a saved version

PU_size <- 1000 # in km2 at the surface
Shape <- "hexagon" # shape of PUs
PUs <- fSpatPlan_Get_PlanningUnits(BndryABNJ, world, PU_size, Shape) # modified PUs with ABNJ only

saveRDS(PUs, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "PUs.rds"))

#' Plot PUs

(ggPU <- fSpatPlan_PlotPUs(PUs, world)) 

#' ## Get the conservation features 
#' 
#' 
fSpatPlan_Get_Polyg <- function(filen, PUs){ # to overlay PUs & GIS vector files (such as filen.shp)

polyg_crs <- st_crs(PUs) 
Polyg1 <- st_read(filen) %>% 
  st_transform(polyg_crs) # apply projection of PUs to filen.shp

Polyg1PUs <- PUs %>% 
  st_centroid() %>% # 50% overlap
  st_within(Polyg1, sparse = FALSE) %>% 
  rowSums() %>%
  as.logical() # TRUE if >50% feature in 1 PU

# Polyg1PUs <- PUs %>% 
#   mutate(locked_in = Polyg1PUs) ##assign new column to PUs (call it name of shp)
# 
return(Polyg1PUs)
}

#' Important Bird Areas

#loading file 
IBA <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp" # no need to st_transform here because done in function when crossing with PUs

#crossing shp & PUs
PUxIBA <- fSpatPlan_Get_Polyg(IBA, PUs)  
PUs <- PUs %>% 
  mutate(IBA = PUxIBA) #add column to PU table with IBAs

#' Important Marine Mammal Areas

#loading file
IMMA <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IMMAs/iucn-imma.shp"

#crossing shp & PUs
PUxIMMA <- fSpatPlan_Get_Polyg(IMMA, PUs)  
PUs <- PUs %>% 
  mutate(IMMA = PUxIMMA)

#' Ecologically & Biologically Significant Areas

NEIO_09<- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp"
PUxNEIO_09 <- fSpatPlan_Get_Polyg(NEIO_09, PUs)  
PUs <- PUs %>% 
  mutate(NEIO_09 = PUxNEIO_09) # Sumatra Upwelling Zone

NWIO_14 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp"
PUxNWIO_14 <- fSpatPlan_Get_Polyg(NWIO_14, PUs)  
PUs <- PUs %>% 
  mutate(NWIO_14 = PUxNWIO_14) # Arabian Oxygen Min Zone

SIO_11 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp"
PUxSIO_11 <- fSpatPlan_Get_Polyg(SIO_11, PUs)  
PUs <- PUs %>% 
  mutate(SIO_11 = PUxSIO_11) # Agulhas Front

SIO_19 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp"
PUxSIO_19 <- fSpatPlan_Get_Polyg(SIO_19, PUs)  
PUs <- PUs %>% 
  mutate(SIO_19 = PUxSIO_19) # Mozambique Channel

SIO_22 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp"
PUxSIO_22 <- fSpatPlan_Get_Polyg(SIO_22, PUs)  
PUs <- PUs %>% 
  mutate(SIO_22 = PUxSIO_22) # Walters Shoals

SIO_23 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp"
PUxSIO_23 <- fSpatPlan_Get_Polyg(SIO_23, PUs)  
PUs <- PUs %>% 
  mutate(SIO_23 = PUxSIO_23) # Coral Seamount And Fracture Zone

SIO_30 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp"
PUx_SIO30 <- fSpatPlan_Get_Polyg(SIO_30, PUs)  
PUs <- PUs %>% 
  mutate(SIO_30 = PUx_SIO30) # Atlantis Seamount (may not be included based on Bndry)

SIO_32 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp"
PUxSIO_32 <- fSpatPlan_Get_Polyg(SIO_32, PUs)  
PUs <- PUs %>% 
  mutate(SIO_32 = PUxSIO_32) # Saya de Malha Bank

SIO_35<- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp"
PUxSIO_35 <- fSpatPlan_Get_Polyg(SIO_35, PUs)  
PUs <- PUs %>% 
  mutate(SIO_35 = PUxSIO_35) # Rusky Knoll

SIO_36 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp"
PUxSIO_36 <- fSpatPlan_Get_Polyg(SIO_36, PUs)  
PUs <- PUs %>% 
  mutate(SIO_36 = PUxSIO_36) # Fool's Flat

SIO_37  <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp"
PUx_SIO37 <- fSpatPlan_Get_Polyg(SIO_37, PUs)  
PUs <- PUs %>% 
  mutate(SIO_37 = PUx_SIO37) # East Broken Ridge

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
VentsActiveSF <- st_as_sf(VentsActive, coords = c(x = "Longitude", y = "Latitude"), crs = 4326) %>% # EPSG:4326 = WGS84
  st_transform(cCRS) # apply Mollweide projection

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

saveRDS(PUxVentsInactive, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "PUxVentsInactive.rds"))
saveRDS(PUxVentsActive, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "PUxVentsActive.rds"))

saveRDS(PUs, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "PUs_cons_feats.rds"))

# Check features size

length(PUs$IMMA[PUs$IMMA==TRUE]) # 1107
length(PUs$IBA[PUs$IBA==TRUE]) # 538
length(PUs$NWIO_14[PUs$NWIO_14==TRUE]) # 1352
length(PUs$NEIO_09[PUs$NEIO_09==TRUE]) # 72 sumatra upwelling
length(PUs$SIO_11[PUs$SIO_11==TRUE]) # 3718 agulhas upwelling
length(PUs$SIO_19[PUs$SIO_19==TRUE]) # 142 mozambique channel
length(PUs$SIO_22[PUs$SIO_22==TRUE]) # 118 walter's shoals
length(PUs$SIO_23[PUs$SIO_23==TRUE]) # 6 fracture zone
length(PUs$SIO_30[PUs$SIO_30==TRUE]) # 0 atlantis seamount
length(PUs$SIO_32[PUs$SIO_32==TRUE]) # 66 saya de malha bank
length(PUs$SIO_35[PUs$SIO_35==TRUE]) # 1 rusky knoll
length(PUs$SIO_36[PUs$SIO_36==TRUE]) # 1 fool's flat
length(PUs$SIO_37[PUs$SIO_37==TRUE]) # 6 east broken ridge
length(PUs$Active_Vents[PUs$Active_Vents==TRUE]) # 50
length(PUs$Inactive_Vents[PUs$Inactive_Vents==TRUE]) # 6
length(PUs$Seamounts[PUs$Seamounts==TRUE]) # 451
length(PUs$Plateaus[PUs$Plateaus==TRUE]) # 1923

#' ## Get the locked-out mining exploration areas and reserved areas

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

## lock out exploration areas and reserved areas

PUs["locked_out_areas"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN # create locked_out column of exploration and reserved mining areas (logical). If =<1 is TRUE, then locked_out is TRUE # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

##' Get the mining-specific cost layer (mineral resources value)

#' Get CFC extent within planning region
# script to get area modified from https://rpubs.com/rural_gis/255550

#load CFC areas polygon shapefile 
CFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")
CFCMoll <- st_transform(CFC, cCRS) # Mollweide proj

#run the intersect function, converting the output to a tibble in the process
#int <- as_tibble(st_intersection(CFCRobin, BndryABNJ))
int <- st_intersection(CFCMoll, BndryABNJ)

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int$CFCarea <- st_area(int)

#plot the layers to visually check result of intersect
plot (CFCMoll$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaCFCtotal <- int %>%
  summarise(areaCFCtot = sum(CFCarea)) #  2.620621e+12 [m^2] with Mollweide but = 2.256596e+12 [m^2] with Robinson

#change data type of areaArable field to numeric (to remove m^2 suffix)
AreaCFCtotal$areaCFCtot <- as.numeric(AreaCFCtotal$areaCFCtot)

## check: load county areas polygon shapefile 

CFC <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp"
PUxCFC <- fSpatPlan_Get_Polyg(CFC, PUs)  
PUs <- PUs %>% 
  mutate(CFC = PUxCFC)

#'area of CFC in PUs < 2607*1000 sqkm but RECHECK
sum(PUxCFC) #'[1] 2607 RECHECK but seems in line with result above

#' PMN value

#load PMN areas polygon shapefile 

PMN <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")
PMNMoll <- st_transform(PMN, cCRS) # Mollweide proj

#run the intersect function, converting the output to a tibble in the process
#int <- as_tibble(st_intersection(CFCRobin, BndryABNJ))
int2 <- st_intersection(PMNMoll, BndryABNJ)

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int2$PMNarea <- st_area(int2)

#plot the layers to visually check result of intersect
plot (PMNMoll$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int2, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaPMNtotal <- int2 %>%
  summarise(areaPMNtot = sum(PMNarea)) # 2.423256e+12 [m^2] with Mollweide but = 1.933646e+12 [m^2] with Robinson

#change data type of areaArable field to numeric (to remove m^2 suffix)
AreaPMNtotal$areaPMNtot <- as.numeric(AreaPMNtotal$areaPMNtot)

## check: load county areas polygon shapefile 

PMN <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp"
PUxPMN <- fSpatPlan_Get_Polyg(PMN, PUs)  
PUs <- PUs %>% 
  mutate(PMN = PUxPMN)

#'area of PMN in PUs < 2422*1000 sqkm but RECHECK
sum(PUxPMN) #'[1] 2422 RECHECK but seems in line with result above

#' Possible valuation of PMN for Indian Ocean = 320 USD/t (CRU, 2019) * 0.0056 t/sqm (Sharma et al., 2011) * 1.933646e+12 sqm (total PMN area in Indian Ocean, from shapefiles from ISA & chunk above) = USD 3.46509363e+12 for all PMN in the Indian Ocean (seems too much)

#' Get mineral shp with USD value into raster

# load CFC areas polygon shapefile 

sfCFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")%>%
 st_transform(cCRS) # Mollweide proj

sfPMN <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")%>%
  st_transform(cCRS) # Mollweide proj

# attribute average cost to mineral type (based on area covered, calculation of area explained with intersection above) + equation to get sqkm value outlined in Msc thesis

sfCFC$value <- 94569.447 # sqkm USD value (might depend based on estimates of metal concentration and bulk density and crust thickness, see Mitzell et al., 2022 iin Sharma, 2022; might depend on metal prices, see Li et al., 2022)

sfPMN$value <- 14307812.057 # sqkm USD value (idem)

# create template raster based on PU to matching resolution, and intersect with sf objects (also because high resolution, e.g. res=100 generates error)

rPUs <- raster::raster(PUs, res=1000) # Formal class RasterLayer, resolution 1 sqkm
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

max(PUs$MiningCost) # USD 14307812352, correct

##' Get the fishing-specific cost layer

# the method below is the same as for shipping and mining - might be better to have consistent method to compare? also, max and mean value make more sense like this, I think, and we are sure to have the sum of values per PU

FishingRaster <- terra::rast("Data/Cost/Cost_Raster_Sum.grd")
FishingRobin <- terra::project(FishingRaster, "ESRI:54009") # transform to Mollweide proj
PUxFishing <- exactextractr::exact_extract(FishingRobin, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(FishingCost = PUxFishing) # add column to PU data frame

max(PUs$Fishing)
#[1] 9523.993 USD/PU/yr
mean(PUs$Fishing)
#[1] 57.31064 USD/PU/yr

#(ggCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

##' Get the shipping-specific cost layer

Shipping <- terra::rast("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Shipping/shipping.tif")
ShippingRobin <- terra::project(Shipping, "ESRI:54009") # transform to Mollweide proj
PUxShipping <- exactextractr::exact_extract(ShippingRobin, PUs, "sum") # attribute sum of value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(ShippingIntensity = PUxShipping) # add column to PU data frame

max(PUs$ShippingIntensity)
#[1] 382940.5
mean(PUs$ShippingIntensity)
#[1] 6336.166

saveRDS(PUs, paste0("~/Documents/MscThesis/GitHub/SpatialPlanning/rds", "PUs_full.rds"))

##' Prioritizr problems_ sector-specific plans

#' FishPlan problem 

p_FishPlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

# or with cbc solver

p_FishPlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_cbc_solver(gap = 0, verbose = TRUE)

p_FishPlanA_Sol <- solve(p_FishPlanA)

#' MinePlan problem (no longhurst, no shipping)

p_MinePlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #0.6 vents otherwise not feasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

# or with cbc

p_MinePlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #0.6 vents otherwise not feasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_cbc_solver(gap = 0, verbose = TRUE)

p_MinePlanA_Sol <- solve(p_MinePlanA, force = TRUE) # need force = TRUE because planning units with very high (> 1e+6) cost values


#' ShipPlan problems (ShipPlan = shipping intensity as cost, no USD values)

p_ShipPlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents = 0.68 otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

# with cbc

p_ShipPlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents = 0.68 otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_cbc_solver(gap = 0, verbose = TRUE)

p_ShipPlanA_Sol <- solve(p_ShipPlanA, force = TRUE)

#' solution dfs of sector-specific plans

df <- p_FishPlanA_Sol%>%
  dplyr::select(geometry, solution_1)

df2 <- p_MinePlanA_Sol%>%
  dplyr::select(geometry, solution_1)

df3 <- p_ShipPlanA_Sol%>%
  dplyr::select(geometry, solution_1)

#' ## Ferrier scores

FerrierFishA <- eval_ferrier_importance(p_FishPlanA, df)
FerrierMineA <- eval_ferrier_importance(p_MinePlanA, df2)
FerrierShipA <- eval_ferrier_importance(p_ShipPlanA, df3)

#' ## Replacement cost scores
# https://www.sciencedirect.com/science/article/pii/S0006320706001820?casa_token=nbc0AGdxJ8IAAAAA:4XoXBy88sLnIzJMfRAq6gz5FkI5UMCNA7oRcxgXAcew0_dpDCJkbdzJI7eIq8WjL-sKaRsFuN-Qi
# A replacement-cost value of zero tells us that there exists an alternative solution with the same properties as the current (best) solution has, i.e. same cost, and same biodiversity value (although probably obtained with representation levels that are different from those in the original optimal selection). A replacement cost larger than zero means that any alternative solution including/excluding the focal site will have either a lower biodiversity value or a larger cost than the optimal one.

ReplacementFishA <- eval_replacement_importance(p_FishPlanA, df, rescale = FALSE)
ReplacementMineA <- eval_replacement_importance(p_MinePlanA, df2, rescale = FALSE, force = TRUE)
ReplacementShipA <- eval_replacement_importance(p_ShipPlanA, df3, rescale = FALSE, force = TRUE)

# Summary fishing

FishBoundSumA <- eval_boundary_summary(p_FishPlanA, df) #88647878
FishCostSumA <- eval_cost_summary(p_FishPlanA, df) #372652
FishFeatRepSumA <- eval_feature_representation_summary(p_FishPlanA, df)
FishPUSumA <- eval_n_summary(p_FishPlanA, df)  #2607
FishTargetSumA <- eval_target_coverage_summary(p_FishPlanA, df)

# Summary mining

MineBoundSumA <- eval_boundary_summary(p_MinePlanA, df2) #75728850
MineCostSumA <- eval_cost_summary(p_MinePlanA, df2) #542316362224
MineFeatRepSumA <- eval_feature_representation_summary(p_MinePlanA, df2)
MinePUSumA <- eval_n_summary(p_MinePlanA, df2) #7274
MineTargetSumA <- eval_target_coverage_summary(p_MinePlanA, df2)

# Summary shipping

ShipBoundSumA <- eval_boundary_summary(p_ShipPlanA, df3) #93680119
ShipCostSumA <- eval_cost_summary(p_ShipPlanA, df3) #1468585
ShipFeatRepSumA <- eval_feature_representation_summary(p_ShipPlanA, df3)
ShipPUSumA <- eval_n_summary(p_ShipPlanA, df3) #5630
ShipTargetSumA <- eval_target_coverage_summary(p_ShipPlanA, df3)

##' Cross-sectoral plan
#' 
#' need to know the relative costs ('budgets') from sector-specific plans as minimum input for i, j, k (can also do without, i.e., make loop from 0-100% but will take a really long time to run and you will have to look through numerous error messages)
#' relative cost = cost(sector_specific_plan_solution)/sum(PUs$sectoral_cost_layer)*100
#' depending on available time and precision needed, change increments
#' 
for (i in seq(0.0, 0.02, by= 0.01)){ #0.00 because min shipping budget to solve shipping-specific problem = 0.6% = 0.006 < 0.01 and we use increments of 0.01
  for (j in seq(0.01, 0.03, by= 0.01)){ #0.01 because min mining budget to solve mining-specific problem = 1.6% = 0.016
    for (k in seq(0.19, 0.21, by= 0.01)){ #0.19 because min fishing budget to solve fishing-specific problem = 19.5% = 0.19
      tryCatch({
        p_cross_plan <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Area") %>% 
          add_min_set_objective() %>%
          add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
          add_linear_constraints(threshold = i*211501215, # i*sum(PUs$ShippingIntensity)
                                 sense = "<=",
                                 data = PUs$ShippingIntensity) %>%
          add_linear_constraints(threshold = j*3.49048e+13, # j*sum(PUs$MiningCost)
                                 sense = "<=",
                                 data = PUs$MiningCost) %>%
          add_linear_constraints(threshold = k*1913029, # k*sum(PUs$Fishing)
                                 sense = "<=",
                                 data = PUs$Fishing) %>%
          add_binary_decisions() %>%
          add_cbc_solver(gap = 0.01, verbose = TRUE)
        
        p_cross_plan_sol <- solve(p_cross_plan, force = TRUE) 
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      print(i)
      print(j)
      print(k)
      
    }
  }
}
# now, check what the smallest values are for i, j and k where you DON'T get an error value. You should find i = 0.01, j = 0.02 and k = 0.2.

# create a problem with those values just to check

p_PlanFinal <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Area") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.67, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
  add_linear_constraints(threshold = 0.01*211501215, # plan should have less than 2% of total ships in selected PUs
                         sense = "<=",
                         data = PUs$ShippingIntensity) %>%
  add_linear_constraints(threshold = 0.02*3.49048e+13, # plan should cost less than 2% of the total mining opportunity cost
                         sense = "<=",
                         data = PUs$MiningCost) %>%
  add_linear_constraints(threshold = 0.20*1913029, # plan should cost less than 20% of the total fishing opportunity cost
                         sense = "<=",
                         data = PUs$Fishing) %>%
  add_locked_out_constraints("locked_out_mining") %>%
  add_binary_decisions() %>%
  add_cbc_solver(gap = 0, verbose = TRUE)

p_PlanFinal_Sol <- solve(p_PlanFinal, force = TRUE) # all good!


##' Sensitivity analysis
##' Sector-specific plans
## Importance of each cons feat for the fishing plan

library(prioritizr)
library(tidyverse)
library(magrittr)

p_FishPlanA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

# no SOI32

p_FishPlan_No_SIO32 <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_SIO32_Sol <- solve(p_FishPlan_No_SIO32)

df_NoSIO32 <- p_FishPlan_No_SIO32_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_NoSIO32 <- eval_boundary_summary(p_FishPlan_No_SIO32, df_NoSIO32) #88314357
FishCost_NoSIO32 <- eval_cost_summary(p_FishPlan_No_SIO32, df_NoSIO32) # 49727  and 372652 with SIO32, and 49727 /sum(PUs$Fishing)*100 = 2.599385 -> when removing SIO32, the reserve for fishing cost 2.6% of the total fishing budget 
FishPU_NoSIO32 <- eval_n_summary(p_FishPlan_No_SIO32, df_NoSIO32)  #2607 compared to 2605

# no EBSA

p_FishPlan_No_EBSA <- problem(PUs, features = c("IBA", "IMMA","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_EBSA_Sol <- solve(p_FishPlan_No_EBSA)

df_No_EBSA <- p_FishPlan_No_EBSA_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_EBSA <- eval_boundary_summary(p_FishPlan_No_SIO32, df_No_EBSA) #58993951
FishCost_No_EBSA <- eval_cost_summary(p_FishPlan_No_SIO32, df_No_EBSA) # 22384
FishPU_No_EBSA <- eval_n_summary(p_FishPlan_No_EBSA, df_No_EBSA)  #1402

## no IBA

p_FishPlan_No_IBA <- problem(PUs, features = c("IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_IBA_Sol <- solve(p_FishPlan_No_IBA)

df_No_IBA <- p_FishPlan_No_IBA_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_IBA <- eval_boundary_summary(p_FishPlan_No_IBA, df_No_IBA) #87774838
FishCost_No_IBA <- eval_cost_summary(p_FishPlan_No_IBA, df_No_IBA) # 369218
FishPU_No_IBA <- eval_n_summary(p_FishPlan_No_IBA, df_No_IBA)  #2545

## no IMMA

p_FishPlan_No_IMMA <- problem(PUs, features = c("IBA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_IMMA_Sol <- solve(p_FishPlan_No_IMMA )

df_No_IMMA  <- p_FishPlan_No_IMMA_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_IMMA <- eval_boundary_summary(p_FishPlan_No_IMMA, df_No_IMMA) #77690737
FishCost_No_IMMA  <- eval_cost_summary(p_FishPlan_No_IMMA, df_No_IMMA) # 370991
FishPU_No_IMMA  <- eval_n_summary(p_FishPlan_No_IMMA, df_No_IMMA)  #2312

## no seamounts

p_FishPlan_No_seam <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_seam_Sol <- solve(p_FishPlan_No_seam )

df_No_seam  <- p_FishPlan_No_seam_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_seam <- eval_boundary_summary(p_FishPlan_No_seam, df_No_seam) #66782644
FishCost_No_seam  <- eval_cost_summary(p_FishPlan_No_seam, df_No_seam) # 369248
FishPU_No_seam  <- eval_n_summary(p_FishPlan_No_seam, df_No_seam)  #2211

## no plateaus

p_FishPlan_No_plat <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Active_Vents", "Inactive_Vents"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_plat_Sol <- solve(p_FishPlan_No_plat )

df_No_plat  <- p_FishPlan_No_plat_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_plat <- eval_boundary_summary(p_FishPlan_No_plat, df_No_plat) #53373144
FishCost_No_plat  <- eval_cost_summary(p_FishPlan_No_plat, df_No_plat) # 363723
FishPU_No_plat  <- eval_n_summary(p_FishPlan_No_seam, df_No_plat)  #1897

## no vents

p_FishPlan_No_vents <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Plateaus"), cost_column = "FishingCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3,0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_FishPlan_No_vents_Sol <- solve(p_FishPlan_No_vents)

df_No_vents  <- p_FishPlan_No_vents_Sol%>%
  dplyr::select(geometry, solution_1)

FishBound_No_vents <- eval_boundary_summary(p_FishPlan_No_vents, df_No_vents) #63143343
FishCost_No_vents  <- eval_cost_summary(p_FishPlan_No_vents, df_No_vents) # 368477.
FishPU_No_vents  <- eval_n_summary(p_FishPlan_No_vents, df_No_vents)  #2180

## Importance of each cons feat for the shipping plan

# no SOI32

p_ShipPlan_No_SIO32 <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas")%>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_SIO32_Sol <- solve(p_ShipPlan_No_SIO32, force = TRUE)

df_NoSIO32 <- p_ShipPlan_No_SIO32_Sol%>% # ! careful, ship and mine dfs have same name - be sure to  run all code
  dplyr::select(geometry, solution_1)

ShipBound_NoSIO32 <- eval_boundary_summary(p_ShipPlan_No_SIO32, df_NoSIO32) 
ShipCost_NoSIO32 <- eval_cost_summary(p_ShipPlan_No_SIO32, df_NoSIO32) # 
ShipPU_NoSIO32 <- eval_n_summary(p_ShipPlan_No_SIO32, df_NoSIO32)  #

# no EBSA

p_ShipPlan_No_EBSA <- problem(PUs, features = c("IBA", "IMMA","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_EBSA_Sol <- solve(p_ShipPlan_No_EBSA, force = TRUE)

df_No_EBSA <- p_ShipPlan_No_EBSA_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_EBSA <- eval_boundary_summary(p_ShipPlan_No_SIO32, df_No_EBSA) #
ShipCost_No_EBSA <- eval_cost_summary(p_ShipPlan_No_SIO32, df_No_EBSA) # 
ShipPU_No_EBSA <- eval_n_summary(p_ShipPlan_No_EBSA, df_No_EBSA)  #

# no IBA

p_ShipPlan_No_IBA <- problem(PUs, features = c("IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_IBA_Sol <- solve(p_ShipPlan_No_IBA, force = TRUE)

df_No_IBA <- p_ShipPlan_No_IBA_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_IBA <- eval_boundary_summary(p_ShipPlan_No_IBA, df_No_IBA) #
ShipCost_No_IBA <- eval_cost_summary(p_ShipPlan_No_IBA, df_No_IBA) # 
ShipPU_No_IBA <- eval_n_summary(p_ShipPlan_No_IBA, df_No_IBA)  #

# no IMMA

p_ShipPlan_No_IMMA <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_IMMA_Sol <- solve(p_ShipPlan_No_IMMA, force = TRUE)

df_No_IMMA  <- p_ShipPlan_No_IMMA_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_IMMA <- eval_boundary_summary(p_ShipPlan_No_IMMA, df_No_IMMA) #
ShipCost_No_IMMA  <- eval_cost_summary(p_ShipPlan_No_IMMA, df_No_IMMA) # 
ShipPU_No_IMMA  <- eval_n_summary(p_ShipPlan_No_IMMA, df_No_IMMA)  #

# no seamounts

p_ShipPlan_No_seam <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_seam_Sol <- solve(p_ShipPlan_No_seam, force = TRUE)

df_No_seam  <- p_ShipPlan_No_seam_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_seam <- eval_boundary_summary(p_ShipPlan_No_seam, df_No_seam) #
ShipCost_No_seam  <- eval_cost_summary(p_ShipPlan_No_seam, df_No_seam) # 
ShipPU_No_seam  <- eval_n_summary(p_ShipPlan_No_seam, df_No_seam)  #

# no plateaus

p_ShipPlan_No_plat <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Active_Vents", "Inactive_Vents"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_plat_Sol <- solve(p_ShipPlan_No_plat, force = TRUE)

df_No_plat  <- p_ShipPlan_No_plat_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_plat <- eval_boundary_summary(p_ShipPlan_No_plat, df_No_plat) #
ShipCost_No_plat  <- eval_cost_summary(p_ShipPlan_No_plat, df_No_plat) #
ShipPU_No_plat  <- eval_n_summary(p_ShipPlan_No_seam, df_No_plat)  #

# no vents

p_ShipPlan_No_vents <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Plateaus"), cost_column = "ShippingIntensity") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_ShipPlan_No_vents_Sol <- solve(p_ShipPlan_No_vents, force = TRUE)

df_No_vents  <- p_ShipPlan_No_vents_Sol%>%
  dplyr::select(geometry, solution_1)

ShipBound_No_vents <- eval_boundary_summary(p_ShipPlan_No_vents, df_No_vents) #
ShipCost_No_vents  <- eval_cost_summary(p_ShipPlan_No_vents, df_No_vents) # 
ShipPU_No_vents  <- eval_n_summary(p_ShipPlan_No_vents, df_No_vents)  #

## Importance of each cons feat for the mining plan

# no SOI32

p_MinePlan_No_SIO32 <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas")%>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_SIO32_Sol <- solve(p_MinePlan_No_SIO32)

df_NoSIO32 <- p_MinePlan_No_SIO32_Sol%>% #!  careful same df name as asbove
  dplyr::select(geometry, solution_1)

MineBound_NoSIO32 <- eval_boundary_summary(p_MinePlan_No_SIO32, df_NoSIO32) #
MineCost_NoSIO32 <- eval_cost_summary(p_MinePlan_No_SIO32, df_NoSIO32) # 
MinePU_NoSIO32 <- eval_n_summary(p_MinePlan_No_SIO32, df_NoSIO32)  #

# no EBSA

p_MinePlan_No_EBSA <- problem(PUs, features = c("IBA", "IMMA","SIO_36","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% # this is a little trick to avoid errors - I'm keeping an EBSA whose size is 1 PU that overlaps with a CFC area, otherwise cost is 0 and the equation cannot be solved properly
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_EBSA_Sol <- solve(p_MinePlan_No_EBSA, force = TRUE)

df_No_EBSA <- p_MinePlan_No_EBSA_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_EBSA <- eval_boundary_summary(p_MinePlan_No_EBSA, df_No_EBSA) #
MineCost_No_EBSA <- eval_cost_summary(p_MinePlan_No_EBSA, df_No_EBSA) # 
MinePU_No_EBSA <- eval_n_summary(p_MinePlan_No_EBSA, df_No_EBSA)  #

# no IBA

p_MinePlan_No_IBA <- problem(PUs, features = c("IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_IBA_Sol <- solve(p_MinePlan_No_IBA, force = TRUE)

df_No_IBA <- p_MinePlan_No_IBA_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_IBA <- eval_boundary_summary(p_MinePlan_No_IBA, df_No_IBA) #
MineCost_No_IBA <- eval_cost_summary(p_MinePlan_No_IBA, df_No_IBA) # 
MinePU_No_IBA <- eval_n_summary(p_MinePlan_No_IBA, df_No_IBA)  #

# no IMMA

p_MinePlan_No_IMMA <- problem(PUs, features = c("IBA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_IMMA_Sol <- solve(p_MinePlan_No_IMMA, force = TRUE)

df_No_IMMA  <- p_MinePlan_No_IMMA_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_IMMA <- eval_boundary_summary(p_MinePlan_No_IMMA, df_No_IMMA) #
MineCost_No_IMMA  <- eval_cost_summary(p_Minelan_No_IMMAA, df_No_IMMA) # 
MinePU_No_IMMA  <- eval_n_summary(p_MinePlan_No_IMMA, df_No_IMMA)  #

# no seamounts

p_MinePlan_No_seam <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_seam_Sol <- solve(p_MinePlan_No_seam, force = TRUE)

df_No_seam  <- p_FishPlan_No_seam_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_seam <- eval_boundary_summary(p_MinePlan_No_seam, df_No_seam) #
MineCost_No_seam  <- eval_cost_summary(p_MinePlan_No_seam, df_No_seam) # 
MinePU_No_seam  <- eval_n_summary(p_MinePlan_No_seam, df_No_seam)  #

# no plateaus

p_MinePlan_No_plat <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_plat_Sol <- solve(p_MinePlan_No_plat, force = TRUE)

df_No_plat  <- p_MinePlan_No_plat_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_plat <- eval_boundary_summary(p_MinePlan_No_plat, df_No_plat) #
MineCost_No_plat  <- eval_cost_summary(p_MinePlan_No_plat, df_No_plat) # 
MinePU_No_plat  <- eval_n_summary(p_MinePlan_No_seam, df_No_plat)  #

# no vents

p_MinePlan_No_vents <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37", "Seamounts", "Plateaus"), cost_column = "MiningCost") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3,0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3)) %>% #inactive vents target = 0.3, active vents  = 0.68, otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out_areas") %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_MinePlan_No_vents_Sol <- solve(p_MinePlan_No_vents, force = TRUE)

df_No_vents  <- p_MinePlan_No_vents_Sol%>%
  dplyr::select(geometry, solution_1)

MineBound_No_vents <- eval_boundary_summary(p_MinePlan_No_vents, df_No_vents) #
MineCost_No_vents  <- eval_cost_summary(p_MinePlan_No_vents, df_No_vents) # 
MinePU_No_vents  <- eval_n_summary(p_MinePlan_No_vents, df_No_vents)  #

##' Sensitivity of cross-sectoral plan to budgets
## Change increments and max budgets as needed

PUs$Area <- 1

df_inc0.01 <- data.frame() #create empty dataframe to add to with each iteration
for (i in seq(0.01, 0.11, by= 0.01)){ #0.01 because min shipping budget to solve cross-sectoral problem = 1%
  for (j in seq(0.02, 0.12, by= 0.01)){ #0.02 for mining
    for (k in seq(0.2, 0.30, by= 0.01)){ #0.2 for fishing
      p_PlanF <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Area") %>% 
        add_min_set_objective() %>%
        add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.68, 0.3)) %>% #inactive vents target = 0.3 otherwise unfeasible
        add_linear_constraints(threshold = i*211501215, # i*sum(PUs$ShippingIntensity)
                               sense = "<=",
                               data = PUs$ShippingIntensity) %>%
        add_linear_constraints(threshold = k*1913029, # k*sum(PUs$FishingCost)
                               sense = "<=",
                               data = PUs$Fishing) %>%
        add_linear_constraints(threshold = j*3.49048e+13, # j*sum(PUs$MiningCost)
                               sense = "<=",
                               data = PUs$MiningCost) %>%
        add_binary_decisions() %>%
        add_cbc_solver(gap = 0, verbose = TRUE)
      
      p_PlanF_sol <- solve(p_PlanF, force = TRUE) 
      
      #extract the information you're interested in
      #  sol <- p_PlanF_Sol$solution_1
      df_sol <- p_PlanF_sol%>%
        dplyr::select(geometry, solution_1)
      
      cost <- eval_cost_summary(p_PlanF, df_sol) # cost = area because cost layer = 1 for each PU
      #    area <- eval_n_summary(p_PlanF, df_sol)
      
      total <- cbind(data.frame(cost), data.frame(i), data.frame(j), data.frame(k)) #
      df_inc0.01 <- rbind(df_inc0.01, data.frame(total)) #add to to dataframe
    }
  }
}

##' Plotting

library(tidyverse)
library(sf)
library(ggplot2)

# PUs_plot

PUs_plot <- ggplot() + 
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  geom_sf(data = PUs,  aes(color="D"), show.legend = "polygon") +
  scale_color_manual(values = c("D" = "lightgrey"),
                     labels = c("ABNJ"),
                     name = "Legend") +
  ggtitle("Planning region") + 
  #coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

PUs_plot

# vents

df_vents <- PUs %>% dplyr::mutate(Legend = case_when((Active_Vents == TRUE & Inactive_Vents == FALSE) ~ "Active vents",
                                                     (Active_Vents == FALSE & Inactive_Vents == TRUE) ~ "Inactive vents",
                                                     (Active_Vents == TRUE & Inactive_Vents == TRUE) ~ "Both",
                                                     (Active_Vents == FALSE & Inactive_Vents == FALSE) ~ "Neither")) %>% 
  dplyr::filter(Legend != "Neither")

Vents_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_vents, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("Active vents" = "red", 
                               "Inactive vents" = "tan2",
                               "Both" = "dodgerblue")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  ggtitle("Hydrothermal vents") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

Vents_plot # correct

# seamounts

df_seamounts <- PUs %>% dplyr::mutate(Legend = case_when((Seamounts == TRUE ) ~ "Seamount",
                                                         (Seamounts == FALSE ) ~ "Absent")) %>% 
  dplyr::filter(Legend != "Absent")

Seamounts_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_seamounts,  aes(fill= Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("Seamount" = "darkblue")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  ggtitle("Seamounts") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))
 
Seamounts_plot

# IMMAs

df_IMMA <- PUs %>% dplyr::mutate(Legend = case_when((IMMA == TRUE) ~ "IMMA",
                                                    (IMMA == FALSE) ~ "Absent"))%>% 
  dplyr::filter(Legend != "Absent")

IMMA_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_IMMA, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("IMMA" = "lightblue")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Important Marine Mammal Areas") 

IMMA_plot

# IBAs

df_IBA <- PUs %>% dplyr::mutate(Legend = case_when((IBA == TRUE) ~ "IBA",
                                                   (IBA == FALSE) ~ "Absent"))%>% 
  dplyr::filter(Legend != "Absent")

IBA_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_IBA, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("IBA" = "cornflowerblue")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  ggtitle("Important Bird Areas") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14)) 

IBA_plot

# plateaus 

df_plateaus <- PUs %>% dplyr::mutate(Legend = case_when((Plateaus == TRUE) ~ "Plateau",
                                                        (Plateaus == FALSE) ~ "Absent"))%>% 
  dplyr::filter(Legend != "Absent")

Plateaus_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_plateaus, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("Plateau" = "cadetblue")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  ggtitle("Plateaus") +
theme(axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14)) 

Plateaus_plot

# EBSAs

df_EBSA <- PUs %>% dplyr::mutate(Legend = case_when((NWIO_14 == TRUE) ~ "Arabian Oxygen Minimum Zone",
                                                    (NEIO_09 == TRUE) ~ "Sumatra Upwelling Zone",
                                                    (SIO_11 == TRUE) ~  "Agulhas Front",
                                                    (SIO_19 == TRUE) ~ "Mozambique Channel",
                                                    (SIO_22 == TRUE) ~ "Walters Shoals",
                                                    (SIO_32 == TRUE) ~ "Saya de Malha Bank",
                                                    (SIO_23 == TRUE & SIO_11 == TRUE) ~ "Coral Seamount Fracture Zone",
                                                    (SIO_35 == TRUE) ~ "Rusky Knoll",
                                                    (SIO_36 == TRUE) ~ "Fools' Flat",
                                                    (SIO_37 == TRUE) ~ "East Broken Ridge",
                                                    (NWIO_14 == FALSE & NEIO_09 == FALSE & SIO_11 == FALSE & SIO_19 == FALSE & SIO_22 == FALSE & SIO_23 == FALSE & SIO_32 == FALSE & SIO_35 == FALSE & SIO_36 == FALSE & SIO_37 == FALSE) ~ "Neither")) %>% 
  #(Active_Vents == FALSE & Inactive_Vents == FALSE) ~ "Neither")) %>% 
  dplyr::filter(Legend != "Neither")

length(SIO_23 == TRUE) #5 -> small on map

EBSA_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_EBSA, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("Arabian Oxygen Minimum Zone" = "dodgerblue4", 
                               "Sumatra Upwelling Zone" = "dodgerblue1",
                               "Agulhas Front" = "deepskyblue1",
                               "Mozambique Channel" = "lightskyblue2",
                               "Walters Shoals" = "paleturquoise",
                               "Saya de Malha Bank" = "palegreen",
                               "Coral Seamount Fracture Zone" = "greenyellow",
                               "Rusky Knoll" = "olivedrab3",
                               "Fools' Flat" = "yellow",
                               "East Broken Ridge" = "darkgoldenrod")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Ecologically and Biologically Significant Areas") 

EBSA_plot

### contractors (mining), see ISA reserved areas and exploration leases

df_contracts <- PUs %>% dplyr::mutate(Legend = case_when((ExplPMN == TRUE) ~ "Nodules, Exploration",
                                                         (ResPMN == TRUE) ~ "Nodules, Reserved",
                                                         (ExplPMS == TRUE) ~ "Sulphides, Exploration",
                                                         (ExplPMN == FALSE & ExplPMS == FALSE & ResPMN == FALSE) ~ "Neither"))%>% 
  dplyr::filter(Legend != "Neither")


Contracts_plot <- ggplot() + 
  geom_sf(data = PUs, color = "lightgrey", show.legend = FALSE) +
  geom_sf(data = df_contracts, aes(fill = Legend), color = NA, show.legend = TRUE, size = 0.1) +
  scale_fill_manual(values = c("Nodules, Exploration" = "goldenrod2",
                               "Nodules, Reserved" = "goldenrod3",
                               "Sulphides, Exploration" = "red2")) +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Locked-out areas") 

Contracts_plot

## cost layer plots

library(RColorBrewer)
library(patchwork)
library(sf)

# Defining palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu"))) #RdYlBu
nam <- parse(text = "log[10] (cost)")
namrc <- parse(text = "log[10] (score)")
sc_fish <- scale_colour_gradientn(name = nam, colours = myPalette(10), na.value = "lightgrey", limits=c(0, 4), aesthetics = c("color","fill"))
sc_fish_rc <- scale_colour_gradientn(name = namrc, colours = myPalette(10), na.value = "lightgrey", limits=c(0, 4), aesthetics = c("color","fill"))

# fishing cost plot

fishing_plot <- ggplot()+
  geom_sf(data = PUs, aes(color = log10(Fishing), fill = log10(Fishing))) +
  sc_fish +
  geom_sf(data = world, size = 0.05, fill = "grey20") +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  ggtitle("Opportunity cost for the fishing sector") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

fishing_plot

# FishPlan plot

p_FishPlanA_Sol_Plot <- ggplot() +
  ggtitle("Optimal reserves for the fishing sector") +
  geom_sf(data = world, fill = "grey20", size = 0.05)  +
  geom_sf(data = p_FishPlanA_Sol, aes(fill = as.factor(solution_1)), size = 0.01) +
  scale_fill_manual(values = c("lightgrey", "#009E73"), labels = c("Not selected", "Selected")) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

p_FishPlanA_Sol_Plot

# fishplot rc scores plot

rc_fish_plot <- ggplot() + 
  sc_fish_rc +
  geom_sf(data = world, fill = "grey20", size = 0.05) +
  geom_sf(data = ReplacementFishA, aes(color=log10(rc), fill = log10(rc))) + 
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  ggtitle("Replacement cost score for fishing-specific reserves") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

rc_fish_plot

#load packages 

# Mining cost plot - plot cost layer of mineral resources

# Defining palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu"))) #RdYlBu
sc_mine <- scale_colour_gradientn(name = nam, colours = myPalette(100), na.value = "lightgrey", limits=c(0, 11), aesthetics = c("color","fill"))
sc_mine_rc <- scale_colour_gradientn(name = namrc, colours = myPalette(100), na.value = "lightgrey", limits=c(0, 11), aesthetics = c("color","fill"))

mining_plot <- ggplot()+
  geom_sf(data = PUs, aes(color = log10(MiningCost), fill = log10(MiningCost))) +
  sc_mine +
  geom_sf(data = world, size = 0.05, fill = "grey20") +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  ggtitle("Opportunity cost for the deep-sea mining sector") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

mining_plot

#' PLotting MinePlan - mining-specific plan

p_MinePlanA_Sol_Plot <- ggplot() +
  ggtitle("Optimal reserves for the deep-sea mining sector") +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE)+
  geom_sf(data = p_MinePlanA_Sol, aes(fill = as.factor(solution_1)), size = 0.01) +
  scale_fill_manual(values = c("lightgrey", "#009E73"), labels = c("Not selected", "Selected")) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

p_MinePlanA_Sol_Plot

# rc mine - plot replacement cost scores of PUs in the mining-specific plan

rc_mine_plot <- ggplot() + 
  sc_mine_rc +
  geom_sf(data = world, fill = "grey20", size = 0.05) +
  geom_sf(data = ReplacementMineA, aes(color=log10(rc), fill = log10(rc))) + 
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  ggtitle("Replacement cost score for mining-specific reserves") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

rc_mine_plot

# shipping cost plot

# Defining palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu"))) #RdYlBu
sc_ship <- scale_colour_gradientn(name = nam, colours = myPalette(100), na.value = "lightgrey", limits=c(0, 5), aesthetics = c("color","fill"))
sc_ship_rc <- scale_colour_gradientn(name = namrc, colours = myPalette(100), na.value = "lightgrey", limits=c(0, 5), aesthetics = c("color","fill"))

# for plotting purposes
#PUs$ShippingIntensity[PUs$ShippingIntensity <  1] <-  1

shipping_plot <- ggplot() +
  geom_sf(data = PUs, aes(color = log10(ShippingIntensity), fill = log10(ShippingIntensity))) +
  sc_ship +
  ggtitle("Opportunity cost for the shipping sector") +
  geom_sf(data = world, size = 0.05, fill = "grey20") +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

shipping_plot

p_ShipPlanA_Sol_Plot <- ggplot() +
  ggtitle("Optimal reserves for the shipping sector") +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE)+
  geom_sf(data = p_ShipPlanA_Sol, aes(fill = as.factor(solution_1)), size = 0.01) +
  scale_fill_manual(values = c("lightgrey", "#009E73"), labels = c("Not selected", "Selected")) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

p_ShipPlanA_Sol_Plot

# rc ship

rc_ship_plot <- ggplot() + 
  sc_ship_rc +
  ggtitle("Replacement cost score for shipping-specific reserves") +
  geom_sf(data = world, fill = "grey20", size = 0.05)+
  geom_sf(data = ReplacementShipA, aes(color=log10(rc), fill = log10(rc))) + 
  sc_ship_rc +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

rc_ship_plot

# constraints

cross_sectoral_plot <- ggplot() +
  geom_sf(data = world, fill = "grey20", size = 0.05, show.legend = FALSE) +
  geom_sf(data = p_PlanFinal_Sol, aes(fill = as.factor(solution_1)), size = 0.01) +
  scale_fill_manual(values = c("lightgrey", "#009E73"), labels = c("Not selected", "Selected")) +
  coord_sf(xlim = c(st_bbox(Bndry)$xmin, st_bbox(Bndry)$xmax), 
           ylim = c(st_bbox(Bndry)$ymin, st_bbox(Bndry)$ymax),
           expand = TRUE) +
  theme_bw() +
  # ggtitle("Optimal reserves for all sectors") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  theme(panel.grid.major=element_line(colour="lightgrey", size = 0.2),
        panel.grid.minor=element_line(colour="lightgrey", size = 0.2)) +
  theme(plot.title = element_text(size=14))

cross_sectoral_plot

# kappa correlation matrix plot

list_plans <- list(spatial_plan1, spatial_plan2) #with the "solution_1" column of each of the spatial plans already in your desired "plan names"
# This function creates Cohen's Kappa Correlation Matrix
create_corrmatrix <- function(list_plans) {
  pacman::p_load(irr)
  
  y = 1
  s_matrix <- list() # empty list
  for(i in 1:length(list_plans)){
    for(j in 1:length(list_plans)){
      kappa_temp <- irr::kappa2(bind_cols(list_plans[[i]], list_plans[[j]]))
      kappa_corrvalue <- kappa_temp$value
      kappa_pvalue <- kappa_temp$p.value
      s_matrix[[y]] <- cbind(colnames(list_plans[[i]]), colnames(list_plans[[j]]), kappa_corrvalue, kappa_pvalue)
      y = y+1
    }
  }
  
  s_matrix_all <- do.call(rbind, s_matrix) %>% 
    as_tibble()
  colnames(s_matrix_all)[1:2] <- c('plan1','plan2')
  
  matrix <- s_matrix_all %>% 
    as_tibble() %>% 
    dplyr::select(-kappa_pvalue) %>% 
    pivot_wider(names_from = plan2, values_from = kappa_corrvalue) %>% 
    as.matrix()
  
  return(matrix)
}

# This plots the Correlation Matrix.
num = length(list_plans)
plot_corrplot <- function(matrix, num) {
  pacman::p_load(corrplot)
  # creating corrplot
  rownames(matrix) <- matrix[,1]
  n <- num + 1 # num represents the number of inputted spatial plans
  matrix_f <- matrix[,2:n]
  class(matrix_f) <- "numeric"
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  plot <- corrplot(matrix_f, method = "shade", col.lim = c(-0.2,1), tl.col = "black", addCoef.col = "black",
                   col=col(200), tl.srt=45)
  return(plot)
}

### apply function to plans of interest

spatial_plan1 <- as.data.frame(p_FishPlanA_Sol$solution_1)
spatial_plan2 <- as.data.frame(p_ShipPlanA_Sol$solution_1)
spatial_plan3 <- as.data.frame(p_MinePlanA_Sol$solution_1)
spatial_plan4 <- as.data.frame(p_PlanFinal_Sol$solution_1)

#spatial_plan1 <- p_FishPlanA_Sol to check
#spatial_plan2 <- p_ShipPlanA_Sol to check

list_plans <- list(spatial_plan1, spatial_plan2, spatial_plan3, spatial_plan4) #with the "solution_1" column of each of the spatial plans already in your desired "plan names"
# This function creates Cohen's Kappa Correlation Matrix
create_corrmatrix <- function(list_plans) {
  pacman::p_load(irr)
  
  y = 1
  s_matrix <- list() # empty list
  for(i in 1:length(list_plans)){
    for(j in 1:length(list_plans)){
      kappa_temp <- irr::kappa2(bind_cols(list_plans[[i]], list_plans[[j]]))
      kappa_corrvalue <- kappa_temp$value
      kappa_pvalue <- kappa_temp$p.value
      s_matrix[[y]] <- cbind(colnames(list_plans[[i]]), colnames(list_plans[[j]]), kappa_corrvalue, kappa_pvalue)
      y = y+1
    }
  }
  
  s_matrix_all <- do.call(rbind, s_matrix) %>% 
    as_tibble()
  colnames(s_matrix_all)[1:2] <- c('plan1','plan2')
  
  write_csv(s_matrix_all, "~/Documents/MscThesis/GitHub/SpatialPlanning/s_matrix_all.csv")
  
  matrix <- s_matrix_all %>% 
    as_tibble() %>% 
    dplyr::select(-kappa_pvalue) %>% 
    pivot_wider(names_from = plan2, values_from = kappa_corrvalue) %>% 
    as.matrix()
  
  return(matrix)
}

matrix_try <- create_corrmatrix(list_plans = list_plans)

# This plots the Correlation Matrix.
num = length(list_plans)
plot_corrplot <- function(matrix, num) {
  pacman::p_load(corrplot)
  # creating corrplot
  rownames(matrix) <- matrix[,1]
  n <- num + 1 # num represents the number of inputted spatial plans
  matrix_f <- matrix[,2:n]
  class(matrix_f) <- "numeric"
  
  colnames(matrix_f) <- c("Fishing-specific plan", "Shipping-specific plan", "Mining-specific plan", "Cross-sectoral plan")
  rownames(matrix_f) <- c("Fishing-specific plan", "Shipping-specific plan", "Mining-specific plan", "Cross-sectoral plan")
  
  #col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  scaleblue <- colorRampPalette(brewer.pal(4, "Blues"))(4)
  plot <- corrplot(matrix_f, method = "shade", type = 'upper', col.lim = c(0,1), bg = "white", tl.col = "black", addCoef.col = "black",
                   col= scaleblue, tl.srt=45, cl.pos="n")
  return(plot)
}

kappa_plot <- plot_corrplot(matrix = matrix_try, num = num)

ggsave("pdfs/kappa.jpeg", width = 10, height = 5, dpi = 600) # ideal dimensions for single figure

# circular barplot of conservation targets reached / exceeded, all sectors

Targets <- c(FishFeatRepSumA$relative_held*100, MineFeatRepSumA$relative_held*100, ShipFeatRepSumA$relative_held*100, FeatRepF$relative_held*100)
Features <- c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents")
Scenario=c( rep('Fishing', 17), rep('Mining', 17), rep('Shipping', 17), rep('Cross-sectoral', 17))

cbp_df <- data.frame(Targets, Scenario, Features)
cbp_df_sliced <- dplyr::slice(cbp_df, -60, -43, -26, -9) %>% #remove NaN rows
  dplyr::mutate(individual = dplyr::case_when(stringr::str_detect(Features, pattern = "IBA") ~ "Important Bird Areas",
                                              str_detect(Features, pattern = "IMMA") ~ "Important Marine Mammal Areas",
                                              str_detect(Features, pattern = "SIO_19") ~ "Mozambique Channel",
                                              str_detect(Features, pattern = "SIO_22") ~ "Walters' Shoals",
                                              str_detect(Features, pattern = "SIO_23") ~ "Coral Seamount and Fracture Zone",
                                              str_detect(Features, pattern = "SIO_32") ~ "Saya de Malha Bank",
                                              #  str_detect(Features, pattern = "SIO_34") ~ "Central Indian Ocean Basin",
                                              str_detect(Features, pattern = "SIO_11") ~ "Agulhas Upwelling Zone",
                                              str_detect(Features, pattern = "NEIO_09") ~ "Sumatra Upwelling Zone",
                                              str_detect(Features, pattern = "NWIO_14") ~ "Arabian Oxygen Minimum Zone",
                                              # str_detect(Features, pattern = "SIO_30") ~ "Atlantis seamount",
                                              str_detect(Features, pattern = "SIO_32") ~ "Saya de Malha Bank",
                                              str_detect(Features, pattern = "SIO_35") ~ "Rusky Knoll",
                                              str_detect(Features, pattern = "SIO_36") ~ "Fools' Flat",
                                              str_detect(Features, pattern = "SIO_37") ~ "East Broken Ridge",
                                              str_detect(Features, pattern = "Seamounts") ~ " Seamounts",
                                              str_detect(Features, pattern = "Inactive_Vents") ~ "Inactive vents",
                                              str_detect(Features, pattern = "Active_Vents") ~ "Active vents",
                                              str_detect(Features, pattern = "Plateaus") ~ "Plateaus"))

#plot 1

cbp_df_sliced_30 <- dplyr::slice(cbp_df_sliced, -1,-6:-12, -15, -17, -22:-28, -31, -33, -38:-44, -47, -49, -54:-60, -63)
target <- as.numeric(30)

colfill <- c("grey85", "seagreen3","dodgerblue2", "cadetblue2")

gg30 <- ggplot(cbp_df_sliced_30, aes(x=individual, y=Targets, fill = Scenario)) + 
  geom_bar(position= 'dodge', stat="identity", width= 0.5) + #labs(title = "Conservation targets reached by sector", x = "Conservation features", y = 'Targets reached (%)') +
  # theme_void() +
  scale_fill_manual(values = colfill) +
  coord_polar(start = 0) +
  ylim(-50,110) +
  geom_errorbar(aes(y = target, ymax = target, ymin = target), 
                color = 'red', 
                linetype = 'dashed', 
                size = 0.4) +
  geom_errorbar(aes(y = 10, ymax = 10, ymin = 10), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 50, ymax = 50, ymin = 50), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 70, ymax = 70, ymin = 70), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 100, ymax = 100, ymin = 100), 
                color = 'grey', 
                linetype = 'solid', 
                size = 0.4) +
  theme(
    axis.line=element_blank(),
    #axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.margin=unit(c(-0.25,0,0,0), "cm"),
    panel.spacing=unit(c(0,0,0,0), "cm"),
    plot.background=element_blank()) +
  
  # Annotate custom scale inside plot
  annotate(
    x = 7.5, 
    y = 15, 
    label = "10", 
    geom = "text", 
    color = "gray"
  ) +
  annotate(
    x = 7.5, 
    y = 55, 
    label = "50", 
    geom = "text", 
    color = "gray"
  ) +
  annotate(
    x = 7.5, 
    y =75, 
    label = "70", 
    geom = "text", 
    color = "gray"
  ) +
  annotate(
    x = 7.5, 
    y =90, 
    label = "100", 
    geom = "text", 
    color = "grey"
  ) +
  annotate(
    x = 7.5, 
    y =35, 
    label = "30", 
    geom = "text", 
    color = "red"
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6)) # to avoid labels going out of panel when too long #https://statisticsglobe.com/wrap-long-axis-labels-ggplot2-plot-into-multiple-lines-r


gg30 +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))#https://statisticsglobe.com/wrap-long-axis-labels-ggplot2-plot-into-multiple-lines-r

## plot 2

cbp_df_sliced_70 <- dplyr::slice(cbp_df_sliced, 1,6:12, 15, 17, 22:28, 31, 33, 38:44, 47, 49, 54:60, 63)

target2 <- as.numeric(70)

#cbp_df_sliced_30 <- dplyr::slice(cbp_df_sliced, -1,-6:-12, -15, -17, -22:-28, -31, -33, -38:-44, -47)
#target <- as.numeric(30)

gg70 <- ggplot(cbp_df_sliced_70, aes(x=individual, y=Targets, fill = Scenario)) + 
  geom_bar(position= 'dodge', stat="identity", width= 0.5) + #labs(title = "Conservation targets reached by sector", x = "Conservation features", y = 'Targets reached (%)') +
  # theme_void() +
  scale_fill_manual(values = colfill) +
  coord_polar(start = 0) +
  ylim(-50,110) +
  geom_errorbar(aes(y = target2, ymax = target2, ymin = target2), 
                color = 'red', 
                linetype = 'dashed', 
                size = 0.4) +
  geom_errorbar(aes(y = 10, ymax = 10, ymin = 10), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 50, ymax = 50, ymin = 50), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 30, ymax = 30, ymin = 30), 
                color = 'grey', 
                linetype = 'dotted', 
                size = 0.4) +
  geom_errorbar(aes(y = 100, ymax = 100, ymin = 100), 
                color = 'grey', 
                linetype = 'solid', 
                size = 0.4) +
  theme(
    axis.line=element_blank(),
    #axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.margin=unit(c(-0.25,0,0,0), "cm"),
    panel.spacing=unit(c(0,0,0,0), "cm"),
    plot.background=element_blank()) +
  
  # Annotate custom scale inside plot
  annotate(
    x = 9.5, 
    y = 15, 
    label = "10", 
    geom = "text", 
    color = "gray"
  ) +
  annotate(
    x = 9.5, 
    y = 55, 
    label = "50", 
    geom = "text", 
    color = "gray"
  ) +
  annotate(
    x = 9.5, 
    y =75, 
    label = "70", 
    geom = "text", 
    color = "red"
  ) +
  annotate(
    x = 9.5, 
    y =90, 
    label = "100", 
    geom = "text", 
    color = "grey"
  ) +
  annotate(
    x = 9.5, 
    y =35, 
    label = "30", 
    geom = "text", 
    color = "grey"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) # to avoid labels going out of panel when too long #https://statisticsglobe.com/wrap-long-axis-labels-ggplot2-plot-into-multiple-lines-r


gg70 +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # to avoid labels going out of panel when too long #https://statisticsglobe.com/wrap-long-axis-labels-ggplot2-plot-into-multiple-lines-r

## sensitivity analysis plots

library(plotly)

# co-vary shipping and mining thresholds: keep constant fish

df_constant_fish001 <- subset(df_inc0.01, subset=!(k > 0.2))

p_const_fish <- ggplot(df_constant_fish001, aes(i, j)) +
  geom_raster(aes(fill=cost)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  scale_fill_distiller(palette = 'RdYlBu') +
  labs(x="Shipping threshold",
       y="Mining threshold",
 title = "Co-varying thresholds for the shipping and mining sectors")

ggplotly(p_const_fish)


# constant mine

df_constant_mine001 <- subset(df_inc0.01, subset=!(j > 0.02))

p_const_mine <- ggplot(df_constant_mine001, aes(i, k)) +
  geom_raster(aes(fill=cost)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  scale_fill_distiller(palette = "RdYlBu") +
  labs(x="Shipping threshold",
       y="Fishing threshold",
       title = "Co-varying thresholds for the shipping and fishing sectors")

ggplotly(p_const_mine)

# constant ship

df_constant_ship001 <- subset(df_inc0.01, subset=!(i > 0.01))
p_const_ship <- ggplot(df_constant_ship001, aes(j, k)) +
  geom_raster(aes(fill=cost)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.7,"cm")) +
  labs(x="Mining threshold",
       y="Fishing threshold",
       title = "Co-varying thresholds for the mining and fishing sectors")

ggplotly(p_const_ship)

##' plotting patchworks

# patchworks

# basic inputs and features

base_plot_4 <- PUs_plot + Contracts_plot + IMMA_plot + IBA_plot + EBSA_plot + Vents_plot + Plateaus_plot + Seamounts_plot +
  plot_layout(nrow = 4, ncol = 2) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold')) & 
  theme(legend.justification = "left")

base_plot_4

ggsave("pdfs/base_plot_4.jpeg", width = 20, height = 15, dpi = 1200)
ggsave("pdfs/base_plot_4.jpg", width = 20, height = 15, dpi = 1200)

# cost layers

costs_full_2 <- fishing_plot + shipping_plot + mining_plot +
  plot_layout(nrow = 2) +
  #plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold')) & 
  theme(legend.justification = "left")

ggsave("pdfs/costs_full.jpeg", width = 15, height = 15, dpi = 1200)
ggsave("pdfs/costs_full.jpg", width = 15, height = 15, dpi = 1200)

# sector-specific solutions

sol_rc_plots <- p_FishPlanA_Sol_Plot + rc_fish_plot + p_ShipPlanA_Sol_Plot + rc_ship_plot + p_MinePlanA_Sol_Plot + rc_mine_plot +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold')) & 
  theme(legend.justification = "left")

ggsave("pdfs/sol_rc_plots.jpeg", width = 15, height = 20, dpi = 1200)
ggsave("pdfs/sol_rc_plots.jpg", width = 15, height = 20, dpi = 1200)

# circular barplot

circbplots <- gg30 + gg70 +
  plot_layout(nrow = 2, guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(legend.position='bottom') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave("pdfs/circbarplots.jpeg", width = 10, height = 12, dpi = 1200)
ggsave("pdfs/circbarplots.jpeg", width = 10, height = 12, dpi = 1200)

# co-varying thresholds plots

p_thresholds <-  p_const_ship + p_const_mine + p_const_fish +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold')) & 
  theme(legend.justification = "left")

ggsave("pdfs/p_thresholds.jpg", width = 15, height = 10, dpi = 1200)
ggsave("pdfs/p_thresholds.jpeg", width = 15, height = 10, dpi = 1200)

# visual abstract: see https://www.urbandemographics.org/post/figures-map-layers-r/

rotate_data <- function(data, x_add = 0, y_add = 0) {
  
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x){ 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geometry = .$geometry * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}

rotate_data_geom <- function(data, x_add = 0, y_add = 0) {
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geom = .$geom * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}

# change parameters as wanted

x = -141.25
color = 'black'
bbox <- st_bbox(PUs)
world_valid <- sf::st_buffer(world, dist = 0) # otherwise not able to crop with bbox
world_valid <- st_crop(world_valid, bbox)
temp1 <- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = df_IBA%>% rotate_data(), aes(fill = Legend), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("IBA" = "red")) +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void() 

temp1b <- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = df_EBSA%>% rotate_data(), aes(fill = Legend), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("Arabian Oxygen Minimum Zone" = "dodgerblue4", 
                               "Sumatra Upwelling Zone" = "dodgerblue1",
                               "Agulhas Front" = "deepskyblue1",
                               "Mozambique Channel" = "lightskyblue2",
                               "Walters Shoals" = "paleturquoise",
                               "Saya de Malha Bank" = "palegreen",
                               "Coral Seamount Fracture Zone" = "greenyellow",
                               "Rusky Knoll" = "olivedrab3",
                               "Fools' Flat" = "yellow",
                               "East Broken Ridge" = "darkgoldenrod")) +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void() 

temp2 <- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = df_seamounts%>% rotate_data(x_add = .1), aes(fill = Legend), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("Seamount" = "blue")) +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void() 

temp3 <- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = df_IMMA%>% rotate_data(x_add = .1), aes(fill = Legend), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("IMMA" = "lightgreen")) +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void() 

temp4 <- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = df_plateaus%>% rotate_data(x_add = .1), aes(fill = Legend), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("Plateau" = "lightblue")) +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void() 

temp5<- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = PUs%>% rotate_data(), aes(color = log10(Fishing), fill = log10(Fishing))) +
  sc_fish +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void()

temp6<- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = PUs%>% rotate_data(), aes(color = log10(ShippingIntensity), fill = log10(ShippingIntensity))) +
  sc_ship +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void()

temp7<- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = PUs%>% rotate_data(), aes(color = log10(MiningCost), fill = log10(MiningCost))) +
  sc_mine +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void()

temp7<- ggplot() +
  
  geom_sf(data = PUs%>% rotate_data(), color = "white", show.legend = FALSE) +
  geom_sf(data = p_FishPlanA_Sol%>% rotate_data(), aes(fill = as.factor(solution_1)), size = 0.01) +
  sc_mine +
  geom_sf(data = world_valid%>% rotate_data(), fill = "grey20", size = 0.05, show.legend = FALSE) +
  theme_void()

##' THAT'S ALL; WELL DONE!


