#' ---
#' title: "Create a Spatial Plan"
#' author: "Jason D. Everett, modified by Lea Fourchault"
#' affiliation: "UQ/CSIRO/UNSW/Tropimundo"
#' date: "Last compiled on March 8th, 2022"
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
#' Script below modified by Lea Fourchault for Cross-sectoral conservation of the Indian Ocean project
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
    st_centroid() %>% #Second, get all the PUs with < 50 % area on land #?
    st_within(Polyg1, sparse = FALSE) %>% 
    rowSums() %>%
    as.logical()
  
  # Polyg1PUs <- PUs %>% 
  #   mutate(locked_in = Polyg1PUs) ##assign new column to PUs (call it name of shp)
  # 
  return(Polyg1PUs)
}

#' ### Set user parameters 

# You can set a region if it is defined in fSpatPlan_Get_PlanningUnits.
# If you add regions, please submit a PR so others can benefit.
# You can also define a region with square boundaries

#' Define region and boundary

cCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(cCRS)
#Region <- c(xmin = 18, xmax = 140, ymin = -39, ymax = 34) 
Region <- c(xmin = 20, xmax = 120, ymin = -60, ymax = 20) 
# Indian Ocean limits: CIA World Factbook (but see: According to the International Hydrographic Organization (IHO) and the United Nations Oceans Atlas, the area from 40 degrees S latitude to 60 degrees S latitude is included in the Indian Ocean. The area that encircles the globe from 60 degrees S latitude to the coast of Antarctica is called The Great Southern Ocean. The Indian Ocean's width extends from 45 degrees E longitude to 110 degrees East longitude.) 
Bndry <- fSpatPlan_Get_Boundary(Region, cCRS) # boundary

#ggplot(data = Bndry) + geom_sf()

#' Intersect bndry and high seas

ABNJ <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/World_High_Seas_v1/High_Seas_v1.shp") 
ABNJRobin <- st_transform(ABNJ, cCRS) # high seas in robin proj, does not work when using pipe
ABNJRobBuff  <- st_buffer(ABNJRobin, 0) # need to create buffer of 0 to avoid Error in CPL_geos_op2(op, x, y) : Evaluation error: TopologyException: Input geom 0 is invalid: Self-intersection at or near point -> overwrite offending geometries
BndryABNJ <- st_intersection(ABNJRobBuff, Bndry)

#' Create PUs

PU_size <- 1000 # in km2, but check resolution for this project?
Shape <- "Hexagon" # shape of PUs
MinDepth <- 0
MaxDepth <- 200 # but check for deep-sea?
PUs <- fSpatPlan_Get_PlanningUnits(BndryABNJ, world, PU_size, Shape) # modified PUs with ABNJ only
#PUs <- fSpatPlan_Get_PlanningUnits(Bndry, world, PU_size, Shape) # original function

#' Plot PUs

(ggPU <- fSpatPlan_PlotPUs(PUs, world)) 

#' ## Get the conservation features 
#' IBAs

#loading file 
IBA <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp"

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
  mutate(NEIO_09 = PUxNEIO_09)

NWIO_14 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp"
PUxNWIO_14 <- fSpatPlan_Get_Polyg(NWIO_14, PUs)  
PUs <- PUs %>% 
  mutate(NWIO_14 = PUxNWIO_14)

SIO_11 <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp"
PUxSIO_11 <- fSpatPlan_Get_Polyg(SIO_11, PUs)  
PUs <- PUs %>% 
  mutate(SIO_11 = PUxSIO_11)

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
PUx_SIO37 <- fSpatPlan_Get_Polyg(SIO_37 , PUs)  
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

#' Plotting features

world <- ne_countries(scale = "medium", returnclass = "sf")
#IndianO <- c(xmin = 17.47941, ymin = -38.7402, xmax = 129.552, ymax = 33.82872)
IndianO <- c(xmin = 20, xmax = 120, ymin = -60, ymax = 20)

IMMA <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IMMAs/iucn-imma.shp")
IBA <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp")
NEIO_09<- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp")
NWIO_14 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp")
SIO_11 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp")
SIO_19 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp")
SIO_22 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp")
SIO_23 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp")
SIO_30 <- st_read  ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp")
SIO_32 <- st_read  ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp")
SIO_35<- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp")
SIO_36 <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp")
SIO_37  <-st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp")
Seamounts <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Seamounts.shp")
Plateaus <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Biodiv_features/BlueHabs/Plateaus.shp")

ggplot() + 
  geom_sf(data = world) +
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
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE)+
  scale_color_manual(values = c("A" = "grey", "B" = "black", "C" = "blue", "D" = "green", "E" = "purple", "F" = "red", "G" = "brown"), 
                     labels = c("EBSAs", "Seamounts", "Plateaus", "IBAs", "IMMAs", "Active vents", "Inactive Vents"),
                     name = "Legend") 

#' Longhurst provinces

longh <- st_read(file.path("Data","LonghurstProvinces","Longhurst_world_v4_2010.shp")) %>% 
  st_make_valid()
longhRob <- st_transform(longh, cCRS)

nr <- st_nearest_feature(PUs, longhRob)

PUs <- PUs %>% 
  mutate(ProvCode = longh$ProvCode[nr],
         ProvDescr = longh$ProvDescr[nr])

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE) +
  geom_sf(data = PUs, aes(fill = ProvCode), colour = "grey50", size = 0.05, alpha = 0.7, show.legend = TRUE) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE)


#' ### Aquamaps 

#' Get aquamaps data for our species and region 
#'TODO Check the overlay of raw v processed species distribution
# if(reprocess){
#   aqua_sf <- fSpatPlan_Get_AquaMaps(PUs, cCRS, MinDepth, MaxDepth, CutOff = CO)
#   
#   saveRDS(aqua_sf, file.path("Output", 
#                              paste(save_name, "PU", paste0(PU_size,"km2"), 
#                                    "AquaMaps_Output.rds", sep = "_"))) # Save rds so you don't have to reprocess everytime.
# } else {
#   aqua_sf <- read_rds(file.path("Output", 
#                                 paste(save_name, "PU", paste0(PU_size,"km2"), 
#                                       "AquaMaps_Output.rds", sep = "_")))
# }

#' ### IUCN 
# TODO
# IUCN <- fSpatPlan_Get_IUCN(PUs, cCRS)

#' ### MICO 
# TODO
# MICO <- fSpatPlan_Get_MICO(PUs, cCRS)

#' ### Seafloor 
# TODO
# SeaFloor <- fSpatPlan_Get_SeaFloor(PUs, cCRS)

#' ### Geomorphic features 
# TODO
# GeoMorph <- fSpatPlan_Get_GeoMorph(PUs, cCRS)

#' ## Intersection of Features with Longhurst 
# TODO

#' ## Setting targets for conservation features

#Cons_feats <- tibble(Features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), Targets = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.7))

#' ## Get locked in areas 
#LockedIn <- fSpatPlan_Get_MPAs(PUs, cCRS)
#(ggMPA <- fSpatPlan_PlotMPAs(LockedIn, world)) # Plot Locked in areas

#' ## Get locked out areas
#' polymetallic nodules / polymetallic sulfides: exploration & reserved areas

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

#' plotting locked out areas

world <- ne_countries(scale = "medium", returnclass = "sf")
IndianO <- c(xmin = 20, xmax = 120, ymin = -60, ymax = 20)
Expl_PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/01_pmn_exploration_areas.shp")
Expl_PMS <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/02_pms_exploration_areas.shp")
Res_PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_contractors/ContractorAreas/04_pmn_reserved_areas.shp")

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = Expl_PMN, aes(color="A"), show.legend = "line") + 
  geom_sf(data = Expl_PMS, aes(color="B"), show.legend = "line") +
  geom_sf(data = Res_PMN, aes(color="C"), show.legend = "line") +
  scale_color_manual(values = c("A" = "red", "B" = "yellow", "C" = "black"),
                     labels = c("Expl_PMN", "Expl_PMS", "Res_PMN"),
                     name = "Legend") +
  ggtitle("Exploration & Reserved areas") + 
  coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

#' ## Get cost 
#' Fishing cost
#FishingCost <- fSpatPlan_Get_Cost(PUs, cCRS) ##make sure cost function is set to fetch .grd in proper directory
##cannot open .grd despite changing pathway in fx

FishingCost <- terra::rast("Data/Cost/Cost_Raster_Sum.grd") %>% 
  terra::as.polygons(trunc = FALSE, dissolve = FALSE, na.rm=FALSE) %>% # Convert to polygon data
  st_as_sf() %>% # Convert to sf
  st_transform(cCRS) %>% # transform to robinson
  st_interpolate_aw(PUs, extensive = FALSE) %>% ## intersect with PUs
  rename(FishingCost = layer)

PUs <- PUs %>% 
  mutate(Fishing = FishingCost$FishingCost)

LogFishing <- log10(FishingCost$FishingCost)
hist(LogFishing)
PUs <- PUs %>% 
  mutate(LogFishing = LogFishing)


#' Plotting log of fishing cost

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = PUs, aes(color=LogFishing)) + 
  ggtitle("Fishing Cost (Log10)") + 
  coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

#(ggCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

#' CFC value
#'script to get area modified from https://rpubs.com/rural_gis/255550

CFC <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp" 
PUxCFC <- fSpatPlan_Get_Polyg(CFC, PUs)  # function returns out as a logical vector
PUs <- PUs %>% 
  mutate(CFC = PUxCFC)

#load CFC areas polygon shapefile 
CFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")
CFCRobin <- st_transform(CFC, cCRS) #robin proj

#run the intersect function, converting the output to a tibble in the process
#int <- as_tibble(st_intersection(CFCRobin, BndryABNJ))
int <- st_intersection(CFCRobin, BndryABNJ)

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int$CFCarea <- st_area(int)

#plot the layers to visually check result of intersect
plot (CFCRobin$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaCFCtotal <- int %>%
  summarise(areaCFCtot = sum(CFCarea)) #2.256596e+12 [m^2], seems a lot? Would be 2256596 sqkm, I don't even have that many PUs

#' area of CFC in PUs < 2753 sqkm = 2.753e+09 sqm because
#' sum(PUxCFC) 
#' [1] 2753

#change data type of areaArable field to numeric (to remove m^2 suffix)
AreaCFCtotal$areaCFCtot <- as.numeric(AreaCFCtotal$areaCFCtot)

#' CFC value for Indian Ocean =  325 m2/g (mean specific-surface area, in Hein et al., 2000; found in Hein et al., 2013) * 2.256596e+12 [m^2]

#' PMN value

#'script to get area modified from https://rpubs.com/rural_gis/255550

#load county areas polygon shapefile 

PMN <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp"
PUxPMN <- fSpatPlan_Get_Polyg(PMN, PUs)  
PUs <- PUs %>% 
  mutate(PMN = PUxPMN)

PMN <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")
PMNRobin <- st_transform(PMN, cCRS) #robin proj

#run the intersect function, converting the output to a tibble in the process
int2 <- st_intersection(PMNRobin, BndryABNJ) #need this to avoid Error in UseMethod("st_area") : 
#int2 <- as_tibble(st_intersection(PMNRobin, BndryABNJ)) #no applicable method for 'st_area' applied to an object of class "c('tbl_df', 'tbl', 'data.frame')"

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int2$PMNarea <- st_area(int2)

#plot the layers to visually check result of intersect
plot (PMNRobin$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int2, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaPMNtotal <- int2 %>%
  summarise(areaPMNtot = sum(PMNarea)) #1.933646e+12 [m^2], also seems a lot?

#'area of PMN in PUs < 2075 sqkm because
#'sum(PUxPMN)
#'[1] 2075

#change data type of areaArable field to numeric (to remove m^2 suffix)
AreaPMNtotal$areaPMNtot <- as.numeric(AreaPMNtotal$areaPMNtot)

#' PMN value for Indian Ocean = 320 USD/t (CRU, 2019) * 0.0056 t/sqm (Sharma et al., 2011) * 1.933646e+12 sqm (total PMN area in Indian Ocean, from shapefiles from ISA & chunk above) = USD 3.46509363e+12 for all PMN in the Indian Ocean (seems too much)

#' Plotting mining resources

world <- ne_countries(scale = "medium", returnclass = "sf")
#IndianO <- c(xmin = 17.47941, ymin = -38.7402, xmax = 129.552, ymax = 33.82872)
CFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")
PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")
ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = CFC, aes(color="A"), show.legend = "line") + 
  geom_sf(data = PMN, aes(color="B"), show.legend = "line") + 
  scale_color_manual(values = c("A" = "grey", "B" = "brown"),
                     labels = c("CFC", "PMN"),
                     name = "Legend") +
  ggtitle("Polymetallic Nodules and Cobalt-rich Ferromanganese Crust areas") + 
  coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)


#' ## Prioritizr problems

PUs$Fishing[PUs$Fishing <  0.0001] <-  0.0001 #to avoid Warning in presolve_check.OptimizationProblem(compile(x)) :
# planning units with very high (> 1e+6) or very low (< 1e-6) non-zero cost values note this may be a false positive

#' IBA & Fishing cost

p_IBA_Fishing <- problem(PUs, features = "IBA", cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(0.9) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE) 

p_IBA_Fishing_Sol <- solve(p_IBA_Fishing)

p_IBA_Fishing_Sol_Plot <- ggplot() +
  ggtitle("IBA vs Fishing, Target = 0.9") +
  geom_sf(data = world, colour ="grey70") +
  geom_sf(data = p_IBA_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_IBA_Fishing_Sol_Plot #looks correct 

#' Active Vents & Fishing cost & target

p_Active_Vents_Fishing <- problem(PUs, features = "Active_Vents", cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(0.7) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_Active_Vents_Fishing_Sol <- solve(p_Active_Vents_Fishing)

p_Active_Vents_Fishing_Sol_Plot <- ggplot() +
  ggtitle("Active Vents vs Fishing, Targets = 0.7") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_Active_Vents_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_Active_Vents_Fishing_Sol_Plot # also looks correct

#' All features & Fishing cost & targets

p_ConsFeat_Fishing <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.3, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.7)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_ConsFeat_Fishing_Sol <- solve(p_ConsFeat_Fishing)
#plot(p_ConsFeat_Fishing_Sol, main = "solution")
#lines(p_ConsFeat_Fishing_Sol[p_ConsFeat_Fishing_Sol$solution_1 > 0.5, ], col = "darkgreen", lwd = 2)

p_ConsFeat_Fishing_Sol_Plot <- ggplot() +
  ggtitle("Cons Feats vs Fishing, 0.3<Targets<0.7") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_ConsFeat_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ConsFeat_Fishing_Sol_Plot # shows that fishing cost is reduced between 40°-60°S, but maybe would be good to reduce target for IBAs/IMMAs or assign Longhurst to distribute better across latitudes

#' All features & Fishing cost & targets & locked-out exploration + reserved areas

PUs["locked_out"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN # create locked_out column (logical). If =<1 is TRUE, then locked_out is TRUE # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

#PUs["locked_out"] <- paste(PUs$ExplPMN, PUs$ExplPMS, PUs$ResPMN)
#PUs %>% mutate(sum = rowSums(across(where(is.logical))))rowSums(Y) > 0
#PUs$LockedOut <- ifelse(PUs$ExplPMN == "TRUE" | PUs$ExplPMS == "TRUE" | PUs$ResPMN == "TRUE")

p_ConsFeat_Fishing_LockedOut <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.3, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.7)) %>%
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_ConsFeat_Fishing_LockedOut_Sol <- solve(p_ConsFeat_Fishing_LockedOut)

p_ConsFeat_Fishing_LockedOut_Sol_Plot <- ggplot() +
  ggtitle("Cons Feats vs Fishing, Mining Contractors Out, 0.3<Targets<0.7") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_ConsFeat_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ConsFeat_Fishing_LockedOut_Sol_Plot


#' All features & Fishing cost & targets & locked-out exploration + reserved areas, change targets to check solver works

length(PUs$Plateaus[PUs$Plateaus==TRUE]) # check count of PUs that have plateaus inside

p_check_targets <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.001, 0.001, 0.03, 0.03, 0.03, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.99, 0.99, 0.99, 0.99)) %>%
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_check_targets_sol <- solve(p_check_targets)

p_check_targets_sol_plot <- ggplot() +
  ggtitle("0.001<Targets<0.9") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_ConsFeat_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_check_targets_sol_plot