#' ---
#' title: "Create a Spatial Plan"
#' author: "Jason D. Everett"
#' affiliation: "UQ/CSIRO/UNSW"
#' date: "Last compiled on `r format(Sys.time(), '%A %d %B %Y')`"
#' output: 
#'   html_document
#' ---
#'

##+ setup, include=FALSE
knitr::opts_chunk$set(warning=FALSE, cache=FALSE, message=FALSE)
# knitr::opts_chunk$set(collapse = TRUE, comment = "", warning = "off")

#' ## Overview
#' This code has been written to simplify the process for running a _prioritizr_ analysis on a given region. It is still a work in progress so feel free to submit pull requests with new features and code improvements.
#' 
#' The code depends on `sf`, `terra`, `tidyverse`, `rnaturalearth`, `prioritizr`, `stars`, `patchwork`.    
#' 
#' To use this code, you will need to download and expand `MME1DATA-Q1215/SpatialPlanning/Data.zip` to the directory `GitHub/SpatialPlanning/Data/`. Note that the download is only 2GB, but the expanded data is 35 GB in size. If you need help subsetting the data to your region due to memory or HD space constraints, contact Jason.
#' 
#' ## Preliminaries 
source("SpatPlan_Extras.R") # Load the extras, including functions and libraries
# source("SpatPlan_Process_AquaMaps.R") # This script reprocesses AquaMaps. WARNING: Lots of time (10s hrs) and memory (10s GB)
# source("SpatPlan_Process_MPAs.R") # Only run if you need to reprocess the MPA data. You will need the WDPA files

if (!file.exists("Data")) {
  stop("The Data folder does not exist at SpatialPlanning/Data. Please download from the RDM and then try again.")
}

reprocess<- TRUE # Do we want to reprocess the PUs/AquaMaps data or use saved versions

#' ### Set user parameters 

# You can set a region if it is defined in fSpatPlan_Get_PlanningUnits.
# If you add regions, please submit a PR so others can benefit.
#Region <- "Australia"
#save_name <- "Australia" # Name used in the saving

# You can also define a region with square boundaries
#Region <- c(xmin = 140, xmax = 178, ymin = -40, ymax = -20)
# save_name <- "SEAustralia" # Name used in the saving

Region <- c(xmin = 18, xmax = 130, ymin = -39, ymax = 34)
Region <- "IO"
save_name <- "IO"
  
PU_size <- 1000 # km2
# Region <- "Global"  
# PU_size <- 2620 # km2 (0.5 deg at equator)
# PU_size <- 10000

Shape <- "Hexagon" # "Shape of PUs
MinDepth <- 0
MaxDepth <- 200
#CO <- 0.5

#' For inverse area targets
#minTarget = 0.2
#maxTarget = 0.8

#' Choose CRS for analysis
cCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # Robinson
# cCRS <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # Robinson: Pacific-centred

#' Get the land boundaries to remove overlap
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(cCRS)

#' ## Boundaries for PUs 
Bndry <- fSpatPlan_Get_Boundary(Region, cCRS)

#' ## Create Planning Units 
if(reprocess){
  PUs <- fSpatPlan_Get_PlanningUnits(Bndry, world, PU_size, Shape)
  saveRDS(PUs, file.path("Output", paste(save_name, "PU", paste0(PU_size,"km2"), "Output.rds", sep = "_")))
} else {
  PUs <- read_rds(file.path("Output", paste(save_name, "PU", paste0(PU_size,"km2"), "Output.rds", sep = "_")))
}


#' ## Get the features 

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

#' ## Get locked in areas 
LockedIn <- fSpatPlan_Get_MPAs(PUs, cCRS)

#' ## Get cost 
Cost <- fSpatPlan_Get_Cost(PUs, cCRS) ##make sure cost function is set to fetch .grd in proper directory
                                      ##cannot open .grd despite changing pathway in fx
#' ## Set up targets 

#' ### Identical fixed targets 
Targets <- 0.2


#' ## Do some plotting 
(ggPU <- fSpatPlan_PlotPUs(PUs, world)) # Plot Planning Units
(ggMPA <- fSpatPlan_PlotMPAs(LockedIn, world)) # Plot Locked in areas
(ggFishingCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

ggplot(data = Bndry) + geom_sf()

####################### get cost polygons (Lea) --- crossing cost shp with PUs
##### general get_polyg fx --- crossing shp with PUs

fSpatPlan_Get_Polyg <- function(filen, PUs){
  
  mine_crs <- st_crs(PUs)
  Mining1 <- st_read(filen) %>% ##should be 'Mine1'?
    st_transform(mine_crs) ### ensure same crs for shp and PUs?
  
  Mining1PUs <- PUs %>% 
    st_centroid() %>% #Second, get all the PUs with < 50 % area on land #comment not applicable here?
    st_within(Mining1, sparse = FALSE) %>% ##'Mine1'?
    rowSums() %>%
    as.logical()
  
  # Mining1PUs <- PUs %>% 
  #   mutate(locked_in = Mining1PUs) ##assign new column to PUs (call it name of shp)
  # 
  return(Mining1PUs)
}

###### loading files

#Cfc resources (need to assign cost?)
filen <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp" 
#Nodule resources (need to assign cost?)
filen2 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp"
#Exploration areas (need to lockout)
filen3 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/01_pmn_exploration_areas.shp"
filen4 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/02_pms_exploration_areas.shp"
filen5 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/03_cfc_exploration_areas.shp"
#Reserved areas (need to lockout)
filen6 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/04_pmn_reserved_areas.shp"
filen7 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/05_cfc_reserved_areas.shp"

#####crossing shp & PUs

## Crust areas crossed

filen <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp" 
out <- fSpatPlan_Get_Polyg(filen, PUs)  # function returns out as a logical vector
PUs <- PUs %>% 
  mutate(CFC = out) # Add logical vector to variable called Mine1

## Nodules crossed

filen2 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp"
out2 <- fSpatPlan_Get_Polyg(filen2, PUs)  
PUs <- PUs %>% 
  mutate(PMN = out2)

## Exploration areas crossed
filen3 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/01_pmn_exploration_areas.shp"
out3 <- fSpatPlan_Get_Polyg(filen3, PUs)  
PUs <- PUs %>% 
  mutate(ExplPMN = out3)

filen4 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/02_pms_exploration_areas.shp"
out4 <- fSpatPlan_Get_Polyg(filen4, PUs)  
PUs <- PUs %>% 
  mutate(ExplPMS = out4)

filen5 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/03_cfc_exploration_areas.shp"
out5 <- fSpatPlan_Get_Polyg(filen5, PUs)  
PUs <- PUs %>% 
  mutate(ExplCFC = out5)

## Reserved areas crossed

filen6 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/04_pmn_reserved_areas.shp"
out6 <- fSpatPlan_Get_Polyg(filen6, PUs)  
PUs <- PUs %>% 
  mutate(ResPMN = out6)

filen7 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Costs/Mining_contractors/ContractorAreas/05_cfc_reserved_areas.shp"
out7 <- fSpatPlan_Get_Polyg(filen7, PUs)  
PUs <- PUs %>% 
  mutate(ResCFC = out7)

####################### get biodiv_features polygons (Lea) --- crossing biodiv shp with PUs
###### loading files

##Important Bird Areas
filen8 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp"
##Important Marine Mammal Areas
filen9 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/IMMAs/iucn-imma.shp"
##EBSAs
filen10 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp"
filen11 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp"
filen12 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp"
filen13 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp"
filen14 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SSIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp"
filen15 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp"
filen16 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp"
filen17 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp"
filen18 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp"
filen19 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp"
filen20 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp"

#####crossing shp & PUs

### Bird Areas crossed

filen8 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/IBAs/MarineIBA_IndianOcean.shp"
out8 <- fSpatPlan_Get_Polyg(filen8, PUs)  
PUs <- PUs %>% 
  mutate(IBA = out8)

### Marine Mammal areas crossed

filen9 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/IMMAs/iucn-imma.shp"
out9 <- fSpatPlan_Get_Polyg(filen9, PUs)  
PUs <- PUs %>% 
  mutate(IMMA = out9)

### EBSAs crossed

filen10 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/NEIO_9_EBSA/NEIO_9_EBSA.shp"
out10 <- fSpatPlan_Get_Polyg(filen10, PUs)  
PUs <- PUs %>% 
  mutate(NEIO_09 = out10)

filen11 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/NWIO_14_EBSA/NWIO_14_EBSA.shp"
out11 <- fSpatPlan_Get_Polyg(filen11, PUs)  
PUs <- PUs %>% 
  mutate(NWIO_14 = out11)

filen12 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_11_EBSA-GIS shapefile/SIO_11_EBSA.shp"
out12 <- fSpatPlan_Get_Polyg(filen12, PUs)  
PUs <- PUs %>% 
  mutate(SIO_11 = out12)

filen13 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_19_EBSA-GIS%20shapefile/SIO_19_EBSA.shp"
out13 <- fSpatPlan_Get_Polyg(filen13, PUs)  
PUs <- PUs %>% 
  mutate(SIO_19 = out13)

filen14 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SSIO_22_EBSA-GIS shapefile/SIO_22_EBSA.shp"
out14 <- fSpatPlan_Get_Polyg(filen14, PUs)  
PUs <- PUs %>% 
  mutate(SIO_22 = out14)

filen15 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_23_EBSA-GIS shapefile/SIO_23_EBSA.shp"
out15 <- fSpatPlan_Get_Polyg(filen15, PUs)  
PUs <- PUs %>% 
  mutate(SIO_23 = out15)

filen16 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_30_EBSA-GIS shapefile/SIO_30_EBSA.shp"
out16 <- fSpatPlan_Get_Polyg(filen16, PUs)  
PUs <- PUs %>% 
  mutate(SIO_30 = out16)

filen17 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_32/SIO_32_EBSA.shp"
out17 <- fSpatPlan_Get_Polyg(filen17, PUs)  
PUs <- PUs %>% 
  mutate(SIO_32 = out17)

filen18 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_35_EBSA-GIS shapefile/SIO_35_EBSA.shp"
out18 <- fSpatPlan_Get_Polyg(filen18, PUs)  
PUs <- PUs %>% 
  mutate(SIO_35 = out18)

filen19 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_36_EBSA-GIS shapefile/SIO_36_EBSA.shp"
out19 <- fSpatPlan_Get_Polyg(filen19, PUs)  
PUs <- PUs %>% 
  mutate(SIO_36 = out19)

filen20 <- "~/Documents/MscThesis/GitHub/SpatialPlanning_copy/Shp/Biodiv_features/EBSAs/SIO_37_EBSA-GIS shapefile/SIO_37_EBSA.shp"
out20 <- fSpatPlan_Get_Polyg(filen20, PUs)  
PUs <- PUs %>% 
  mutate(SIO_37 = out20)
