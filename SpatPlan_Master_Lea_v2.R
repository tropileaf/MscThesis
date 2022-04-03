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

fSpatPlan_Get_Polyg_Binary <- function(filen, PUs){
  
  polyg_crs <- st_crs(PUs)
  Polyg1 <- st_read(filen) %>% 
    st_transform(polyg_crs) 
  
  Polyg1PUs <- PUs %>% 
    st_centroid() %>% #Second, get all the PUs with < 50 % area on land #?
    st_within(Polyg1, sparse = FALSE) %>% 
    rowSums()
  
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
  st_transform(cCRS) #

#Region <- c(xmin = 18, xmax = 140, ymin = -39, ymax = 34) 
Region <- c(xmin = 18, xmax = 120, ymin = -45, ymax = 34) 
# Indian Ocean limits: CIA World Factbook (but see: According to the International Hydrographic Organization (IHO) and the United Nations Oceans Atlas, the area from 40 degrees S latitude to 60 degrees S latitude is included in the Indian Ocean. The area that encircles the globe from 60 degrees S latitude to the coast of Antarctica is called The Great Southern Ocean. The Indian Ocean's width extends from 45 degrees E longitude to 110 degrees East longitude.) 
Bndry <- fSpatPlan_Get_Boundary(Region, cCRS) # boundary

#ggplot(data = BndryABNJ) + geom_sf()

#' Intersect bndry and high seas

ABNJ <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/World_High_Seas_v1/High_Seas_v1.shp") 
ABNJRobin <- st_transform(ABNJ, cCRS) # high seas in robin proj, does not work when using pipe
ABNJRobBuff  <- st_buffer(ABNJRobin, 0) # need to create buffer of 0 to avoid Error in CPL_geos_op2(op, x, y) : Evaluation error: TopologyException: Input geom 0 is invalid: Self-intersection at or near point -> overwrite offending geometries
#BndryABNJ <- st_transform(ABNJRobBuff, cCRS)
BndryABNJ <- st_intersection(ABNJRobBuff, Bndry)

#BndryABNJ <- st_buffer(ABNJ, 0)%>%
#  st_transform(ABNJ, cCRS) #Error in s2_geography_from_wkb(x, oriented = oriented, check = check) : 
#Evaluation error: Found 1 feature with invalid spherical geometry.
#[1] Loop 36 is not valid: Edge 0 crosses edge 1108.


#' Create PUs

PU_size <- 1000 # in km2, but check resolution for this project?
Shape <- "Hexagon" # shape of PUs
#MinDepth <- 0
#MaxDepth <- 200 # but check for deep-sea?
PUs <- fSpatPlan_Get_PlanningUnits(BndryABNJ, world, PU_size, Shape) # modified PUs with ABNJ only
#PUsAll <- fSpatPlan_Get_PlanningUnits(Bndry, world, PU_size, Shape) # original function

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

#' Plotting features

#world <- ne_countries(scale = "medium", returnclass = "sf")
#IndianO <- c(xmin = 17.47941, ymin = -38.7402, xmax = 129.552, ymax = 33.82872)
#IndianO <- c(xmin = 18, xmax = 110, ymin = -45, ymax = 34) %>% # Convert to polygon data
# st_as_sf(coords = c("long", "lat"), crs = 4326) %>% # Convert to sf
#st_transform(cCRS)


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
  #coord_sf(xlim = c(18, 130), ylim = c(34, -60), expand = FALSE)+
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)+
  scale_color_manual(values = c("A" = "grey", "B" = "black", "C" = "blue", "D" = "green", "E" = "purple", "F" = "red", "G" = "brown"), 
                     labels = c("EBSAs", "Seamounts", "Plateaus", "IBAs", "IMMAs", "Active vents", "Inactive Vents"),
                     name = "Legend") 

#ggplot() + #just to check coord_sf
# geom_sf(data = PUs, fill = "lightsteelblue2", color = "grey64", size = 0.05, show.legend = FALSE) +
#geom_sf(data = world, colour = "grey20", fill = "grey20", alpha = 0.9, size = 0.1, show.legend = FALSE) +
#coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
#theme_bw() +
#labs(subtitle = "Planning Units")

#cCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
# st_transform(cCRS)
#IndianO <- c(xmin = 18, xmax = 130, ymin = -60, ymax = 22)

#ggplot() + 
# geom_sf(data = world) +
#coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE)
#coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)


#' Longhurst provinces

longh <- st_read(file.path("Data","LonghurstProvinces","Longhurst_world_v4_2010.shp")) %>% 
  st_make_valid()
longhRob <- st_transform(longh, cCRS)
longhRob %>% dplyr::select(ProvCode) %>% head(2) # geometry is sticky, as should be
nr <- st_nearest_feature(PUs, longhRob)
PUs <- PUs %>% 
  mutate(ProvCode = longh$ProvCode[nr],
         ProvDescr = longh$ProvDescr[nr]) # add columns with ProvCode and ProvDescr to PUs df. ProvCode gives the province code for each PU, ie., character such as "ANTA"
#plot

ggplot() + 
  geom_sf(data = world, color = "grey20", fill = "grey20", size = 0.1, show.legend = FALSE) +
  geom_sf(data = PUs, aes(fill = ProvCode), colour = "grey50", size = 0.05, alpha = 0.7, show.legend = TRUE) +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE)
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

# provinces intersection with PUs
# creating a new column for logical intersection between each longhurst province type and PUs
# provinces = WARM, SSTC, SPSG, SANT, NPSW, ISSG, ANTA

WARMprov <- PUs[PUs$ProvCode == "WARM",]
PUxWARM <- st_contains(PUs,WARMprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(WARM = PUxWARM) # adding column of single longhurst province to PUs
length(PUs$WARM[PUs$WARM==TRUE]) # check area size: 2288 PUs = 2288*1000sqkm

SSTCprov <- PUs[PUs$ProvCode == "SSTC",]
PUxSSTC <- st_contains(PUs,SSTCprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(SSTC = PUxSSTC)
length(PUs$SSTC[PUs$SSTC==TRUE]) # area = 3729 PUs

SPSGprov <- PUs[PUs$ProvCode == "SPSG",]
PUxSPSG <- st_contains(PUs,SPSGprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(SPSG = PUxSPSG)
length(PUs$SPSG[PUs$SPSG==TRUE]) # 12295 PUs

SANTprov <- PUs[PUs$ProvCode == "SANT",]
PUxSANT <- st_contains(PUs,SANTprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(SANT = PUxSANT)
length(PUs$SANT[PUs$SANT==TRUE]) # 6140 PUs

NPSWprov <- PUs[PUs$ProvCode == "NPSW",]
PUxNPSW <- st_contains(PUs,NPSWprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(NPSW = PUxNPSW)
length(PUs$NPSW[PUs$NPSW==TRUE]) # 1238 PUs

NEWZprov <- PUs[PUs$ProvCode == "NEWZ",]
PUxNEWZ <- st_contains(PUs,NEWZprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(NEWZ = PUxNEWZ)
length(PUs$NEWZ[PUs$NEWZ==TRUE])  # 881 PUs

ISSGprov <- PUs[PUs$ProvCode == "ISSG",]
PUxISSG <- st_contains(PUs,ISSGprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(ISSG = PUxISSG)
length(PUs$ISSG[PUs$ISSG==TRUE]) # 1617 PUs

ANTAprov <- PUs[PUs$ProvCode == "ANTA",]
PUxANTA <- st_contains(PUs,ANTAprov, sparse = FALSE) %>%
  rowSums() %>% 
  as.logical()
PUs <- PUs %>% 
  mutate(ANTA = PUxANTA)
length(PUs$ANTA[PUs$ANTA==TRUE]) # 7006 PUs, also correct based on map

#ANTA <- PUs["ANTA", "ProvCode"]
#WARMprov <- PUs[PUs$ProvCode == "WARM",]
#PUxWARM <- fSpatPlan_Get_Polyg(ANTA, PUs)
#nr <- st_nearest_feature(PUs, ANTA)
#longh <- st_read(file.path("Data","LonghurstProvinces","Longhurst_world_v4_2010.shp"))
#ANTAcode <- filter(longh, ProvCode == "WARM")
#ANTAprov <- filter(PUs, ProvCode == "ANTA")
#ANTAprov <- PUs[PUs$ProvCode %in% c("E06000001", "E06000020"), ]
#ANTAsf <- st_as_sf(ANTAprov, coords = c(x = "Longitude", y = "Latitude"), crs = 4326) %>% 
# st_transform(cCRS)
#PUxANTA <- fSpatPlan_Get_Polyg(ANTAsf, PUs)  
#PUs <- PUs %>% 
# mutate(ANTA = ANTA)

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
  geom_sf(data = world) +
  geom_sf(data = Expl_PMN, aes(color="A"), show.legend = "line") + 
  geom_sf(data = Expl_PMS, aes(color="B"), show.legend = "line") +
  geom_sf(data = Res_PMN, aes(color="C"), show.legend = "line") +
  scale_color_manual(values = c("A" = "red", "B" = "yellow", "C" = "black"),
                     labels = c("Expl_PMN", "Expl_PMS", "Res_PMN"),
                     name = "Legend") +
  ggtitle("Exploration & Reserved areas") + 
  #coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

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


hist(FishingCost$FishingCost)
hist(PUs$Fishing) # same as FishingCost$FishingCost, correct

max(FishingCost$FishingCost) # USD 5769.704 per sqkm per yr #5077.817 when using 29480 PUs (cut off IO at 45°S)
max(PUs$Fishing) #USD 5769.704 per sqkm per yr, same, correct
mean(PUs$Fishing) # USD 21.0153 per sqkm per yr #28.34043
min(PUs$Fishing) #USD 0 per sqkm per yr  #0.2267226
median(PUs$Fishing) # USD 8.480085 per sqkm per yr #11.04459

Log10_Fishing <- log10(FishingCost$FishingCost)
hist(Log10_Fishing)
PUs <- PUs %>% 
  mutate(Log10_Fishing = Log10_Fishing)


#' Plotting log of fishing cost

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = PUs, aes(color=Log10_Fishing)) + 
  ggtitle("Fishing Cost (Log10)") + 
  #coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$Fishing)) +
  geom_sf(data = world) +
  ggtitle("Fishing Cost") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)
#coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

#(ggCost <- fSpatPlan_PlotCost(Cost, world)) # Plot cost

#' CFC value
#'script to get area modified from https://rpubs.com/rural_gis/255550

CFC <- "~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp" 
PUxCFC <- fSpatPlan_Get_Polyg(CFC, PUs)  # function returns out as a logical vector
PUs <- PUs %>% 
  mutate(CFC = PUxCFC)

#load CFC areas polygon shapefile 
CFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")%>%
  st_transform(cCRS) #robin proj

#run the intersect function, converting the output to a tibble in the process
#int <- as_tibble(st_intersection(CFCRobin, BndryABNJ))
int <- st_intersection(CFC, BndryABNJ)

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int$CFCarea <- st_area(int)

#plot the layers to visually check result of intersect
plot (CFC$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaCFCtotal <- int %>%
  summarise(areaCFCtot = sum(CFCarea)) #2.256596e+12 [m^2], seems a lot? Would be 2256596 sqkm

#' area of CFC in PUs < 2753*1000 sqkm = 2.753e+12 sqm because
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


PMN <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")%>%
  st_transform(cCRS) #robin proj

#run the intersect function, converting the output to a tibble in the process
int2 <- st_intersection(PMN, BndryABNJ) #need this to avoid Error in UseMethod("st_area") : 
#int2 <- as_tibble(st_intersection(PMNRobin, BndryABNJ)) #no applicable method for 'st_area' applied to an object of class "c('tbl_df', 'tbl', 'data.frame')"

#add in an area count column to the tibble (area of each CFC poly in intersect layer)
int2$PMNarea <- st_area(int2)

#plot the layers to visually check result of intersect
plot (PMN$geometry, col='green')
plot(BndryABNJ$geometry, add=T)
plot(int2, col='red', add=T) #correct

#group data by county area and calculate the total arable land area per county
#output as new tibble
AreaPMNtotal <- int2 %>%
  summarise(areaPMNtot = sum(PMNarea)) #1.933646e+12 [m^2], also seems a lot? but should be okay

#'area of PMN in PUs < 2075*1000 sqkm because
#'sum(PUxPMN)
#'[1] 2075

#change data type of areaArable field to numeric (to remove m^2 suffix)
AreaPMNtotal$areaPMNtot <- as.numeric(AreaPMNtotal$areaPMNtot)

#' PMN value for Indian Ocean = 320 USD/t (CRU, 2019) * 0.0056 t/sqm (Sharma et al., 2011) * 1.933646e+12 sqm (total PMN area in Indian Ocean, from shapefiles from ISA & chunk above) = USD 3.46509363e+12 for all PMN in the Indian Ocean (seems too much)

#' Plotting mining resources

#world <- ne_countries(scale = "medium", returnclass = "sf")%>% 
#st_transform(cCRS)
#IndianO <- c(xmin = 17.47941, ymin = -38.7402, xmax = 129.552, ymax = 33.82872)
#CFC <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Crust/Crust_areas_Hein2013.shp")%>% 
# st_transform(cCRS)
#PMN <- st_read ("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")%>% 
# st_transform(cCRS)

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = CFC, aes(color="A"), show.legend = "line") + 
  geom_sf(data = PMN, aes(color="B"), show.legend = "line") + 
  scale_color_manual(values = c("A" = "grey", "B" = "brown"),
                     labels = c("CFC", "PMN"),
                     name = "Legend") +
  ggtitle("Polymetallic Nodules and Cobalt-rich Ferromanganese Crust areas") + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)
# coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

#' assigning mining value to PUs df

#PUxPMNbinary <- as.numeric(PMN) #assigns NAs to all
#PUxPMNbinary <- if(TRUE) 1 else 0
#PUs <- PUs %>% 
# mutate(PMNbinary = PUxPMNbinary)

#PUxCFCbinary <- as.numeric(CFC)  
#PUs <- PUs %>% 
# mutate(CFCbinary = PUxCFCbinary)

#MiningCost <- case_when(
# df$Mining_layer1 > 0.5 ~ 14307812057,
#  df$Mining_layer2 > 0.5 ~ 94569447,
# TRUE ~ 0
#)

#https://dplyr.tidyverse.org/reference/case_when.html 
#https://www.sharpsightlabs.com/blog/case-when-r/ 
#https://dplyr.tidyverse.org/reference/if_else.html

#PUs$CFCcost <- case_when(
# PUs$CFC = FALSE ~ 0,
#TRUE ~ 94569447)

#CFC_val <- 94569447 #in USD per PU
#CFCxPU_col <- PUs$CFC

#MiningCost <- function(Mining_layer1, Mining_layer2) {
# case_when(
#  Mining_layer1 == "TRUE" ~ 14307812057,
# if CFCxPU_col == "TRUE" ~ 94569447,
#TRUE ~ 0
# )
#}

#PMNvalue <- function(PUs) {

# if(PUs$PMN == TRUE) {PUs$PMN[TRUE]=200}
#if(PUxPMN == FALSE) {PUs$PMN[FALSE]=0}

#}

#PMNvalue <- function(x) {

# if(PUs$PMN == TRUE) {PUs$PMN[TRUE]=200}
# if(PUxPMN == FALSE) {print("x was false!")}

#}
#iris$Sepal.Length[3]=999
#PUs <- PUs %>% 
# mutate(WARM = PUxWARM)

#df <- tibble(PU_ID = c(1, 2, 3, 4), Mining_layer1 = c(TRUE, FALSE, TRUE, FALSE), Mining_layer2 = c(FALSE, FALSE, TRUE, FALSE))
#MiningCost <- case_when(
# df$Mining_layer1 > 0.5 ~ 14307812057,
#df$Mining_layer2 > 0.5 ~ 94569447,
#TRUE ~ 0
#)

df <- tibble(PU_ID = c(1, 2, 3, 4), Mining_layer1 = c(TRUE, FALSE, TRUE, FALSE), Mining_layer2 = c(FALSE, TRUE, FALSE, FALSE))
MiningCost <- case_when(
  df$Mining_layer1 == "TRUE" ~ 14307812057,
  df$Mining_layer2 == "TRUE" ~ 94569447,
  TRUE ~ 0
)

df <- df %>% 
  mutate(MiningCostLayer = MiningCost)

#MiningCost <- case_when(
# df$Mining_layer1 == "TRUE" ~ 14307812057,
#df$Mining_layer2 == "TRUE" ~ 94569447,
#df$Mining_layer1 == "TRUE" & df$Mining_layer2 == "TRUE"  ~ 14307812057 + 94569447
#TRUE ~ 0
#)

#df <- df %>% 
## mutate(MiningCostLayer = MiningCost)

MiningCost <- case_when(
  PUs$PMN == "TRUE" ~ 14307812057, #value USD per PU
  PUs$CFC == "TRUE" ~ 94569447,
  TRUE ~ 0
)

PUs <- PUs %>% 
  mutate(MiningCostLayer = MiningCost)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$MiningCostLayer)) +
  geom_sf(data = world) +
  ggtitle("Mining Cost") 

sum(PUs$MiningCostLayer) # 2.917483e+13 looks like it worked
hist(PUs$MiningCostLayer) # bins together 0s and CFC values, issues when changing bin or defining xlim

length(PUs$MiningCostLayer[PUs$MiningCostLayer==94569447]) # = 2736, makes sense
length(PUs$MiningCostLayer[PUs$MiningCostLayer==14307812057]) # = 2021, also makes sense

#' Per sqm because values too high to solve in prioritizr

MiningCostSQM <- case_when(
  PUs$PMN == "TRUE" ~ 14.307812057, # value per sqm
  PUs$CFC == "TRUE" ~ 0.094569447, # value per sqm
  TRUE ~ 0
)

PUs <- PUs %>% 
  mutate(MiningCostSQM = MiningCostSQM)

ggplot() + 
  geom_sf(data = PUs, aes(color=PUs$MiningSQM)) +
  geom_sf(data = world) +
  ggtitle("Mining Cost") 

length(PUs$MiningCostSQM[PUs$MiningCostSQM==0.094569447]) # = 2736, makes sense #2245 when IO cut off at 45°S
length(PUs$MiningCostSQM[PUs$MiningCostSQM==14.307812057]) # = 2021, also makes sense #2078 when IO cut off at 45°S
max(MiningCostSQM) #14.30781
min(MiningCostSQM) # 0
mean(MiningCostSQM) # 0.7329807 #1.015738 when IO cut off at 45°S
median(MiningCostSQM) # 0

MiningCostSQKM <- case_when(
  PUs$PMN == "TRUE" ~ 14307812.057, # value per sqm
  PUs$CFC == "TRUE" ~ 0094569.447, # value per sqm
  TRUE ~ 0
)

PUs <- PUs %>% 
  mutate(MiningCostSQKM = MiningCostSQKM)

# OR https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile

#library(rgdal)
## read "/path/to/files/filename.shp"
#shp <- st_read("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/Nodules_Hein2013.shp")

## add new attribute data (just the numbers 1 to the number of objects)
#PMNwithCost <- 14307812057:nrow(shp)

## write out to a new shapefile
#writeOGR(shp, "~/Users/leafourchault/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Mining_resources/Nodules/", "filename2")  

#' Ant's exercise on fishing and mining: we want to know how many years we need to fish to reach the mining value

#FishingRaster <- terra::rast("Data/Cost/Cost_Raster_Sum.grd")

# fishing raster = USD/sqm/yr; max value : 2,609,303 
# mining = USD/PU

FishingPerPU <- PUs$Fishing * 1000  # * 1000 * 1000 because in sqkm and 1PU = 1000 sqkm
PUs <- PUs %>% 
  mutate(FishingPerPU = FishingPerPU)

max(PUs$FishingPerPU) # USD 5077817 per PU per yr # 5.769704e+12
mean(PUs$FishingPerPU) # USD 28340.43 per PU per yr ## 21015301585
min(PUs$FishingPerPU) # USD 226.7226 per PU per yr #1e-04 
median(PUs$FishingPerPU) # USD 11044.59 per PU per yr #8480084871 


#NumberOfYears  <- PUs$FishingPerPU / PUs$MiningCostLayer
#PUs <- PUs %>% 
# mutate(NumberOfYears = NumberOfYears)

#hist(PUs$NumberOfYears) #NaN for all
#max(PUs$NumberOfYears)
#min(PUs$NumberOfYears) 
#mean(PUs$NumberOfYears) 
#min(PUs$NumberOfYears) 
#median(PUs$NumberOfYears) 

#NumberOfYears2  <- PUs$Fishing / PUs$MiningCostLayer # keeping fishing /sqm instead of /PU
#PUs <- PUs %>% 
# mutate(NumberOfYears2 = NumberOfYears2)

#hist(PUs$NumberOfYears2) #NaN for all again
#max(PUs$NumberOfYears2)
#min(PUs$NumberOfYears2) 
#mean(PUs$NumberOfYears2) 
#median(PUs$NumberOfYears2) 

PUs$MiningCostLayer[PUs$MiningCostLayer <  0.0001] <-  0.0001 # otherwise cannot divide, because many 0s
PUs$FishingPerPU[PUs$FishingPerPU <  0.0001] <-  0.0001 # same as for mining cost layer

#NumberOfYears  <- PUs$FishingPerPU / PUs$MiningCostLayer # fishing value per PU
#PUs <- PUs %>% 
# mutate(NumberOfYears = NumberOfYears)

#hist(PUs$NumberOfYears) #mostly 0 to 1e+08, goes up to 6e+16
#max(PUs$NumberOfYears) #5.769704e+16 = 57697040 billion years, does not make sense
#min(PUs$NumberOfYears) #1.057424e-12
#mean(PUs$NumberOfYears) #1.886501e+14
#median(PUs$NumberOfYears) #6.646313e+13

PUs$MiningCostSQM[PUs$MiningCostSQM <  0.0001] <-  0.0001 # otherwise cannot divide, because many 0s
PUs$Fishing[PUs$Fishing <  0.0001] <-  0.0001 # same as for mining cost layer

PUs$MiningCostSQKM[PUs$MiningCostSQKM <  0.0001] <-  0.0001 # otherwise cannot divide, because many 0s

#NumberOfYears2  <- PUs$Fishing / PUs$MiningCostLayer # keeping fishing /sqm instead of /PU
#PUs <- PUs %>% 
# mutate(NumberOfYears2 = NumberOfYears2)

#hist(PUs$NumberOfYears2) # 0 to 6e+07
#max(PUs$NumberOfYears2) #57697037
#min(PUs$NumberOfYears2) #0
#mean(PUs$NumberOfYears2) #188650.1
#median(PUs$NumberOfYears2) #66463.13

#NumberOfYears3  <- PUs$Fishing / PUs$MiningCostSQM # keeping fishing /sqm instead of /PU
#PUs <- PUs %>% 
# mutate(NumberOfYears3 = NumberOfYears3)

#hist(PUs$NumberOfYears3) # 0 to 6e+07
#max(PUs$NumberOfYears3) #57697037
#min(PUs$NumberOfYears3) #0
#mean(PUs$NumberOfYears3) #188650.1
#median(PUs$NumberOfYears3) #66463.13

#' odd but 

#max(PUs$MiningCostLayer) #14307812057
#max(PUs$MiningCostSQM) #14.30781

#'

NumberOfYearsSQM  <- PUs$MiningCostSQM / PUs$Fishing # keeping fishing /sqm instead of /PU
PUs <- PUs %>% 
  mutate(NumberOfYearsSQM = NumberOfYearsSQM)

hist(PUs$NumberOfYearsSQM) # 
max(PUs$NumberOfYearsSQM) # 945.6945
min(PUs$NumberOfYearsSQM) # 1.733191e-08
mean(PUs$NumberOfYearsSQM) # 3.301491
median(PUs$NumberOfYearsSQM) # 1.504594e-05

NumberOfYears5  <- PUs$MiningCostLayer / PUs$FishingPerPU # fishing per PU and mining per PU
PUs <- PUs %>% 
  mutate(NumberOfYears5 = NumberOfYears5)

hist(PUs$NumberOfYears5) # 
max(PUs$NumberOfYears5) # 945694470000
min(PUs$NumberOfYears5) # 1.733191e-17
mean(PUs$NumberOfYears5) # 2280900162
median(PUs$NumberOfYears5) #1.504594e-14  but should be the same as NumberOfYearsSQM since both layers are per PU?

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = PUs, aes(color=NumberOfYearsSQM)) + 
  ggtitle("Number of Fishing Years to Reach Mining Value") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim)
#coord_sf(xlim = c(IndianO[["xmin"]], IndianO[["xmax"]]), ylim = c(IndianO[["ymin"]], IndianO[["ymax"]]), expand = FALSE)

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = PUs, aes(color=MiningCostSQM)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Mining cost sqm")

#' when looking only at areas with mining resources

intFishCFC <- st_intersection(CFC, FishingCost)
NumberOfYears_intFishCFC  <- 94569.447 / intFishCFC$FishingCost # mining value per sqkm, fishing cost units as created by jase (units unsure, probably sqkm, need to ask)
intFishCFC <- intFishCFC %>% 
  mutate(NumberOfYears_intFishCFC = NumberOfYears_intFishCFC)

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = intFishCFC, aes(color=NumberOfYears_intFishCFC)) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Number of Fishing Years to Reach Mining Value in CFC Areas")

intFishPMN <- st_intersection(PMN, FishingCost)
NumberOfYears_intFishPMN  <- 14307812.057 / intFishPMN$FishingCost # mining value per sqkm, fishing value assumed to be in USD per sqkm - but ask Jase
intFishPMN <- intFishPMN %>% 
  mutate(NumberOfYears_intFishPMN = NumberOfYears_intFishPMN)

ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = intFishPMN, aes(color=NumberOfYears_intFishPMN)) + 
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  ggtitle("Number of Fishing Years to Reach Mining Value in PMN Areas")

#' Shipping penalty

Shipping <- terra::rast("~/Documents/MscThesis/GitHub/SpatialPlanning/Shp/Costs/Shipping/shipping.tif")
#View(Shipping)
#plot(Shipping, main = 'Shipping routes')
ShippingRobin <- terra::project(Shipping, "ESRI:54030") #transform to robin proj
#plot(ShippingRobin, main = 'Shipping routes')
PUxShipping <- exactextractr::exact_extract(ShippingRobin, PUs, "mean") # attribute mean value of pixels in each PU to PU
PUs <- PUs %>% 
  mutate(ShippingIntensity = PUxShipping) # add column to PU data frame

max(PUs$ShippingIntensity)
#[1] 335.2159 # 345.664 when cut off at 45°S although this is probably because I changed from 20°N to 34°N
mean(PUs$ShippingIntensity)
#[1] 3.75951  #5.175909 when cut off at 45°S

#ShippingPlusOne <- Shipping$shipping +1
#ShippingLog <- log10(ShippingPlusOne)
#hist(log10(ShippingPlusOne))  #Warning: [hist] a sample of0% of the cells was used (of which 21% was NA); meaning +1 didn't work?

#ShippingRobin <- terra::project(Shipping, "cCRS") #Error: [project] cannot get output boundaries

#ShippingRobin <- project(Shipping, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 #+datum=WGS84 +units=m +no_defs") # no error but does not assign projection

#ShippingRobin <- project(Shipping, "ESRI:54030", method = "near") #transform to robin proj # error about method "near"



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
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #coord_sf(xlim = c(18, 120), ylim = c(34, -45), expand = FALSE) +
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
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  # coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_Active_Vents_Fishing_Sol_Plot # also looks correct

#' Cons features & Fishing cost & targets & locked-out exploration + reserved areas, change targets to check solver works

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
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  # coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_check_targets_sol_plot

#' All features excluding longhurst vs fishing cost

p_NoLonghurst <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>% #inactive ventstarget = 0.3 otherwise unfeasible
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_NoLonghurst_sol <- solve(p_NoLonghurst)

p_NoLonghurst_sol_plot <- ggplot() +
  ggtitle("No Longhurst vs Fishing, Mining Out") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_NoLonghurst_sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_NoLonghurst_sol_plot


#' All features including  longhurst vs Fishing cost 

p_ConsFeat_Fishing <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents", "WARM", "SSTC", "SPSG", "SANT", "NPSW", "NEWZ", "ISSG", "ANTA"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_ConsFeat_Fishing_Sol <- solve(p_ConsFeat_Fishing)
#plot(p_ConsFeat_Fishing_Sol, main = "solution")
#lines(p_ConsFeat_Fishing_Sol[p_ConsFeat_Fishing_Sol$solution_1 > 0.5, ], col = "darkgreen", lwd = 2)

p_ConsFeat_Fishing_Sol_Plot <- ggplot() +
  ggtitle("Cons Feats including Longhurst vs Fishing, 0.1<Targets<0.7") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_ConsFeat_Fishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_ConsFeat_Fishing_Sol_Plot # shows that fishing cost is reduced between 40°-60°S, but maybe would be good to reduce target for IBAs/IMMAs or assign Longhurst to distribute better across latitudes #lgh assigned

#' All features & Fishing cost & targets & locked-out exploration + reserved areas

PUs["locked_out"] = PUs$ExplPMN | PUs$ExplPMS | PUs$ResPMN # create locked_out column (logical). If =<1 is TRUE, then locked_out is TRUE # https://stackoverflow.com/questions/54507486/merging-two-true-false-dataframe-columns-keeping-only-true

#PUs["locked_out"] <- paste(PUs$ExplPMN, PUs$ExplPMS, PUs$ResPMN)
#PUs %>% mutate(sum = rowSums(across(where(is.logical))))rowSums(Y) > 0
#PUs$LockedOut <- ifelse(PUs$ExplPMN == "TRUE" | PUs$ExplPMS == "TRUE" | PUs$ResPMN == "TRUE")

p_NoLonghurst_LockedOut_Fishing_Shipping <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_linear_penalties(1, "ShippingIntensity") %>% 
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_NoLonghurst_LockedOut_Fishing_Shipping_Sol <- solve(p_NoLonghurst_LockedOut_Fishing_Shipping) # when adding longhurst provinces to normal locked_out problem: Error in .local(a, b = b, ...) : 
#no solution found (e.g. due to problem infeasibility or time limits); error solved when reducing targets for inactive vents to 0.3

p_NoLonghurst_LockedOut_Fishing_Shipping_Sol_Plot <- ggplot() +
  ggtitle("No Longhurst vs Fishing, Mining Contractors Out, Shipping = 1") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_NoLonghurst_LockedOut_Fishing_Shipping_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_NoLonghurst_LockedOut_Fishing_Shipping_Sol_Plot

#' Trial shipping penalty

p_IBA_ShippingFishing <- problem(PUs, features = "IBA", cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(0.9) %>%
  add_linear_penalties(0, "ShippingIntensity") %>% 
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_IBA_ShippingFishing_Sol <- solve(p_IBA_ShippingFishing)

p_IBA_ShippingFishing_Sol_Plot <- ggplot() +
  ggtitle("IBA vs Fishing and Shipping, Penalty = 0") +
  geom_sf(data = world, colour ="grey70") +
  geom_sf(data = p_IBA_ShippingFishing_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  # coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_IBA_ShippingFishing_Sol_Plot # correct

p_IBA_ShippingFishing2 <- problem(PUs, features = "IBA", cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(0.9) %>%
  add_linear_penalties(200, "ShippingIntensity") %>% 
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_IBA_ShippingFishing_Sol2 <- solve(p_IBA_ShippingFishing2)

p_IBA_ShippingFishing_Sol_Plot2 <- ggplot() +
  ggtitle("IBA vs Fishing and Shipping, Penalty = 200") +
  geom_sf(data = world, colour ="grey70") +
  geom_sf(data = p_IBA_ShippingFishing_Sol2, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_IBA_ShippingFishing_Sol_Plot2 # also correct

#' all cons feats, fishing and shipping + contractors locked out

p_FeatsFishShipLockOut <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents", "WARM", "SSTC", "SPSG", "SANT", "NPSW", "NEWZ", "ISSG", "ANTA"), cost_column = "Fishing") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)) %>%
  add_linear_penalties(200, "ShippingIntensity") %>% 
  add_binary_decisions() %>%
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_FeatsFishShipLockOut_Sol <- solve(p_FeatsFishShipLockOut)

p_FeatsFishShipLockOut_Sol_Plot <- ggplot() +
  ggtitle("p_FeatsFishShipLockOut, Penalty = 200") +
  geom_sf(data = world, colour ="grey70") +
  geom_sf(data = p_FeatsFishShipLockOut_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #  coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_FeatsFishShipLockOut_Sol_Plot # correct

#' mining try

#' usual issues with plotting. also, was selecting all PUs except the ones with CFC or PMN on them, probably because cost = 0

PUs$MiningCostSQM[PUs$MiningCostSQM <  0.0001] <-  0.0001

p_IBA_MiningSQM <- problem(PUs, features = c("IBA"), cost_column = "MiningCostSQM") %>% 
  add_min_set_objective() %>%
  add_relative_targets(0.9) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, verbose = FALSE)

p_IBA_MiningSQM_Sol <- solve(p_IBA_MiningSQM) # when adding longhurst provinces to normal locked_out problem: Error in .local(a, b = b, ...) : 
#no solution found (e.g. due to problem infeasibility or time limits); error solved when reducing targets for inactive vents to 0.3

p_IBA_MiningSQM_Sol_Plot <- ggplot() +
  ggtitle("IBA vs Mining") +
  #geom_sf(data = world, colour ="grey80") +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  geom_sf(data = p_IBA_MiningSQM_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_IBA_MiningSQM_Sol_Plot


#' with mining cost layer + locked out contractors, shipping penalty, no longhurst
PUs$MiningCostLayer[PUs$MiningCostLayer <  0.0001] <-  0.0001

p_NoLonghurst_LockedOut_Mining_Shipping <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCostLayer") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_linear_penalties(1, "ShippingIntensity") %>% 
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_NoLonghurst_LockedOut_Mining_Shipping_Sol <- solve(p_NoLonghurst_LockedOut_Mining_Shipping) # error because planning units with very high (> 1e+6) cost values

p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol_Plot <- ggplot() +
  ggtitle("No Longhurst vs Mining, Mining Contractors Out, Shipping = 1") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol_Plot


#' with mining cost layer + locked out contractors, shipping penalty, no longhurst

PUs$MiningCostSQM[PUs$MiningCostSQM <  0.0001] <-  0.0001

p_NoLonghurst_LockedOut_MiningSQM_Shipping <- problem(PUs, features = c("IBA", "IMMA", "NEIO_09", "NWIO_14", "SIO_11", "SIO_19", "SIO_22", "SIO_23", "SIO_30", "SIO_32", "SIO_35", "SIO_36", "SIO_37","Seamounts", "Plateaus", "Active_Vents", "Inactive_Vents"), cost_column = "MiningCostSQM") %>% 
  add_min_set_objective() %>%
  add_relative_targets(c(0.7, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.7, 0.3)) %>%
  add_binary_decisions() %>%
  add_linear_penalties(1, "ShippingIntensity") %>% 
  add_locked_out_constraints(locked_out = "locked_out") %>%
  add_gurobi_solver(gap = 0.1, verbose = FALSE)

p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol <- solve(p_NoLonghurst_LockedOut_MiningSQM_Shipping) # error because planning units with very high (> 1e+6) cost values

p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol_Plot <- ggplot() +
  ggtitle("No Longhurst vs Mining, Mining Contractors Out, Shipping = 1") +
  geom_sf(data = world, colour ="grey80") +
  geom_sf(data = p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol, aes(fill = as.factor(solution_1)), size = 0.05) +
  scale_fill_manual(values = c("NA", "red")) +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
  #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) +
  theme_bw() +
  #scale_colour_gradientn(colours = colorspace::diverge_hcl(7)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol_Plot

#' ## Comparing solutions

# based on Jase's functions:

solA <- p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol
solB <- p_NoLonghurst_LockedOut_Fishing_Shipping_Sol

fCompareSolutions <- function(solA, solB){ #modified function from Jase, run every paragraph separetely to avoid error with: gg  not found
  
  #Select only the ID column and the solution column; change the name of the columns in 'cellID' and 'ctrl_sol'  
  ctrl2 <- solA %>%
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, solA_sol = solution_1)
  
  rmv2 <- solB %>%
    dplyr::select(geometry, solution_1) %>%
    rename(geometry = geometry, solB_sol = solution_1)
  
  #Create an object with the two solutions
  soln <- ctrl2
  soln$solB_sol <-  rmv2$solB_sol
  
  #Combine the result and call them differently if a PUs in present in both the solution, or only in one of them
  soln <- mutate(soln, Combined = as.numeric(solA_sol + solB_sol)) %>%
    mutate(Compare = case_when(Combined == 2 ~ "Same",
                               solA_sol == 1 & solB_sol == 0 ~ "Removed (-)",
                               solA_sol == 0 & solB_sol == 1 ~ "Added (+)",
                               Combined == 0 ~ "Not Selected"),
           Compare = factor(Compare, levels = c("Added (+)", "Same", "Removed (-)", "Not Selected"))) %>%
    filter(!is.na(Compare))
  
  # label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%\nCost: ",round((sum(rmv$solution_1*rmv$cost) - sum(ctrl$solution_1*ctrl$cost))/ sum(ctrl$solution_1*ctrl$cost)*100, 1), "%\n")
  # label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%")
  
  #Plot the results
  gg <- ggplot() +
    ggtitle("Mining versus Fishing, Shipping penalty = 1") +
    geom_sf(data = world, colour ="grey0") +
    geom_sf(data = soln, aes(fill = Compare), colour = NA, size = 0.0001) +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim) +
    #coord_sf(xlim = c(18, 130), ylim = c(22, -60), expand = FALSE) + #put this line here to avoid issues when plotting
    theme_bw() +
    theme(axis.title = element_blank(), text = element_text(size = 20), plot.title = element_text(size = 12)) +
    scale_fill_manual(values = c("Added (+)" = "Red", "Same" = "ivory3", "Removed (-)" = "Blue", "Not Selected" = "ivory4"), drop = FALSE) +
    return(gg)
}

gg

miningVSfishing <- fCompareSolutions(solA, solB)

#' summary statistics Mining plan

SolutionMining <- as_tibble(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol$solution_1)
FeaturesMiningShippingSummary <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping, SolutionMining)

MinBoundSum <- eval_boundary_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinConSum <- eval_connectivity_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinCostSum <- eval_cost_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinFeatRep <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinPUSum <- eval_n_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinTargetSum <- eval_target_coverage_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)
MinFeatAbd <- feature_abundances(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol)

#' summary statistics Fishing plan

SolutionFishing <- as_tibble(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol$solution_1)
FeaturesFishingShippingSummary <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_Fishing_Shipping, SolutionFishing)

FishBoundSum <- eval_boundary_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishConSum <- eval_connectivity_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishCostSum <- eval_cost_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishFeatRep <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishPUSum <- eval_n_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishTargSum <- eval_target_coverage_summary(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)
FishFeatAbd <- feature_abundances(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol)

#' for 3 solutions

#fCompareSolutions <- function(ctrl, rmv, aaa){

#Select only the ID column and the solution column; change the name of the columns in 'cellID' and 'ctrl_sol'  
# ctrl2 <- ctrl %>%
# dplyr::select(ID, solution_1) %>%
#  rename(cellID = ID, ctrl_sol = solution_1)

#  rmv2 <- rmv %>%
#   dplyr::select(ID, solution_1) %>%
#  rename(cellID = ID, rmv_sol = solution_1)

#rmv3 <- aaa %>%
# dplyr::select(ID, solution_1) %>%
#rename(cellID = ID, rmv_sol3 = solution_1)

#Create an object with the three solutions
#  soln <- ctrl2
# soln$rmv_sol <- rmv2$rmv_sol
#soln$rmv_sol3 <- rmv3$rmv_sol3

#Combine the result and call them differently if a PUs in present in both the solution, or only in one of them
#  soln <- mutate(soln, Combined = as.numeric(ctrl_sol + rmv_sol + rmv_sol3)) %>%
#   mutate(Compare = case_when(Combined == 3 ~ "Same",
#                             ctrl_sol == 1 & rmv_sol == 1 & rmv_sol3 == 0 ~ "No species & Total area",
#                            ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 1 ~ "Total area & Fraction of the area",
#                           ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 1 ~ "No species & Fraction of the area",
#                          ctrl_sol == 1 & rmv_sol == 0 & rmv_sol3 == 0 ~ "No species",
#                         ctrl_sol == 0 & rmv_sol == 1 & rmv_sol3 == 0 ~ "Total area",
#                       ctrl_sol == 0 & rmv_sol == 0 & rmv_sol3 == 1 ~ "Fraction of the area",
#                      Combined == 0 ~ "Not Selected"),
# Compare = factor(Compare, levels = c("Same",
#                                     "No species & Total area",
#                                    "Total area & Fraction of the area",
#                                   "No species & Fraction of the area",
#                                  "No species",
#                                 "Total area",
#                                "Fraction of the area", 
#                               "Not Selected"))) %>%
#  filter(!is.na(Compare))
#
# label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%\nCost: ",round((sum(rmv$solution_1*rmv$cost) - sum(ctrl$solution_1*ctrl$cost))/ sum(ctrl$solution_1*ctrl$cost)*100, 1), "%\n")
# label_loc$label <- paste0("Area: ",round((sum(rmv$solution_1) - sum(ctrl$solution_1))/ sum(ctrl$solution_1)*100, 1), "%")

#Plot the results
#gg <- ggplot() +
# geom_sf(data = soln, aes(fill = Compare), colour = "black", size = 0.0001) +
#theme_bw() +
#theme(axis.title = element_blank(), text = element_text(size = 20), plot.title = element_text(size = 12)) +
#scale_fill_manual(values = c("Same" = "red",
#                            "No species & Total area" = "magenta",
#                           "Total area & Fraction of the area" = "blue",
#                          "No species & Fraction of the area" = "darkorchid",
#                         "No species" = "limegreen",
#                        "Total area" = "orange",
#                       "Fraction of the area" = "lightsalmon1",
#                      "Not Selected" = "white"))

#  return(gg)
#}


#' evaluate feature representation

#SolutionMining <- as_tibble(p_NoLonghurst_LockedOut_MiningSQM_Shipping_Sol$solution_1)
#FeaturesMiningShippingSummary <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_MiningSQM_Shipping, SolutionMining)

#SolutionFishing <- as_tibble(p_NoLonghurst_LockedOut_Fishing_Shipping_Sol$solution_1)
#FeaturesFishingShippingSummary <- eval_feature_representation_summary(p_NoLonghurst_LockedOut_Fishing_Shipping, SolutionFishing)

#' circular plot


library(tidyverse)

# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value=sample( seq(10,100), 60, replace=T)
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p

create_circularbplot <- function(df, target) {
  # colors
  scenario.legend_color <- c('Mining and Shipping' = 'tan3', 'Fishing and Shipping' = 'darkseagreeen', 
                             #'SSP585' = 'salmon4', 'uninformed' = 'grey30', 'NA' = NA)
                             scenario.legend_list <- c('Mining', 'Fishing', 'SSP585', 'uninformed', ' ')
                             species.legend_color <- c('Important Bird Areas' = 'cornflowerblue', 'Important Marine Mammal Areas' = 'cadetblue4', 
                                                       'Sumatra Upwelling Zone' = 'skyblue', 'Arabian Oxygen Minimum Zone' = 'paleturquoise3',
                                                       'Mozambique Channel' = 'olivedrab', 'Agulhas Upwelling Zone' = 'color', 'Walters Shoals' = 'darkolivegreen3', 
                                                       'Saya de Malha' = 'darkseagreen', 'Central Indian Ocean Basin' = 'limegreen', 
                                                       'Fracture Zone' = 'springgreen1', 'Atlantis Seamount' = 'color', 'Rusky Knoll' = 'color', 'Fools Flat' = 'color', 'East Broken Ridge' = 'color', 'Seamounts' = 'color', 'Inactive Vents' = 'color', 'Active Vents' = 'color', 'Plateaus' = 'color')
                             group.legend_color <- c('Epipelagic' = 'chartreuse4', 'Mesopelagic' = 'royalblue', 'Bathypelagic' = 'color', 'Abyssopelagic' = 'color')
                             
                             
                             # manipulating df # mutate columns of solution_1 for fishing and solution_1 for mining to df?  df: 1 column = features, 1 column = targets, 1 column = fishing or shipping (plan), 1 column ocean layer
                             target <- as.numeric(target) # 0.3, 0.7 -  need to double
                             target <- 0.3
                             df <- merge(FeaturesFishingShippingSummary, FeaturesMiningShippingSummary, by="feature")
                             
                             df.manip <- df %>%  # create df with all features in 1 column, respective targets reached across each  scenario in matching columns # check prioritizr function with problem and solution objects # mutate fishing and mining : use cbind to join two dataframes
                               dplyr::mutate(scenario = case_when(str_detect(plan, pattern = 'Fishing') ~ 'Fishing and Shipping', #fishing or mining
                                                                  str_detect(plan, pattern = 'Mining') ~ 'Mining and Shipping',
                                                                  #str_detect(plan, pattern = 'SSP585') ~ 'SSP585',
                                                                  #str_detect(plan, pattern = 'uninformed') ~ 'uninformed'),
                                                                  group = case_when(str_detect(features, pattern = "IBA") ~ "Epipelagic",
                                                                                    str_detect(features, pattern = "IMMA") ~ "Epipelagic",
                                                                                    str_detect(features, pattern = "SIO_19") ~ "Epipelagic",
                                                                                    str_detect(features, pattern = "SIO_22") ~ "Epipelagic",
                                                                                    str_detect(features, pattern = "SIO_32") ~ "Epipelagic",
                                                                                    str_detect(features, pattern = "SIO_34") ~ "Epipelagic", # or meso
                                                                                    str_detect(features, pattern = "SIO_11") ~ "Mesopelagic", # or epi
                                                                                    str_detect(features, pattern = "NEIO_09") ~ "Mesopelagic", # or epi
                                                                                    str_detect(features, pattern = "NWIO_14") ~ "Mesopelagic",
                                                                                    str_detect(features, pattern = "SIO_30") ~ "Mesopelagic",
                                                                                    str_detect(features, pattern = "SIO_35") ~ "Mesopelagic",
                                                                                    str_detect(features, pattern = "SIO_36") ~ "Mesopelagic",
                                                                                    str_detect(features, pattern = "SIO_37") ~ "Mesopelagic",
                                                                                    str_detect(features, pattern = "Seamounts") ~ "Bathypelagic",
                                                                                    str_detect(features, pattern = "Inactive_Vents") ~ "Bathypelagic",
                                                                                    str_detect(features, pattern = "Active_Vents") ~ "Bathypelagic",
                                                                                    str_detect(features, pattern = "Plateaus") ~ "Bathypelagic",
                                                                                    str_detect(features, pattern = "SIO_23") ~ "Abyssopelagic")) %>% 
                                               dplyr::mutate(individual = case_when(str_detect(features, pattern = "IBA") ~ "Important Bird Area",
                                                                                    str_detect(features, pattern = "IMMA") ~ "Important Marine Mammal Area",
                                                                                    str_detect(features, pattern = "SIO_19") ~ "Mozambique Channel",
                                                                                    str_detect(features, pattern = "SIO_22") ~ "Walters Shoals",
                                                                                    str_detect(features, pattern = "SIO_32") ~ "Saya de Malha",
                                                                                    str_detect(features, pattern = "SIO_34") ~ "Central Indian Ocean Basin",
                                                                                    str_detect(features, pattern = "SIO_11") ~ "Agulhas Upwelling Zone",
                                                                                    str_detect(features, pattern = "NEIO_09") ~ "Sumatra Upwelling Zone",
                                                                                    str_detect(features, pattern = "NWIO_14") ~ "Arabian Oxygen Minimum Zone",
                                                                                    str_detect(features, pattern = "SIO_30") ~ "Atlantis Seamount",
                                                                                    str_detect(features, pattern = "SIO_35") ~ "Rusky Knoll",
                                                                                    str_detect(features, pattern = "SIO_36") ~ "Fools Flat",
                                                                                    str_detect(features, pattern = "SIO_37") ~ "East Broken Ridge",
                                                                                    str_detect(features, pattern = "Seamounts") ~ "Seamounts",
                                                                                    str_detect(features, pattern = "Inactive_vents") ~ "Inactive Vents",
                                                                                    str_detect(features, pattern = "Active Vents") ~ "Active Vents",
                                                                                    str_detect(features, pattern = "SIO_23") ~ "Fracture Zone",
                                                                                    str_detect(features, pattern = "Plateaus") ~ "Plateaus")) %>% 
                                               dplyr::rename(value = representation)
                                             
                                             # adding NA rows for each feature
                                             ALB <- data.frame(features = 'ALB', plan = NA, value = 0, scenario = NA, group = 'commercial', individual = 'Albacore tuna')
                                             Caretta <- data.frame(features = 'Caretta_caretta_IUCN', plan = NA, value = 0, scenario = NA, group = 'bycatch')
                                             Chelonia <- data.frame(features = 'Chelonia_mydas_IUCN', plan = NA, value = 0, scenario = NA, group = 'bycatch')
                                             Dermochelys <- data.frame(features  = 'Dermochelys_coriacea_IUCN', plan = NA, value = 0, scenario = NA, group = 'bycatch')
                                             Eretmochelys <- data.frame(features  = 'Eretmochelys_imbricata_IUCN', plan = NA, value = 0, scenario = NA, group = 'bycatch')
                                             Lepidochelys <- data.frame(features  = 'Lepidochelys_olivacea_IUCN', plan = NA, value = 0, scenario = NA, group = 'bycatch')
                                             SKP <- data.frame(features  = 'SKP', plan = NA, value = 0, scenario = NA, group = 'commercial')
                                             SWO <- data.frame(features  = 'SWO', plan = NA, value = 0, scenario = NA, group = 'commercial')
                                             YFT <- data.frame(features  = 'YFT', plan = NA, value = 0, scenario = NA, group = 'commercial')
                                             
                                             
                                             data <- df.manip %>%  #for blank bar
                                               bind_rows(ALB) %>% 
                                               bind_rows(Caretta) %>% 
                                               bind_rows(Chelonia) %>% 
                                               bind_rows(Dermochelys) %>% 
                                               bind_rows(Eretmochelys) %>% 
                                               bind_rows(Lepidochelys) %>% 
                                               bind_rows(SKP) %>% 
                                               bind_rows(SWO) %>% 
                                               bind_rows(YFT) %>% 
                                               group_by(group) %>% 
                                               arrange(features)
                                             
                                             # creating plot
                                             # Set a number of 'empty bar' to add at the end of each group
                                             empty_bar <- 1
                                             to_add <- data.frame(matrix(NA, empty_bar*length(unique(data$group)), ncol(data)) ) #data=full df + NAs, group = ocean layer
                                             colnames(to_add) <- colnames(data)
                                             to_add$group <- rep(levels(as.factor(data$group)), each=empty_bar)
                                             data <- rbind(data, to_add)
                                             data <- data %>% arrange(group)
                                             data$id <- seq(1, nrow(data)) 
                                             
                                             # prepare a data frame for base lines
                                             base_data <- data %>% 
                                               group_by(group) %>% 
                                               dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
                                               rowwise() %>% 
                                               dplyr::mutate(title = mean(c(start, end)))
                                             
                                             ##########################
                                             ## Data/lines of species ##
                                             ##########################
                                             species_data <- data %>% 
                                               group_by(individual) %>% 
                                               dplyr::summarize(start = min(id), end = max(id)) %>% 
                                               dplyr::mutate(title = mean(c(start, end)))
                                             
                                             species_data[1,3] <- 30
                                             
                                             # Get the name and the y position of each label for the species
                                             label_sp <- data %>% 
                                               group_by(individual) %>% 
                                               dplyr::summarize(individual = unique(individual, na.rm = TRUE),
                                                                id = mean(id, na.rm = TRUE))
                                             number_of_bar <- nrow(label_sp)
                                             angle <- 360 - 90 * (label_sp$id - 2.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
                                             label_sp$hjust <- ifelse( angle < -90, 1, 0)
                                             label_sp$angle <- ifelse(angle < -90, angle+180, angle)
                                             
                                             # for the percentage lines
                                             grid_data <- base_data
                                             grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1.5
                                             grid_data$start <- grid_data$end - 1
                                             grid_data <- grid_data[-1,]
                                             
                                             p <- ggplot(data, aes(x = as.factor(id), y = value, fill = scenario)) + 
                                               
                                               # plotting the bars
                                               geom_bar(aes(x = as.factor(id), y = value, fill = scenario), 
                                                        stat = "identity", 
                                                        position = 'dodge') +
                                               
                                               # defining colors of the bars
                                               scale_fill_manual(name = "Solution",
                                                                 values = scenario.legend_color,
                                                                 labels = scenario.legend_list) +
                                               
                                               # Add text showing the value of each 100/75/50/25 lines
                                               geom_segment(data = grid_data, 
                                                            aes(x = end, y = 10, xend = start, yend = 10), 
                                                            colour = "grey50", 
                                                            alpha = 1, 
                                                            size = 0.5 , 
                                                            inherit.aes = FALSE ) +
                                               geom_segment(data = grid_data, 
                                                            aes(x = end, y = 20, xend = start, yend = 20), 
                                                            colour = "grey50", 
                                                            alpha = 1,
                                                            size = 0.5,
                                                            inherit.aes = FALSE ) +
                                               geom_segment(data = grid_data, 
                                                            aes(x = end, y = 30, xend = start, yend = 30), 
                                                            colour = "grey50", 
                                                            alpha = 1,
                                                            size = 0.5,
                                                            inherit.aes = FALSE ) +
                                               geom_segment(data = grid_data, 
                                                            aes(x = end, y = 40, xend = start, yend = 40), 
                                                            colour = "grey50", 
                                                            alpha = 1,
                                                            size = 0.5,
                                                            inherit.aes = FALSE ) +
                                               geom_segment(data = grid_data, 
                                                            aes(x = end, y = 50, xend = start, yend = 50), 
                                                            colour = "grey50", 
                                                            alpha = 1,
                                                            size = 0.5,
                                                            inherit.aes = FALSE ) +
                                               annotate("text", x = rep(max(data$id),5), 
                                                        y = c(10, 20, 30, 40, 50), 
                                                        label = c('10','20','30','40','50'), 
                                                        color = "grey50", 
                                                        size=4, 
                                                        angle = -5, 
                                                        fontface = "bold", 
                                                        hjust=0.5) +
                                               
                                               # setting limitations of actual plot
                                               ylim(-50,55) +
                                               theme_minimal() +
                                               coord_polar() + 
                                               
                                               # Add base line information (commercial + bycatch labels)
                                               geom_segment(data = base_data, 
                                                            aes(x = start, y = -5, xend = end, yend = -5), 
                                                            colour = group.legend_color,
                                                            alpha = 0.8, 
                                                            size = 0.6, 
                                                            inherit.aes = FALSE, 
                                                            show.legend = FALSE)  +
                                               geom_text(data = base_data, 
                                                         aes(x = title, y = -8, label = group), 
                                                         color = group.legend_color,
                                                         hjust = c(1,0), 
                                                         alpha = 0.8, 
                                                         size = 4, 
                                                         fontface = "bold", 
                                                         inherit.aes = FALSE, 
                                                         show.legend = FALSE) +
                                               
                                               # Adding the lines for the species
                                               geom_segment(data = species_data, 
                                                            aes(x = start, y = 55, xend = end, yend = 55, color = individual), 
                                                            alpha = 1, 
                                                            size = 5, 
                                                            inherit.aes = FALSE)  +
                                               
                                               # Defining colors of these lines
                                               scale_color_manual(name = "Species",
                                                                  values = species.legend_color) +
                                               
                                               # Adding the target % for each species # need 2 error bars - 0.3, 0.7
                                               geom_errorbar(aes(y = target, ymax = target, ymin = target), 
                                                             color = 'red', 
                                                             linetype = 'dashed', 
                                                             size = 0.8) +
                                               
                                               theme(
                                                 legend.position = "bottom",
                                                 axis.text = element_blank(),
                                                 axis.title = element_blank(),
                                                 panel.grid = element_blank(),
                                                 plot.margin = unit(rep(0.5,4), "cm") 
                                               )
                                             
                                             return(p)
}  

#' kappa correlation plot

fcreate_kappacorrplot <- function(sol1, sol2, sol3, sol4, dir) {
  
  library(irr)
  library(corrplot)
  
  Species <- sol1 %>% 
    as_tibble() %>% 
    dplyr::select(solution_1) %>% 
    dplyr::rename(Species = solution_1)
  Species_Fisheries <- sol2 %>% 
    as_tibble() %>% 
    dplyr::select(solution_1) %>% 
    dplyr::rename(Species_Fisheries = solution_1)
  Species_CoastalProtection <- sol3 %>% 
    as_tibble() %>% 
    dplyr::select(solution_1) %>% 
    dplyr::rename(Species_CoastalProtection = solution_1)
  Species_Carbon <- sol4 %>% 
    as_tibble() %>% 
    dplyr::select(solution_1) %>% 
    dplyr::rename(Species_Carbon = solution_1)
  
  s_list <- list(Species, Species_Fisheries, Species_CoastalProtection, Species_Carbon)
  y = 1
  s_matrix <- list()
  for(i in 1:4){
    for(j in 1:4){
      kappa_temp <- irr::kappa2(bind_cols(s_list[[i]], s_list[[j]]))
      kappa_corrvalue <- kappa_temp$value
      kappa_pvalue <- kappa_temp$p.value
      s_matrix[[y]] <- cbind(colnames(s_list[[i]]), colnames(s_list[[j]]), kappa_corrvalue, kappa_pvalue)
      y = y+1
    }
  }
  s_matrix_all <- do.call(rbind, s_matrix) %>% 
    as_tibble()
  colnames(s_matrix_all)[1:2] <- c('plan1','plan2')
  
  matrix_final1 <- s_matrix_all %>% 
    as_tibble() %>% 
    dplyr::select(-kappa_pvalue) %>% 
    pivot_wider(names_from = plan2, values_from = kappa_corrvalue) %>% 
    as.matrix()
  
  matrix_final2 <- s_matrix_all %>% 
    as_tibble()
  
  write_csv(matrix_final2, paste0(dir,"kappa_matrix.csv"))
  
  # creating corrplot
  rownames(matrix_final1) <- matrix_final1[,1]
  n <- 4 + 1 # 4 is the number of inputted scenarios
  matrix_final2 <- matrix_final1[,2:n]
  class(matrix_final2) <- "numeric"
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  plot <- corrplot(matrix_final2, method = "shade", cl.lim = c(-0.02,1), tl.col = "black", addCoef.col = "black",
                   col=col(200), tl.srt=45)
  return(plot)
}