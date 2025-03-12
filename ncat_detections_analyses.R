#---------------------------------------------------------------------------------------------------#
#R pipeline: telemetry data - NCAT project
#Prepared by Lisette Delgado, email: m.lisette.delgado@dal.ca
#Based on OTN's telemetry workshop: https://ocean-tracking-network.github.io/otn-workshop-base/
#---------------------------------------------------------------------------------------------------#


#--------------------------------- Packages installations ------------------------------------------#
#Once the packages are installed, no need to run this section again
#Most of the packages are at the CRAN repository, so can be installed easily with install.package()
#A couple of packages need a url address, see below:
#Installation devtools package
install.packages("remotes")
library(remotes)
install.packages("devtools")
devtools:::install_github("gearslaboratory/gdalUtils")
#Installation of glatos package
install_url("https://gitlab.oceantrack.org/GreatLakes/glatos/-/archive/master/glatos-master.zip",
            build_opts = c("--no-resave-data", "--no-manual"))
install.packages("pathroutr", repos = "https://jmlondon.r-universe.dev")


#--------------------------------- Library call ---------------------------------------------------#
#Run at the start of every run
library(glatos)
library(actel)
library(ggplot2)
library(gganimate)
library(gifski)
library(spdplyr)
library(tidyverse)
library(dplyr)
library(marmap)
library(stringr)
library(readxl)
library(lubridate)
library(ggmap)
library(plotly)
library(viridis)
library(pathroutr)
library(sf)
library(ggspatial)
library(raster)
library(sfnetworks)
library(dygraphs)
library(geosphere)
library(reshape2)
library(spData)
library(oce)
library(ggrepel)
library(gridExtra)
#library(tidyr) do not run this library, it interferes with another library's function


#--------------------------------- Set working directory ---------------------------------------#
#***Change the path to directory accordingly
setwd("~/Cod/OTN")


#--------------------------------- Base maps ---------------------------------------------------#
#Get base map from NOAA
#Can change the resolution "res" if prefer a higher or lower resolution
#***Change coordinates accordingly to deployment sites
#These maps focuses on south Labrador and Newfoundland
bathyData <- getNOAA.bathy(-63,-46,56,46, res=3, keep=T)
bathyData_2 <- getNOAA.bathy(-60,-47,55.5,44.5, res=3, keep=T)

#Zoom out to include detections in the arctic and south US, to see all the outliers
bathyData_arctic <- getNOAA.bathy(-77,-46,65,34, res=5, keep=T)

#Zoom in to the stations
bathyData_stations <- getNOAA.bathy(-55.7,-48,54.7,47.5, res=1, keep=T)

#Plot base map
#Can add contour is prefer with the function "geom_contour"
base_noaa <- autoplot(bathyData, geom=c("r")) +
  scale_fill_etopo() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill="none") #to remove depth legend

base_noaa_2 <- autoplot(bathyData_2, geom=c("r")) +
  scale_fill_etopo() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill="none") #to remove depth legend

#Plot map including a wider region
base_noaa_arctic <- autoplot(bathyData_arctic, geom=c("r")) +
  scale_fill_etopo() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill="none") #to remove depth legend

#Plot zooming in the ncat stations
base_noaa_stations <- autoplot(bathyData_stations, geom=c("r")) +
  scale_fill_etopo() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill="none")

#--------------------------------- Tagging metadata -----------------------------------------------------#
#Read the tag metadata, in this case, it is an excel file
#***Change the file name accordingly
tag_metadata <- readxl::read_excel('ncat_tagging_metadata_2022-07.xlsx')

#Separate date to get an additional column with the year
tag_metadata_2 <- tag_metadata %>% separate(utc_release_date_time, sep = "-", into = c("year","month","day"))

#Add the catalognmumber column "NCAT + animal_id"
tag_metadata_2$catalognumber <- tag_metadata_2$animal_id
tag_metadata_2$catalognumber <- paste("NCAT-",tag_metadata_2$catalognumber,sep ="")

#Add releases from 2022 (at this time not included in the NCAT tagging metadata)
tag_metadata_2022 <- readxl::read_excel('NCAT_metadata_tagging_2022_04_lynx.xls', 
                                   sheet='Tag Metadata', skip=4)

#Edit file to have the same column names as NCAT tagging metadata
tag_metadata_2022_2 <- tag_metadata_2022 %>% separate(UTC_RELEASE_DATE_TIME, sep = "-", into = c("year","month","day"), remove = FALSE) #get a new column for the year
tag_metadata_2022_2 <- tag_metadata_2022_2[-1,] #remove sample data
colnames(tag_metadata_2022_2)[5] <- 'tag_serial_number'
colnames(tag_metadata_2022_2)[6] <- 'tag_id_code...6'#change column name, to be the same as the previous metadata file
colnames(tag_metadata_2022_2)[23] <- 'length'
colnames(tag_metadata_2022_2)[37] <- 'release_latitude' #change column name
colnames(tag_metadata_2022_2)[39] <- 'release_longitude'#change column name
colnames(tag_metadata_2022_2)[42] <- 'release_day'

tag_metadata_2022_2_subset <- tag_metadata_2022_2[c("tag_serial_number","tag_id_code...6", "length", "release_day","release_latitude","release_longitude","year")] #subset columns that are relevant

#Release date in the correct format
tag_metadata_2022_2_subset$release_day <- format(as.Date(tag_metadata_2022_2_subset$release_day),"%Y-%m-%d")

#Add the catalognumber column
tag_metadata_2022_2_subset$catalognumber <- tag_metadata_2022_2_subset$tag_serial_number
tag_metadata_2022_2_subset$catalognumber <- paste("NCAT-",tag_metadata_2022_2_subset$catalognumber,sep ="")
tag_metadata_2022_2_subset$catalognumber <- paste(tag_metadata_2022_2_subset$catalognumber,"-",tag_metadata_2022_2_subset$release_day,sep ="")

#The same but for releases from 2023 (at this time not included in the NCAT tagging metadata)
#For the genomics paper, these releases won't be included as there is no detection of these cod at the time we were preparing the MS
#tag_metadata_2023 <- readxl::read_excel('ncat_metadata_tagging_2023_QC.xlsx', 
#                                        sheet='Tag Metadata', skip=4)
#Edit file to have the same column names as ncat tagging metadata
#tag_metadata_2023_2 <- tag_metadata_2023 %>% mutate(UTC_RELEASE_DATE_TIME = as.Date(UTC_RELEASE_DATE_TIME))
#tag_metadata_2023_3 <- tag_metadata_2023_2 %>% separate(UTC_RELEASE_DATE_TIME, sep = "-", into = c("year","month","day","hour"), remove = FALSE) #get a new column for the year
#tag_metadata_2023_3 <- tag_metadata_2023_3[-1,] #remove sample data
#colnames(tag_metadata_2023_3)[5] <- 'tag_serial_number'
#colnames(tag_metadata_2023_3)[6] <- 'tag_id_code...6'#change column name, to be the same as the previous metadata file
#colnames(tag_metadata_2023_3)[22] <- 'length'
#colnames(tag_metadata_2023_3)[35] <- 'release_latitude' #change column name
#colnames(tag_metadata_2023_3)[36] <- 'release_longitude'#change column name
#colnames(tag_metadata_2023_3)[37] <- 'release_day'

#tag_metadata_2023_3_subset <- tag_metadata_2023_3[c("tag_serial_number","tag_id_code...6", "length", "release_day","release_latitude","release_longitude","year")] #subset columns that are relevant

#Release date in the correct format
#tag_metadata_2023_3_subset$release_day <- format(as.Date(tag_metadata_2023_3_subset$release_day),"%Y-%m-%d")

#Add the catalognumber column
#tag_metadata_2023_3_subset$catalognumber <- tag_metadata_2023_3_subset$tag_serial_number
#tag_metadata_2023_3_subset$catalognumber <- paste("NCAT-",tag_metadata_2023_3_subset$catalognumber,sep ="")
#tag_metadata_2023_3_subset$catalognumber <- paste(tag_metadata_2023_3_subset$catalognumber,"-",tag_metadata_2023_3_subset$release_day,sep ="")

#Include NLCOD releases of samples that were sequenced
#NLCOD data sent by Emilie
tag_metadata_NLCOD <- read.csv("qualified_detections/sequenced_NLCOD_transmitters.csv")

#Edit file to have the same column names as ncat tagging metadata
tag_metadata_NLCOD$Date <- as.POSIXct(tag_metadata_NLCOD$Date, format = "%m/%d/%Y", tz = "UTC")
tag_metadata_NLCOD_2 <- tag_metadata_NLCOD %>% separate(Date, sep = "-", into = c("year","month","day"), remove = FALSE) #get a new column for the year
columns_to_keep <- c("Transmitter_ID","ID.code","Date","year","Lat..dec.","Long..dec.","Length","SERIAL_N")
tag_metadata_NLCOD_3 <- tag_metadata_NLCOD_2[, names(tag_metadata_NLCOD_2) %in% columns_to_keep] #keep only certain columns
tag_metadata_NLCOD_3$Project <- "NLCOD"
colnames(tag_metadata_NLCOD_3)[3] <- 'tag_id_code...6'
colnames(tag_metadata_NLCOD_3)[8] <- 'length'
colnames(tag_metadata_NLCOD_3)[6] <- 'release_latitude' #change column name
colnames(tag_metadata_NLCOD_3)[7] <- 'release_longitude'#change column name
colnames(tag_metadata_NLCOD_3)[4] <- 'release_day'

#Make sure is the same format (DFO longitudes don not have negative sign)
tag_metadata_NLCOD_3$release_longitude <- -abs(tag_metadata_NLCOD_3$release_longitude)

#Create a catalognumber 
tag_metadata_NLCOD_3$catalognumber <- paste(tag_metadata_NLCOD_3$Project,tag_metadata_NLCOD_3$SERIAL_N,tag_metadata_NLCOD_3$release_day, sep = "-")

#Delete SERIAL-N
tag_metadata_NLCOD_4 <- tag_metadata_NLCOD_3[, !names(tag_metadata_NLCOD_3) %in% "SERIAL_N"]

#Use this function to merge data frames with different columns (same columns are the named the same)
rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

#Merged both files with the rbind.all.columns function
tag_metadata_w_2022 <- rbind.all.columns(tag_metadata_2022_2_subset, tag_metadata_2)
tag_metadata_w_NLCOD <- rbind.all.columns(tag_metadata_NLCOD_4, tag_metadata_w_2022)

#Make sure there are no duplicate rows
tag_metadata_w_NLCOD <- tag_metadata_w_NLCOD %>% distinct()

#Select only relevant columns
tag_metadata_w_NLCOD_subset <- tag_metadata_w_NLCOD[c("catalognumber","release_day","release_latitude","release_longitude","year","length")] #subset columns that are relevant

#Convert to sf file
tag_metadata_w_NLCOD_subset_2 <- st_as_sf(tag_metadata_w_NLCOD_subset, coords = c("release_longitude","release_latitude"), crs=4269, remove = FALSE)
tag_metadata_w_NLCOD_subset_3 <- st_transform(tag_metadata_w_NLCOD_subset_2, 4269)

#Add a column that states if was an inshore or offshore release
#First download a shapefile of the Canadian provinces from https://www.arcgis.com/home/item.html?id=dcbcdf86939548af81efbd2d732336db
canada_shapefile = st_read("Canada/Canada.shp")

#Rename Atlantic provinces where detections can be found
canada_shapefile_rename <- canada_shapefile %>%
  mutate(NAME = case_when(NAME %in% c("Newfoundland and Labrador","Quebec","New Brunswick", "Nova Scotia", "Prince Edward Island") ~ "Atlantic Canada",
                          TRUE ~ NAME))
#Subset to obtain a shapefile of the Atlantic coast
atlantic_coast <- subset(canada_shapefile_rename ,NAME == "Atlantic Canada")

#Turn it into a rough metre-length system by projection
atlantic_coast <- st_transform(atlantic_coast, 4269)

#Make a 46.3km buffer area (~25 nautical miles)
atlantic_coast_buff <- st_buffer(atlantic_coast, 46300)

#Intersect detection points
inshore <- lengths(st_intersects(tag_metadata_w_NLCOD_subset_3, atlantic_coast_buff)) > 0

#Create a vector based on the logical list
location_deploy <- ifelse(inshore == "TRUE", "inshore", "offshore")

#Add the new column to the data frame
tag_metadata_w_NLCOD_final <- cbind(tag_metadata_w_NLCOD_subset_3, location_deploy)

#Remove geometry column, change spatial data frame to tibble
tag_metadata_w_NLCOD_final_2 <- tag_metadata_w_NLCOD_final  %>%  as_tibble()  %>%  dplyr::select(-geometry)

#Subset only relevant columns to be use later
catalognumber_by_release_site <- tag_metadata_w_NLCOD_final_2 %>%
  subset(select=c("catalognumber","location_deploy"))
  
#Subset length by catalognumber
catalognumber_by_length <- tag_metadata_w_NLCOD_final_2 %>%
  subset(select=c("catalognumber","length"))

#Plot release sites by year
map_release_sites <- base_noaa_2 +
  geom_point(data=tag_metadata_w_NLCOD_final_2, aes(x=release_longitude, y=release_latitude, color=year), size=4) +
  scale_color_brewer(palette="YlOrRd") +
  theme(legend.text=element_text(size=13))

#Save map
ggsave(plot = map_release_sites, filename = "map_release_sites_by_year.jpeg", units="in", width=10, height=8)

#Count releases by year
tag_metadata_w_NLCOD_count <- tag_metadata_w_NLCOD_final_2 %>%
  group_by(year) %>%
  summarize(count =n())

#Save table
write.csv(tag_metadata_w_NLCOD_count, "tagging_count_per_year.csv")

#List of releases by year
#Example for 2019
list_releases_2019 <- tag_metadata_w_NLCOD_final_2 %>%
  filter(year == "2019")
list_catalognumber_releases_2019 <- list_releases_2019$catalognumber
#Save file
write.csv(list_catalognumber_releases_2019, "list_catalognumber_releases_2019.csv")

#--------------------------------- Stations sites---------------------------------------------------#
#Import deployment data
#***Change the file name accordingly
deploy_metadata <- read_csv('deployment_locations.csv')
#Remove any duplicate
deploy_metadata_no_dup <- deploy_metadata %>% distinct()
#Add a project column to differentiate from the DFO stations
deploy_metadata_no_dup$project <- 'NCAT'

#Plot deployment stations
map_deployments_ncat <- base_noaa_stations + #using the base map obtained from NOAA
  geom_point(data=deploy_metadata_no_dup, aes(x=Deployed_long, y=Deployed_lat), size=0.6)

#Save map
ggsave(plot = map_deployments_ncat, filename = "map_stations_ncat.jpeg", units="in", width=10, height=8)

#Add inshore DFO stations data
#***change file name accordingly
inshore_metadata <- read_csv('qualified_detections/receivers_2019-2024.csv') #most recent file

#Add receiver column
inshore_metadata <- inshore_metadata %>%
  mutate(receiver = paste(Receiver_Type,Receiver_ID, sep = "-"))

#Change column names
colnames(inshore_metadata)[14] <- "Deployed_lat"
colnames(inshore_metadata)[15] <- "Deployed_long"
colnames(inshore_metadata)[17] <- "Station"

#Keep only relevant (3) columns
inshore_metadata_subset <- subset(inshore_metadata, select = c("Station","Deployed_lat","Deployed_long"))
#Add a project column
inshore_metadata_subset$project <- 'DFO'

#Combine data frames (ncat stations plus inshore DFO stations)
stations_deployments <- rbind(deploy_metadata_no_dup,inshore_metadata_subset)

#Count unique stations
unique_stations_deployments <- unique(stations_deployments$Station)

#Count stations by project
stations_deployments_count_by_project <- stations_deployments %>%
  group_by(project) %>%
  summarize(count =n())

#Print stations by projects NCAT or DFO
stations_deployments_count_by_project

#Plot deployment stations
map_deployments <- base_noaa_2 +
  geom_point(data=stations_deployments, aes(x=Deployed_long, y=Deployed_lat), size=3)

#Save map
ggsave(plot = map_deployments, filename = "map_stations_including_dfoinshore.jpeg", units="in", width=10, height=8)

#new map: deployment + releases
map_deploy <- base_noaa_2 +
  geom_point(data=deploy_metadata_no_dup, aes(x=Deployed_long, y=Deployed_lat), size=2, color="#eeba30")

map_deploy_releases <- map_deploy +
  geom_point(data=tag_metadata_w_NLCOD_final_2, aes(x=release_longitude, y=release_latitude, color=year), size=4, shape=15) +
  scale_color_manual(values = c("#fa7e1e","#d62976","#962fbf","#30afc1"))
#Save map
ggsave(plot = map_deploy_releases, filename = "map_deployments_and_releases.pdf", units="in", width=10, height=8)

#--------------------------------- Detection data - preparation ---------------------------------------------------#
#Read all csv files
#***Add new data as data become available
proj_dets_2019 <- read.csv("ncat_matched_detections_2019.csv")
proj_dets_2020 <- read.csv("ncat_matched_detections_2020.csv")
proj_dets_2021 <- read.csv("ncat_matched_detections_2021.csv")
proj_dets_2022 <- read.csv("ncat_matched_detections_2022.csv")
proj_dets_2023 <- read.csv("ncat_matched_detections_2023.csv")
proj_dets_2024 <- read.csv("ncat_matched_detections_2024.csv")

#Combine detections from multiple years
#***Add more years as data become available
proj_dets <- rbind(proj_dets_2019,proj_dets_2020,proj_dets_2021,proj_dets_2022,proj_dets_2023,proj_dets_2024)

#Add additional inshore detections from NLCOD project (data sent by email, not in the OTN website)
#Merged both files with the rbind.all.columns function
inshore_dets <- read.csv("DFO_inshore_raw/inshore_detections_to_add.csv")
proj_dets_2 <- rbind.all.columns(proj_dets , inshore_dets)

#Remove duplicates, any duplicate entry
proj_dets_no_dup <- proj_dets_2 %>% distinct()

#Make sure the collection date is format correctly
proj_dets_no_dup %>% mutate(datecollected=ymd_hms(datecollected))

#Formatting and adding columns to help the following analyses
#Renaming some columns in the detection extract files
actel_dets_no_dup <- proj_dets_no_dup %>% dplyr::mutate(receiver_sn = as.integer(receiver),
                                          detection_timestamp_utc = datecollected,
                                          transmitter_codespace = extractCodeSpaces(tagname),
                                          transmitter_id = extractSignals(tagname))

#Add a column year-month
actel_dets_no_dup_2 <- actel_dets_no_dup %>% unite(year_month, c("yearcollected","monthcollected"),remove = FALSE)

#Add a column year-month-day
actel_dets_no_dup_3 <- actel_dets_no_dup_2 %>% unite(year_month_day, c("yearcollected","monthcollected","daycollected"),remove = FALSE) %>%
  mutate(year_month_day=ymd(year_month_day))

#To add season information to the data frame
#First create a function for the seasons based on the month number
get_season <- function(date) {
  month <- month(date)
  day <- day(date)
  
  if ((month == 3 && day >= 20) || month %in% c(4, 5) || (month == 6 && day < 21)) {
    return("Spring")
  } else if ((month == 6 && day >= 21) || month %in% c(7, 8) || (month == 9 && day < 22)) {
    return("Summer")
  } else if ((month == 9 && day >= 22) || month %in% c(10, 11) || (month == 12 && day < 21)) {
    return("Fall")
  } else {
    return("Winter")
  }
}

#Create the "seasons" column
actel_dets_no_dup_3$seasons = sapply(actel_dets_no_dup_3$year_month_day, get_season)

#Add the seasons column
actel_dets_no_dup_4 <- actel_dets_no_dup_3 %>% unite(season_year, c("seasons","yearcollected"),remove = FALSE)

#Add a nafo division column
#Download shapefile from: https://www.nafo.int/Data/GIS
#Open shapefile
divisions <- st_read("Divisions/NAFO_Divisions_SHP/NAFO_Divisions_2021_poly_clipped.shp")
st_crs(divisions)
#Create a column with a "point" longitude + latitude, this will delete those columns so first make a duplicate
actel_dets_no_dup_4$longitude1 = actel_dets_no_dup_4$longitude
actel_dets_no_dup_4$latitude1 = actel_dets_no_dup_4$latitude
actel_dets_no_dup_5 <- st_as_sf(actel_dets_no_dup_4, coords = c("longitude1","latitude1"), crs=4269)

#Make sure are in the same format
actel_dets_no_dup_6 <- st_transform(actel_dets_no_dup_5, 4269)

#Use st_join function to assign the division to each detection
actel_dets_no_dup_7 <- st_join(actel_dets_no_dup_6, divisions['Label'], join = st_intersects)

#***Location of Hant's Harbour appears outside the nafo divisions, thus is show as NA, in this case, replace NA to 3L:
actel_dets_no_dup_7$Label <- actel_dets_no_dup_7$Label %>% replace_na('3L')

#Change label header to nafo
actel_dets_no_dup_7 <- actel_dets_no_dup_7 %>%
  rename(nafo=Label)

#Add a column of inshore or offshore location
#First download a shapefile of the Canadian provinces from https://www.arcgis.com/home/item.html?id=dcbcdf86939548af81efbd2d732336db
canada_shapefile = st_read("Canada/Canada.shp")

#Rename Atlantic provinces where detections can be found
canada_shapefile_rename <- canada_shapefile %>%
  mutate(NAME = case_when(NAME %in% c("Newfoundland and Labrador","Quebec","New Brunswick", "Nova Scotia", "Prince Edward Island") ~ "Atlantic Canada",
                          TRUE ~ NAME))
#Subset to obtain a shapefile of the Atlantic coast
atlantic_coast <- subset(canada_shapefile_rename ,NAME == "Atlantic Canada")

#Turn it into a rough metre-length system by projection
atlantic_coast = st_transform(atlantic_coast, 4269)

#Make a 46.3km buffer area (~25 nautical miles)
atlantic_coast_buff = st_buffer(atlantic_coast, 46300)

#Intersect detection points
inshore <- lengths(st_intersects(actel_dets_no_dup_7, atlantic_coast_buff)) > 0

#Create a vector based on the logical list
location <- ifelse(inshore == "TRUE", "inshore", "offshore")

#Add the new column to the data frame
actel_dets_no_dup_8 <- cbind(actel_dets_no_dup_7, location)

#Remove geometry column, change spatial data frame to tibble
actel_dets_no_dup_9 <- actel_dets_no_dup_8  %>%  as_tibble()  %>%  dplyr::select(-geometry)

#Add deployment location "catalognumber_by_release_site"
actel_dets_no_dup_9$location_deploy = NA
actel_dets_no_dup_10<- dplyr::inner_join(actel_dets_no_dup_9, catalognumber_by_release_site, by="catalognumber")
actel_dets_no_dup_11 <- subset(actel_dets_no_dup_10, select=-location_deploy.x) 
colnames(actel_dets_no_dup_11)[45] = "locations_deploy"

#Add length based on "catalognumber_by_length"
actel_dets_no_dup_11$length = NA
actel_dets_no_dup_12 <- dplyr::inner_join(actel_dets_no_dup_11, catalognumber_by_length,by="catalognumber")
actel_dets_no_dup_13 <- subset(actel_dets_no_dup_12, select=-length.x)
colnames(actel_dets_no_dup_13)[46] = "length"
#Add a category based on length
actel_dets_no_dup_13$length_cat = NA
actel_dets_no_dup_13$length_cat <- ifelse(actel_dets_no_dup_13$length > 0.83, "large", ifelse(actel_dets_no_dup_13$length > 0.67, "medium", "small"))

#This is the final dataframe the will be used for further analyses
proj_dets_final <- actel_dets_no_dup_13

#Filter the release information
proj_dets_final_wo_release <- proj_dets_final %>% filter(receiver != 'release')

#Filter only release information
proj_dets_final_only_release <- proj_dets_final %>% filter(receiver == 'release')

#--------------------------------- Overall detection analyses and outliers identification ---------------------------------------------------#
#Detections from
print(min(proj_dets_final_wo_release$year_month_day))
#to
print(max(proj_dets_final_wo_release$year_month_day))
#Total number of detections
print(nrow(proj_dets_final_wo_release))

#from/to NLCOD
NLCOD_only <- proj_dets_final_wo_release %>% filter(detectedby == "NLCOD")
print(min(NLCOD_only$year_month_day))
print(max(NLCOD_only$year_month_day))

#Plot all detections, use the base_noaa_arctic
all_dets_plot <- base_noaa_arctic + geom_point(data = proj_dets_final_wo_release,
                                               mapping = aes(x = longitude, y = latitude), size =2)

#Save map of all detections
ggsave(plot = all_dets_plot, filename = "map_detections_all.jpeg", units="in", width=10, height=8)

#Count of detections by different projects (NCAT, ASF, etc)
proj_dets_count_by_projects <- proj_dets_final_wo_release %>% count(detectedby)

#Save file
write.csv(proj_dets_count_by_projects,"detections_count_by_project.csv")

#Check which individual was found too far north
individual_north <- proj_dets_final %>% filter(latitude >= 60)
individual_north$tagname
#and too far south
individual_south <- proj_dets_final %>% filter(latitude <= 45)
individual_south$tagname

#Check this individuals movement
#Filter first
individual_outliers <- proj_dets_final %>% 
  filter(grepl("A69-9001-10884|A69-9001-10933", tagname))
#Then plot 
abacus_plot_outliers <- individual_outliers %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))
#Then visualize
abacus_plot_outliers

#The individuals had not other important detections (i.e., multiple or in different months), so they will be eliminated
proj_dets_final_wo_outlier <- proj_dets_final %>% filter(!grepl("A69-9001-10884|A69-9001-10933", tagname))
print(nrow(proj_dets_final_wo_outlier))
#Remove outlier(s) in file without releases information
proj_dets_final_wo_release_wo_outlier <- proj_dets_final_wo_release %>% filter(!grepl("A69-9001-10884|A69-9001-10933", tagname))
print(nrow(proj_dets_final_wo_release_wo_outlier))

#Detect other possible outliers...individuals with unusual number of detections
abacus_plot_individuals_detection_days <- proj_dets_final_wo_outlier %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#Visualize
abacus_plot_individuals_detection_days

#Save abacus plot ***change file name accordingly
htmlwidgets::saveWidget(
  widget = abacus_plot_individuals_detection_days, #the plotly object
  file = "abacus_plot_by_indv_all.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)

#The previous abacus plot showed 6 individuals with multiple consecutive detections in the same site, a little suspicious, thus will look closer
proj_dets_possible_outlier <- proj_dets_final_wo_outlier %>% 
  filter(grepl("A69-9001-11275|A69-9001-11191|A69-9001-10355|A69-9001-10339|A69-9001-10360|A69-9001-11360", tagname))

#Plot again
abacus_plot_individuals_possible_outlier <- proj_dets_possible_outlier %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))
#Visualize
abacus_plot_individuals_possible_outlier

#Count detections of these 2 outliers
outlier2 <- proj_dets_final_wo_outlier %>% 
  filter(grepl("A69-9001-11191", tagname))
print(nrow(outlier2))

outlier3 <- proj_dets_final_wo_outlier %>% 
  filter(grepl("A69-9001-10355", tagname))
print(nrow(outlier3))

#Remove 2 more outliers from the main data_frames
proj_dets_final_wo_outlier <- proj_dets_final_wo_outlier %>% 
  filter(!grepl("A69-9001-11191|A69-9001-10355", tagname))

proj_dets_final_wo_release_wo_outlier <- proj_dets_final_wo_release_wo_outlier %>% 
  filter(!grepl("A69-9001-11191|A69-9001-10355", tagname))

#Plot NLCOD
dets_NLCOD_only <- proj_dets_final_wo_release %>% filter(collectioncode == "NLCOD")

#Plot again
abacus_plot_individuals_dets_NLCOD <- dets_NLCOD_only %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#Plot all detections without the outliers
all_dets_wo_outliers_plot <- base_noaa_2 + geom_point(data = proj_dets_final_wo_release_wo_outlier,
                                               mapping = aes(x = longitude, y = latitude), size =2)

#Save map of all detections
ggsave(plot = all_dets_wo_outliers_plot, filename = "map_detections_all_wo_outliers.jpeg", units="in", width=10, height=8)

#Number of detections left after removing outliers
print(nrow(proj_dets_final_wo_release_wo_outlier))

#Number of individuals detected overall after removing outliers
unique_tagnames <- proj_dets_final_wo_release_wo_outlier %>%
  group_by(tagname) %>% 
  summarize(count =n())

print(nrow(unique_tagnames))

#List of unique tagnames
list_unique_tagnames <- unique_tagnames$tagname

#Count of detections by different projects (NCAT, ASF, etc)
proj_dets_count_by_projects_wo_outliers <- proj_dets_final_wo_release_wo_outlier %>% count(detectedby)

#Save file
write.csv(proj_dets_count_by_projects_wo_outliers,"detections_count_by_project_wo_outliers.csv")

#How many individuals were detected only the day of release
#Count of detections by individual
dets_count_by_tagname <- proj_dets_final_wo_release_wo_outlier %>% 
  count(tagname)

#Number of individuals that have been detected at least once after release, including the same date of release
print(nrow(dets_count_by_tagname))

#Individuals detected only the day of release
#Summary of detections by tagname and day, use data frame that includes releases
dets_summary_count_by_tagname_day  <- proj_dets_final_wo_outlier  %>% 
  group_by(tagname, year_month_day) %>% 
  summarize(count =n())

#Total number of days detected by each individual
number_days_dets_by_indv <- dets_summary_count_by_tagname_day %>%
  group_by(tagname)  %>% 
  summarize(count =n())

#Filter out individuals only found in 1 day, which would be the release day
number_days_dets_by_indv_filtered <- number_days_dets_by_indv %>% filter(count > 1)

#Get number of individuals detected at least once after the date of release
print(nrow(number_days_dets_by_indv_filtered))

#Get list of tagnames
tagnames_list_filtered <- number_days_dets_by_indv_filtered$tagname

#Keep detections of only individuals detected the day after release
proj_dets_final_filtered <- proj_dets_final_wo_outlier %>% 
  filter(grepl(paste(tagnames_list_filtered, collapse = "|"), tagname))

proj_dets_final_wo_release_filtered <- proj_dets_final_wo_release_wo_outlier %>% 
  filter(grepl(paste(tagnames_list_filtered, collapse = "|"), tagname))

#Plot
all_dets_filtered_plot <- base_noaa_2 + geom_point(data = proj_dets_final_wo_release_filtered,
                                                      mapping = aes(x = longitude, y = latitude), size =2)

#Save map of all detections
ggsave(plot = all_dets_filtered_plot, filename = "map_detections_all_filtered.jpeg", units="in", width=10, height=8)

#Number of detections left after removing outliers
print(nrow(proj_dets_final_wo_release_filtered))

#Number of individuals detected overall after removing outliers
unique_tagnames <- proj_dets_final_wo_release_filtered %>%
  group_by(tagname) %>% 
  summarize(count =n())

print(nrow(unique_tagnames))

#List of unique tagnames
list_unique_tagnames <- unique_tagnames$tagname

#Count of detections by different projects (NCAT, ASF, etc)
proj_dets_count_by_projects_filtered <- proj_dets_final_wo_release_filtered %>% count(detectedby)

#Save file
write.csv(proj_dets_count_by_projects_filtered,"detections_count_by_project_filtered.csv")

#--------------------------------- Overall detections by year and month ------------------------------------------------------------#
#Bar plot of all detections by month per year
plot_detections_by_month_year <- proj_dets_final_wo_release_filtered %>% 
  mutate(year_month=ym(year_month)) %>%
  group_by(year_month) %>% #can group by station, species etc.
  summarize(count =n()) %>% #how many dets per year_month
  ggplot(aes(x = (month(year_month) %>% as.factor()), 
             y = count, 
             fill = (year(year_month) %>% as.factor())
  )
  )+ 
  geom_col(position = position_dodge2(preserve = "single"))+ #bars same width
  xlab("Month")+
  ylab("Total Detection Count")+
  ggtitle('Detections by Month')+ #title
  labs(fill = "Year")+
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

#Save plot
ggsave(plot = plot_detections_by_month_year, filename = "Detections_by_month_year.jpeg", units="in", width=10, height=8)

#Bar plot but by projects and stations
proj_dets_final_wo_release_wo_outlier_ncat <- proj_dets_final_wo_release_filtered %>%
  filter(grepl("NCAT", detectedby))
proj_dets_final_wo_release_wo_outlier_dfo <- proj_dets_final_wo_release_filtered %>%
  filter(grepl("NLCOD", detectedby))

#Plot NCAT stations
plot_detections_by_month_year_ncat <- proj_dets_final_wo_release_wo_outlier_ncat %>% 
  mutate(year_month=ym(year_month)) %>%
  group_by(year_month) %>% #can group by station, species etc.
  summarize(count =n()) %>% #how many dets per year_month
  ggplot(aes(x = (month(year_month) %>% as.factor()), 
             y = count, 
             fill = (year(year_month) %>% as.factor())
  )
  )+ 
  geom_col(position = position_dodge2(preserve = "single"))+ #bars same width
  xlab("Month")+
  ylab("Total Detection Count")+
  ggtitle('Detections by Month')+ #title
  labs(fill = "Year")+
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

#Save plot
ggsave(plot = plot_detections_by_month_year_ncat, filename = "Detections_by_month_year_ncat_stations.pdf", units="in", width=10, height=8)

#Now by DFO stations
plot_detections_by_month_year_dfo <- proj_dets_final_wo_release_wo_outlier_dfo %>% 
  mutate(year_month=ym(year_month)) %>%
  group_by(year_month) %>% #can group by station, species etc.
  summarize(count =n()) %>% #how many dets per year_month
  ggplot(aes(x = (month(year_month) %>% as.factor()), 
             y = count, 
             fill = (year(year_month) %>% as.factor())
  )
  )+ 
  geom_col(position = position_dodge2(preserve = "single"))+ #bars same width
  xlab("Month")+
  ylab("Total Detection Count")+
  ggtitle('Detections by Month')+ #title
  labs(fill = "Year")+
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

#Save plot
ggsave(plot = plot_detections_by_month_year_dfo, filename = "Detections_by_month_year_dfo_stations.pdf", units="in", width=10, height=8)

#Check those dfo detections in march-may
detects_march_may_dfo <- proj_dets_final_wo_release_wo_outlier_dfo %>%
  filter(!grepl("3|4|5", monthcollected))

detects_march_may_dfo_names <- detects_march_may_dfo %>%
  distinct(tagname)
#Plot
abacus_plot_dfo_march_may <- detects_march_may_dfo %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#--------------------------------- Overall detections by NAFO and season ---------------------------------------------------#

#Plot of detections based on nafo division
plot_detections_by_nafo <- proj_dets_final_wo_release_filtered%>% 
  group_by(nafo) %>% #can group by station, species etc.
  summarize(count =n()) %>% #how many dets per division
  ggplot(aes(x = nafo, 
             y = count))+ 
  geom_bar(stat = "identity", position = "dodge2")+ 
  xlab("NAFO division")+
  ylab("Total Detection Count")+
  ggtitle('Detections by NAFO division') #title

#Save plot
ggsave(plot = plot_detections_by_nafo, filename = "Detections_by_nafo.jpeg", units="in", width=10, height=8)

#Count detections by nafo division
count_detections_by_nafo <- proj_dets_final_wo_release_filtered %>%
  group_by(nafo) %>%
  summarize(count =n())

#Save file
write.csv(count_detections_by_nafo,"detections_count_by_nafo.csv")

#Count of individuals detected by month
count_indv_by_month <- proj_dets_final_wo_release_filtered %>%
  group_by(monthcollected, tagname) %>%
  summarize(count =n())

count_indv_by_month_2 <- count_indv_by_month %>%
  group_by(monthcollected) %>%
  summarize(count =n())

#Save file
write.csv(count_indv_by_month_2,"number_indv_detected_by_month.csv")

#Count of individuals detected by inshore/offshore
count_indv_by_shore <- proj_dets_final_filtered %>%
  group_by(location, tagname) %>%
  summarize(count =n())

count_indv_by_shore_2 <- count_indv_by_shore %>%
  group_by(location) %>%
  summarize(count =n())

#Save file
write.csv(count_indv_by_shore_2,"number_indv_detected_by_shore.csv")

#How many individuals were detected inshore and offshore
#Count the individuals based on their locations
result_shore <- proj_dets_final_filtered %>%
  group_by(tagname) %>%
  summarise(
    inshore_count = sum(location == "inshore"),
    offshore_count = sum(location == "offshore")
  ) %>%
  mutate(
    inshore_only = as.numeric(inshore_count >= 1 & offshore_count == 0),
    offshore_only = as.numeric(inshore_count == 0 & offshore_count >= 1),
    in_both = as.numeric(inshore_count >= 1 & offshore_count >= 1)
  )

# Get the counts of individuals in each category
count_inshore_only <- sum(result_shore$inshore_only)
count_offshore_only <- sum(result_shore$offshore_only)
count_in_both <- sum(result_shore$in_both)

write.csv(result_shore,"resutls_by_inshore_offshore.csv")

#Count of individuals detected by nafo division
count_indv_by_nafo <- proj_dets_final_wo_release_filtered %>%
  group_by(nafo, tagname) %>%
  summarize(count =n())

count_indv_by_nafo_2 <- count_indv_by_nafo %>%
  group_by(nafo) %>%
  summarize(count =n())

#Save file
write.csv(count_indv_by_nafo_2,"number_indv_detected_by_nafo.csv")

#Run in a loop to obtain detections by season
#First create a list of all year_month
uniquemonths = proj_dets_final_wo_release_wo_outlier$year_month
uniquemonths_list = unique(uniquemonths)

#Loop
for (m in uniquemonths_list){
  dets_by_month <- proj_dets_final_wo_release_wo_outlier %>% 
    filter(year_month==(m))
  
  map_dets_by_month <- 
    base_noaa + 
    ylab("Latitude") +
    xlab("Longitude") +
    geom_point(data = dets_by_month, #filtering for recent deployments
               aes(x = longitude, y = latitude), #specify the data
               shape = 19, size = 4) #lots of aesthetic options here!
  #save
  ggsave(plot = map_dets_by_month, filename = paste((m),".jpeg",sep = ""), units="in", width=15, height=8) 
}

#By season
dets_by_season <- proj_dets_final_wo_release_filtered

#Plot
map_dets_by_seasons <- 
  base_noaa_2 + 
  ylab("Latitude") +
  xlab("Longitude") +
  geom_point(data = dets_by_season, #filtering for recent deployments
             aes(x = longitude, y = latitude, color=seasons),
             position = position_jitter(h=0.08,w=0.08),
             shape = 19, alpha = 0.5, size = 3) +
  scale_color_manual(values=c("orange","darkgreen","yellow","blue"))+
  ggtitle('All detections') #***change year

#Save map
ggsave(plot = map_dets_by_seasons, filename = "detection_by_seasons.jpeg",units="in", width=10, height=6)

#By season and year
dets_by_season_by_year <- proj_dets_final_wo_release_filtered %>% 
  filter(yearcollected == "2023") #change year

list_seasons <- dets_by_season_by_year$seasons %>% unique()

#Plot
map_dets_by_season_by_year <- 
  base_noaa_2 + 
  ylab("Latitude") +
  xlab("Longitude") +
  geom_point(data = dets_by_season_by_year, #filtering for recent deployments
             aes(x = longitude, y = latitude, color=season_year),
             position = position_jitter(h=0.08,w=0.08),
             shape = 19, alpha = 0.5, size = 3) +
  scale_color_manual(values=c("orange","darkgreen","yellow","blue"))+
  ggtitle('Detections 2023') #***change year

#Save map
#***Change the year in the name
ggsave(plot = map_dets_by_season_by_year, filename = "detection_by_season_2023.jpeg",units="in", width=10, height=6)

#Plots by season
#By season and year
dets_by_season <- proj_dets_final_wo_release_filtered %>% 
  filter(seasons == "Fall") #change season

#Plot
map_dets_by_season <- 
  base_noaa_2 + 
  ylab("Latitude") +
  xlab("Longitude") +
  geom_point(data = dets_by_season, #filtering for recent deployments
             aes(x = longitude, y = latitude, color=seasons),
             position = position_jitter(h=0.08,w=0.08),
             shape = 19, alpha = 0.5, size = 1) +
  scale_color_manual(values=c("orange")) + #***change color
  ggtitle('Detections Fall') #***change year

#Save map
#***Change the year in the name
ggsave(plot = map_dets_by_season, filename = "detection_by_season_fall.pdf",units="in", width=10, height=6)

#Put all maps in one figure
combined_plot_seasons <- grid.arrange(map_dets_by_season_winter + labs(title = "b) Winter"), 
                              map_dets_by_season_spring + labs(title = "c) Spring"), 
                              map_dets_by_season_summer + labs(title = "d) Summer"), 
                              map_dets_by_season_fall + labs(title = "e) Fall"),
                              nrow = 2, widths = c(1, 1), heights = c(1, 1))

#Save the combined plot as an image file
ggsave("combined_plots_2.png", combined_plot_seasons, width = 10, height = 8)

#Find patterns by nafo and month
#Plot by each division by year with a loop
uniquenafo <- proj_dets_final_wo_release_wo_outlier$nafo
uniquenafo_list <- unique(uniquenafo)

for (n in uniquenafo_list){
  dets_by_nafo <- proj_dets_final_wo_release_wo_outlier  %>% 
    filter(nafo==(n)) #***Change division
  
  #Barplot
  plot_dets_by_nafo_month_year <- dets_by_nafo %>% 
    mutate(year_month=ym(year_month)) %>% #make datetime
    group_by(year_month) %>% #can group by station, species etc.
    summarize(count =n()) %>% #how many dets per year_month
    ggplot(aes(x = (month(year_month) %>% as.factor()), 
               y = count, 
               fill = (year(year_month) %>% as.factor())
    )
    )+ 
    geom_col(position = position_dodge2(preserve = "single"))+
    scale_fill_manual("legend", values = c("2019" = "coral2", "2020" = "green3", "2021" = "turquoise3", "2022" = "mediumpurple2", "2023"= "pink"))+
    xlab("Month")+
    ylab("Total Detection Count")+
    ggtitle(n)+ 
    labs(fill = "Year")+
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12))
  #Save plot
  #***Change division name
  ggsave(plot = plot_dets_by_nafo_month_year, filename = paste((n),"_detections_by_month.jpeg",sep = ""), units="in", width=10, height=8)
}

#--------------------------------- Efficiency of NCAT stations ---------------------------------------------------#
#Summary of detections by station and by project
count_detections_by_stations <- proj_dets_final_wo_release_wo_outlier  %>% 
  group_by(detectedby, station) %>% 
  summarize(count =n())

#Save file
write.csv(count_detections_by_stations,"detections_count_by_stations_by_project.csv")

#Summary of detections by station
NCAT_detections <- proj_dets_final_wo_release_wo_outlier  %>% 
  filter(detectedby == "NCAT")

count_detections_by_NCAT_stations <- NCAT_detections  %>% 
  group_by(station) %>% 
  summarize(count =n())

#Save file
write.csv(count_detections_by_NCAT_stations,"detections_count_by_NCAT_stations.csv")

#Combine count with NCAT stations data
stations_deployments$count <- count_detections_by_NCAT_stations$count[match(stations_deployments$Station,count_detections_by_NCAT_stations$station)]
#Filter only NCAT stations
stations_deployments_ncat_by_number_detections <- stations_deployments %>% filter(project == 'NCAT')

#Assign categories
stations_deployments_ncat_by_number_detections <- within(stations_deployments_ncat_by_number_detections, {   
  frequency_cat <- NA # need to initialize variable
  frequency_cat[count < 1] <- "zero"
  frequency_cat[count >= 5 & count < 100] <- "<100"
  frequency_cat[count >= 100 & count < 500] <- "<500"
  frequency_cat[count >= 500 & count < 1000] <- "<1000"
  frequency_cat[count >= 1000 & count < 2000] <- "<2000"
  frequency_cat[count >= 2000 & count < 4000] <- "<4000"
  frequency_cat[count >= 4000 & count < 8000] <- "<8000"
} )

#Order
stations_deployments_ncat_by_number_detections$frequency_cat <- factor(stations_deployments_ncat_by_number_detections$frequency_cat, levels = c("<8000", "<4000", "<2000", "<1000", "<500", "<100", "NA"))

#Plot
map_stations_by_number_detections <- base_noaa_stations + #using the base map obtained from NOAA
  geom_point(data=stations_deployments_ncat_by_number_detections, aes(x=Deployed_long, y=Deployed_lat, color=frequency_cat)) +
  scale_color_manual(values = c("#960000", "#9D290B", "#B05216", "#C47B21", "#D8A32B", "#EBCC36", "gray")) +
  guides(color=guide_legend("Number of total detections"))

#Save plot
ggsave(plot = map_stations_by_number_detections, filename = "map_deployments_count_number_detections.jpeg", units="in", width=10, height=8)

#Summary of detections by stations and tagname
count_detections_by_station_and_tagname  <- NCAT_detections  %>% 
  group_by(station, tagname) %>% 
  summarize(count =n())

#Count of individuals by station
count_detections_by_station_and_tagname_2 <- count_detections_by_station_and_tagname %>%
  group_by(station) %>%
  summarize(count=n())

#Save
write.csv(count_detections_by_station_and_tagname_2,"count_indv_detect_by_NCAT_stations.csv")

#Combine count with NCAT stations data
stations_deployments$count <- count_detections_by_station_and_tagname_2$count[match(stations_deployments$Station,count_detections_by_station_and_tagname_2$station)]
#Filter only NCAT stations
stations_deployments_ncat <- stations_deployments %>% filter(project == 'NCAT')

write.csv(stations_deployments_ncat,"indiv_count_by_ncat_stations.csv")

#Assign categories
stations_deployments_ncat <- within(stations_deployments_ncat, {   
  frequency_cat <- NA # need to initialize variable
  frequency_cat[count < 1] <- "NA"
  frequency_cat[count >= 1 & count < 5] <- "<5"
  frequency_cat[count >= 5 & count < 10] <- "<10"
  frequency_cat[count >= 10 & count < 20] <- "<20"
  frequency_cat[count >= 20 & count < 40] <- "<40"
  frequency_cat[count >= 40 & count < 60] <- "<60"
  frequency_cat[count >= 60 & count < 80] <- "<80"
  frequency_cat[count >= 80 & count < 100] <- "<100"
} )

#Order
stations_deployments_ncat$frequency_cat <- factor(stations_deployments_ncat$frequency_cat, levels = c("<100", "<80", "<60", "<40", "<20", "<10", "<5", "zero"))

stations_deployments_ncat_2 <-  stations_deployments_ncat %>%
  filter((count >= 20 & count < 200))

#Plot
map_deployments_count_indv <- base_noaa_stations + #using the base map obtained from NOAA
  geom_point(data=stations_deployments_ncat, aes(x=Deployed_long, y=Deployed_lat, colour=frequency_cat), size=1) +
  scale_color_manual(values = c("#960000", "#9D290B", "#B05216", "#C47B21", "#D8A32B", "#EBCC36", "#FFF541")) +
  ggrepel::geom_label_repel(data = stations_deployments_ncat_2,
                            mapping = aes(x = Deployed_long, y = Deployed_lat, label = count), min.segment.length = unit(0, 'lines') ) +
  guides(color=guide_legend("Number of cod detected"))

#Save plot
ggsave(plot = map_deployments_count_indv, filename = "map_deployments_count_ncatindv.jpeg", units="in", width=10, height=8)

#Get map of efficiency by month
#Obtain a list of different year_months
list_year_month <- proj_dets_final_wo_release_wo_outlier$year_month %>%
  unique()
#Get the dates in format so it can be sorted in chronological order
list_year_month_sorted_final <- c("2019_7","2019_8","2019_9","2019_10","2019_11","2019_12","2020_1","2020_2","2020_3","2020_4","2020_5","2020_6","2020_7","2020_8","2020_9","2020_10","2020_11","2020_12","2021_1","2021_2","2021_3","2021_4","2021_5","2021_6","2021_7","2021_8","2021_9","2021_10","2021_11","2021_12","2022_1","2022_2","2022_3","2022_4","2022_5","2022_6","2022_7","2022_8","2022_9","2022_10","2022_11","2022_12","2023_1","2023_2","2023_3","2023_4","2023_5","2023_6","2023_7","2023_8","2023_9","2023_10")
#remove 2019, no NCAT detections
list_year_month_sorted_final_2 <- list_year_month_sorted_final[!grepl("^2019", list_year_month_sorted_final)]

#Create a loop to plot detection of individuals by month
for (m in list_year_month_sorted_final_2) {
  NCAT_detecions_filter_by_month <- NCAT_detections  %>% 
    filter(year_month==(m))

  NCAT_count_detections_by_station_and_tagname  <- NCAT_detecions_filter_by_month %>% 
    group_by(station, tagname) %>% 
    summarize(count =n())

  NCAT_count_detections_by_station_and_tagname_2 <- NCAT_count_detections_by_station_and_tagname %>%
    group_by(station) %>%
    summarize(count=n())

#Combine count with NCAT stations data
  stations_deployments$count <- NCAT_count_detections_by_station_and_tagname_2$count[match(stations_deployments$Station,NCAT_count_detections_by_station_and_tagname_2$station)]
#Filter only NCAT stations
  stations_deployments_ncat <- stations_deployments %>% filter(project == 'NCAT')

#Filter stations with at least 1 indv detection
  stations_deployments_ncat_2 <- stations_deployments_ncat %>% filter(count != 'NA')

#Plot
  map_deployments_count_indv <- base_noaa_stations + #using the base map obtained from NOAA
    geom_point(data=stations_deployments_ncat_2, aes(x=Deployed_long, y=Deployed_lat)) +
    ggrepel::geom_label_repel(data = stations_deployments_ncat_2,
                            mapping = aes(x = Deployed_long, y = Deployed_lat, label = count), size = 10, min.segment.length = unit(0, 'lines'), max.overlaps = Inf )
#Save plot
  ggsave(plot = map_deployments_count_indv, filename = paste((m),"_map_deployments_count_ncatindv.jpeg", sep = ""), units="in", width=10, height=8)
}

#Get list of seasons
list_seasons_sorted_final <- c("Spring_2020","Summer_2020","Fall_2020","Winter_2021","Spring_2021","Summer_2021","Fall_2021","Winter_2022","Spring_2022","Summer_2022","Fall_2022","Winter_2023","Spring_2023","Summer_2023","Fall_2023")

#Create a loop to plot detection of individuals by seasons
for (m in list_seasons_sorted_final) {
  NCAT_detecions_filter_by_month <- NCAT_detections  %>% 
    filter(season_year==(m))
  
  NCAT_count_detections_by_station_and_tagname  <- NCAT_detecions_filter_by_month %>% 
    group_by(station, tagname) %>% 
    summarize(count =n())
  
  NCAT_count_detections_by_station_and_tagname_2 <- NCAT_count_detections_by_station_and_tagname %>%
    group_by(station) %>%
    summarize(count=n())
  
  #Combine count with NCAT stations data
  stations_deployments$count <- NCAT_count_detections_by_station_and_tagname_2$count[match(stations_deployments$Station,NCAT_count_detections_by_station_and_tagname_2$station)]
  #Filter only NCAT stations
  stations_deployments_ncat <- stations_deployments %>% filter(project == 'NCAT')
  
  #Filter stations with at least 1 indv detection
  stations_deployments_ncat_2 <- stations_deployments_ncat %>% filter(count != 'NA')
  
  #Plot
  map_deployments_count_indv <- base_noaa_stations + #using the base map obtained from NOAA
    geom_point(data=stations_deployments_ncat_2, aes(x=Deployed_long, y=Deployed_lat)) +
    ggrepel::geom_label_repel(data = stations_deployments_ncat_2,
                              mapping = aes(x = Deployed_long, y = Deployed_lat, label = count), min.segment.length = unit(0, 'lines'), max.overlaps = Inf )
  #Save plot
  ggsave(plot = map_deployments_count_indv, filename = paste((m),"_map_deployments_count_ncatindv.jpeg", sep = ""), units="in", width=10, height=8)
}

#--------------------------------- Analysis of individual survival ---------------------------------------------------#
#Count of detections by individual
dets_count_by_tagname <- proj_dets_final_wo_release_wo_outlier %>% 
  count(tagname)

#Number of individuals that have been detected at least once after release, including the same date of release
print(nrow(dets_count_by_tagname))

#Summary of detections by tagname and day, use data frame that includes releases
dets_summary_count_by_tagname_day  <- proj_dets_final_wo_outlier  %>% 
  group_by(tagname, year_month_day) %>% 
  summarize(count =n())

#Total number of days detected by each individual
number_days_dets_by_indv <- dets_summary_count_by_tagname_day %>%
  group_by(tagname)  %>% 
  summarize(count =n())

#Filter out individuals only found in 1 day, which would be the release day
number_days_dets_by_indv_filtered <- number_days_dets_by_indv %>% filter(count > 1)

#Get number of individuals detected at least once after the date of release
print(nrow(number_days_dets_by_indv_filtered))

#Get list of tagnames
tagnames_list_filtered <- number_days_dets_by_indv_filtered$tagname

#Save file
write.csv(number_days_dets_by_indv,"days_detect_by_indv.csv")

number_days_dets_by_indv_filtered_2 <- number_days_dets_by_indv_filtered %>% arrange(count)
#Get the mean
mean_count <- mean(number_days_dets_by_indv_filtered_2$count)
#Get the median
median_count <- median(number_days_dets_by_indv_filtered_2$count)

#Plot
plot_distirbution <- ggplot(data = number_days_dets_by_indv_filtered_2, aes(x = reorder(tagname, count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = paste("Mean:", round(mean_count, 2))), x = 0.1, y = max(number_days_dets_by_indv_filtered_2$count), hjust = 0, vjust = 1, color = "red") +
  geom_text(aes(label = paste("Median:", round(median_count, 2))), x = 0.1, y = max(number_days_dets_by_indv_filtered_2$count), hjust = 0, vjust = 2, color = "blue") +
  labs(title = "Distribution of Individual Counts", x = "Individuals", y = "Count") +
  theme_minimal()

#Save
ggsave(plot = plot_distirbution, filename = "plot_number_days_dets_by_indv.pdf", units="in", width=10, height=8)

##Survival after release
#Filter from data frame only individuals detected only after release date
dets_after_release_final_filtered <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_list_filtered) 

#Sort the dataframe by individual and detection_date
dets_after_release_final_filtered_2 <- dets_after_release_final_filtered %>% arrange(tagname, datecollected)

dets_by_indv_number_days_more_than_one_all <- dets_after_release_final_filtered_2 %>% 
  group_by(tagname) %>%
  mutate(days_between = c(NA, diff(datecollected))) %>%
  ungroup()

#Save csv file
write.csv(dets_by_indv_number_days_more_than_one_all,"day_indiv_detected_w_release_more_one_day_detection_all.csv")

#Obtain a summary of different days that were detected
dets_by_indv_number_days_more_than_one <- dets_after_release_final_filtered %>% 
  group_by(tagname) %>%
  summarise(num_detections = length(datecollected),
            start = min(datecollected),
            end = max(datecollected),
            det_days=length(unique(as.Date(datecollected))))

#Add number of days from release to last detection
dets_by_indv_number_days_more_than_one$date_diff <- as.numeric(as.Date(as.character(dets_by_indv_number_days_more_than_one$end)) - as.Date(as.character(dets_by_indv_number_days_more_than_one$start)))

#Save csv file
write.csv(dets_by_indv_number_days_more_than_one,"day_indiv_detected_w_release_more_one_day_detection.csv")

#Filter according last day seen since release:
detected_2days <- dets_by_indv_number_days_more_than_one %>% filter(date_diff <=2)
detected_week <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=3) %>% filter(date_diff <=7)
detected_2weeks <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=8) %>% filter(date_diff <=14)
detected_month <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=15) %>% filter(date_diff <=30)
detected_2months <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=31) %>% filter(date_diff <=60)
detected_6months <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=61) %>% filter(date_diff <=180)
detected_more6months <- dets_by_indv_number_days_more_than_one %>% filter(date_diff >=181)
detected_less1month <- dets_by_indv_number_days_more_than_one %>% filter(date_diff <=31)

#Get list of individuals from each category
tagnames_2days <- detected_2days$tagname
tagnames_week <- detected_week$tagname
tagnames_2weeks <- detected_2weeks$tagname
tagnames_month <- detected_month$tagname
tagnames_2months <- detected_2months$tagname
tagnames_6months <- detected_6months$tagname
tagnames_more6months <- detected_more6months$tagname

#Get full data by category
days_detected_2days <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_2days)
days_detected_week <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_week)
days_detected_2weeks <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_2weeks)
days_detected_month <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_month)
days_detected_2months <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_2months)
days_detected_6months <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_6months)
days_detected_more6months <- proj_dets_final_wo_outlier %>% filter(tagname %in% tagnames_more6months)

#Abacus plot showing detections 
#***change day_detected_xxx accordingly
abacus_plot_individuals_detection_2days <- days_detected_2days %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#Save abacus plot ***change file name accordingly
htmlwidgets::saveWidget(
  widget = abacus_plot_individuals_detection_2days, #the plotly object
  file = "abacus_plot_by_indv_2days.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)

##First and last detection per individual
#Create a new data product, det_days, that give you the unique dates that an animal was seen by a station
dets_by_station_by_indv_number_days <- proj_dets_final_wo_release_wo_outlier %>% 
  group_by(station) %>%
  summarise(num_detections = length(datecollected),
            start = min(datecollected),
            end = max(datecollected),
            uniqueIDs = unique(tagname), 
            det_days=length(unique(as.Date(datecollected))))

write.csv(dets_by_station_by_indv_number_days,"day_indiv_detected_by_station.csv")

#Including release
dets_by_indv_number_days <- proj_dets_final_wo_outlier %>% 
  group_by(tagname) %>%
  summarise(num_detections = length(year_month_day),
            start = min(year_month_day),
            end = max(year_month_day),
            det_days=length(unique(as.Date(datecollected))),
            days_between=as.numeric(end-start))
#Save
write.csv(dets_by_indv_number_days,"day_indiv_detected_w_release.csv")

#Plot
plot_all <- dets_by_indv_number_days %>%
  ggplot(aes(x=tagname))+
  geom_linerange(aes(ymin=start, ymax=end, x=tagname ),
                 size=1.5) +
  geom_point(aes(y=start, colour = "#CB5416")) +
  geom_point(aes(y=end, colour = "#267266")) +
  coord_flip() +
  ggtitle(label = "First and last detection from all releases") +
  ylab("year") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "none")

#save plot
ggsave(plot = plot_all, filename = "plot_first_last_detec_all.pdf", units="in", width=20, height=20)


#Get percentage of days in each division
percentage_days_in_each_nafo <- days_in_each_nafo %>%
  group_by(tagname) %>%
  mutate(percent = count/sum(count)) %>% 
  dplyr::select(-count) %>%
  spread(nafo, percent)

write.csv(percentage_days_in_each_nafo, "percentage_in_each_nafo_by_indv.csv")

#Keep only individuals detected in more than one division
detections_by_indv_by_nafo_multiple <- subset(detections_by_indv_by_nafo, duplicated(tagname) | duplicated(tagname, fromLast=TRUE))

#count how many indv were detected in more than 1 NAFO division
indv_detect_in_multiple_nafo <- detections_by_indv_by_nafo_multiple %>% count(tagname)

#Save
write.csv(indv_detect_in_multiple_nafo, "indv_detect_in_multiple_nafo.csv")

#Get a list of individuals found in multiple nafo divisions
indv_detect_in_multiple_nafo_list <- indv_detect_in_multiple_nafo$tagname

#Plot detection counts by nafo division of each indiv
plot_detection_indiv_by_nafo <- detections_by_indv_by_nafo_multiple %>% #how many dets per division
  ggplot(aes(y = fct_rev(tagname), 
             x = count, fill=nafo), color=nafo)+ 
  geom_bar(stat = "identity")+ 
  xlab("Total Detection Count")+
  ylab("Individuals")+
  guides(fill=guide_legend(title="NAFO division")) +
  ggtitle('Detections of individuals by NAFO division') #title

#Save plot
ggsave(plot = plot_detection_indiv_by_nafo, filename = "Barplot_detections_individuals_by_multiple_nafo.jpeg", units="in", width=10, height=8)

#Analyses inshore-offshore
#Counts of each individual inshore vs offshore
detections_by_indv_by_inshore_offshore <- proj_dets_final_wo_release_wo_outlier %>% 
  group_by(tagname, year_month_day, location) %>% #can group by station, species etc.
  summarize(count =n()) 

#Keep only individuals detected inshore (DFO)
days_inshore_offshore <- detections_by_indv_by_inshore_offshore %>%
  group_by(tagname, location) %>%
  summarize(count =n()) 

#Get percentage of days in each division
percentage_days_inshore_offshore <- days_inshore_offshore %>%
  group_by(tagname) %>%
  mutate(percent = count/sum(count)) %>% 
  dplyr::select(-count) %>%
  spread(location, percent)

write.csv(percentage_days_inshore_offshore, "percentage_inshore_offshore.csv")

#Plot percentage of days inshore vs offshore
percentage_days_inshore_offshore <- reshape2::melt(percentage_days_inshore_offshore, id.vars="tagname")

plot_percentage_days_inshore_offshore <- percentage_days_inshore_offshore %>%
  ggplot(aes(y = fct_rev(tagname), 
             x = value, fill=variable), color=variable) + 
  geom_bar(stat = "identity") + 
  xlab("Percentage of total days detected") +
  ylab("Individuals")+
  guides(fill=guide_legend(title="inshore/offshore")) +
  scale_fill_discrete(labels=c('Inshore', 'Offshore')) +
  ggtitle('Percentage inshore and offshore') #title

#Save plot
ggsave(plot = plot_percentage_days_inshore_offshore, filename = "Percentange of days spend inshore offshore.jpeg", units="in", width=10, height=8)

##Abacus plot individuals in multiple nafo divisions
abacus_plot_individuals_multiple_division <- ind_multiple_nafo %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, marker=list(color = ~latitude, colorscale="Viridis", showscale=TRUE)) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#Save plot
htmlwidgets::saveWidget(
  widget = abacus_plot_individuals_multiple_division, #the plotly object
  file = "abacus_plot_by_indv_in_multiple_division.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)

#Color by station
abacus_plot_individuals_multiple_days_by_nafo <- ind_multiple_nafo %>%
  plot_ly(x = ~detection_timestamp_utc, y = ~tagname, type = "scatter",  mode = "markers",text = ~station, color = ~nafo) %>%
  layout(
    xaxis = list(
      type = 'date',
      tickformat = "%B<br>%Y"
    ))

#Save plot
htmlwidgets::saveWidget(
  widget = abacus_plot_individuals_multiple_days_by_nafo, #the plotly object
  file = "abacus_plot_by_indv_in_multiple_days_by_nafo.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)

#--------------------------------- Analysis of individual movement ---------------------------------------------------#
##Plot each individual movement
#Example of individual with tagnanme: A69-9001-11393, as in glatos tutorial
#See below to do this in a loop!!
actel_dets_plot_indv <- proj_dets_final_wo_outlier %>% 
  filter(tagname=="A69-9001-11393")

actel_dets_plot_indv_map <- 
  base_noaa + 
  ylab("Latitude") +
  xlab("Longitude") +
  geom_point(data = actel_dets_plot_indv,
             aes(x = longitude, y = latitude, color=year_month), #color by year_month
             shape = 19, size = 4) #lots of aesthetic options here!

#Save map
ggsave(plot = actel_dets_plot_indv_map, filename = "A69-9001-11393.pdf", units="in", width=15, height=8) 

#For loop to get maps for each individual using our uniquenames_list
for (i in tagnames_list_filtered) {
  actel_dets_plot_indv <- proj_dets_final_wo_outlier %>% 
    filter(tagname==(i))
  
  actel_dets_plot_indv_map <- 
    base_noaa_2 + 
    ylab("Latitude") +
    xlab("Longitude") +
    geom_point(data = actel_dets_plot_indv, 
               aes(x = longitude, y = latitude, color=year_month),
               shape = 19, size = 4)
  
  #Save map
  ggsave(plot = actel_dets_plot_indv_map, filename = paste((i),".jpeg",sep = ""), units="in", width=15, height=8) 
  
}

#Try with lines connecting points
for (t in tagnames_list_filtered) {
  actel_dets_plot_indv <- proj_dets_final_wo_outlier %>% 
    filter(tagname==(t))
  
  actel_dets_plot_indv_map <- 
    base_noaa + 
    ylab("Latitude") +
    xlab("Longitude") +
    geom_point(data = actel_dets_plot_indv,
               aes(x = longitude, y = latitude, color=tagname), #color by year_month
               shape = 19, size = 4) +
    geom_path(data=actel_dets_plot_indv, aes(x = longitude, y = latitude, color=tagname))
  
  #Save map
  ggsave(plot = actel_dets_plot_indv_map, filename = paste((t),"_trace.jpeg",sep = ""), units="in", width=15, height=8) 
}

##Distance between points
ind_distribution_more_1_day <-  filter(proj_dets_final_wo_outlier, tagname %in% dets_by_indv_1_day_list,) #more than 1 day list

#Subset 4 columns: name, date, long, lat
ind_distribution <- subset(ind_distribution_more_1_day, select = c("tagname","year_month_day","latitude","longitude"))
#Make sure there are no duplicates, no same day
ind_distribution <- ind_distribution %>% distinct()
#Make sure lat and lon are numeric
ind_distribution$latitude <- as.numeric(ind_distribution$latitude)  
ind_distribution$longitude<- as.numeric(ind_distribution$longitude)

#Analysis by individual distance
for (t in dets_by_indv_1_day_list) {
  dets_by_indv <- ind_distribution %>% filter(tagname==(t)) #***change individual name 
  #Order dates
  dets_by_indv_order <- dets_by_indv %>% arrange(ymd(dets_by_indv$year_month_day))
  #Subset long and lat so dist Haversine function works
  dets_by_indv_subset <- subset(dets_by_indv_order, select = c("longitude","latitude"))
  #Distance matrix
  #pairwise_dist <- distm(dets_by_indv_subset, fun=distGeo)
  distance <- distHaversine(dets_by_indv_subset)
  #Add a space for the empty row
  distance <- append(distance,'NA',after=0)
  #Combine distance to the individual data frame
  dets_by_indv_order_w_distance <- cbind(dets_by_indv_order, distance) #***change individual name
  #Save dataframe
  write.csv(dets_by_indv_order_w_distance, gsub('.csv','_distance.csv',t))
}

distance_files <- dir(path = "~/Cod/OTN/distance/", full.names = TRUE, recursive = TRUE)
all_distance_files <- ldply(as.list(files), read.csv)
write.csv(all_distance_files, "all_distance_files.csv")

##Map with animation using glatos tutorial
#Create detections event variable
#When using data frame all are considered potentially false, when opening a csv file with all data, only 2% is considered false.
detection_events <- 
  read_otn_detections("merged.csv") %>% # a combine csv file with all detections
  false_detections(tf = 3600) %>%  #find false detections
  filter(passed_filter != FALSE) %>% 
  detection_events(location_col = 'station', time_sep=86400) #one day in seconds
#2.36% potentially false

plot_data <- detection_events %>% 
  dplyr::select(animal_id, mean_longitude,mean_latitude, first_detection)

#NCAT-1282963-2019-07-09, NCAT-1287714-2019-07-23
#Example with one fish based on animal_id name
one_fish <- plot_data[plot_data$animal_id == "NCAT-1287714-2019-07-23",] 

#left = min(plot_data$mean_longitude), 
#bottom = min(plot_data$mean_latitude), 
#right = max(plot_data$mean_longitude), 
#top = max(plot_data$mean_latitude)),
register_stadiamaps(key = "83e7ad95-4d0d-4d0b-8158-b081ac359672", write = TRUE)
#Get a new base map 
basemap <- get_stadiamap(
  bbox = c(left = -60, 
           bottom = 46, 
           right = -46, 
           top = 55),
  maptype = "stamen_terrain_background",
  crop = FALSE, 
  zoom = 6)

ggmap(basemap)

#Plot static map
act.plot <-
  ggmap(basemap) +
  geom_point(data = one_fish, aes(x = mean_longitude, y = mean_latitude, group = animal_id, color = animal_id), size = 2) +
  geom_path(data = one_fish, aes(x = mean_longitude, y = mean_latitude, group = animal_id, color = animal_id)) +
  labs(title = "ACT animation",
       x = "Longitude", y = "Latitude", color = "Tag ID")

ggplotly(act.plot)

#Animate!
act.plot <-
  act.plot +
  labs(subtitle = 'Date: {format(frame_along, "%d %b %Y")}') +
  transition_reveal(first_detection) +
  shadow_mark(past = TRUE, future = FALSE) +
  theme(plot.subtitle = element_text(size = 15, face = "bold"))

gganimate::animate(act.plot, duration = 12)

#Save animation as a gif
anim_save("NCAT-1287714-2019-07-23.gif", animation = act.plot)

##Plot multiple individuals in one animation
#Get a list of individuals animal_id you want to combine
multiple_fish <-  read_lines("catalog_ids.txt")
multiple_fish_list <- as.list(multiple_fish)

for (m in multiple_fish_list) {
  multiple_fish_plot <- plot_data[plot_data$animal_id == (m),] 
  
  #plot static map
  act.plot_multiple <-
    ggmap(basemap) + #make sure you have get basemap
    geom_point(data = multiple_fish_plot, aes(x = mean_longitude, y = mean_latitude, group = animal_id, color = animal_id), size = 2) +
    geom_path(data = multiple_fish_plot, aes(x = mean_longitude, y = mean_latitude, group = animal_id, color = animal_id)) +
    labs(title = "ACT animation", x = "Longitude", y = "Latitude", color = "Tag ID")
  
  ggplotly(act.plot_multiple)
  
  #animate!
  act.plot_multiple <- act.plot_multiple +
    labs(subtitle = 'Date: {format(frame_along, "%d %b %Y")}') +
    transition_reveal(first_detection) +
    shadow_mark(past = TRUE, future = FALSE) +
    theme(legend.position = "none")
  
  gganimate::animate(act.plot_multiple)
}

#save
anim_save(animation=act.plot_multiple, filename="multiple_individuals_narrow.gif")