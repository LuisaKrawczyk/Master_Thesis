# Download weather data ERA5

# I. Preparation----
rm(list = ls())
packages <- c("ecmwfr","raster","rlist", "purrr","ncdf4")
lapply(packages, require, character.only = TRUE)

library(ncdf4)
library(ecmwfr)
install.packages("ecmwfr")
library(raster)
install.packages("raster")
library(rlist)
install.packages("rlist")
library(purrr)
#install.packages("purrr")


# Enter user ID. Get User ID and API key registering at CDS
# @Luisa- du musst dich bei CDS 
# (https://cds.climate.copernicus.eu/user/register?destination=%2F%23!%2Fhome) registrieren und dann
# hier deine User-ID und deinen API-key angeben.

user_id <- "95450"
API_key <- "c1059400-14f8-4a4b-85c4-f537010f64d8"
#local_storage <- "../data/weather_era5/" # hier deinen Speicherort angeben
local_storage <- "E:/Wetterdaten/" # hier deinen Speicherort angeben


# set a key to the keychain
wf_set_key(user = user_id, key = API_key, service = "cds")

# variable selection
# @Luisa du kannst dir die verfÃ¼gbaren Daten unter folgenden Links angucken: 
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
# Ich denke, wir kÃ¶nnen die Auswahl hier aber einfach beibehalten.
variables <- c('10m_u_component_of_wind', '10m_v_component_of_wind', 
               '2m_temperature', 'relative_humidity','total_precipitation',
               'total_cloud_cover',  'surface_pressure')
nicht_mehr <- c("ozone_mass_mixing_ratio",'downward_uv_radiation_at_the_surface', 
                'vertical_velocity' )

# II. Download function----
fx <- function(yr, i){
  var <- variables[i]
  file_on_disk <- paste0(yr,"_",variables[i],".nc")
  request <- list(
    "product_type" = "reanalysis","variable" = var,  "year" = yr,
    "month" = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
    "day" = paste(c("01","02","03", "04", "05", "06", "07", "08", "09", 10:31)), 
    "time" = c('00:00', '01:00', '02:00',
               '03:00', '04:00', '05:00',
               '06:00', '07:00', '08:00',
               '09:00', '10:00', '11:00',
               '12:00', '13:00', '14:00',
               '15:00', '16:00', '17:00',
               '18:00', '19:00', '20:00',
               '21:00', '22:00', '23:00'),
    "format" = "netcdf",
    "area" = "49/-125/23/-60", # North, West, South, East 
                           # @Luisa - hier die Koordinaten für die USA eintragen
    "target" = file_on_disk)
  
  if (var %in% c( 'ozone_mass_mixing_ratio', 'relative_humidity',
                  'specific_humidity', 'vertical_velocity')) {
    request <- list.append(request,"dataset_short_name" = "reanalysis-era5-pressure-levels",
                           "pressure_level" = "1000")
  } else {
    request <- list.append(request,"dataset_short_name" = "reanalysis-era5-single-levels")       
  }
  
  wf_request(user     = user_id,        # user ID (for authentification)
             request  = request,        # the request
             transfer = TRUE,           # download the file
             path     = local_storage)  # set download path
  
  
  file_from_disk <- raster::brick(paste0(local_storage,file_on_disk))
}

# III. Download----

# @Luisa - die Verbindung zum Server lÃ¤uft (zumindest bei mir) nicht immer stabil. Es lohnt sich das
# im Blick zu behalten und den Download ggf. neu zu starten. Wenn die Jahresdateien zu groß sind,
# sollten wir kleinteiliger (je Monat) downloaden.

fx("2019",4)
fx("2019",5)
fx("2019",6)
fx("2019",7)

for(i in c(1:length(variables))){
map(dates, fx, i = 1)
}

for(i in c(1:length(variables))){
  map(yr, fx, i = 1)
}

for(i in c(1:length(variables))){
  fx("2018", variables[i])
}
for(i in c(1:length(variables))){
  fx("2021", variables[i])
}
for(i in c(4:length(variables))){
  map("2021", fx, i)
}
fx("2019",1)
fx("2019",2)
for(i in c(1:length(variables))){
  fx("2020", variables[i])
}
