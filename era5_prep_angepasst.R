#===================================================================================================
# Description: This script files performs a series of tasks related to weather data. 
# 0- Preparation
# 1- Prepare population raster data   (@Luisa - diesen Teil nicht ausf√ºhren!)
# 2- Create transformation "P" matrix (@Luisa - diesen Teil nicht ausf√ºhren!)
# 3- Split files by month (~ 2.6h)
# 4- Create temperature, precipitation, and heat, exposure bins-(~ ca. 4h)
# 5- Calculate heat index for every day (~ ca. 15min)
# 6- Aggregate to county level (~ ca. 15min)
# 7- Generate final dataset
#===================================================================================================

# 0. Preparation------------------------------------------------------------------------------------

# Clean up workspace
rm(list=ls())
gc()
getwd()
setwd("E:/Luisa Master Thesis/")

# Clean and load packages
wants <- c("rgdal","sp","raster","ncdf4","data.table","dplyr","sf","fasterize",
           "dplyr","plyr","purrr", "exactextractr", "Matrix", "lubridate","furrr")
needs <- wants[!(wants %in% installed.packages()[,"Package"])]
if(length(needs)) install.packages(needs)
lapply(wants, function(i) require(i, character.only=TRUE))
rm(needs,wants)

# raster template 
weather_ras  <- raster("./data/weather_era5/era5/2m_temperature/2018_2m_temperature.nc")

mask <- weather_ras
mask[] <- 1
plot(mask)
id <- mask
id[] <- 1:length(id)
plot(id)
names(id) <- names(mask) <- "cellid"

# 1. Prepare population raster data (nicht durchf√ºhren!)--------------------------------------------

# load nc raster weather raster file and translate into polygon
weather_ras  <- raster("./data/weather_era5/era5_land/2m_temperature/2020_2m_temperature.nc",
                       level = 1) 

mask <- weather_ras
mask[] <- 1
plot(mask)
id <- mask
id[] <- 1:length(id)
plot(id)
names(id) <- names(mask) <- "cellid"

grid <- st_as_sf(rasterToPolygons(id))
grid <- st_set_crs(grid, 4326)
grid <- st_transform(grid , st_crs(4326))

# function extracts the population counts per era5-polygon
ras      <- raster(paste0("./data/population/gpw_v4_population_count_rev11_2020_30_sec.tif"))
ras      <- crop(ras, extent(weather_ras))
crs(ras) <- CRS("+init=epsg:4326")
tmp      <- exact_extract(ras, grid, fun = "sum")
pop      <- as.data.frame(tmp) 

colnames(pop) <- paste0("pop")
pop <- bind_cols(grid, pop)

k <- pop %>% mutate(pop = ifelse(pop > 10000, 10000, pop))
plot(k[, "pop"], border = NA)

save(pop, file = "./data/population/population_era5_grid.RData")

# 2. Create transformation "P" matrix to aggregate raster to zip areas- (nicht durchf√ºhren!)--------

# load population data and raster template
load("./data/population/population_era5_grid.RData")
r  <- raster("./data/weather_era5/era5_land/2m_temperature/2020_2m_temperature.nc",  level = 1) 

# load county shape file 
load("./data/county_shape/us_counties.RData")

# add grid cell id for cells over land
mask <- r
mask[] <- 1 #0
mask[!is.na(r[])] <- 1
plot(mask)
id <- mask
id[] <- 1:length(id)
id[mask[]==0] <- NA
plot(id)
names(id) <- names(mask) <- "cellid"

# rasterize population data
pop       <- fasterize(pop, r, field = "pop")
pop       <- list(pop)
pop[[2]]  <- id
pop       <- brick(pop)
#pop [is.na(pop [])] <- 0 

# change crs
st_crs(shape)
crs(pop) <- CRS("+init=epsg:4326")

# extract grid cells per zip area
county_pop  <- exact_extract(pop, shape, .parallel = TRUE)

# generate weight by coverage fraction and population
fx <- function(x){
  tmp           <- county_pop[[x]]         
  tmp[,c(1)]    <- tmp[,c(1)]*tmp[,3]
  tmp           <- as.data.frame(tmp[,c(1)])
  tmp           <- t(t(tmp) / colSums(tmp, na.rm = TRUE))
  tmp           <- as.data.frame(tmp) %>% mutate(ID = county_pop[[x]][,2])
  colnames(tmp) <- c("pop", "cellid")
  tmp$county    <- shape$NAME_2[[x]]
  tmp$GID_2     <- shape$GID_2[[x]]
  tmp$state     <- shape$NAME_1[[x]]
  tmp$GID_1     <- shape$GID_1[[x]]
  tmp$county_id <- x
  tmp
}

county_pop <- bind_rows(map(c(1:length(county_pop)), fx))

# create projection matrix

# Say you have a
# A <- P %*% G

# where:
# G is a matrix with all the climate data, obtained as G <- g[], where g is a raster stack
# M has a dimension n x t, where n (190008) is the number of grid cells in the climate data, 
#   and t (24*365*2) the number of time periods or raster layers
# A is a matrix with N (3117) rows (number of zip) and t (24*365*2) columns (number of layers in g)
# P is our weight/projection matrix we need to construct
# P dimensions? n (190008) rows (number of grid cells) by N (3117) columns (number of county areas)

# Let's start
g <- stack(r) 

# Create P matrix
P <- sparseMatrix(j=county_pop$county_id,
                  i=county_pop$cellid,
                  x=county_pop$pop,
                  dims=c(ncell(g), length(unique(county_pop$GID_2))))

colnames(P) <- as.character(shape$GID_2) 

P_matrix <- P

# Double check that all columns sum to 1
print(unique(colSums(P_matrix)))


# Test
G <- g[]
A <- t(P_matrix) %*% G
dim(A)
head(A)
A <- as.data.frame(as.matrix(A))
A$GID_2 = row.names(A)

# Plot 
shp <- left_join(shape, A, by = "GID_2")
plot(shp[, "X2.metre.temperature"])
plot(r)

# Save projection matrix to disk
save(P_matrix, file = "./data/population/P_matrix.RData")




# 3. Split files by months--------------------------------------------------------------------------
# Luisa: ".nc" files are transformed to RMD files in this step 

f_split <- function(f){
tmp_full <- nc_open(paste0("./data/weather_era5/era5/", f))
time <- as_datetime(c(tmp_full$var[[1]]$dim[[3]]$val*60*60), origin="1900-01-01",  tz = "UTC")
tmp_full <- ncvar_get(tmp_full, tmp_full$var[[1]]$name)
fm <- function(m){
  idx <- which(substr(time, 6, 7) == m)
  tmp <- t(brick(tmp_full[, , idx], xmn=extent(mask)[3] , xmx=extent(mask)[4] , 
                  ymn=extent(mask)[1] , ymx=extent(mask)[2] , crs=CRS("+init=epsg:4326")))
  names(tmp) <- time[which(substr(time, 6, 7) == m)]
  crs(tmp) <- "+proj=longlat +datum=WGS84"
#saveRDS(tmp, paste0("./data/weather_era5/month/", substr(f, 1, nchar(f)-3), "_", m, ".RDS"))
writeRaster(tmp, paste0("./data/weather_era5/month/era5/", substr(f, 1, nchar(f)-3), "_", m, ".grd"))
print(paste0(m,"---", f)) }
map(c("01","02","03","04", "05", "06", "07", "08", "09", 10:12), fm)
}

path = "./data/weather_era5/era5/"
fileslist <- list.files(path, recursive = TRUE)
#files <- fileslist[1:27]

for(i in fileslist[1]){f_split(i)}
#files = c("era5/2m_temperature/2020_2m_temperature.nc")
# Use this when list.files param recursive is not working
# actual_files <- c()
#for(var in fileslist){actual_files = append(actual_files,paste0(var,"/",list.files(paste0(path,var))))}


# 4. Calculate max/min temperature, averages for every day-----------------------------------------

func <- function(m, var){
  
  # load first day of following month
  f <- list.files(paste0("./data/weather_era5/month/era5/", var , "/"))
  f <- f[grepl(".grd", f, fixed =T)]
  n <- length(f)
#  raster_brick <-readRDS(paste0("./data/weather_era5/month/era5/", var , "/", f[[m]] ))
  raster_brick <- stack(paste0("./data/weather_era5/month/era5/", var , "/", f[[m]] )) 
  if(m < n){
  next_day <- stack(paste0("./data/weather_era5/month/era5/", var , "/", f[[m+1]]))
  next_day <- next_day[[c(1:24)]]
  
  raster_brick <- stack(raster_brick, next_day)
  }
  
  # adjust time zones
  time  <- substr(names(raster_brick), 2, nchar(names(raster_brick)))
  time  <- as_datetime(time , tz = "UTC")
  
  # extract maximum per day
  f_max <- function(y, zone){
    rb         <- raster_brick
    names(rb)  <- as_datetime(time , tz = zone)
    t          <- names(rb)[grepl(y, names(rb),  fixed = TRUE)]
    tmp        <- subset(rb, t)
    tmp        <- calc(tmp, function(x){max(x, na.rm = TRUE)})
    names(tmp) <- y
    print(y)
    tmp
  }
  # extract sum per day
  f_sum <- function(y, zone){
    rb         <- raster_brick
    names(rb)  <- as_datetime(time , tz = zone)
    t       <- names(rb)[grepl(y,names(rb),  fixed = TRUE)]
    tmp     <- subset(rb, t)
    tmp     <- calc(tmp, function(x){sum(x, na.rm = TRUE)})
    names(tmp) <- y
    print(y)
    tmp
  }
  # extract mean per day
  f_avg <- function(y, zone){
    rb         <- raster_brick
    names(rb)  <- as_datetime(time , tz = zone)
    t          <- names(rb)[grepl(y,names(rb),  fixed = TRUE)]
    tmp        <- subset(rb, t)
    tmp        <- calc(tmp, function(x){mean(x, na.rm = TRUE)})
    names(tmp) <- y
    print(y)
    tmp
  }
  # extract duration per day during daylight ?
  f_dur <- function(y, zone){
    rb              <- raster_brick
    names(rb)       <- as_datetime(time , tz = zone)
    t               <- names(rb)[grepl(y,names(rb),  fixed = TRUE)]
    tmp             <- subset(raster_brick, t)
    tmp[tmp > 0]    <- 1
    tmp             <- calc(tmp, function(x){sum(x, na.rm = TRUE)})
    names(tmp) <- y
    print(y)
    tmp
  }
  # temp bin
  f_bin <- function(y, zone){
    rb         <- raster_brick
    names(rb)  <- as_datetime(time , tz = zone)
    t          <- names(rb)[grepl(y, names(rb),  fixed = TRUE)]
    tmp        <- subset(rb, t)
    bins     <- seq(-10, 40, 5)        # bins to that will be saved
    binsbase <- seq(-15, 50, 5)        # base bins over which computation are done
    nbins    <- length(binsbase) - 1L  # no of bins
    
    # compute bin counts
    ras <- as.array(tmp)
    M <- apply(ras, 3, FUN = function(x){findInterval(x, binsbase, rightmost.closed=TRUE)})
    M <- t(apply(M, 1, tabulate, nbins = nbins))
    
    # aggregate extreme bins
    right  <- rowSums(M[, which(match(binsbase, max(bins)) != "NA"):ncol(M)])
    left   <- rowSums(M[, 1:which(match(binsbase, min(bins)) != "NA")])
    mid    <- M[,(which(match(binsbase, min(bins)) != "NA")+1):(which(match(binsbase, max(bins)) != "NA")-1)]
    
    # export
    mat <- cbind(left,mid,right)
    colnames(mat) <- bins
    
    # Store bins in a raster file
    mout <- lapply(split(seq_len(ncol(mat)),(seq_len(ncol(mat))-1) %/%1 +1),
                   function(i) matrix(mat[,i], nrow=dim(mask)[[1]], byrow = FALSE))
    
    mout  <- array(unlist(mout),  dim = c(105, 261, length(mout)))
    
    r <- brick(mout, xmn=extent(mask)[1] , xmx=extent(mask)[2] , 
               ymn=extent(mask)[3] , ymx=extent(mask)[4] , crs=CRS("+init=epsg:4326"))
    
    names(r) <-   paste0(y,"_", gsub("-","m",paste("bin_",bins,sep="")))
    
    print(y)
    r
  }
  
  if(var %in% c("2m_temperature")){
    raster_brick <- raster_brick -273.15
    fz <- function(z){
#    day_max  <- map(unique(substr(names(raster_brick), 2, 11)), f_max, zone = z)
    }
 #   day_max <- map(c("US/Pacific", "US/Mountain", "US/Central", "US/Eastern"), fz)
    fz <- function(z){
      dm <- brick(day_max[[z]])
      dm <- dm[[c(1:(dim(dm)[[3]]-1))]]
    }
#    day_max <- map(c(1:4), fz)
 #   saveRDS(day_max, paste0("./data/weather_era5/day/", var,"/max_",   f[[m]]))
    
    ####
    fz <- function(z){
      day_bin  <- map(unique(substr(names(raster_brick), 2, 11)), f_bin, zone = z)
    }
    day_bin <- map(c("US/Pacific", "US/Mountain", "US/Central", "US/Eastern"), fz)
    
    fz <- function(z){
      dm <- brick(day_bin[[z]])
      dm <- dm[[c(1:(dim(dm)[[3]]-11))]]
    }
    day_bin <- map(c(1:4), fz)
    saveRDS(day_bin, paste0("./data/weather_era5/day/", var,"/bin_",   f[[m]]))
  }
  if(var %in% c("total_precipitation")){
    # sum of precipitation
    fz <- function(z){
      day_sum  <- map(unique(substr(names(raster_brick), 2, 11)), f_sum, zone = z)
    }
    day_sum <- map(c("US/Pacific", "US/Mountain", "US/Central", "US/Eastern"), fz)
    fz <- function(z){
      dm <- brick(day_sum[[z]])
      dm <- dm[[c(1:(dim(dm)[[3]]-1))]]
    }
    day_sum <- map(c(1:4), fz)
    saveRDS(day_sum, paste0("./data/weather_era5/day/", var,"/sum_",  f[[m]]))
    
    # precipitation hours
    fz <- function(z){
      day_dur  <- map(unique(substr(names(raster_brick), 2, 11)), f_dur, zone = z)
    }
    day_dur <- map(c("US/Pacific", "US/Mountain", "US/Central", "US/Eastern"), fz)
    
    fz <- function(z){
      dm <- brick(day_dur[[z]])
      dm <- dm[[c(1:(dim(dm)[[3]]-1))]]
    }
    
    day_dur <- map(c(1:4), fz)
    saveRDS(day_dur, paste0("./data/weather_era5/day/", var,"/hours_",  f[[m]]))
  }
  if(var %in% c('10m_u_component_of_wind','10m_v_component_of_wind', 'relative_humidity',
                'surface_pressure', 'total_cloud_cover')){
    fz <- function(z){
      day_mean  <- map(unique(substr(names(raster_brick), 2, 11)), f_avg , zone = z)
    }
    day_mean  <- map(c("US/Pacific", "US/Mountain", "US/Central", "US/Eastern"), fz)
    
    fz <- function(z){
      dm <- brick(day_mean[[z]])
      dm <- dm[[c(1:(dim(dm)[[3]]-1))]]
    }
    day_mean  <- map(c(1:4), fz)
    saveRDS(day_mean, paste0("./data/weather_era5/day/", var,"/mean_",  substr(f[[m]], 1, (nchar(f[[m]])-3)), "RDS"))
  }
}

n <- length(list.files(paste0("./data/weather_era5/month/era5/2m_temperature/")))

for(i in c('total_cloud_cover','total_precipitation')){
  map(c(1:36), func, var = i)
}

# st¸ckchenweise: 
i="2m_temperature"
map(c(1:36),func,var=i)
func(27,i)
func(29,i)


map(c(7:10), func, var = i)
map(c(13:24), func, var = i)

# output f¸r temperature: 
# max_Monatsfile und bin_Monatsfile

# 5. Create heat index-----------------------------------------------------------------------------

f <- list.files("./data/weather_era5/day/2m_temperature/")
f <- f[which(substr(f, 1, 3) == "max")]

fz <- function(z, zone){
 max <- readRDS(paste0("./data/weather_era5/day/2m_temperature/", z))  
 max <- max[[zone]]
 }
 
day_max_pac <- stack(map(f, fz, zone = 1))
day_max_mou <- stack(map(f, fz, zone = 2))
day_max_cen <- stack(map(f, fz, zone = 3))
day_max_eas <- stack(map(f, fz, zone = 4))

fh <- function(h){
  day_hwi <- h
  day_hwi[[1]]  <- day_hwi[[1]]-day_hwi[[1]]
  for(i in c(2:dim(h)[[3]])){
    tmp             <- h[[i]]
    tmp[tmp <  32]  <- NA                    
    tmp[tmp >= 32]  <- tmp[tmp >= 32]-32        #h+
    tmp[is.na(tmp) == TRUE]  <- -10^10          #h-
    tmp          <- day_hwi[[i-1]] + tmp
    tmp[tmp < 0] <- 0
    day_hwi[[i]]  <- tmp
    print(i)
  }
names(day_hwi) <- names(h)
day_hwi 
}

hwi <- map(list(day_max_pac, day_max_mou, day_max_cen, day_max_eas), fh) 
saveRDS(hwi,"./data/weather_era5/day/2m_temperature/hwi.RDS")


# 6. Aggregate daily raster variables to the county level -----------------------------------------

# Loop over variables: ~ 35 mins
load("./data/population/P_matrix.RData")

# find time zone of every county
# load county shape file 
load("./data/county_shape/us_counties.RData")
tz <- st_read("./data/county_shape/ne_10m_time_zones.shp")

centroid <- st_centroid(shape)
idx <- unlist(st_intersects(centroid, tz))
shape$tz <- tz$zone[idx]

tz_pacific <- shape %>% filter(tz == "-8") %>% select(GID_2); tz_pacific$geometry <- NULL
tz_mountain <-shape %>% filter(tz == "-7") %>% select(GID_2); tz_mountain$geometry <- NULL
tz_central <- shape %>% filter(tz == "-6") %>% select(GID_2); tz_central$geometry <- NULL
tz_eastern <- shape %>% filter(tz == "-5") %>% select(GID_2); tz_eastern$geometry <- NULL

fx <- function(x){
  # Read data to a matrix
  R <- readRDS(paste0("./data/weather_era5/day/", x))
  fz <- function(z){
    r <- R[[z]];    G <- r[];    P <- P_matrix
    # Aggregate to county level
    A  <- t(P) %*% G;    ct <- as.character(rownames(A))
    A  <- as.data.frame(as.matrix(A))
  }
  tmp1 <- fz(1); tmp1 <- tmp1 %>% filter(rownames(tmp1) %in% tz_pacific$GID_2)
  tmp2 <- fz(2); tmp2 <- tmp2 %>% filter(rownames(tmp2) %in% tz_mountain$GID_2)
  tmp3 <- fz(3); tmp3 <- tmp3 %>% filter(rownames(tmp3) %in% tz_central$GID_2)
  tmp4 <- fz(4); tmp4 <- tmp4 %>% filter(rownames(tmp4) %in% tz_eastern$GID_2)
  A <- bind_rows(tmp1, tmp2, tmp3, tmp4)
  saveRDS(A, paste0("./data/weather_era5/county/", x))
}

files <- list.files("./data/weather_era5/day/", recursive = TRUE)  
#map(files, fx)

# st¸ckchenweise:
x=files[109]
# noch zu computen: 
map(files[73:108],fx)
# map(files[114:145],fx) done (max_temperature)
x=files[288]

# Temperatur 73:
files <- files[grepl("max_", files, fixed =T)]
map(files,fx)

df <- readRDS("./data/weather_era5/county/2m_temperature/max_2018_2m_temperature_01.grd")

df19 <- stack("./data/weather_era5/month/era5/total_cloud_cover/2018_total_cloud_cover_05.grd")
plot(df19[[1]])

# 7. Combine all data -----------------------------------------------------------------------------

fx <- function(x, v){
  tmp       <- readRDS(paste0("./data/weather_era5/county/", x))
  tmp$GID_2 <- row.names(tmp)
 
   if(v == "2m_temperature"){
    vname = "2m_temperature_max"
  }else if(v %in% c("total_precipitation_sum", "total_precipitation_hours", "heat_index")){
    vname = paste0(v)
  }else{
    vname =paste0(v, "_mean") 
  }
  tmp1       <- reshape(tmp, varying=c(1:(ncol(tmp)-1)),
                       direction="long", idvar=c("GID_2"),
                       v.names=vname, timevar="day")
  dt <- colnames(tmp)[c(1:(length(colnames(tmp))-1))]
  tmp1 <- tmp1 %>% mutate(day = substr(dt, 2, nchar(dt))[day])
}

# monthly data
files <- list.files("./data/weather_era5/county/", recursive = TRUE)  
files <- files[!grepl("bin", files,  fixed = TRUE)]
files <- files[!grepl("hwi", files,  fixed = TRUE)]

for(i in c('2m_temperature', 'total_precipitation_sum', 'total_precipitation_hours',
'10m_u_component_of_wind', '10m_v_component_of_wind', 
'relative_humidity','surface_pressure', 
'total_cloud_cover')){
  
  if(i == 'total_precipitation_hours'){
    file <- files[grepl("hour", files,  fixed = TRUE)]
  }else if(i == 'total_precipitation_sum'){
    file <- files[grepl("sum", files,  fixed = TRUE)]
  }else{
  file <- files[grepl(i, files,  fixed = TRUE)]
  }
  
  df   <- bind_rows(map(file, fx, v = i))
  saveRDS(df, paste0("./data/weather_era5/county/final/", i, ".RDS"))
}


# hwi
files <- list.files("./data/weather_era5/county/", recursive = TRUE)  
files <- files[grepl("hwi", files,  fixed = TRUE)]

hwi <- bind_rows(map(files, fx, v = "heat_index"))
saveRDS(hwi, paste0("./data/weather_era5/county/final/heat_index.RDS"))

#binned temperature
files <- list.files("./data/weather_era5/county/", recursive = TRUE)  
files <- files[grepl("bin", files,  fixed = TRUE)]

fx <- function(x){
  tmp <- readRDS(paste0("./data/weather_era5/county/", x))
  bin <- unique(substr(colnames(tmp), 13, nchar(colnames(tmp))))
  fg <- function(g){
  tmp1 <- tmp[, colnames(tmp)[grepl(bin[[g]], colnames(tmp),  fixed = TRUE)]]
  tmp1$GID_2 <- row.names(tmp1)
  tmp2       <- reshape(tmp1, varying=c(1:(ncol(tmp1)-1)),
                        direction="long", idvar=c("GID_2"),
                        v.names= paste0("temperature_", bin[[g]]), timevar="day")
  dt <- colnames(tmp1)[c(1:(length(colnames(tmp1))-1))]
  tmp2 <- tmp2 %>% mutate(day = substr(dt, 2, 11)[day])
  }
  
  tmp <- map(c(1:length(bin)), fg)
  tmp <- join_all(tmp, by=c("GID_2", "day"), type='left')
}

df   <- bind_rows(map(files, fx))
saveRDS(df, "./data/weather_era5/county/final/binned_temperature.RDS")

# bis hier ohne Probleme gelaufen

# combine all
files <- list.files("./data/weather_era5/county/final/", recursive = TRUE)  
fl <- function(l){
  tmp <- readRDS(paste0("./data/weather_era5/county/final/", l))
}
data <- map(files, fl)
weather_data <- join_all(data, by=c("GID_2", "day"), type='left')

save(weather_data, file = "./data/weather_era5/weather_data.RData")


write.csv(weather_data, file="./data/weather_era5/weather_data.csv")


