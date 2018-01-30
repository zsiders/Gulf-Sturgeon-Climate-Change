setwd("/Volumes/HDD/Users/Zach/Documents/Projects/Gulf_Sturgeon/Gulf_Sturgeon_GitHub/Gulf-Sturgeon-Climate-Change")

###-----------------------------------------------------
#		READING IN CLIMATE DATA
###-----------------------------------------------------

library(ncdf4)

null_f <- nc_open(paste0(getwd(),"/climate_CM3/Rp 2.6 - Low/", list.files("./climate_CM3/Rp 2.6 - Low/")[1]))
latvec <- ncvar_get(null_f, varid='lat_bnds')
lonvec <- ncvar_get(null_f, varid='lon_bnds')


pr_26 <- vector(mode='list', length=length(list.files("./climate_CM3/Rp 2.6 - Low/")))

for(i in 1:length(pr_26)) pr_26[[i]] <- ncvar_get(nc_open(paste0(getwd(),"/climate_CM3/Rp 2.6 - Low/", list.files("./climate_CM3/Rp 2.6 - Low/")[i])),varid="pr")

pr_45 <- vector(mode='list', length=length(list.files("./climate_CM3/Rp 4.5 - Mid/")))

for(i in 1:length(pr_45)) pr_45[[i]] <- ncvar_get(nc_open(paste0(getwd(),"/climate_CM3/Rp 4.5 - Mid/", list.files("./climate_CM3/Rp 4.5 - Mid/")[i])),varid="pr")

pr_60 <- vector(mode='list', length=length(list.files("./climate_CM3/Rp 6.0 - MidHigh/")))

for(i in 1:length(pr_60)) pr_60[[i]] <- ncvar_get(nc_open(paste0(getwd(),"/climate_CM3/Rp 6.0 - MidHigh/", list.files("./climate_CM3/Rp 6.0 - MidHigh/")[i])),varid="pr")

pr_85 <- vector(mode='list', length=length(list.files("./climate_CM3/Rp 8.5 - High/")))

for(i in 1:length(pr_85)) pr_85[[i]] <- ncvar_get(nc_open(paste0(getwd(),"/climate_CM3/Rp 8.5 - High/", list.files("./climate_CM3/Rp 8.5 - High/")[i])),varid="pr")

###-----------------------------------------------------
#		TRANSFORMING PRECIPITATION DATA
###-----------------------------------------------------

stat_coord <- read.table('station_GPS_2.txt',sep='\t',header=TRUE)
sw_statcoor <- stat_coord[stat_coord$drainage=='sw',]
sw_statcoor <- sw_statcoor[-1,]
ap_statcoor <- stat_coord[stat_coord$drainage=='ap',]

wtmat<-function(lat,lon,latvec,lonvec){
  #FINDS THE LAT LON BOUNDS
  latt<-latvec[which.min(abs(latvec-lat))]
  if(latt > lat){
	latdn <- latt-2
	latup <- latt
  } else {
	latup <- latt+2
	latdn <- latt
  }
  lont<-lonvec[which.min(abs(lonvec-lon))]
  if(lont > lon){
	londn <- lont-2.5
	lonup <- lont
  } else {
	lonup <- lont+2.5
	londn <- lont
  }  
  latl <- c(latdn,latup)
  lonl <- c(londn,lonup)
  #CALCULATES A WEIGHT MATRIX
  wt <- matrix(NA,2,2)
  latwt <- c(2-(lat-latdn),2-(latup-lat))
  lonwt <- c(2.5-(lon-londn),2.5-(lonup-lon))
  latwt <- latwt/2
  lonwt <- lonwt/2.5
  wt[1,1] <- latwt[1]+lonwt[1]
  wt[1,2] <- latwt[2]+lonwt[1]
  wt[2,1] <- latwt[1]+lonwt[2]
  wt[2,2] <- latwt[2]+lonwt[2]
  wt <- wt/2
  list(latl,lonl,wt,latwt,lonwt)
}

#weight storage lists
apwt <- list(NA,nrow(ap_statcoor))
swwt <- list(NA,nrow(ap_statcoor))
apwtm <- array(NA,dim=c(2,2,nrow(ap_statcoor)))
swwtm <- array(NA,dim=c(2,2,nrow(sw_statcoor)))
#for loop running wtmat function
for(i in 1:nrow(ap_statcoor))
{
	apwt[[i]] <- wtmat(ap_statcoor[i,3],abs(ap_statcoor[i,4])+180,latvec[1,],lonvec[1,])
	apwtm[,,i] <- apwt[[i]][[3]]
}
for (i in 1:nrow(sw_statcoor)){
	swwt[[i]] <- wtmat(sw_statcoor[i,3],abs(sw_statcoor[i,4])+180,latvec[1,],lonvec[1,])
	swwtm[,,i] <- swwt[[i]][[3]]
}

#gets the coodinates from the NC file of climate data by extracting them from the wtmat values
lonap2 <- matrix(NA,nrow(ap_statcoor),2)
latap2 <- matrix(NA,nrow(ap_statcoor),2)
for (i in 1:nrow(ap_statcoor)){
	lonap2[i,] <- apwt[[i]][[2]]
	latap2[i,] <- apwt[[i]][[1]]
}
lonsw2 <- matrix(NA,nrow(sw_statcoor),2)
latsw2 <- matrix(NA,nrow(sw_statcoor),2)
for (i in 1:nrow(sw_statcoor)){
	lonsw2[i,] <- swwt[[i]][[2]]
	latsw2[i,] <- swwt[[i]][[1]]
}

summat<- function(var,lonvec,latvec,lon,lat,station,year,wt,month)
{
	if(length(month) > 1){
		tmpb <- rep(NA,length(month))
		for(i in 1:length(month))
		{
			tmp <- matrix(NA,2,2) #makes matrix
			m=month[i]
			W=station
			Q=year
			tmp[1,1] <- var[which(lonvec==lon[W,1]), which(latvec==lat[W,1]), eval(m+12*(Q-1))]
			tmp[1,2] <- var[which(lonvec==lon[W,1]), which(latvec==lat[W,2]), eval(m+12*(Q-1))]
			tmp[2,1] <- var[which(lonvec==lon[W,2]), which(latvec==lat[W,1]), eval(m+12*(Q-1))]
			tmp[2,2] <- var[which(lonvec==lon[W,2]), which(latvec==lat[W,2]), eval(m+12*(Q-1))]
			tmp2 <- tmp*wt[[W]][[3]]
			tmpb[i] <- sum(tmp2)/2
		}
		if(mean(tmpb)<1){
			#0.039701*3600*24*30 is 1kg=1liter=0.0393701inches of water per m^2 per s times 3600 seconds per hour times 24 hours per day times 30,31,30 days per month
			tmpa <- mean(tmpb*c(30,31,30)*0.039701*3600*24)
		}else{
			tmpa <- mean(tmpb)*1.8+32
		}
	}else{
		tmpa <- NA
		tmp <- matrix(NA,2,2) #makes matrix
		m=month
		W=station
		Q=year
		tmp[1,1] <- var[which(lonvec==lon[W,1]), which(latvec==lat[W,1]), eval(m+12*(Q-1))]
		tmp[1,2] <- var[which(lonvec==lon[W,1]), which(latvec==lat[W,2]), eval(m+12*(Q-1))]
		tmp[2,1] <- var[which(lonvec==lon[W,2]), which(latvec==lat[W,1]), eval(m+12*(Q-1))]
		tmp[2,2] <- var[which(lonvec==lon[W,2]), which(latvec==lat[W,2]), eval(m+12*(Q-1))]
		tmp2 <-tmp*wt[[W]][[3]]
		tmpb <- sum(tmp2)/2
		if(tmpb < 1){
			tmpa <- tmpb*0.039701*3600*24*30
		}else{
			tmpa <- tmpb*1.8+32
		}
	}
	return(tmpa)
}
slider <- seq(0,100-5,by=5)
var <- list(pr_26, pr_45, pr_60, pr_85)
apcarr <- array(NA,dim=c(95,nrow(ap_statcoor),length(var)),dimnames=list(seq(2006,2100), 1:nrow(ap_statcoor), paste0('rp',c('2.6','4.5','6.0','8.5'))))
swcarr <- array(NA,dim=c(95,nrow(sw_statcoor),length(var)),dimnames=list(seq(2006,2100), 1:nrow(sw_statcoor), paste0('rp',c('2.6','4.5','6.0','8.5'))))

for (g in 1:length(var)){
	for (i in 1:nrow(ap_statcoor)){
		for (j in 1:5){
			for(h in 1:length(pr_26))
			{
				apcarr[j+slider[h],i,g] <- summat(var[[g]][[h]], lonvec[1,], ceiling(latvec[1,]), lonap2, latap2, i, j, apwt, 9:11)
			}
		}
	}
}

for (g in 1:length(var)){
	for (i in 1:nrow(sw_statcoor)){
		for (j in 1:5){
			for(h in 1:length(pr_26))
			{
				swcarr[j+slider[h],i,g] <- summat(var[[g]][[h]], lonvec[1,], ceiling(latvec[1,]), lonsw2, latsw2, i, j, swwt, 9:11)
			}
		}
	}
}











