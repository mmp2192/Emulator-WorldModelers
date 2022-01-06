#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly = TRUE)

#### required libraries

library(R.matlab)
library(RANN)

today = as.numeric(format(Sys.Date(),"%Y"))
allyears = c(1983:2016)

lat = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/EastAfrica_AgGRIDlat.mat")
lat = lat$sublat
lon = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/EastAfrica_AgGRIDlon.mat")
lon = lon$sublon
lats=lat[,1]
lons=lon[1,]

dir = getwd()

thiscrop = args[1]
southlat = as.numeric(args[2])
northlat = as.numeric(args[3])
westlon = as.numeric(args[4])
eastlon = as.numeric(args[5])
tempchange = as.numeric(args[6])
precipchange = as.numeric(args[7])

southjj = nn2(lats,southlat)
southjj = southjj$nn.idx[1]
northjj = nn2(lats,northlat)
northjj = northjj$nn.idx[1]
westii = nn2(lons,westlon)
westii = westii$nn.idx[1]
eastii = nn2(lons,eastlon)
eastii = eastii$nn.idx[1]

gridlatidx = c(southjj:northjj)
gridlonidx = c(westii:eastii)
gridlats = lats[gridlatidx]
gridlons = lons[gridlonidx]

#### import global emulator variable indices and coefficients

regioncoeffs = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_EastAfrica_coefficients.mat", sep = ""))
regioncoeffs = regioncoeffs$regioncoeffs
regionvars = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_EastAfrica_variables.mat", sep = ""))
regionvars = regionvars$regionvars

### import growing season dates

regionpday = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_pday_EastAfrica.mat", sep = ""))
regionpday = regionpday$regionpday
regiongslength = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_gslength_EastAfrica.mat", sep = ""))
regiongslength = regiongslength$regiongslength

### import mean CHIRTS/CHIRPS emulated yield anomaly
### and mean CHIRTS/CHIRPS year

meanyieldanom = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_EastAfrica_mean_emulated_yield_anomaly.mat", sep = ""))
meanyieldanom = meanyieldanom$meanyieldanom

meanyieldanom = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_EastAfrica_mean_emulated_yield_anomaly.mat", sep = ""))
meanyieldanom = meanyieldanom$meanyieldanom

tasmax = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/CHIRTS_CHIRPS/chirts_05deg_daily_tasmax_mean_EastAfrica.mat", sep=""))
tasmax = tasmax$regiontmaxmean
tasmin = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/CHIRTS_CHIRPS/chirts_05deg_daily_tasmin_mean_EastAfrica.mat", sep=""))
tasmin = tasmin$regiontminmean
pr = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/CHIRTS_CHIRPS/chirps_05deg_daily_pr_mean_EastAfrica.mat", sep=""))
pr = pr$regionprmean

tasmax = tasmax+tempchange
tasmin = tasmin+tempchange
pr = pr*((precipchange/100)+1)
tas = (tasmax+tasmin)/2

#allcliminput = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/EastAfrica/Emulator_Parameters/rf", thiscrop, "_EastAfrica_climate_inputs.mat", sep = ""))
#allcliminput = allcliminput$allcliminput

output = array(data = NA, dim=c(length(gridlats),length(gridlons)))

for (ii in c(1:length(gridlons))){
	subii = gridlonidx[ii]
	for (jj in c(1:length(gridlats))){
		subjj = gridlatidx[jj]
		thislat = lat[subjj,subii]
		thislon = lon[subjj,subii]
		if (thislon>180){
	        thislon = thislon-180
		}
		if ((!is.na(regionpday[subjj,subii]))&&(!is.na(regiongslength[subjj,subii]))){
			
			### isolate growing season dates for this pixel
			pday = round(regionpday[subjj,subii])
			gslength = round(regiongslength[subjj,subii])
			hday = pday+gslength
			if (hday>365){
				hday=hday-365
			}
			if (hday<pday){
				seasonidx = c(pday:365, 1:hday)
			}else{
				seasonidx = c(pday:hday)
			}
			
			#### calculate planting window and anthesis time periods
			# planting window: 11-day window around planting date
			plantingidx = c((pday-5):(pday+5))
			plantingidx[plantingidx>365] = plantingidx[plantingidx>365]-365
			plantingidx[plantingidx<1] = plantingidx[plantingidx<1]+365
			# before anthesis: planting date to beginning of anthesis
			beforeidx = c(pday:(round(pday+(gslength/2))-5))
			beforeidx[beforeidx>365] = beforeidx[beforeidx>365]-365
			beforeidx[beforeidx<1] = beforeidx[beforeidx<1]+365
			# during anthesis: 11-day window around middle of growing season
			duringidx = c((round(pday+(gslength/2))-5):(round(pday+(gslength/2))+5))
			duringidx[duringidx>365] = duringidx[duringidx>365]-365
			duringidx[duringidx<1] = duringidx[duringidx<1]+365
			# after anthesis: end of anthesis to end of growing season
			afteridx = c((round(pday+(gslength/2))+5):(round(pday+gslength)))
			afteridx[afteridx>365] = afteridx[afteridx>365]-365
			afteridx[afteridx<1] = afteridx[afteridx<1]+365
			
	    	#### isolate emulator input variables and coefficients for this pixel
		    coeffs = regioncoeffs[subjj,subii,]
		    vars = regionvars[subjj,subii,]
		    
   			#### isolate pixel from CHIRTS/CHIRPS data
   			thistasmax = tasmax[subjj,subii,]
   			thistasmin = tasmin[subjj,subii,]
   			thispr = pr[subjj,subii,]
   			thistas = tas[subjj,subii,]
		    if (sum(!is.na(thispr))==365){
		    	
		    inputs = array(data=NA, dim=c(1,40))
		    
			# mean growing season temperature
			if (is.element(1,vars)){
				inputs[1,1] = mean(thistas[seasonidx])
			}
			# mean growing season precipitation
			if (is.element(2,vars)){
				inputs[1,2] = mean(thispr[seasonidx])
			}
			# number of growing season days where tasmax exceeds 30C
			if (is.element(3,vars)){
				inputs[1,3] = length(which((thistasmax[seasonidx])>30))
			}
			# number of growing season days where tasmax exceeds 35C
			if (is.element(4,vars)){
				inputs[1,4] = length(which((thistasmax[seasonidx])>35))
			}
			# number of growing season days where tasmin is less than 0C
			if (is.element(5,vars)){
				inputs[1,5] = length(which((thistasmin[seasonidx])<0))
			}
			# number of growing season days where tasmin is less than 5C
			if (is.element(6,vars)){
				inputs[1,6] = length(which((thistasmin[seasonidx])<5))
			}
			# number of growing season days where pr exceeds 1mm
			if (is.element(7,vars)){
				inputs[1,7] = length(which((thispr[seasonidx])>1))
			}
			# maximum consecutive growing season days where pr<0.01mm
			if (is.element(8,vars)){
				subpr = thispr[seasonidx]
				streak=numeric(0)
				if (subpr[1]<0.01){
					streak[1] = 1
				} else {
					streak[1] = 0
				}
				for (dd in c(2:length(subpr))){
					if (subpr[dd]<0.01){
						streak[dd] = streak[dd-1] + 1
						streak[dd-1] = 0
					} else {
						streak[dd] = 0
					}
				}
				inputs[1,8] = max(streak)
			}
			# mean temperature during planting window
			if (is.element(9,vars)){
				inputs[1,9] = mean(thistas[plantingidx])
			}
			# mean temperature before anthesis
			if (is.element(10,vars)){
				inputs[1,10] = mean(thistas[beforeidx])
			}
			# mean temperature during anthesis
			if (is.element(11,vars)){
				inputs[1,11] = mean(thistas[duringidx])
			}
			# mean temperature after anthesis
			if (is.element(12,vars)){
				inputs[1,12] = mean(thistas[afteridx])
			}
			# mean precipitation during planting window
			if (is.element(13,vars)){
				inputs[1,13] = mean(thispr[plantingidx])
			}
			# mean precipitation before anthesis
			if (is.element(14,vars)){
				inputs[1,14] = mean(thispr[beforeidx])
			}
			# mean precipitation during anthesis
			if (is.element(15,vars)){
				inputs[1,15] = mean(thispr[duringidx])
			}
			# mean precipitation after anthesis
			if (is.element(16,vars)){
				inputs[1,16] = mean(thispr[afteridx])
			}
			# number of planting window days where tasmax exceeds 30C
			if (is.element(17,vars)){
				inputs[1,17] = length(which((thistasmax[plantingidx])>30))
			}
			# number of days before anthesis where tasmax exceeds 30C
			if (is.element(18,vars)){
				inputs[1,18] = length(which((thistasmax[beforeidx])>30))
			}
			# number of days during anthesis where tasmax exceeds 30C
			if (is.element(19,vars)){
				inputs[1,19] = length(which((thistasmax[duringidx])>30))
			}
			# number of days after anthesis where tasmax exceeds 30C
			if (is.element(20,vars)){
				inputs[1,20] = length(which((thistasmax[afteridx])>30))
			}
			# number of planting window days where tasmax exceeds 35C
			if (is.element(21,vars)){
				inputs[1,21] = length(which((thistasmax[plantingidx])>35))
			}
			# number of days before anthesis where tasmax exceeds 35C
			if (is.element(22,vars)){
				inputs[1,22] = length(which((thistasmax[beforeidx])>35))
			}
			# number of days during anthesis where tasmax exceeds 35C
			if (is.element(23,vars)){
				inputs[1,23] = length(which((thistasmax[duringidx])>35))
			}
			# number of days after anthesis where tasmax exceeds 35C
			if (is.element(24,vars)){
				inputs[1,24] = length(which((thistasmax[afteridx])>35))
			}
			# number of planting window days where tasmin is less than 0C
			if (is.element(25,vars)){
				inputs[1,25] = length(which((thistasmin[plantingidx])<0))
			}
			# number of days before anthesis where tasmin is less than 0C
			if (is.element(26,vars)){
				inputs[1,26] = length(which((thistasmin[beforeidx])<0))
			}
			# number of days during anthesis where tasmin is less than 0C
			if (is.element(27,vars)){
				inputs[1,27] = length(which((thistasmin[duringidx])<0))
			}
			# number of days after anthesis where tasmin is less than 0C
			if (is.element(28,vars)){
				inputs[1,28] = length(which((thistasmin[afteridx])<0))
			}
			# number of planting window days where tasmin is less than 5C
			if (is.element(29,vars)){
				inputs[1,29] = length(which((thistasmin[plantingidx])<5))
			}
			# number of days before anthesis where tasmin is less than 5C
			if (is.element(30,vars)){
				inputs[1,30] = length(which((thistasmin[beforeidx])<5))
			}
			# number of days during anthesis where tasmin is less than 5C
			if (is.element(31,vars)){
				inputs[1,31] = length(which((thistasmin[duringidx])<5))
			}
			# number of days after anthesis where tasmin is less than 5C
			if (is.element(32,vars)){
				inputs[1,32] = length(which((thistasmin[afteridx])<5))
			}
			# number of planting window days where pr exceeds 1mm
			if (is.element(33,vars)){
				inputs[1,33] = length(which((thispr[plantingidx])>1))
			}
			# number of days before anthesis where pr exceeds 1mm
			if (is.element(34,vars)){
				inputs[1,34] = length(which((thispr[beforeidx])>1))
			}
			# number of days during anthesis where pr exceeds 1mm
			if (is.element(35,vars)){
				inputs[1,35] = length(which((thispr[duringidx])>1))
			}
			# number of days after anthesis where pr exceeds 1mm
			if (is.element(36,vars)){
				inputs[1,36] = length(which((thispr[afteridx])>1))
			}
			# maximum consecutive planting window days where pr<0.01mm
			if (is.element(37,vars)){
				subpr = thispr[plantingidx]
				streak=numeric(0)
				if (subpr[1]<0.01){
					streak[1] = 1
				} else {
					streak[1] = 0
				}
				for (dd in c(2:length(subpr))){
					if (subpr[dd]<0.01){
						streak[dd] = streak[dd-1] + 1
						streak[dd-1] = 0
					} else {
						streak[dd] = 0
					}
				}
				inputs[1,37] = max(streak)
			}
			# maximum consecutive days before anthesis where pr<0.01mm
			if (is.element(38,vars)){
				subpr = thispr[beforeidx]
				streak=numeric(0)
				if (subpr[1]<0.01){
					streak[1] = 1
				} else {
					streak[1] = 0
				}
				for (dd in c(2:length(subpr))){
					if (subpr[dd]<0.01){
						streak[dd] = streak[dd-1] + 1
						streak[dd-1] = 0
					} else {
						streak[dd] = 0
					}
				}
				inputs[1,38] = max(streak)
			}
			# maximum consecutive days during anthesis where pr<0.01mm
			if (is.element(39,vars)){
				subpr = thispr[duringidx]
				streak=numeric(0)
				if (subpr[1]<0.01){
					streak[1] = 1
				} else {
					streak[1] = 0
				}
				for (dd in c(2:length(subpr))){
					if (subpr[dd]<0.01){
						streak[dd] = streak[dd-1] + 1
						streak[dd-1] = 0
					} else {
						streak[dd] = 0
					}
				}
				inputs[1,39] = max(streak)
			}
			# maximum consecutive days after anthesis where pr<0.01mm
			if (is.element(40,vars)){
				subpr = thispr[afteridx]
				streak=numeric(0)
				if (subpr[1]<0.01){
					streak[1] = 1
				} else {
					streak[1] = 0
				}
				for (dd in c(2:length(subpr))){
					if (subpr[dd]<0.01){
						streak[dd] = streak[dd-1] + 1
						streak[dd-1] = 0
					} else {
						streak[dd] = 0
					}
				}
				inputs[1,40] = max(streak)
			}
			
			v = inputs[vars]
			
			emuyield = coeffs[1] + (coeffs[2]*v[1]) + (coeffs[3]*v[2]) + (coeffs[4]*v[3]) + (coeffs[5]*v[4]) + (coeffs[6]*v[5]) + (coeffs[7]*v[1]*v[2]) + (coeffs[8]*v[1]*v[3]) + (coeffs[9]*v[1]*v[4]) + (coeffs[10]*v[1]*v[5]) + (coeffs[11]*v[2]*v[3]) + (coeffs[12]*v[2]*v[4]) + (coeffs[13]*v[2]*v[5]) + (coeffs[14]*v[3]*v[4]) + (coeffs[15]*v[3]*v[5]) + (coeffs[16]*v[4]*v[5]) + (coeffs[17]*v[1]*v[1]) + (coeffs[18]*v[2]*v[2]) + (coeffs[19]*v[3]*v[3]) + (coeffs[20]*v[4]*v[4]) + (coeffs[21]*v[5]*v[5])
			
			meanemuyield = meanyieldanom[subjj,subii]
			output[jj,ii] = emuyield/meanemuyield
		    }
		}
	}
}

output[which(output>100)] = 100
output[which(output<(-100))] = -100

ncells = length(gridlats)*length(gridlons)

df = data.frame(latitude = NA, longitude = NA, yield_anomaly = NA, time = NA)

counter = 1

for (ii in c(1:length(gridlons))){
	subii = gridlonidx[ii]
	for (jj in c(1:length(gridlats))){
		subjj = gridlatidx[jj]
		thislat = lat[subjj,subii]
		thislon = lon[subjj,subii]
		df[counter,1] = thislat
		df[counter,2] = thislon
		df[counter,3] = output[jj,ii]
		df[counter,4] = today
		counter = counter+1
	}
}

write.csv(df, paste(dir, "/output/output.csv", sep=""), row.names=F)
