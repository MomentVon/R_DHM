###### input data and some other Parame data ##########
###### period Parame ####
infoStarDay = as.Date("1989-1-1") ##*##
infoEndDay = as.Date("1993-12-31") ##*##
infoPeriod = 24 ##*## 24,12,8,6,4,3,2,1h<-day/h/min/s
UHPeriodN = 30 ##*## how many period in a all abfluss process
infoSubDay = as.integer(24 / infoPeriod) ##how many period hat in one day
DateDay = seq(infoStarDay,infoEndDay,1)
infoDayN = length(DateDay)
infoPeriodN = infoDayN * infoSubDay
NDay = toNDayofYear(DateDay)
NMonth = as.integer(substr(DateDay,6,7)) 
Date4ET = data.frame(NDay, NMonth)
# NDayTem2 = array(rep(NDayTem,subDay), c(dayN,subDay))
# NDay = array(t(NDayTem2), periodN)


###### Paramter Grid(commen) ##########
infoGridN = 311 #the number of grids
# infoGridRowN = 21 #the rows number of FLOWDRIC 
# infoGridColN = 35 #the clows number of FLOWDRIC
infoGridLength  = 55000 * 10^3 ##mm #the length every Grid
infoGridArea = infoGridLength^2 ##mm2 #the area every Grid
infoRiverGridN = 95 #one demesion
infoEstuaryN =4 #one demesion


###### param allprocess #########
infoZoneN = 6
infoSoilParamN = 7 ##*## ##
infoLandusParamN = 41 ##*##  ##id, rate, high, root depth, rarc, rmin, LAI X 12, rough X 12, DIS X 12
paramDepth23 = 400
paramDepth4 = 300


###### ET #########
infoWindH = 10 ##*## in which high measure the Wind speed, because the geoDAta don't have the infomation, so in there input it. If the Geo data have the infomation, it is not neccessary.


###### Interception ########
infoCoefficientFreeThroughfall = 0.9  ##*##
# CanopyStorageCapacity #s
# CoefficientFreeThroughfall  #, p
# TrunkPartitioningCoefficient # pt
# NecessaryToSaturateCanopy # Pg (mm)
# RainfallDuringSaturation # R (mmh- ‘)
# Evapotranspiration # E (mm h- ’ ),


###### Runoff ###########
paramSoilMoistureCapacityB = 0.21  ##*## ##b is the soil moisture capacity shape parameter, which is a measure of the spatial variability of the soil moisture capacity
paramInfiltrationRateB = 0.7  ##*## ##b is the soil infiltrtionrate shape parameter,


###### base flow #######
paramExponentARNOBase = 1.3
paramSoilMoistureVolumeARNOBaseThresholdRadio = 0.1
infoDrainageLossMax = 20 ##*## #mm
infoDrainageLossMin = 1 ##*## #mm


###### ground water ###########
paramClappHornbergerB = 1.5  ##*##
paramSoilPorosityToRootZone = 1.2 ##*## ##root zone is lecker than fest ston, so Porosity is biger than SoilPorosity
paramSoilPorosityToRootZone = 1.1  ##*## #in root zone soil is loker than other zone
paramClappHornbergerB = 1.1
infoLiquidDensity = 1.0 ##*## # Unit: 1000 kg/m3
infoTensionHeadMax = 1500
# BubblingPressure #ys is the air-entry tension head


###### rout ############
UHUnitTranslate = 35##*## ##need to caculate




