###### read other data ##########
setwd("F:/DisHydroModell/DistributeHydrologyModell/OtherData")
###### read Grid Info include river, estuary, hydrostation infomation ######
GridID = read.table("ct_id5.txt", skip = 6)
infoGridRowN = dim(GridID)[1]  #the rows number of FLOWDRIC 
infoGridColN = dim(GridID)[2]   #the clows number of FLOWDRIC
GridAllIDVector = array(as.matrix(GridID),c(infoGridRowN * infoGridColN))
GridGridID = matrix(seq(1,infoGridN,1))

filePathLocation = "ct_location.txt"
GridLocation = read.table(filePathLocation, header = T)

FlowDirection = read.table("ct_direction5.txt", skip = 6)

RiverGridID = read.table("ct_river_id.txt") #one demesion
TranslateMatrixRiverInterflow = array(0.0, c(infoGridN))
TranslateMatrixRiverInterflow[as.matrix(RiverGridID)] = 1

EstuaryID = read.table("ct_estuary_id.txt") #one demesion

HydroStationID = read.table("ct_hydrostation_id.txt") #one wert

###### read DEM #########
filePathDEM = "ct_dem025.txt"
GridDEM = fctClassify(filePathDEM,20,0,5) ##m
GridEvalution = tapply(GridDEM[,2] * GridDEM[,3], GridDEM[,1], sum) * 1000  ##mm
id = as.data.frame(GridAllIDVector)
names(id) = "id"
GridDEM_m_Matrix = join(id,as.data.frame(cbind(id = seq(1,infoGridN,1),GridEvalution / 1000)))
GridDEM4UH = array(GridDEM_m_Matrix[,2], dim = c(infoGridRowN,infoGridColN))
GridDEM4UH[is.na(GridDEM4UH)] = constNull


###### make translate matrix #############
Grid2River = fctGrid2AimGrid (GridGridID, RiverGridID, GridDEM. =  GridDEM4UH)
River2Estuary = fctGrid2AimGrid(RiverGridID, EstuaryID, GridDEM. =  GridDEM4UH)
Estuary2Station = fctGrid2AimGrid(EstuaryID, HydroStationID, GridDEM. =  GridDEM4UH)

TranslateMatrixGrid2River = fctMakeGridTranslateMatrix(Grid2River, RiverGridID)
TranslateMatrixRiver2Estuary = fctMakeGridTranslateMatrix(River2Estuary, EstuaryID)
TranslateMatrixEstuary2Station = fctMakeGridTranslateMatrix(Estuary2Station, HydroStationID)
##############

###### read LanduseData ###########
filePathLanduse = "ct_landuse025.txt"
GridLanduse = fctClassify(filePathLanduse,20,0,5)
filePathLanduseLib = "landuse_lib_qinguha1ji.txt"
LanduseLib = read.table(filePathLanduseLib, header = T, row.names = 1)
LAI = as.data.frame(cbind(id = LanduseLib$code, LanduseLib[,6:17]))
RoughLength = as.data.frame(cbind(id = LanduseLib$code, LanduseLib[,18:29]))
DisplaceHeight = as.data.frame(cbind(id = LanduseLib$code, LanduseLib[,30:41]))

FrameGridLanduse = as.data.frame(GridLanduse)
names(FrameGridLanduse) = c("GridID", "code", "Rate")
GridLanduseParamTem = join(FrameGridLanduse, LanduseLib)[,-2]
GridLanduseParam = array(0.0, c(infoGridN, infoLandusParamN))
GridLanduseParam[,1] = 1:infoGridN
for (i in 2:infoLandusParamN) {
  GridLanduseParam[,i] = tapply(GridLanduseParamTem[,2] * GridLanduseParamTem[,1 + i], GridLanduseParamTem[,1], sum)
}
colnames(GridLanduseParam) = names(LanduseLib)
GridLanduseParam = as.data.frame(GridLanduseParam)

###### read SoilData ###########
filePathSoil = "ct_soil025.txt"
GridSoil = fctClassify(filePathSoil,20,0,5)
filePathSoilInterpolation = "soil_interpolation.txt"
SoilMatricPotential = read.table(filePathSoilInterpolation, row.names = 1)
filePathMU2Class = "MU2Class.txt"
MU2Class = read.table(filePathMU2Class, header = T)
filePathSoilLib = "soil_lib.txt"
SoilLib = read.table(filePathSoilLib, header = T, row.names = 1)
infoSoilLibN = dim(SoilLib)[1]
MU2Class[which(MU2Class[,3] == 0),3] = MU2Class[which(MU2Class[,3] == 0),2]

FrameGridSoil = as.data.frame(GridSoil)
names(FrameGridSoil) = c("GridID", "MU_GLOBAL", "Rate")
GridSoilClass = join(FrameGridSoil, MU2Class)

FrameGridTopSoil = GridSoilClass[,c(1,3,4)]
names(FrameGridTopSoil) = c("GridID",  "Rate", "Code")
TopSoilTem = join(FrameGridTopSoil, SoilLib)

FrameGridSubSoil = GridSoilClass[,c(1,3,5)]
names(FrameGridSubSoil) = c("GridID",  "Rate", "Code")
SubSoilTem = join(FrameGridSubSoil, SoilLib)

TopSoilParam = SubSoilParam = array(0.0, c(infoGridN, infoSoilParamN))
TopSoilParam[,1] = SubSoilParam[,1] = 1:infoGridN
for (i in 2:infoSoilParamN) {
  TopSoilParam[,i] = tapply(TopSoilTem[,2] * TopSoilTem[,2 + i], TopSoilTem[,1], sum)
  SubSoilParam[,i] = tapply(SubSoilTem[,2] * SubSoilTem[,2 + i], SubSoilTem[,1], sum)
}

colnames(TopSoilParam) = paste("T_",names(SoilLib),sep = "")
colnames(SubSoilParam) = paste("S_",names(SoilLib),sep = "")
GridSoilParam = as.data.frame(cbind(TopSoilParam, SubSoilParam[,2:infoSoilParamN]))
###### zone Paramter ###########
ZoneDepth = data.frame(Depth0 = GridLanduseParam$SL_mm, Depth1 = GridLanduseParam$root_depth_mm , Depth23 = paramDepth23, Depth4 = paramDepth4)
# ZoneDepth$Depth0 = ZoneDepth$Depth0 * GridLanduseParam[,5 + NMonth]
GridLanduseParam2Join = as.data.frame(cbind(seq(1,12,1), t(GridLanduseParam)[6:17,]))
names(GridLanduseParam2Join) = c("id", GridGridID)
NMonth2Join = as.data.frame(NMonth)
names(NMonth2Join) = c("id")
JoinTem1 = t(join(NMonth2Join, GridLanduseParam2Join)[,-1])


########################################
###### read MetrologyData #########
infoMetroFieldN = 8  ##*## #c(2,6,7,8,9,10,11,13) 1.lat	2.Temp	3.Tmax	4.Tmin	5.U	6.sun	7.humi	8.Rainfall
setwd("F:/DisHydroModell/DistributeHydrologyModell/MetroData_Cuntan")
StationLocationTable = read.table("StationLocation.txt", header = T)
StationLocation = StationLocationTable[,c(5,6)]
infoStationN = length(StationLocationTable[,1])
StationData = array(0.0, c(infoPeriodN, infoStationN, infoMetroFieldN))
for (i in 1:infoStationN) {
  Excel = read.xlsx(paste(StationLocationTable$id[i],".xlsx",sep = ""))[6942:(6941 + infoPeriodN),c(2,6,7,8,9,10,11,13)]
  StationData[,i,]= as.matrix(Excel)
}
###### fist deal and check Geo data #########
StationDealData = StationData
which(StationDealData == -32744) ## 32744 32766 -32766 -32744
StationDealData[which(StationDealData == 32700)] = 0.4
StationDealData[which(StationDealData >= 32000)] = StationDealData[which(StationDealData >= 32000)] - 32000
StationDealData[which(StationDealData >= 31000)] = StationDealData[which(StationDealData >= 31000)] - 31000
StationDealData[which(StationDealData >= 30000)] = StationDealData[which(StationDealData >= 30000)] - 30000

StationDealData[,,c(2,3,4,5,6,8)] = StationDealData[,,c(2,3,4,5,6,8)] * 0.1

GridMetroData = Station2GridIDW(StationLocation, GridLocation, StationDealData)

