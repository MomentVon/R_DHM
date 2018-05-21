############### make a VIC M_VIC #############
M_VICRET = M_VICEvapC = M_VICTransp = M_VICEvapS = M_VICInterception = M_VICBaseflow = M_VICRunoff = M_VICInfiltration = M_VICRoutSurFlow = M_VICRoutBasFlow = array(0.0, c(infoPeriodN, infoGridN))
VICZoneDepth = data.frame(Depth0 = rep(10,infoGridN),Depth1 = rep(200,infoGridN),Depth2 = rep(300,infoGridN),Depth3 = rep(500,infoGridN))
M_VICMoistureVolume = array(0.0, c(infoPeriodN, infoGridN, 4))
M_VICMoistureVolume[1,,1:2] = 0.1 * as.matrix(VICZoneDepth[,1:2])

############ groundwater paramters #######################
SoilParam4GroundWaterTem = GridSoilParam[,c(2,3,4,6,7)]
names(SoilParam4GroundWaterTem) = c("Porosity", "FieldCapacity", "WiltingPoint", "SaturatedSoilSuctionHead", "SaturatedHydraulicConductivity")
VICSoilParam4GroundWater = cbind(VICZoneDepth, SoilParam4GroundWaterTem)

a = Sys.time()
for (M_VIC_i in 2:infoPeriodN) {
  ### evaptranspiration ############
  print(M_VIC_i)
  RET = ReferenceET_PenmanMonteith(NDay[M_VIC_i], GridEvalution / 1000.0, GridLocation[,2], GridMetroData[M_VIC_i,,3], GridMetroData[M_VIC_i,,4], GridMetroData[M_VIC_i,,2], GridMetroData[M_VIC_i,,5], infoWindH, GridMetroData[M_VIC_i,,6], GridMetroData[M_VIC_i,,7])
  M_VICRET[M_VIC_i,] = RET
  AerodynamicResistance = fctAerodynamicResistance(GridLanduseParam[,29 + NMonth[M_VIC_i]], GridLanduseParam[,17 + NMonth[M_VIC_i]], GridMetroData[M_VIC_i,,5])
  
  ET = ET_VIC(RET, GridMetroData[M_VIC_i,,8], M_VICMoistureVolume[(M_VIC_i-1),,1], VICZoneDepth$Depth0, M_VICMoistureVolume[(M_VIC_i-1),,2] + M_VICMoistureVolume[(M_VIC_i-1),,3], (VICZoneDepth$Depth1 + VICZoneDepth$Depth2) * GridSoilParam$T_Porosity_, AerodynamicResistance, GridLanduseParam$rarc, GridLanduseParam$rmin)
  M_VICEvapC[M_VIC_i,] = ET$EvapC
  M_VICTransp[M_VIC_i,] = ET$Transp
  M_VICEvapS[M_VIC_i,] = ET$EvapS
  
  ### interception ###############
  Interception = INTERCEPTION_Gash(ModelMoistureVolume[(M_VIC_i-1),,1], ZoneDepth$Depth0, GridMetroData[M_VIC_i,,8], ET$EvapC)
  M_VICInterception[M_VIC_i,] = Interception
  
  
  ### baseflow #########
  Baseflow = BASEFLOW_ARNO (M_VICMoistureVolume[(M_VIC_i-1),,4], VICZoneDepth$Depth3 * GridSoilParam$S_Porosity_)
  M_VICBaseflow[M_VIC_i,] = Baseflow
  ### runoff ############
  
  InfiltrationRateMax = fctInfiltrationGreenAmpt(GridSoilParam$T_SaturatedHydraulicConductivity_mm_day, GridSoilParam$T_WettingFrontSoilSuctionHead_mm, M_VICMoistureVolume[(M_VIC_i-1),,2] / VICZoneDepth$Depth1, GridSoilParam$T_Porosity_, M_VICMoistureVolume[(M_VIC_i-1),,2])
  
  Runoff = RUNOFF_VIC(GridMetroData[M_VIC_i,,8] - Interception, (paramSoilMoistureCapacityB + 1) * (VICZoneDepth$Depth1 + VICZoneDepth$Depth2) * GridSoilParam$T_Porosity_, maxSVector(0.0,M_VICMoistureVolume[(M_VIC_i-1),,2] + M_VICMoistureVolume[(M_VIC_i-1),,3]), InfiltrationRateMax)
  
  M_VICRunoff[M_VIC_i,] = Runoff$Runoff
  M_VICInfiltration[M_VIC_i,] = Runoff$Infiltration
  
  ### ground water replan / interflow ##########
  GroundWaterIn = as.data.frame(as.matrix(M_VICMoistureVolume[M_VIC_i - 1,,]))
  colnames(GroundWaterIn) = c("Volum0", "Volum1", "Volum2","Volum3")
  Groundwater = GROUNDWATER_VIC(ET, Interception, Runoff$Infiltration, Baseflow, GroundWaterIn, VICSoilParam4GroundWater)
  
  M_VICMoistureVolume[M_VIC_i,,] = as.matrix(Groundwater)
  M_VICMoistureVolume[M_VIC_i,,][is.na(M_VICMoistureVolume[M_VIC_i,,])] = 0.0
  
  ### route ###############
  
  M_VICRoutSurFlow[M_VIC_i,] = Runoff$Runoff
  M_VICRoutBasFlow[M_VIC_i,] = Baseflow
}
b = Sys.time()
b-a
M_VICUHAll = fctUHALLMake(fctMakeUH_Shipeng, fctMakeUH_Shipeng, fctMakeUH_Shipeng, fctMakeUH_Shipeng)
M_VICStationFlowQ = CONFLUENCE(M_VICRoutSurFlow, M_VICRoutBasFlow, M_VICUHAll)
plot(200:730,M_VICStationFlowQ[200:730])
write.table(M_VICStationFlowQ,"try1.txt")
RoutSurFlow[1:30,]
testn = 1826
mean(GridMetroData[2:testn,,8])
mean(M_VICRET[2:testn,])
mean(M_VICEvapC[2:testn,] + M_VICTransp[2:testn,] + M_VICEvapS[2:testn,])
mean(M_VICEvapC[2:testn,])
mean(M_VICTransp[2:testn,])
mean(M_VICEvapS[2:testn,])
mean(M_VICInterception[2:testn,])
mean(M_VICRunoff[2:testn,]) 
mean(M_VICInfiltration[2:testn,]) 
mean(M_VICMoistureVolume[2:testn, , 1])
mean(M_VICMoistureVolume[2:testn, , 2])
mean(M_VICMoistureVolume[2:testn, , 3])
mean(M_VICMoistureVolume[2:testn, , 4])
mean(M_VICRoutSurFlow[2:testn,])
mean(M_VICRoutBasFlow[2:testn,])

mean(M_VICStationFlowQ[2:testn])
dim(M_VICMoistureVolume)


write.table(M_VICStationFlowQ,"try1989-1932.txt",row.names = F)








