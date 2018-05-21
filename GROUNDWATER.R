###### Groung water ###########
fctHydraulicConductivity <- function(VolumetricSoilMoistureContent, SoilPorosity, SaturatedHydraulicConductivity, paramClappHornbergerB. = paramClappHornbergerB){
  HydraulicConductivity = SaturatedHydraulicConductivity * (VolumetricSoilMoistureContent / SoilPorosity)^(2 * paramClappHornbergerB. + 3)
  HydraulicConductivity[is.nan(HydraulicConductivity)] = 0.0 
  return(HydraulicConductivity)
}
fctHydraulicDiffusivity <- function(VolumetricSoilMoistureContent, SoilPorosity, SaturatedSoilSuctionHead, SaturatedHydraulicConductivity, paramClappHornbergerB. = paramClappHornbergerB){
  return(paramClappHornbergerB. * SaturatedSoilSuctionHead * SaturatedHydraulicConductivity / SoilPorosity * (VolumetricSoilMoistureContent / SoilPorosity)^(paramClappHornbergerB. + 2))
}

###### caculate the paramters #############
fctTensionHeadInterpolation <- function(GridID, VolumetricSoilMoistureContent, SoilMatricPotential){
  MatricPotential = array(0.0, c(infoSoilLibN))
  TensionHead = array(0.0, c(infoSoilLibN))
  infoInterpolationN = (dim(SoilMatricPotential)[2] - 1) / 2
  for (i in 1:infoSoilLibN) {
    MatricPotential[i] = as.numeric(interpolation(VolumetricSoilMoistureContent, SoilMatricPotential[i,(2 : (1 + infoInterpolationN))], SoilMatricPotential[i,((2 + infoInterpolationN) :( 1 + 2 * infoInterpolationN))]))
  }
  TensionHead = MatricPotential / (9.8 * infoLiquidDensity)
  return(TensionHead)
}

fctTensionHead <- function(Volum, Depth, SoilPorosity, SoilMatricPotential, GridSoilClass){
  Volum = as.matrix(Volum)
  Depth = as.matrix(Depth)
  HeadTensionTEM = array(0.0, c(infoGridN * infoSoilLibN,3))
  for (i in 1:infoGridN) {
    HeadTensionTEM[(1 + infoSoilLibN * (i - 1)) : (i * infoSoilLibN),1] = i
    HeadTensionTEM[(1 + infoSoilLibN * (i - 1)) : (i * infoSoilLibN),2] = SoilMatricPotential[,1]
    HeadTensionTEM[(1 + infoSoilLibN * (i - 1)) : (i * infoSoilLibN),3] = ifelse(Depth[i] == 0.0, 0.0, fctTensionHeadInterpolation(i, Volum[i] / Depth[i], SoilMatricPotential))
  }
  FrameGridSoilClass = as.data.frame(GridSoilClass)
  FrameHeadTensionTEM = as.data.frame(HeadTensionTEM)
  names(FrameGridSoilClass) = c("GridID",  "Rate", "Code")
  names(FrameHeadTensionTEM) = c("GridID",  "Code", "TensionHead")
  HeadTensionTem2 = join(FrameGridSoilClass, FrameHeadTensionTEM, by = c("GridID", "Code"))
  HeadTension = tapply(HeadTensionTem2[,2] * HeadTensionTem2[,4], HeadTensionTem2[,1],sum)
  HeadTension[is.na(HeadTension)] = infoTensionHeadMax
  return(HeadTension)
}


fctChangeVolum <- function(BaseDepth, HeadHight, HeadLow, HydraulicConductivity){
  ChangeVolum = HydraulicConductivity * (HeadHight) * (HeadHight - HeadLow) / infoGridArea
  return(ChangeVolum)
}

fctChangeHead <- function(ChangeVolum, SoilPorosity){
  return(ChangeVolum / SoilPorosity)
}

###### ground water bewegen ############
fctWaterHeadBalance <- function(Volum, Head,  ZoneDepth, BaseDepth, Porosity, SaturatedHydraulicConductivity){
  GridID2Join = as.data.frame(GridAllIDVector)
  names(GridID2Join) = "id"
  GridPara2Join = as.data.frame(cbind(seq(1,infoGridN,1),Volum, Head,  ZoneDepth, BaseDepth, Porosity, SaturatedHydraulicConductivity))
  names(GridPara2Join) = c("id", "Volum", "Head",  "ZoneDepth", "BaseDepth", "Porosity", "SaturatedHydraulicConductivity")
  WertVector = join(GridID2Join, GridPara2Join)
  WertMatrix = array(as.matrix(WertVector),c(infoGridRowN, infoGridColN,7))
  WertMatrix[which(WertMatrix == constNull)] = NA
  for (i in 1:(infoGridRowN - 1)){
    for (j in 1:(infoGridColN - 1)){
      ChangeMatrix = JudgeMatrix = WertMatrix[i : (i+1), j : (j+1),]
      if(length(is.na(JudgeMatrix)) >= 3) next
      HighHead = max(JudgeMatrix[,,3] + JudgeMatrix[,,5])
      LowHead = min(JudgeMatrix[,,3] + JudgeMatrix[,,5])
      HighLocal = which((JudgeMatrix[,,3] + JudgeMatrix[,,5]) == HighHead)[2]
      LowLocal = which((JudgeMatrix[,,3] + JudgeMatrix[,,5]) == LowHead)[2]
      MaxVolum = JudgeMatrix[,,2][HighLocal]
      MinVolum = JudgeMatrix[,,2][LowLocal]
      MaxDepth = JudgeMatrix[,,4][HighLocal]
      MinDepth = JudgeMatrix[,,4][LowLocal]
      MaxBaseDepth = JudgeMatrix[,,5][HighLocal]
      MinBaseDepth = JudgeMatrix[,,5][LowLocal]
      MaxPorosity = JudgeMatrix[,,6][HighLocal]
      MinPorosity = JudgeMatrix[,,6][LowLocal]
      MaxSaturatedHydraulicConductivity = JudgeMatrix[,,7][HighLocal]
      
      HydraulicConductivity = fctHydraulicConductivity (MaxVolum / MaxDepth, MaxPorosity, MaxSaturatedHydraulicConductivity)
      ChangeVolum = fctChangeVolum(MaxBaseDepth - MinBaseDepth, HighHead, LowHead, HydraulicConductivity)
      MaxChangeHead = fctChangeHead(ChangeVolum, MaxPorosity)
      MinChangeHead = fctChangeHead(ChangeVolum, MinPorosity)
      
      ChangeMatrix[,,2][HighLocal] = ChangeMatrix[,,2][HighLocal] - ChangeVolum
      ChangeMatrix[,,2][LowLocal] = ChangeMatrix[,,2][LowLocal] + ChangeVolum
      ChangeMatrix[,,3][HighLocal] = ChangeMatrix[,,3][HighLocal] -  MaxChangeHead
      ChangeMatrix[,,3][LowLocal] = ChangeMatrix[,,3][LowLocal] + MaxChangeHead
      WertMatrix[i : (i+2), j : (j+2),] = ChangeMatrix
    }
  }
  WertVectorTem = as.data.frame(array(WertMatrix,c(infoGridRowN * infoGridColN,7)))
  names(WertVectorTem) = c("id", "Volum", "Head",  "ZoneDepth", "BaseDepth", "Porosity", "SaturatedHydraulicConductivity")
  
  GridID2Join2 = as.data.frame(seq(1,infoGridN,1))
  names(GridID2Join2) = "id"
  WertOut = join(GridID2Join2, WertVectorTem)
  return(WertOut)
}

###### ground water function ############
###### ground water VIC function ############
GROUNDWATER_VIC <- function(Evapotranspiration, Interception, Infiltration, BaseFlow, GroundWaterIn, GridSoilParame){ ## Infiltration = Infiltration - Evapotranspiration
  Infiltration0 = Infiltration
  HydraulicConductivity1_2 = fctHydraulicConductivity(GroundWaterIn$Volum1 / GridSoilParame$Depth1, GridSoilParame$Porosity, GridSoilParame$SaturatedHydraulicConductivity)
  HydraulicConductivity2_3 = fctHydraulicConductivity(GroundWaterIn$Volum2 / GridSoilParame$Depth2, GridSoilParame$Porosity, GridSoilParame$SaturatedHydraulicConductivity)
  HydraulicDiffusivity1_2 = fctHydraulicDiffusivity(GroundWaterIn$Volum1 / GridSoilParame$Depth1, GridSoilParame$SaturatedSoilSuctionHead, GridSoilParame$Porosity, GridSoilParame$SaturatedHydraulicConductivity)  # Unit ist depth: mm ##must translate to a time phase
  HydraulicDiffusivity2_3 = fctHydraulicDiffusivity(GroundWaterIn$Volum2 / GridSoilParame$Depth2, GridSoilParame$SaturatedSoilSuctionHead, GridSoilParame$Porosity, GridSoilParame$SaturatedHydraulicConductivity)  # Unit ist depth: mm ##must translate to a time phase
  InterFlowFlux1 = HydraulicConductivity1_2 + HydraulicDiffusivity1_2 
  InterFlowFlux2 = HydraulicConductivity2_3 + HydraulicDiffusivity2_3 
  InterFlowFlux1[is.na(InterFlowFlux1)] = 0.0
  InterFlowFlux2[is.na(InterFlowFlux2)] = 0.0
  InterFlowFlux1 = minVector(InterFlowFlux1, GroundWaterIn$Volum1)
  InterFlowFlux2 = minVector(InterFlowFlux2, GroundWaterIn$Volum2)
  GroundWaterOut = GroundWaterIn
  GroundWaterOut$Volum0 = maxSVector(0.0, GroundWaterIn$Volum0 + as.matrix(Interception) - as.matrix(Evapotranspiration$EvapC))
  GroundWaterOut$Volum1 = maxSVector(0.0, GroundWaterIn$Volum1 - as.matrix(Evapotranspiration$EvapS + Evapotranspiration$Transp) + Infiltration0 - InterFlowFlux1)
  GroundWaterOut$Volum2 = maxSVector(0.0, GroundWaterIn$Volum2 + InterFlowFlux1 - InterFlowFlux2)
  GroundWaterOut$Volum3 = maxSVector(0.0, GroundWaterIn$Volum3 + InterFlowFlux2 - BaseFlow)
  return(GroundWaterOut)
}
###### ground water LK3L function ############
GROUNDWATER_LK3L <- function(Evapotranspiration, Infiltration, BaseFlow, GroundWaterIn, GridSoilParame){
  InterFlowFlux0 = minVector(Infiltration, Evapotranspiration)
  InterFlowFlux1 = maxVector(0.0, Infiltration - InterFlowFlux0)
  InterFlowFlux2 = minVector(Infiltration1, BaseFlow)
  GroundWaterOut = GroundWaterIn
  GroundWaterOut$Volum0 = GroundWaterIn$Volum0 + Infiltration0 - Evapotranspiration
  GroundWaterOut$Volum1 = minVector(GroundWaterIn$Volum1 + InterFlowFlux1 - InterFlowFlux2, GridSoilParame$Depth1 * GridSoilParame$Porosity)
  GroundWaterOut$Volum2 = GroundWaterIn$Volum2 + InterFlowFlux2 - BaseFlow
  GroundWaterOut$Volum3 = maxSVector(0.0, (GroundWaterIn$Volum1 + InterFlowFlux1 - InterFlowFlux2) - GridSoilParame$Depth1 * GridSoilParame$Porosity)
  Head = fctVolumToHead(GroundWaterOut$Volum1, GridSoilParame$Depth1, GridSoilParame$Porosity, GridSoilParame$Type)
  Wert = fctWaterHeadBalance(GroundWaterOut$Volum1, Head)
  GroundWaterOut$Volum1 = Wert$Volum
  return(GroundWaterOut)
}
###### ground water LK5Z function ############
GROUNDWATER_LK5Z <- function(Evapotranspiration, Interception, Infiltration, BaseFlow, GroundWaterIn, GridSoilParame, SoilMatricPotential, GridSoilClass.){
  GroundWaterOut = GroundWaterIn
  InterFlowFlux0 = as.matrix(Interception)
  GroundWaterOut$Volum0 = maxSVector(0.0, GroundWaterIn$Volum0 + InterFlowFlux0 - as.matrix(Evapotranspiration$EvapC))
  InterFlowFlux1 = Infiltration + GroundWaterIn$Volum0 + InterFlowFlux0 - (as.matrix(Evapotranspiration$Transp) + as.matrix(Evapotranspiration$EvapS)) 
  GroundWaterOut$Volum1 = minVector(GroundWaterIn$Volum1 + InterFlowFlux1, GridSoilParame$Depth1 * GridSoilParame$Porosity * paramSoilPorosityToRootZone)
  InterFlowFlux2 = maxSVector(0.0, GroundWaterIn$Volum1 + InterFlowFlux1 - GridSoilParame$Porosity * GridSoilParame$Depth1 * paramSoilPorosityToRootZone) 
  GroundWaterOut$Volum1 = maxVector(GridSoilParame$WiltingPoint * GridSoilParame$Depth1, GroundWaterOut$Volum1 - InterFlowFlux2)
  
  InterFlowFlux4 = minVector(GroundWaterIn$Volum3, GridSoilParame$Depth4 * GridSoilParame$Porosity_Sub - GroundWaterIn$Volum4)
  GroundWaterOut$Volum4 = GroundWaterIn$Volum4 + InterFlowFlux4 - BaseFlow
  # GroundWaterOut$Volum4 = minVector(GridSoilParame$Depth4 * GridSoilParame$Porosity_Sub, GroundWaterIn$Volum4 + InterFlowFlux4 - BaseFlow)
  # InterFlowFlux5_4 = maxSVector(0.0, GroundWaterIn$Volum4 + InterFlowFlux4 - BaseFlow - (GridSoilParame$Depth4 * GridSoilParame$Porosity_Sub))
  
  Volum23 = minVector(GridSoilParame$Porosity * GridSoilParame$Depth23, GroundWaterIn$Volum3 + GroundWaterIn$Volum2 - InterFlowFlux4 + InterFlowFlux2)
  InterFlowFlux5_23 = maxSVector(0.0, (GroundWaterIn$Volum3 + GroundWaterIn$Volum2 - InterFlowFlux4 + InterFlowFlux2) - (GridSoilParame$Porosity * GridSoilParame$Depth23))
  Depth2TermVolum_1_2 = GridSoilParame$Porosity * GridSoilParame$Depth23 - Volum23
  
  Depth2Judge_2_3 = GridSoilParame$Depth23 * 0.5 * (GridSoilParame$FieldCapacity + GridSoilParame$Porosity)
  Depth2_2 = 2 * Depth2TermVolum_1_2 / (GridSoilParame$Porosity - GridSoilParame$FieldCapacity)
  GroundWaterOut$Depth2 = Depth2_2
  GroundWaterOut$Depth2[which(InterFlowFlux5_23 > 0)] = 0.0
  GroundWaterOut$Depth2[which(Volum23 < Depth2Judge_2_3)] = GridSoilParame$Depth23[which(Volum23 < Depth2Judge_2_3)]
  GroundWaterOut$Depth3 = GridSoilParame$Depth23 - GroundWaterOut$Depth2
  
  GroundWaterOut$Volum2 = 0.5 * (GridSoilParame$Porosity + GridSoilParame$FieldCapacity) * GroundWaterOut$Depth2
  GroundWaterOut$Volum2 = minVector(GroundWaterOut$Volum2, Volum23)
  GroundWaterOut$Volum3 = Volum23 - GroundWaterOut$Volum2
  
  TensionHead = fctTensionHead(GroundWaterOut$Volum2, GroundWaterOut$Depth2, GridSoilParame$Porosity, SoilMatricPotential, GridSoilClass.)
  BalanceChanged = fctWaterHeadBalance(GroundWaterOut$Volum3, maxSVector(0.0, GroundWaterOut$Depth3 - 0.8 * as.matrix(TensionHead)),  GroundWaterOut$Depth3, GridEvalution - GridSoilParame$Depth1 - GridSoilParame$Depth23 - GridSoilParame$Depth4, GridSoilParame$Porosity, GridSoilParame$SaturatedHydraulicConductivity)
  GroundWaterOut$Volum3 = BalanceChanged$Volum
  RiverGet = GroundWaterOut$Volum3 * as.matrix(TranslateMatrixRiverInterflow)
  GroundWaterOut$Volum3 = GroundWaterOut$Volum3 - RiverGet
  GroundWaterOut$Volum5 = as.matrix(RiverGet) + as.matrix(InterFlowFlux5_23) ##as.matrix(InterFlowFlux5_4) + 
  return(GroundWaterOut)
}




