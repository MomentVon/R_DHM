###### runoff ###############
fctMoistureCapacity <- function(SoilMoistureVolume, SoilMoistureCapacityMax){
  SoilMoistureCapacity = SoilMoistureCapacityMax * (1 - (1 - SoilMoistureVolume / SoilMoistureCapacityMax * (paramSoilMoistureCapacityB + 1))^(1 / (paramSoilMoistureCapacityB + 1)))  #A is the fraction of an area for which the soil moisture capacity is less than or equal to i
  SoilMoistureCapacity[is.na(SoilMoistureCapacity)] = SoilMoistureCapacityMax[is.na(SoilMoistureCapacity)]
  return(SoilMoistureCapacity)
}

fctSaturatedArea <- function(SoilMoistureCapacity, SoilMoistureCapacityMax){
  SaturatedArea = 1.0 - (1 - SoilMoistureCapacity / SoilMoistureCapacityMax)^paramSoilMoistureCapacityB
  return(SaturatedArea)
} 


###### InfiltrationGreenAmpt ##############
fctInfiltrationGreenAmpt <- function(HydraulicConductivity, WettingFrontSoilSuction, SoilMoistureContent, EffectivePorosity, SoilMoistureVolume){
  InfiltrationRat = HydraulicConductivity * (1 + WettingFrontSoilSuction * (EffectivePorosity - SoilMoistureContent) / SoilMoistureVolume)
  return(maxSVector(0.0, InfiltrationRat))  ## P > Ks
}

###### InfiltrationSER ################
fctSoilInfiltrationSER <- function(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureCapacity){
  TEMMin = SoilMoistureCapacityMax / (paramSoilMoistureCapacityB + 1) * ((1 - SoilMoistureCapacity / SoilMoistureCapacityMax)^(paramSoilMoistureCapacityB + 1) - (1 - (SoilMoistureCapacity + PrecipitationHoch) / SoilMoistureCapacityMax)^(paramSoilMoistureCapacityB + 1))
  TEMMax = SoilMoistureCapacityMax / (paramSoilMoistureCapacityB + 1) * (1 - SoilMoistureCapacity / SoilMoistureCapacityMax)^(paramSoilMoistureCapacityB + 1)
  TEMDiff = PrecipitationHoch - (SoilMoistureCapacityMax - SoilMoistureCapacity) 
  TEM = TEMMin
  TEM[which(TEMDiff > 0.0)] = TEMMax[which(TEMDiff > 0.0)]
  return(TEM)
}
###### InfiltrationOIER ######################
fctSoilInfiltrationOIER <- function(PrecipitationHoch, InfiltrationRateMax){
  TEMMin = InfiltrationRateMax / (paramInfiltrationRateB + 1) * (1 - (1 - PrecipitationHoch / (InfiltrationRateMax + 0.0000000001))^(paramInfiltrationRateB + 1))
  TEMMax = InfiltrationRateMax / (paramInfiltrationRateB + 1)
  TEMDiff = PrecipitationHoch - InfiltrationRateMax 
  TEM = TEMMin
  TEM[which(TEMDiff > 0.0)] = TEMMax[which(TEMDiff > 0.0)]
  return(TEM)
}

###### SaturationExcessRunoff, return a list mit $Runoff and $Infiltration #########
RUNOFF_SER <- function(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureVolum){
  SoilMoistureCapacity = fctMoistureCapacity(SoilMoistureVolum, SoilMoistureCapacityMax)
  SoilInfiltrationSER = fctSoilInfiltrationSER(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureCapacity)
  SaturationExcessRunoff = PrecipitationHoch - SoilInfiltrationSER
  return(list(Runoff = SaturationExcessRunoff, Infiltration = SoilInfiltrationSER))
}

###### OverInfiltrationExcessRunoff, return a list mit $Runoff and $Infiltration #########
RUNOFF_OIER <- function(PrecipitationHoch, InfiltrationRateMax){
  SoilInfiltrationOIER = fctSoilInfiltrationOIER(PrecipitationHoch, InfiltrationRateMax)
  OverInfiltrationExcessRunoff = maxSVector(0.0,PrecipitationHoch - SoilInfiltrationOIER)
  return(list(Runoff = OverInfiltrationExcessRunoff, Infiltration = SoilInfiltrationOIER))
}

###### VIC, return a list mit $Runoff and $Infiltration #########
RUNOFF_VIC <- function(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureVolume, InfiltrationRateMax){
  SoilMoistureCapacity = fctMoistureCapacity(SoilMoistureVolume, SoilMoistureCapacityMax)
  fctProcess <- function(rate, PrecipitationHoch. = PrecipitationHoch[i], SoilMoistureCapacityMax. = SoilMoistureCapacityMax[i], SoilMoistureCapacity. = SoilMoistureCapacity[i], InfiltrationRateMax. = InfiltrationRateMax[i]){
    SoilInfiltrationSER = fctSoilInfiltrationSER(rate * PrecipitationHoch., SoilMoistureCapacityMax., SoilMoistureCapacity.)
    Pr_RunoffSER = rate * PrecipitationHoch. - SoilInfiltrationSER
    SaturatedArea = fctSaturatedArea(SoilMoistureCapacity. + rate * PrecipitationHoch., SoilMoistureCapacityMax.)
    SoilInfiltrationOIER = fctSoilInfiltrationOIER(PrecipitationHoch. - Pr_RunoffSER, InfiltrationRateMax.)
    Pr_RunoffOIER = PrecipitationHoch. - Pr_RunoffSER - SoilInfiltrationOIER
    returnTem = Pr_RunoffOIER + rate * PrecipitationHoch.
    returnTem[is.na(returnTem)] = PrecipitationHoch.
    return(returnTem)
  }
  nVector = length(PrecipitationHoch)
  RunoffVIC = SoilInfiltrationVIC = SoilInfiltrationVICSER = SoilInfiltrationVICOIER = RatSER = array(0.0, dim = c(nVector,1))
  for (i in 1:nVector) {
    if(PrecipitationHoch[i] == 0.0) {
      SoilInfiltrationVIC[i] = 0.0
      RunoffVIC[i] = 0.0
    }
    else {
      RatSER[i] = eindim618funk(0, 1, PrecipitationHoch[i], processfk = fctProcess, umkriesn = 5)
      SoilInfiltrationVICSER[i] = fctSoilInfiltrationSER(RatSER[i] * PrecipitationHoch[i], SoilMoistureCapacityMax[i], SoilMoistureCapacity[i])
      SoilInfiltrationVICOIER[i] = fctSoilInfiltrationOIER((1 - RatSER[i]) * PrecipitationHoch[i], InfiltrationRateMax[i])
      SoilInfiltrationVIC[i] = SoilInfiltrationVICSER[i] + SoilInfiltrationVICOIER[i]
      RunoffVIC[i] = PrecipitationHoch[i] - SoilInfiltrationVIC[i]
    }
  }
  SoilInfiltrationVIC[is.na(SoilInfiltrationVIC)] = fctSoilInfiltrationSER(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureCapacity)[is.na(SoilInfiltrationVIC)]
  RunoffVIC = as.matrix(PrecipitationHoch) - SoilInfiltrationVIC
  return(list(Runoff = RunoffVIC, Infiltration = SoilInfiltrationVIC, RatSER = RatSER))
}
###### VitcalMisingExcessRunoff, return a list mit $Runoff and $Infiltration #########
RUNOFF_VM <- function(PrecipitationHoch, SoilMoistureCapacityMax, SoilMoistureCapacity, InfiltrationRateMax){
  SoilInfiltrationOIER = fctSoilInfiltrationOIER(PrecipitationHoch, InfiltrationRateMax)
  OverInfiltrationExcessRunoff = PrecipitationHoch - SoilInfiltrationOIER
  SoilInfiltrationSER = fctSoilInfiltrationSER(SoilInfiltrationOIER, SoilMoistureCapacityMax, SoilMoistureCapacity)
  SaturationExcessRunoff = SoilInfiltrationOIER - SoilInfiltrationSER
  return(list(Runoff = OverInfiltrationExcessRunoff + SaturationExcessRunoff, Infiltration = SoilInfiltrationSER))
}






