###### baseflow ########
BASEFLOW_ARNO <- function(SoilMoistureVolume, SoilMoistureVolumeMax, infoDrainageLossMax. = infoDrainageLossMax, infoDrainageLossMin. = infoDrainageLossMin, paramExponentARNOBase. = paramExponentARNOBase, paramSoilMoistureVolumeARNOBaseThresholdRadio. = paramSoilMoistureVolumeARNOBaseThresholdRadio){
  SoilMoistureVolumeARNOBaseThreshold = paramSoilMoistureVolumeARNOBaseThresholdRadio. * SoilMoistureVolumeMax
  TEMMin = infoDrainageLossMin. * SoilMoistureVolume / SoilMoistureVolumeMax
  TEMMax = TEMMin + (infoDrainageLossMax. - infoDrainageLossMin.) * ((SoilMoistureVolume - SoilMoistureVolumeARNOBaseThreshold) / (SoilMoistureVolumeMax - SoilMoistureVolumeARNOBaseThreshold))^paramExponentARNOBase.
  TEMDiff = SoilMoistureVolume - SoilMoistureVolumeARNOBaseThreshold
  TEM = TEMMin
  TEM[which(TEMDiff > 0.0)] = TEMMax[which(TEMDiff > 0.0)]
  return(minVector(SoilMoistureVolume, TEM))
}

