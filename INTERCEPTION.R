###### Interception ########
INTERCEPTION_Gash <- function(Volum, CanopyStorageCapacity, RainfallDuringSaturation, Evaporation, infoCoefficientFreeThroughfall. = infoCoefficientFreeThroughfall){
  NecessaryToSaturateCanopy = (-1 * RainfallDuringSaturation * CanopyStorageCapacity / Evaporation) * log(1 - Evaporation / (RainfallDuringSaturation * (1 - infoCoefficientFreeThroughfall.)))
  NecessaryToSaturateCanopy[which(RainfallDuringSaturation == 0.0)] = 0.0
  NecessaryToSaturateCanopy[is.nan(NecessaryToSaturateCanopy)] = 0.0
  NecessaryToSaturateCanopy[is.infinite(NecessaryToSaturateCanopy)] = 0.0
  return(minVector(RainfallDuringSaturation, minVector(CanopyStorageCapacity, NecessaryToSaturateCanopy)))
}





