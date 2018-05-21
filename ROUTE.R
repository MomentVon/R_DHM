###### Route ######
fctUHALLMake <- function(fctMakeUH4G2RScf, fctMakeUH4G2RBas, fctMakeUH4R2E, fctMakeUH4E2S){
  UHGrid2RiverScf = fctMakeUH4G2RScf(Grid2River)
  UHGrid2RiverBas = fctMakeUH4G2RBas(Grid2River)
  UHRiver2Estuary = fctMakeUH4R2E(River2Estuary)
  UHEstuary2Station = fctMakeUH4E2S(Estuary2Station)
  return(list(UHGrid2RiverScf = UHGrid2RiverScf, UHGrid2RiverBas = UHGrid2RiverBas, UHRiver2Estuary = UHRiver2Estuary, UHEstuary2Station = UHEstuary2Station))
}


CONFLUENCE <- function(SfcFlowAllGrid, BasFlowAllGrid, UHAll){
  CutGridSfc = fctCutGridFlow(SfcFlowAllGrid, GridGridID, RiverGridID)
  ScfRiverFlow = fctUHConfluence(CutGridSfc$FlowOtherGrid, UHAll$UHGrid2RiverScf) %*% TranslateMatrixGrid2River + CutGridSfc$FlowAimGrid
  
  CutGridBas = fctCutGridFlow(BasFlowAllGrid, GridGridID, RiverGridID)
  BasRiverFlow = fctUHConfluence(CutGridBas$FlowOtherGrid, UHAll$UHGrid2RiverBas) %*% TranslateMatrixGrid2River + CutGridBas$FlowAimGrid
  
  RiverFlow = ScfRiverFlow + BasRiverFlow
  
  CutRiver = fctCutGridFlow(RiverFlow, RiverGridID, EstuaryID)
  EstuaryFlow = fctUHConfluence(CutRiver$FlowOtherGrid, UHAll$UHRiver2Estuary) %*% TranslateMatrixRiver2Estuary + CutRiver$FlowAimGrid
  
  CutEstuary = fctCutGridFlow(EstuaryFlow, EstuaryID, HydroStationID)
  StationFlow = fctUHConfluence(CutEstuary$FlowOtherGrid, UHAll$UHEstuary2Station) %*% TranslateMatrixEstuary2Station + CutEstuary$FlowAimGrid
  
  StationFlowQ = StationFlow * UHUnitTranslate
  return(StationFlowQ)
}




