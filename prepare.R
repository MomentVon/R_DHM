###### Prepra rout Infomation ####
fctGIUH_Nash <- function(t){
  ut = 1.0 / (coefficientStorage * gamma(paramShape)) * (t / coefficientStorage)^(paramShape - 1) * exp(-1.0 * t / coefficientStorage)
  return(ut)
}

fctGIUH_Shipeng <- function(t){
  ut = t / (2 * paramRiverheadNumber * (paramStreamLength / WaveVelocity)^2) * exp(-1.0 * WaveVelocity^2 * t^2 / (4 * paramRiverheadNumber * paramStreamLength^2))
  return(ut)
}

fctIUH_RuiXF<-function(t){
  ut = RiverLength / (4 * pi * AttenuationCoeffcient * t^3)^0.5 * exp(-1.0 * (WaveVelocity * t - RiverLength)^2/(4 * AttenuationCoeffcient * t)) #km2/h
  return(ut)
}

###### find every Grid which Grid will go (confluence way) prepare for UH make######
fctGrid2AimGridFind<-function(l,m, AimGridID, GridID. = GridID, FlowDirection. = FlowDirection){
  n = 0
  j = l
  k = m
  while((GridID.[j,k]!=constNull)){
    f=FlowDirection.[j,k]
    if(f==1) {
      j=j-1
      k=k
      n=n+1
    }
    if(f==2) {
      j=j-1
      k=k+1
      n=n+1.42
    }
    if(f==3) {
      k=k+1
      j=j
      n=n+1
    }
    if(f==4) {
      j=j+1
      k=k+1
      n=n+1.42
    }
    if(f==5) {
      j=j+1
      k=k
      n=n+1
    }
    if(f==6) {
      j=j+1
      k=k-1
      n=n+1.42
    }
    if(f==7) {
      k=k-1
      j=j
      n=n+1
    }
    if(f==8) {
      j=j-1
      k=k-1
      n=n+1.42
    }
    j=j
    k=k
    if(GridID.[j,k] %in% as.matrix(AimGridID))  break
  }
  
  tem=array(0,dim = c(3))
  tem[1]=j
  tem[2]=k
  tem[3]=n
  return(tem)
}

###### return a Matrix wirh:1. GridID, 2.AimGridID, 3.length(Unit whith not m/km, but ist 1Unit = Grid length), 4.DiffElevation(Unit is m) ###########
fctGrid2AimGrid <- function(OriginalGridID, AimGridID, GridID. = GridID, GridDEM. = GridDEM){
  Grid2AimGrid=matrix(0, dim(OriginalGridID)[1] - dim(AimGridID)[1], 4)
  i=1
  for(m in 1:infoGridColN){
    for(l in 1:infoGridRowN){
      judge=GridID.[l,m]
      if((judge %in% setdiff(as.matrix(OriginalGridID), as.matrix(AimGridID))) == F) next
      Grid2AimGrid[i,1]=GridID.[l,m]
      demHigh=GridDEM.[l,m]
      tem=fctGrid2AimGridFind(l,m, AimGridID)
      Grid2AimGrid[i,2]=GridID.[tem[1],tem[2]]
      demLow=GridDEM.[tem[1],tem[2]]
      Grid2AimGrid[i,3]=tem[3]
      Grid2AimGrid[i,4]=max(0.00000001, (demHigh - demLow))
      i=i+1
    }
  }
  return(Grid2AimGrid)
}

###### make the TranslateMatrix, which will mit Acculate from GridList to AimGridList. #######
fctMakeGridTranslateMatrix <- function(Grid2AimGrid, AimID){
  GridN = dim(Grid2AimGrid)[1]
  AimGridN = dim(AimID)[1]
  TranslateMatrix = matrix(0.0, GridN, AimGridN)
  for(k in 1:GridN) {
    TranslateMatrix[k, match(Grid2AimGrid[k,2],as.matrix(AimID))] = 1
  }
  return(TranslateMatrix)
}

###### UHParame include 5.RiverVelocity, 6.StreamLength, 7.StreamLengthNextLowOder, 8.StreamBifurcation, 9.StreamBifurcationNextHighOder, 10.StreamArea, 11.StreamAreaNextLowOder ############## 
fctMakeUH_Nash <- function(Grid2AimGridPlusUHParam, UHPeriodN. = UHPeriodN){
  GridN = dim(Grid2AimGridPlusUHParam)[1]
  UH = matrix(0.0, UHPeriodN., GridN)
  for(j in 1:GridN){
    RiverLength = Grid2AimGridPlusUHParam[j,3]
    RiverVelocity = Grid2AimGridPlusUHParam[j,5]
    StreamLength = Grid2AimGridPlusUHParam[j,6]
    StreamLengthNextLowOder = Grid2AimGridPlusUHParam[j,7]
    StreamBifurcation = Grid2AimGridPlusUHParam[j,8]
    StreamBifurcationNextHighOder = Grid2AimGridPlusUHParam[j,9]
    StreamArea = Grid2AimGridPlusUHParam[j,10]
    StreamAreaNextLowOder = Grid2AimGridPlusUHParam[j,11]
    ratioLength = StreamLength / StreamLengthNextLowOder
    ratioBifurcation = StreamBifurcation / StreamBifurcationNextHighOder
    ratioArea = StreamArea / StreamAreaNextLowOder
    coefficientStorage <<- 3.29 * ratioBifurcation^0.78 * ratioArea^-0.78 * ratioLength^0.07
    paramShape <<- 0.7 * RiverLength / RiverVelocity * ratioBifurcation^-0.48 * ratioArea^-0.48 * ratioLength^-0.48
    bondh = bondl = 0.0 
    for(i in 1:UHPeriodN.){
      bondh = bondl+1
      Term = integrate(fctGIUH_Nash, bondl, bondh)$value
      UH[i,j] = round(Term,7)
      bondl = bondh
    }
  }
  return(UH)
}
###### UHParame include 5.coefficientStorage, 6.paramShape #######
fctMakeUH_Nash_simple <- function(Grid2AimGridPlusUHParam, UHPeriodN. = UHPeriodN){
  GridN = dim(Grid2AimGridPlusUHParam)[1]
  UH = matrix(0.0, UHPeriodN., GridN)
  for(j in 1:GridN){
    coefficientStorage <<- Grid2AimGridPlusUHParam[j,5]
    paramShape <<- Grid2AimGridPlusUHParam[j,6]
    bondh = bondl = 0.0 
    for(i in 1:UHPeriodN.){
      bondh = bondl+1
      Term = integrate(fctGIUH_Nash, bondl, bondh)$value
      UH[i,j] = round(Term,7)
      bondl = bondh
    }
  }
  return(UH)
}

###### UHParame include 5.average paramStreamLength(内链平均长度), 6.paramRiverheadNumber(外链个数), 7.WaveVelocity ###########
fctMakeUH_Shipeng <- function(Grid2AimGridPlusUHParam, UHPeriodN. = UHPeriodN){
  GridN = dim(Grid2AimGridPlusUHParam)[1]
  UH = matrix(0.0, UHPeriodN., GridN)
  for(j in 1:GridN){
    paramStreamLength <<- Grid2AimGridPlusUHParam[j,3]
    paramRiverheadNumber <<- 1
    WaveVelocity <<- 2
    
    bondh = bondl = 0.0 
    for(i in 1:UHPeriodN.){
      bondh = bondl+1
      Term = integrate(fctGIUH_Shipeng,bondl,bondh)$value
      UH[i,j] = round(Term,7)
      bondl = bondh
    }
  }
  return(UH)
}

###### UHParame include 5.Discharge, 6.RiverWidth, 7.RiverVelocity ##############
fctMakeUH_RuiXF <- function(Grid2AimGridPlusUHParam, UHPeriodN. = UHPeriodN){
  GridN = dim(Grid2AimGridPlusUHParam)[1]
  UH = matrix(0.0, UHPeriodN., GridN)
  for(j in 1:GridN){
    RiverLength <<- 1.0 * (Grid2AimGridPlusUHParam[j,3]) * infoGridlength #km
    WaterSurfaceGradient = 1.0*(Grid2AimGridPlusUHParam[j,4] / RiverLength)/1000.0 #km/km
    Discharge = Grid2AimGridPlusUHParam[j,5] #m3/s
    RiverWidth = Grid2AimGridPlusUHParam[j,6] #m
    RiverVelocity = Grid2AimGridPlusUHParam[j,7] #m/s
    AttenuationCoeffcient <<- Discharge / (2 * WaterSurfaceGradient * RiverWidth)
    WaveVelocity <<- 1.6 * RiverVelocity #Manning is 5/3 but checy is 1.5, so i take 1.6
    
    bondh = bondl = 0.0 
    for(i in 1:UHPeriodN.){
      bondh = bondl + 1
      Term = integrate(fctIUH_RuiXF, bondl, bondh)$value
      UHs[i,j] = round(Term,7)
      bondl = bondh
    }
  }
  return(UH)
}

###### cut GridFlow from AllGridFlow to AimGridFlow and OtherGridFlow, from 1 cut to 2. Get list with $OtherGridFlow and $AimGridFlow###########
fctCutGridFlow <- function(OriginalGridFlow, GridID, AimID){
  FlowTem = as.data.frame(cbind(GridID, t(OriginalGridFlow)))
  names(FlowTem) = c("id",as.character(seq(1,infoPeriodN,1)))
  DiffID = as.data.frame(setdiff(as.matrix(GridID), as.matrix(AimID)))
  names(DiffID) = "id"
  DFAimID = as.data.frame(AimID)
  names(DFAimID) = "id"
  
  JoinTem = join(DiffID, FlowTem)
  FlowOtherGrid=t(as.matrix(JoinTem[,-1])) 
  
  JoinTem = join(DFAimID, FlowTem)
  FlowAimGrid=t(as.matrix(JoinTem[,-1]))    #the surfaceflow from VIC with subday just only other grids
  
  return(list(FlowOtherGrid = FlowOtherGrid, FlowAimGrid = FlowAimGrid))
}
###### caculate the Runoff mit selbst UH, but not have to the AimGrid, just know how many Confluence have to its AimGrid######
fctUHConfluence <- function(Flow, UH){
  GridN = dim(Flow)[2]
  FlowResult = matrix(0.0,infoPeriodN,GridN)    #the result of baseflow confluence
  for(k in 1:GridN){
    FlowResult[,k] = discrete_convolution(UH[,k], Flow[,k])[1:infoPeriodN]
  }
  return(FlowResult)
}















