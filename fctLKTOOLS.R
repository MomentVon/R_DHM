Station2GridIDW <- function(StationLocation, GridLocation, StationData){
  fctlengthtem <- function(a,b){
    return((a-b)^2)
  }
  
  periodN = dim(StationData)[1]
  stationN = dim(StationData)[2]
  fieldN = dim(StationData)[3]
  gridN = dim(GridLocation)[1]
  
  GridData = array(0.0,dim = c(periodN, gridN, fieldN))
  LengthG2STemLong = outer(GridLocation$longitude, StationLocation$longitude, fctlengthtem)
  LengthG2STemLati = outer(GridLocation$latitude, StationLocation$latitude, fctlengthtem)
  LengthG2STem = LengthG2STemLong + LengthG2STemLati
  LengthG2S = LengthG2STem^0.5
  LengthG2Ssamllst5S = array(0.0, c(gridN, stationN))
  for (i in 1:gridN) {
    LengthG2Ssamllst5S[i,order(LengthG2S[i,])[1:5]] = LengthG2S[i,order(LengthG2S[i,])[1:5]]
  }
  LengthSum = rowSums(LengthG2Ssamllst5S)
  WeightG2S5S = array(0.0, c(gridN, stationN))
  for (i in 1:gridN) {
    WeightG2S5S[i,] = LengthG2Ssamllst5S[i,] / LengthSum[i]
  }
  
  for (i in 1:fieldN) {
    GridData[,,i] = StationData[,,i] %*% t(WeightG2S5S)
  }
  
  return(GridData)
}

toNDayofYear <- function(DateDay){
  dn = length(DateDay)
  NDay = array(-1, dim = c(dn))
  for (i in 1:dn) {
    NDay[i] = as.integer(DateDay[i] - as.Date(paste(substr(DateDay[i],1,4),1,1,sep = "-"))) + 1
  }
  return(NDay)
}

discrete_convolution <- function(VectorSmall,VectorBig){
  lengthSmall = length(VectorSmall)
  lengthBig = length(VectorBig)
  Tem = matrix(0.0,(lengthBig + lengthSmall -1), lengthSmall)
  for (i in 1:lengthSmall) {
    Tem[i:(lengthBig-1+i),i] = VectorBig * VectorSmall[i]
  }
  Result = rowSums(Tem)
  return(Result)
}

interpolation <- function(bktprm,BKTPRMLST,ZLPRMLST){
  
  n=length(BKTPRMLST)
  if(bktprm < BKTPRMLST[1]){
    zp=(bktprm - BKTPRMLST[1])*(ZLPRMLST[1] - ZLPRMLST[2])/(BKTPRMLST[1] - BKTPRMLST[2]+0.0000000001) + ZLPRMLST[1]
    
  }
  if(bktprm > BKTPRMLST[n]){
    zp=(bktprm - BKTPRMLST[n-1])*(ZLPRMLST[n-1] - ZLPRMLST[n])/(BKTPRMLST[n-1] - BKTPRMLST[n]+0.0000000001) + ZLPRMLST[n-1]
    
  }
  if(bktprm >= BKTPRMLST[1]) {
    for(i in 1:n){
      if(bktprm < BKTPRMLST[i]){
        zp=(bktprm - BKTPRMLST[i-1])*(ZLPRMLST[i-1] - ZLPRMLST[i])/(BKTPRMLST[i-1] - BKTPRMLST[i]+0.000000001) + ZLPRMLST[i-1]
        break()
      }
    }
  }
  return(zp)
}
########## reclassify ######
fctClassify <- function(filePath, classifyN, leftDropN = 0, upDropN = 0){
  MetaData = read.table(filePath, nrows = 6)
  valueNoData = MetaData[6,2]
  colN = MetaData[1,2]
  rowN = MetaData[2,2]
  OriginalWert = read.table(filePath, skip = 6 + upDropN, na.strings = valueNoData)[,(1 + leftDropN):colN]
  GridWert = array(valueNoData,c(1,3))
  for (i in 1:infoGridColN) {
    for (j in 1:infoGridRowN) {
      if(GridID[j,i] != valueNoData){
        LanduseTem = as.matrix(OriginalWert[(1 + (j - 1) * classifyN):(j * classifyN), (1 + (i - 1) * classifyN):(i * classifyN)])
        LanduseArtTem = table(unlist(LanduseTem))
        lengthTem = length(LanduseArtTem)
        MatrixTem = array(GridID[j,i], c(lengthTem,3))
        MatrixTem[,2] = as.integer(dimnames(LanduseArtTem)[[1]])
        RateTem = as.integer(LanduseArtTem)
        rateTem2 = sum(RateTem)
        MatrixTem[,3] = RateTem / rateTem2
        GridWert = rbind(GridWert,MatrixTem)
      }
      else next
    }
  }
  return(GridWert[-1,])
}
