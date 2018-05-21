###### other paramter caculate ###############
fctMoistureCapacity <- function(SoilMoistureVolume, SoilMoistureCapacityMax){
  SoilMoistureCapacity = SoilMoistureCapacityMax * (1 - (1 - SoilMoistureVolume * (paramSoilMoistureCapacityB + 1) / SoilMoistureCapacityMax)^(1 / (paramSoilMoistureCapacityB + 1)))  #A is the fraction of an area for which the soil moisture capacity is less than or equal to i
  SoilMoistureCapacity[is.na(SoilMoistureCapacity)] = (SoilMoistureVolume / SoilMoistureCapacityMax)[is.na(SoilMoistureCapacity)]
  return(SoilMoistureCapacity)
}

fctSaturatedArea <- function(SoilMoistureCapacity, SoilMoistureCapacityMax){
  SaturatedArea = 1.0 - (1 - SoilMoistureCapacity / SoilMoistureCapacityMax)^paramSoilMoistureCapacityB
  return(SaturatedArea)
} 
fctAerodynamicResistance <- function(VegetationDisplacementHeight, VegetationRoughLength, WindSpeed){
  AerodynamicResistance = 4.72 * (log(VegetationDisplacementHeight / VegetationRoughLength))^2 / (1 + 0.54 * WindSpeed)
  return(AerodynamicResistance)
}
###### ET_PenmanMonteith UNIT:[date, m, o, oC, oC, oC, m/s, m, W/m2, %(0-100), %(0-100)] ################
ReferenceET_PenmanMonteith <- function(JDay, Elevation, Latitude, Tmax, Tmin, Tmean,  Wind, WindH, SunHour, RelativeHumidity){
  # DailySolarRadiationMJ = 0.0864 * DailySolarRadiationW  #Equation 6  ##1W m-2 = 0.0864 MJ m-2
  Wind2 = Wind * 4.87 / log(67.8 * WindH - 5.42)  #Equation 7
  Delt = 4098 * (0.6108 * exp(17.27 * Tmean / (Tmean + 237.3))) / (Tmean + 237.3)^2  #Equation 9
  Pressure = 101.3 * ((293 - 0.0065 * Elevation) / 293)^5.26  #Equation 10
  Gama = 0.000665 * Pressure  #Equation 11  ##γ = psychrometric constant, kPa °C-1;P = atmospheric pressure, kPa, [Eq. 10];λ = latent heat of vaporization, 2.45, MJ kg-1;cp = specific heat at constant pressure, 1.013 10-3, MJ kg-1°C-1;μ = ratio molecular weight of water vapour/dry air = 0.622.
  DeltTerm = Delt / (Delt + Gama * (1 + 0.34 * Wind2))  #Equation 12
  PressureTerm = Gama / (Delt + Gama * (1 + 0.34 * Wind2))  #Equation 13
  TemperatureTerm = (900 / (Tmean + 273)) * Wind2  #Equation 14
  ETemperature <- function(Temperature){
    return(0.6108 * exp(17.27 * Temperature / (Temperature + 237.3)))
  }  #Equation 15
  SaturationVaporPressure = (ETemperature(Tmax) + ETemperature(Tmin)) / 2  #Equation 18
  # ActualVaporPressure = (ETemperature(Tmin) * MaximumRelativeHumidity / 100 + ETemperature(Tmax) * MinimumRelativeHumidity / 100) /2  #Equation 19
  ActualVaporPressure = RelativeHumidity / 100.0 * (ETemperature(Tmax) + ETemperature(Tmin)) / 2
  EarthSun = 1 + 0.033 * cos(2 * pi * JDay / 365.25)  #Equation 23
  SolarDeclination = 0.409 * sin(2* pi *JDay / 365.25 - 1.39)  #Equation 24
  LatitudeRad = pi * Latitude / 180 #Equation 25
  SunsetHourAngle = acos(-1 * tan(LatitudeRad) * tan(SolarDeclination))  #Equation 26
  ExtraterrestrialRadiation = 24 * 60 / pi * 0.082 * EarthSun * (SunsetHourAngle * sin(LatitudeRad) * sin(SolarDeclination) + cos(LatitudeRad) * cos(SolarDeclination) * sin(SunsetHourAngle))  #Equation 27 ##Gsc = solar constant = 0.0820 MJ m-2 min-1;
  ClearSkySolarRadiation = ExtraterrestrialRadiation * (0.75 + 2*10^-5 * Elevation)  #Equation 28
  PossibleSunshineHours = 24 / pi * SunsetHourAngle
  DailySolarRadiationMJ = ExtraterrestrialRadiation * (0.14233 + 0.658667 * (SunHour / PossibleSunshineHours) + 0.028  * (77.6667 / 100) - 0.11533* (EarthSun / PossibleSunshineHours) * (77.6667 / 100)) ### in zhaona: Solar radiation estimation using sunshine hour and air pollution index in China
  NetShortwaveRadiation = (1 - 0.23) * DailySolarRadiationMJ  #Equation 29  ##α = albedo or canopy reflection coefficient, which is 0.23 for the hypothetical grass reference crop, dimensionless
  NetOutgoingLongwaveRadiation = 4.903*10^-9 * (((Tmax + 273.16)^4 + (Tmax + 273.16)^4) / 2) * (0.34 - 0.14 * ActualVaporPressure^0.5) * (1.35 * DailySolarRadiationMJ / ClearSkySolarRadiation - 0.35)  #Equation 30  ##σ = Stefan-Boltzmann constant [4.903 10-9 MJ K-4 m-2 day-1],
  NetRadiation = NetShortwaveRadiation - NetOutgoingLongwaveRadiation  #Equation 31
  EquivalentNetRadiation = 0.408 * NetRadiation  #Equation 32
  ETRadiation = DeltTerm * EquivalentNetRadiation #Equation 33
  ETWind = PressureTerm * TemperatureTerm * (SaturationVaporPressure - ActualVaporPressure)  #Equation 34
  ReferenceEvapotranspiration = ETRadiation + ETWind  #Equation 35
  return(ReferenceEvapotranspiration)
}

###### ET_Hargreaves UNIT:[date, m, o, oC, oC, oC] ###########
ReferenceET_Hargreaves <- function(JDay, Latitude, Tmax, Tmin, Tmean){
  EarthSun = 1 + 0.033 * cos(2 * pi * JDay / 365.25)  #Equation 23
  SolarDeclination = 0.409 * sin(2* pi *JDay / 365.25 - 1.39)  #Equation 24
  LatitudeRad = pi * Latitude / 180 #Equation 25
  SunsetHourAngle = acos(-1 * tan(LatitudeRad) * tan(SolarDeclination))  #Equation 26
  ExtraterrestrialRadiation = 24 * 60 / pi * 0.082 * EarthSun * (SunsetHourAngle * sin(LatitudeRad) * sin(SolarDeclination) + cos(LatitudeRad) * cos(SolarDeclination) * sin(SunsetHourAngle))  #Equation 27 ##Gsc = solar constant = 0.0820 MJ m-2 min-1;
  ExtraterrestrialRadiationMM = ExtraterrestrialRadiation / 2.45  #B.25 ##(Ra in mm d-1 = Ra in MJ m-2 d-1 / 2.45).
  ReferenceEvapotranspiration = 0.0023 * (Tmax - Tmin)^0.5 * (Tmean + 17.8) * ExtraterrestrialRadiationMM 
  return(ReferenceEvapotranspiration)
}
###### VIC ###############
ET_VIC <- function(ReferenceEvapotranspiration, PrecipitationHoch, MoistureVolume, MoistureCapacityMax, MoistureVolume1, MoistureCapacityMax1, AerodynamicResistance, ArchitecturalResistance, StomatalResistance){
  EvaporationCnopyMax = (MoistureVolume / MoistureCapacityMax)^0.6667 * (AerodynamicResistance / (AerodynamicResistance + ArchitecturalResistance)) * ReferenceEvapotranspiration
  Rate = minSVector(1, (MoistureVolume + PrecipitationHoch) / (EvaporationCnopyMax + 0.0000000000001))
  Rate = maxSVector(0.0,minSVector(1.0,Rate))
  EvaporationCnopy = minVector(PrecipitationHoch, Rate * EvaporationCnopyMax)
  Rate1 = minSVector(1, (MoistureVolume1 + PrecipitationHoch) / (MoistureCapacityMax1 + 0.0000000000001))
  Rate2 = AerodynamicResistance / ( AerodynamicResistance + ArchitecturalResistance + StomatalResistance)
  Rate1 = maxSVector(0.0,minSVector(1.0,Rate1))
  Rate2 = maxSVector(0.0,minSVector(1.0,Rate2))
  Transpiration = ((1.0 - Rate1) * Rate2 + Rate1 * Rate2 * (1 - (MoistureVolume1 / MoistureCapacityMax1)^0.667)) * ReferenceEvapotranspiration
  Transpiration[is.na(Transpiration)] = 0.0
  Transpiration = minVector(MoistureVolume1, Transpiration)
  MoistureCapacity = fctMoistureCapacity(MoistureVolume1, MoistureCapacityMax1 * (paramSoilMoistureCapacityB + 1))
  SaturatedArea = fctSaturatedArea(MoistureCapacity, MoistureCapacityMax1 * (paramSoilMoistureCapacityB + 1))
  # fctRate3 <- function(A){
  #   return(1 + (1 - (1 - AS)^(1 / paramSoilMoistureCapacityB)) / (1 - (1 - A)^(1 / paramSoilMoistureCapacityB)))
  # }
  # Rate3 = array(0.0, length(SaturatedArea))
  # for (i in 1:length(SaturatedArea)) {
  #   AS <<- SaturatedArea[i]
  #   Rate3[i] = AS + integrate(fctRate3,min(0.001,SaturatedArea[i]),SaturatedArea[i])$value
  # }
  # 
  # Rate3 = SaturatedArea * (1.0 + log(SaturatedArea))
  Rate3 = maxSVector(0.0,minSVector(1.0,SaturatedArea + (1 - SaturatedArea) * 0.3))
  # print(c("rate3",mean(SaturatedArea),mean(Rate3)))
  EvaporationSoil = Rate3 * ReferenceEvapotranspiration
  EvaporationSoil = minVector(MoistureVolume1 - Transpiration,EvaporationSoil)
  return(list(EvapC = EvaporationCnopy, Transp = Transpiration, EvapS = EvaporationSoil))
}


