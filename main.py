from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo, GCGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import writeToNetCDF

# CMAQ Example
CMAQ_file = "/Volumes/DataStorage/DataFusionData/CMAQ/CCTM_ACONC_combined_EQUATES_2010_2019.nc"
obs_file = "/Volumes/DataStorage/DataFusionData/CMAQ/PM25_obs_2010_2019_CONUS_ug.csv"
data_fusion_output = "/Volumes/DataStorage/DataFusionData/fused/CCTM_ACONC_fused_EQUATES_2010_2019_PM.nc"
CMAQ_pollutant = "PM25_AVG"
obs_pollutant = "PM25"

geo = CMAQGridInfo(CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)
OBSList = extractOBS(obs_pollutant, obs_file, geo)

datafusion_results = dataFusion(CMAQList, OBSList, geo)
writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fusion_output)

# # GC Example
# GC_file = "/Volumes/DataStorage/DataFusionData/fused/GC_Combine_2015_O3_fused.nc4"
# obs_file = "/Volumes/DataStorage/DataFusionData/GC/NO2_obs_2010_2019_CONUS_GC.csv"
# data_fusion_output = "/Volumes/DataStorage/DataFusionData/fused/GC_Combine_2015_fused.nc4"
# CMAQ_pollutant = "NO2_AVG"
# obs_pollutant = "NO2"

# geo = GCGridInfo(GC_file)
# CMAQList = extractCMAQ(CMAQ_pollutant, GC_file, geo)
# OBSList = extractOBS(obs_pollutant, obs_file, geo)

# datafusion_results = dataFusion(CMAQList, OBSList, geo)
# writeToNetCDF(datafusion_results, CMAQ_pollutant, GC_file, data_fusion_output)