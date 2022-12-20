from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import writeToNetCDF


CMAQ_file = "./data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
obs_file = "./data/obs_2016_2020.csv"
data_fusion_output = "./results/CCTM.ACONC.combined.FIRE.hires4.2016_2020_fused.nc"
CMAQ_pollutant = "PM25_TOT_AVG"
obs_pollutant = "PM25"

geo = CMAQGridInfo(CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file)
OBSList = extractOBS(obs_pollutant, obs_file, geo)

datafusion_results = dataFusion(CMAQList, OBSList, geo)
writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fusion_output)