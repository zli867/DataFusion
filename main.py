from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo, GCGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import writeToNetCDF


if __name__ == '__main__':
    # # CMAQ Example
    # CMAQ_file = "/home/zli867/DFManuscript/DataFusion/data/fused/CCTM_ACONC_fused_EQUATES_2010_2019_PM_O3.nc"
    # obs_file = "/home/zli867/DataFusion/data/CMAQ/NO2_obs_2010_2019_CONUS_ppb.csv"
    # data_fusion_output = "/home/zli867/DFManuscript/DataFusion/data/fused/CCTM_ACONC_fused_EQUATES_2010_2019_PM_O3_NO2.nc"
    # CMAQ_pollutant = "NO2_AVG"
    # obs_pollutant = "NO2"

    # geo = CMAQGridInfo(CMAQ_file)
    # CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)
    # OBSList = extractOBS(obs_pollutant, obs_file, geo)
    # datafusion_results = dataFusion(CMAQList, OBSList, geo, 'default', 10)
    # writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fusion_output)
    
    # GC Example
    GC_file = "/Volumes/DataStorage/DataFusionData/GC/GC_Combine_2015.nc4"
    obs_file = "/Volumes/DataStorage/DataFusionData/GC/O3_obs_2010_2019_CONUS_GC.csv"
    data_fusion_output = "/Volumes/DataStorage/DataFusionData/test/GC_fused.nc4"
    GC_pollutant = "O3_MDA8"
    obs_pollutant = "O3"

    geo = GCGridInfo(GC_file) # the only change is here
    CMAQList = extractCMAQ(GC_pollutant, GC_file, geo)
    OBSList = extractOBS(obs_pollutant, obs_file, geo)
    datafusion_results = dataFusion(CMAQList, OBSList, geo, 'default')
    writeToNetCDF(datafusion_results, GC_pollutant, GC_file, data_fusion_output)