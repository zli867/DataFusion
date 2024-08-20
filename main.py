from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo, GCGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import writeToNetCDF
import tracemalloc
import time

# # CMAQ Example
# CMAQ_file = "/Volumes/DataStorage/DataFusionData/CMAQ/CCTM_ACONC_combined_EQUATES_2010_2019.nc"
# obs_file = "/Volumes/DataStorage/DataFusionData/CMAQ/PM25_obs_2010_2019_CONUS_ug.csv"
# data_fusion_output = "/Volumes/DataStorage/DataFusionData/fused/CCTM_ACONC_fused_EQUATES_2010_2019_PM.nc"
# CMAQ_pollutant = "PM25_AVG"
# obs_pollutant = "PM25"

# geo = CMAQGridInfo(CMAQ_file)
# CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)
# OBSList = extractOBS(obs_pollutant, obs_file, geo)
# datafusion_results = dataFusion(CMAQList, OBSList, geo)
# writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fusion_output)
if __name__ == '__main__':
    # GC Example
    GC_file = "/Volumes/DataStorage/DataFusionData/GC/GC_Combine_2015.nc4"
    obs_file = "/Volumes/DataStorage/DataFusionData/GC/O3_obs_2010_2019_CONUS_GC.csv"
    data_fusion_output = "/Volumes/DataStorage/DataFusionData/fused/GC_Combine_2015_fused_serial_o3.nc4"
    CMAQ_pollutant = "O3_MDA8"
    obs_pollutant = "O3"

    geo = GCGridInfo(GC_file)
    CMAQList = extractCMAQ(CMAQ_pollutant, GC_file, geo)
    OBSList = extractOBS(obs_pollutant, obs_file, geo)
    tracemalloc.start()
    start = time.time()
    datafusion_results = dataFusion(CMAQList, OBSList, geo)
    end = time.time()
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()
    print("The time of execution of above program is :", (end-start) * 10**3, "ms")
    writeToNetCDF(datafusion_results, CMAQ_pollutant, GC_file, data_fusion_output)