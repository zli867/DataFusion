from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo
from DataExtraction.ExtractOBS import extractOBS
from Evaluation.StatisticalMetrics import stats_metrics


original_CMAQ_file = "./data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
obs_file = "./data/obs_2016_2020.csv"
fused_CMAQ_file = "./results/CCTM.ACONC.combined.FIRE.hires4.2016_2020_fused.nc"
CMAQ_pollutant = "PM25_TOT_AVG"
obs_pollutant = "PM25"

geo = CMAQGridInfo(original_CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, original_CMAQ_file)
OBSList = extractOBS(obs_pollutant, obs_file, geo)
fused_CMAQ_List = extractCMAQ(CMAQ_pollutant, fused_CMAQ_file)

cmaq_performance = stats_metrics(CMAQList[0], OBSList[0], geo)
fused_cmaq_performance = stats_metrics(fused_CMAQ_List[0], OBSList[0], geo)
print(cmaq_performance)
print(fused_cmaq_performance)

