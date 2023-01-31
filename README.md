## Authorship
* The model is based on the paper: 
Method for Fusing Observational Data and Chemical 
Transport Model Simulations To Estimate Spatiotemporally 
Resolved Ambient Air Pollution: https://pubs.acs.org/doi/full/10.1021/acs.est.5b05134
  
* The code is implemented by Zongrun Li (zli867@gatech.edu)

## Folder Structure
* DataExtraction: includes functions to extract observation or CMAQ
    * ExtractCMAQ.py: functions to extract CMAQ
    * ExtractOBS.py: functions to extract observations
    
* DataFusion
    * CMAQatOBS.py: extracts CMAQ values at observation locations
    * CorrelationParameter.py: generates correlation parameters R<sub>coll</sub>, R<sub>1</sub> and R<sub>2</sub>
    * DataFusionDriver.py: does all the routines for data fusion.
    * KrigingOBS.py: kriging interpolation for observations.
    * adjustCMAQ.py: functions for generating parameters for FC<sub>1</sub>.
    * dataFusion.py: functions for generating parameters for FC<sub>2</sub>. Includes helper functions for 
      FC<sub>1</sub>, FC<sub>2</sub> and FC<sub>opt</sub>
    
* DataWriter
    * DataWriter.py: defines a function which writes data fusion results to NetCDF. Define functions here 
      for other format outputs.

* Evaluation
  * StatisticalMetrics.py: defines some statistical metrics for evaluating CMAQ or data fused CMAQ results performance. 
    Current statistical metrics includes R<sup>2</sup> and RMSE.
* DataFusionPerformance.py: an example code for evaluating original CMAQ performance or data fused results performance.

* main.py: an example code for doing data fusion (PM<sub>2.5</sub> is included as an example, output format is NetCDF).
* FusedBurnImpacts: Calculate fire impacts based on the following formula: $$\frac{CMAQ_{fire} - CMAQ_{no\ fire}}{CMAQ_{fire}} \cdot CMAQ_{fused\ fire}$$
* util: provides utilities for plot figures or preprocessing CMAQ results or observations.
    * CombineNetCDF.py: provides a function to combine multiple daily CMAQ outputs to one NetCDF file. 
    * CombineObs.py: provides a function to combine multiple daily observations to one csv file. 
    * GeoProcess.py: provides a function to plot US states.
    
* Visualized_Data_Fusion.ipynb: provides part of visualization for data fusion process.

* results: a folder used to store data fusion results.

* data
    * geo: includes a US shapefile for plotting figures.
    * You can put input data here.
* datafusion.yml: conda environment for this project.

## Build Up
* Set up conda environment 
```
conda env create -f datafusion.yml
```

## Run the Data Fusion Code
### Input data
#### observation data:
Download EPA observation data at: https://aqs.epa.gov/aqsweb/airdata/download_files.html.
Processed the observation data and generate a daily observation data in csv format. 

| Time | obs_pollutant |siteCode|Latitude|Longitude|
| ------ | ----| ------ | ----| ------ |
| YYYY-MM-DD | observation pollutants values| unique for different sites | site location| site location |

* The csv files at least includes the information mentioned above. Notice that the header ```obs_pollutant``` will
be used in ```main.py``` to extract specific observation pollutants.
* Sort the data by time (oldest first). The time standard should match the CMAQ (eg: if CMAQ uses LST, observation 
  should be LST)

#### CMAQ data:
* Use ```hr2day``` to process CMAQ hourly data to daily.
* Combine all the CMAQ data you want to do the data fusion (sort the data by time, oldest first). You can use utilities 
  I provided

### Data fusion
You need to revise following variables to run data fusion:
```
CMAQ_file = "./data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
obs_file = "./data/obs_2016_2020.csv"
data_fusion_output = "./results/CCTM.ACONC.combined.FIRE.hires4.2016_2020_fused_new_env.nc"
CMAQ_pollutant = "PM25_TOT_AVG"
obs_pollutant = "PM25"
```

* CMAQ_file: CMAQ daily output
* obs_file: csv files which follows the standard I mentioned before.
* data_fusion_output: the data fusion output location and filename.
* CMAQ_pollutant: pollutant variable name you want to fused in CMAQ.
* obs_pollutant: pollutant variable name you want to fused in observation file.

