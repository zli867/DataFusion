## Authorship
* The model is used to fuse observational data and chemical 
transport model simulations
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
    * DataWriter.py: defines a function that writes data fusion results to NetCDF. Define functions here for other format outputs.

* Evaluation
  * StatisticalMetrics.py: defines statistical metrics for evaluating CMAQ or data fused CMAQ results performance.

* Tools
  * CombineCMAQ.py: a utility that combines different years of CMAQ outputs.
  * CombineObs.py: a utility that combines different observation CSV files.
  * Hr2DayGC.py: a utility to calculate daily average, 1-hr maximum, MDA8, etc. (standard local time) concentration criteria from hourly GC outputs (UTC-0).

* DataFusion_Evaluation.ipynb: utility for evaluating original CMAQ performance or data fused results performance.

* Visualized_Data_Fusion.ipynb: provides part of visualization for the data fusion process.

* main.py: an example code for doing data fusion (PM<sub>2.5</sub> is included as an example. Output format is NetCDF).

* FusedBurnImpacts.py: Calculate fire impacts based on the following formula: $$\frac{CMAQ_{fire} - CMAQ_{no\ fire}}{CMAQ_{fire}} \cdot CMAQ_{fused\ fire}$$

* util: provides utilities for plot figures or preprocessing CMAQ results or observations.
    * GeoProcess.py: provides a function to plot US states.
    * DailyMetrics.py: convert hourly data to daily for different criteria, e.g., daily average, 1-hr maximum, MDA8.
    * VisualizationFunc.py: visualization helper functions.
    * FormatConverter.py: a format converter for the TFLAG variable in CMAQ.

* DataWithholdingCMAQ.py: an example for doing 5-fold cross-validations under site-wise and random withholding strategies on CMAQ simulations.

* DataWithholdingGC.py: an example for doing 5-fold cross-validations under site-wise and random withholding strategies on GC simulations.

* data
    * geo: includes a US shapefile for plotting figures.
    * timezone: includes a shapefile for the timezone region and its timezone.
* datafusion.yml: conda environment for this project (tested on MacOS).
* requirements.txt: conda environment for this project (tested on Linux).

## Build Up
* Set up conda environment 
```
conda env create -f environment.yml
```
or
```
pip install -r requirements.txt
```
## Run the Data Fusion Code
### Input data
#### observation data:
Download EPA observation data at: https://aqs.epa.gov/aqsweb/airdata/download_files.html.
Processed the observation data and generate a daily observation data in csv format. 

| Time | obs_pollutant |siteCode|Latitude|Longitude|
| ------ | ----| ------ | ----| ------ |
| YYYY-MM-DD | observation pollutants values| unique for different sites | site location| site location |

* The CSV files at least include the information mentioned above. Notice that the header ```obs_pollutant``` will
be used in ```main.py``` to extract specific observation pollutants.
* Sort the data by time (oldest first). The time standard should match the CMAQ (e.g.: if CMAQ uses LST, observation 
  should be LST)

#### CMAQ data:
* Use ```hr2day``` to process CMAQ hourly data to daily.
* Combine all the CMAQ data you want to do the data fusion (sort the data by time, oldest first). You can use utilities 
  I provided

### Data fusion
You need to revise the following variables to run data fusion:
```
CMAQ_file = "./data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
obs_file = "./data/obs_2016_2020.csv"
data_fusion_output = "./results/CCTM.ACONC.combined.FIRE.hires4.2016_2020_fused_new_env.nc"
CMAQ_pollutant = "PM25_TOT_AVG"
obs_pollutant = "PM25"
```

* CMAQ_file: CMAQ daily output
* obs_file: CSV files which follows the standard I mentioned before.
* data_fusion_output: the data fusion output location and filename.
* CMAQ_pollutant: pollutant variable name you want to fuse in CMAQ.
* obs_pollutant: pollutant variable name you want to fuse in the observation file.

### Other Chemical Transport Model
1. Create the combined netCDF format daily simulation outputs.
2. Implement a function in DataExtraction/ExtractCMAQ.py. Examples of GEOS-Chem (GCGridInfo), CMAQ (CMAQGridInfo), and WRF-Chem (WRFGridInfo) are provided. The return value is a `dict` which should include the following geographic information of your CTM simulations: 
* crs: coordinate reference system which is used to convert (lat, lon) $\Leftrightarrow$ (Y, X)
* X: the X coordinates of each grid cell's center. It is used to calculate distance.
* Y: the Y coordinates of each grid cell's center. It is used to calculate distance.
* X_bdry: minimum and maximum of the domain's X.
* Y_bdry: minimum and maximum of the domain's Y.
* Lat: latitude of each grid cell's center.
* Lon: longitude of each grid cell's center.
* Lat_bdry: minimum and maximum of the domain's latitude.
* Lon_bdry: minimum and maximum of the domain's longitude.
* time: timestamp for the simulation time dimension.

3. In main.py, replace 
```
geo = CMAQGridInfo(CMAQ_file)
```
with the function you implemented. Then, you are good to run the code.
<br>
**Notice**: For global model (for example, GC). You need a projection to convert latitude and longitude to X and Y since the model needs to calculate distance.
