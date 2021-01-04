# Regionalization: Spatial-Functional Mixture Models

## Data
The folder **data** provides necessary datasets in the analysis of China or BTH. The file **China_PM25_daily.csv** includes city-level daily PM2.5 concentrations of entire China (338 cities) for two years (Jan 1st, 2015 - Dec 31st, 2016). The file **BTH_PM25_monthly.csv** includes station-level monthly PM2.5 concentrations from 73 stations in the BTH region for four years (Jan 2013 - Dec2016). The file **city_location_China.csv** and **station_location_BTH.csv** provide the topographic information involving longitude and latitude for China and BTH respectively. The file **China_elevation_1km.tif** contains China's elevation data with 1km resolution, which is used for the neighborhood selection in the analysis of China. The file **Neighbor338_neighbor.Rdata** contains the selected neighbors of 338 cities of China, which was resulted by the script **GetNeighbors.R** in the folder **RealData_China** and was saved for the convenience in the analysis of China. 
The folder **China_boundary_map** includes necessary files for the plot of China's boundary. The folder **BTH_boundary_map** includes necessary files for the plot of BTH's boundary, while the file **BTH_map.Rdata** contains the map data for BTH.

## R Code

The folder **code** provides necessary R code to implement the spatial-functional mixture model.  The R script **functions.R** contains some basic functions used in the programs. Three different parts of programs are included.

### Simulation
The folder **Simulation** implements the simulation studies in the manuscript.
The **homo** folder provides code for the homoscedastic case in simulation. The **Simulation_homo.R** is the main function for regionalization and reproduction of *Table 1* and *Table 2*. The associated script **SimulateMRF.R** simulates the Markov random fields.
The associated scripts **homo_iid_multi.R**, **homo_iid_mrf.R**, **homo_spatial_multi.R**, **homo_spatial_mrf.R** perform the EM algorithms of "FMM", "FMM-MRF", "SFMM", "SFMM-MRF" in the manuscript, respectively.
The **hetero** folder provides code for the heteroscedastic case in simulation. The **Simulation_hetero.R** is the main function for regionalization and reproduction of *Table 3*, *Table 4*,  and *Figure 3*. The associated script **SimulateHeteroData.R** simulates the heteroscedastic data.

### RealData_China
The folder **RealData_China** contains R functions for the real data analysis of China. The R script **Regionalization_China.R** is the main function for regionalization and reproduction of *Figure 1*, *Figure 4*, *Figure 5*, *Figure 6*, and *Figure 9* in the manuscript. The associated script **regional_China.R** includes the main EM algorithm for regionalization. The associated script **GetNeighbors.R** performs neighborhood selection by Markov random field combining geographic information.

### RealData_BTH
The folder **RealData_BTH** contains R functions for the real data analysis of BTH region.  The R script **Regionalization_BTH.R** is the main function for regionalization and the reproduction of *Figure 2*, *Figure 7*, *Figure 8*, and *Figure 10*, while the associated script **regional_BTH.R** includes the EM algorithm.