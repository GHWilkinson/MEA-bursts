# MEA-bursts
An R script using Spike List csv files to analyse Axion MEA data.

This script generates data about Spike activity and Network bursting for multi-electrode data from Axion Maestro Pro recordings of cyto-view MEA plates.
The script will prompt you to upload a Spike List csv file that can be generated during recordings from Axis Navigator and a Plate Map csv file (examples provided).

The following R packages are required:
tidyverse
dplyr
pracma
ggforce
purrr

The script will generate 2 csv files containing spike metrics and synchronised burst metrics respectively. Barcode and recording name will automatically be added to the naming.

Spike metrics provide measures on the general excitability of cultures such as spike rate.

Network activity is analysed using array–wide spike detection rate (ASDR) with a bin width of 200 ms. Synchronized bursts are detected from binned data by a 3-step process: The start of a synchronized burst is detected by a spike count of at least 40% of the maximum ASDR. If the subsequent bin(s) also contained a spike count above the threshold it is included in the synchronized burst. The end of a synchronized burst is determined by a period of 400 ms or more without a spike count about the threshold.

The script will also generate tiff files of ASDR and raster plots for each well that contain barcode, recording name, well and treatment/cell line information in the name.

The script was written by G Wilkinson and J Haddon.
