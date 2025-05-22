# Repository for fitting the CSOM to audio data collected from Tentsmuir Forest, Fife, Scotland UK and processed with BirdNET

Data Folder -> Raw:
-Includes folders for each site
-Within site folder, each contains the CSV ouput from birdnet where 1 file is equivalent to 1hr recording processed
-x.rds - list of format to fit CSOM which contains scores

format dat.R
- Function to format data into that required by the CSOM as not all 3-sec clips will be analysed (to generate x.rds)

Fit null CSOM.R 
- Code modified from the accompanying material of Rhinehart et al 2022 to fit the CSOM to the data - requires the x list to be run using format dat.R

