Genome-wide analysis reveals the regulatory role of position-specific DNA G-quadruplex structure on gene expression
---------------------------------------------------------------------------------------
This resource provides the R code and processed data to reproduce key results described in Dan H, et al. Genome-wide analysis reveals the regulatory role of position-specific DNA G-

quadruplex structure on gene expression

### Getting started
**1.** Clone Github repository. 
**1.** Clone Github repository and get the data from google drive. 
```
https://github.com/Hd0909/G4_regulation.git
```
download the rdata from the  following link
https://drive.google.com/file/d/1AIHGz3qQJlqTqoGSNNqeFJbF7VIdyKVa/view?usp=sharing


All the data required to generate the figures are in the ./data folder 
The bash script/Rscript  used to download or preprocess the raw data are also in the ./data folder


**2.** Run the following script to get the figures and supplementary figures
```
cd ./Rscript

## Figure 1 and supplementary Figure1
Rscript all_figure1.R
# Figure 2 and supplementary Figure2
Rscript all_figure2.R
# Figure 3,4 and supplementary Figure3-5
Rscript all_figure34.R
# Figure 5,6 and supplementary Figure6-11
Rscript all_figure56.R
# supplementary Figure12-13
Rscript all_figure12_13.R
```
### Contact
danhuang2018dana@gmail.com
