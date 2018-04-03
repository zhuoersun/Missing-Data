# Missing-Data
This is a list of R codes necessary for simulations and the real data analysis for the paper of "Semiparametric Estimation in Regression with Missing Covariates Using Single-Index Models". 

Simulation
----------------

* Please put all .R files under the same directory;

* *AWEE_MAR.R*, *CC_MAR.R*, *full_MAR.R*, *ip_aug_MAR.R* and *ip_noaug_MAR.R* are files of different functions 
   for different methods, including generating data, point estimation, SE & coverages calculation;

* Please run "Table_x_n=..." to get the corresponding results in the tables. The functions in the above
   files will be sourced;

* Detailed explanations are given as comments in the .R files.


YSS Data Analysis
----------------

* Please run *yss.R* to get the results in Table 7 of the paper;

* The user guide and the codebook of the survey study are provided;

* To get the original data set, please send a request through https://uwaterloo.ca/canadian-student-tobacco-alcohol-drugs-survey.
