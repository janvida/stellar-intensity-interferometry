This is a README file on how to run rndmsiistat.c


This code takes a set of parameters defined by a pilot file. 
The default pilot file is default.pilot



To make a new pilot file, make sure that all of the parameters in default.pilot are specified. The order doesnt matter but the * on each line of the parameters matters.



Flag:
-h shows help
-p specifies that a pilot file other than default.pilot will be used



How to run the code with default pilot file:
> ./rndmsiistat 

How to run the code with specified pilot file:
> ./rndmsiistat -p pilotfilename.pilot

How to run code and save results:
> ./rndmsiistat > filename.txt




The standard output of the code appears as follows:


0 1.071980 0.045280 1.072995 0.045250
SNR Sim 3.062610 +/- 0.309370 HBT 3.162278
10 0.919922 0.047095 0.918491 0.047066
SNR Sim 2.944581 +/- 0.297448 HBT 3.162278
20 0.543098 0.042855 0.545318 0.042884
SNR Sim 3.235897 +/- 0.326875 HBT 3.162278
30 0.160235 0.043545 0.163912 0.043511



The first number is the parameter that is changing (probably the baseline separating the two telescopes)
The second column is |g|^2 for the photon correlation
The third column is the standard deviation of the mean of |g|^2 for the photon correlation.
The fourth column is |g|^2 for the signal correlation
The fifth column is the standard deviation of the mean of |g|^2 for the signal correlation.
Then the line below shows the signal to noise ratio (and +/- error) calculated by the simulation and the predicted value from RHB.


If the results are saved to a file, then the file looks like the above data except without the lines showing the SNR.
