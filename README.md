# MCCC

**M**arkov Chain Monte **C**arlo for **C**orrelation **C**oefficient

## Purpose
The purpose of this program is to examine the correlation coefficient of a given linear fitting, with the consideration of data error in both coordinates.

## Languages
Markov Chain Monte Carlo is iterated by [WinBUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/), and manipulated by [Matlab](https://www.mathworks.com/).

## Codes
`matbugs.m` is an Matlab interface to WinBUGS, [details](https://github.com/matbugs/matbugs). 

`Correlation.txt` defines the correlation coefficient model for WinBUGS.

`Regression.txt` defines the linear regression model for WinBUGS.

`EisoDeltaTCorr.m` gives an example of computing correlation coefficient for Matlab.

`EisoDeltaTReg.m` gives an example of linear regression for Matlab.

Other files are the data and results.