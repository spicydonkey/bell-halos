# Bell test in scattering halos

## Experiment
1. Entangled source
2. Local operation
3. Detection
4. Correlation

## 1. Entangled source
[] user-friendly and fast monitor of halo mode occupancy


## 2. Local operation: Rotation
### Characterisation
Automated scanning of Raman amplitude (using AOM1 for state rotation) using code in dir_watch.

[x] load parameter log file
[x] organise data directory according to shots with same R-amp
[] analysis codebase
    [x] categorise data by mf
    [x] sort data by raman amplitude (easier to plot, fit, etc.)
    [] smooth halo sph density before further processing incl population, theta,...
    [] bad zones: pad with NaN (at density eval)
    	[x] poles (use inZone.m)
    	[] dark spots
    	[x] --spontaneous halo-- general region selector: inZone.m
    [] scat mode histogram of theta
    	[x] assuming mf=0 rotate like 1
    [] processing the 0's?
    [] verify each halo capture
        [] return distribution (nn) - for further processing
        [] graphical
    [] data for higher raman amp
[] general theta (no modulo pi) by fitting Rabi flopping to population ratio

## 3. Detection

## 4. Correlation
[] halo fit
    [] to unit sphere
    [] fast g2 amplitude algorithm for BB-pairs
    

## MISC
[] plotFlatMap