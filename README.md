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
    	[x] ~~dark spots~~: `holes` in mF=1 halos came from bad halo capture at an initial radial culling stage
    	[x] --spontaneous halo-- general region selector: inZone.m
        [] implement it!
    [x] scat mode histogram of theta
    	[x] assuming mf=0 rotate like 1
    [] processing the 0's?
    [] verify each halo capture
        [x] return distribution (nn) - for further processing
        [x] graphical
        [] 2D gaussian smoothing
        [] print out
    [x] data for higher raman amp
[] general theta (no modulo pi) by fitting Rabi flopping to population ratio
    [] Rabi oscillation model: quite tricky - there are difficulties from AOM amp --> diff efficiency relation
    [] interpolation: a smoothing fit
    [] how does this affect dTheta histogram?
[] LOOP pulse after varying delay from SOURCE
    [] evolution of Pop/Theta histogram
[] better halo capture
    [x] liberal inclusion of points for input to halo capture
    [x] u-sph fitted halo needs filtering
        [] TEST IT!
    [] verbose output
	[] output a figure at each stage

% TODO
% [x] improve algorithm
% [x] 1. capture BEC: get BEC positions: jitter, indexible bec counts
% [x] 2. 1st stage halo capture
%   [x] approx to halo center
%   [x] capture rough halo counts and center shot-wise
%   [x] 2.1. liberal filtering of counts: radial, BEC, thermal
        [x] radial - kdR_hicap is currently hardcoded into function
        [x] elev filter
            ~[] need to define field "elev_max" for each configs.halo element~
            - configs.halo{I}.zcap isn't used anymore
            [x] translate zcap --> elev_max if undefined
        [x] thermal fraction around poles/BEC
        [] filter stages need to be handled better

% [x] 2.2. ellipsoid mapping --> unit sphere
        [x] packaged mapping into a function
% [x] 3. clean halo
%   [x] 3.1. radial filter
%   [x] 3.2. pole filter
%

## 3. Detection

## 4. Correlation
[] halo fit
    [] to unit sphere
    [] fast g2 amplitude algorithm for BB-pairs
    

## MISC
[x] plotFlatMap
    see plotFlatMapWrappedRad
    [x] rad
    [x] wrapping
