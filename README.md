# Bell test in scattering halos

## Experiment
1. Entangled source
2. Local operation
3. Detection
4. Correlation

## 1. Entangled source
- [] user-friendly and fast monitor of halo mode occupancy


## 2. Local operation: Rotation
### Characterisation
Automated scanning of Raman amplitude (using AOM1 for state rotation) using code in dir_watch.

- [x] load parameter log file
- [x] organise data directory according to shots with same R-amp
- [] analysis codebase
    - [x] categorise data by mf
    - [x] sort data by raman amplitude (easier to plot, fit, etc.)
    - [] smooth halo sph density before further processing incl population, theta,...
    - [] bad zones: pad with NaN (at density eval)
    	- [x] poles (use inZone.m)
    	- [x] ~~dark spots~~: `holes` in mF=1 halos came from bad halo capture at an initial radial culling stage
    	- [x] --spontaneous halo-- general region selector: inZone.m
        - [x] implement it!
    - [x] scat mode histogram of theta
    	- [x] assuming mf=0 rotate like 1
    - [] processing the 0's?
    - [] verify each halo capture
        - [x] return distribution (nn) - for further processing
        - [x] graphical
        - [] 2D gaussian smoothing
        - [] print out
    - [x] data for higher raman amp
- [] general theta (no modulo pi) by fitting Rabi flopping to population ratio
    - [] Rabi oscillation model: quite tricky - there are difficulties from AOM amp --> diff efficiency relation
    - [] interpolation: a smoothing fit
    - [] how does this affect dTheta histogram?
- [] LOOP pulse after varying delay from SOURCE
    - [] evolution of Pop/Theta histogram
- [] better halo capture
    - [x] liberal inclusion of points for input to halo capture
    - [x] u-sph fitted halo needs filtering
        - [x] TEST IT!
    - [x] verbose output
	- [x] output a figure at each stage
		- [x] bec + thermal
		- [x] first stage halo capture
		- [x] ellipsoid fit
		- [x] cleaned halo

### Data
#### 2017/11/14
- [ ] 3 ms
- [ ] 3.5 ms
- [ ] 4 ms
- [ ] 3 ms - hi power


## 3. Detection

## 4. Correlation
- [] halo fit
    - [] to unit sphere
    - [] fast g2 amplitude algorithm for BB-pairs
- [] analyse bell data 

## MISC
- [x] plotFlatMap
    - see plotFlatMapWrappedRad
    - [x] rad
    - [x] wrapping
- generalise processes by letting functions take array (shot) inputs
    - complete experiment is equivalent to a large shot
    - [] capture_bec
    - [] haloZoneCount
    - [] halo_zone_density: the lat-lon zone-counting function

