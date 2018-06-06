Entangled source in uniform and constant B-field
2018.04.10


c:\remote\settings201810Apr143056.xml

## mini-exp 1: gradual trap SO and B-dynamics

c:\remote\settings201810Apr152535.xml (T=0.5 ms)

exponential switch-off can reduce ringing, but not the total fluctuating transient duration.
above setting works well in reducing overshoot compared to abrupt SO.

c:\remote\settings201810Apr153430.xml (T=0.8 ms)

'''2ms transient B-field after trap SO'''


## setting up source pulse 2ms after trap SO
* [x] record orig (current) waveplate configs


## Raman pulse 2ms after trap SO and alignment
* MultiPulse2 @ 22.46081 = 1.97(1) ms ~ 2 ms after trap SO

c:\remote\settings201810Apr171948.xml

c:\remote\settings201810Apr180533.xml

* [x] Done: seems to have >90% transfer
	* even at fairly high number
	* [ ] quick scan against low N - BEC

	
---

% 2. XZ-uniform bias (nuller after trap SO)
UB_src=1.583e6;      % Zeeman splitting [Hz] - nuller bias
T_src=4.45e-6;             % Tpi ~ 8.7 us; Tpi2 ~ 4.4 us
K_src=0.2;              % max=0.5 (AOM saturation)

Gs_src=3;

phi1=0;
phi2=0;

% evaluate field freqs
dF_Raman_src=(UB_src-abs(Ek));         % [Hz]    Raman detuning
f1_src=f0_AOM-dF_Raman_src/2;            % [Hz]    Raman L1
f2_src=f0_AOM+dF_Raman_src/2;            % [Hz]    Raman L2

---

* took 19 shots 
	* needs longer push to completely separate two halos


## Test +2ms config against longer delays
### 3ms
* MultiPulse2 @ 22.46181 (3ms)
	* B field would have changed
	* alignment should be OK

* seems OK --> would be a more stab'd field if mode occ is OK
	
	
## Global beam alignment
* displaced +Z by ~0.3 mm
	* most likely aligned to the new halo - kick is upwards (source gen'd in 2-3ms after trap SO)


## scattered number
### 3ms
c:\remote\settings201811Apr100114.xml


### 2ms
c:\remote\settings201811Apr101303.xml



## Sequence program for global rotation
* SRC @ ~ms delay after trap SO
* global rotation ~ 1ms after SRC
	* ASAP but need to form scattered pairs

### SRC @ 2ms after TSO
c:\remote\settings201811Apr104855.xml

5ms RF-SWITCH 0.4ms after SRC (Raman trig)
SHUTTERS OPEN

#### TEST
* 2pi @ +1ms SRC
* pretty good - but global beam seems to oscillate off-center slowly?
c:\remote\settings201811Apr112248.xml (END)


### SRC @ 3ms after TSO
c:\remote\settings201811Apr112856.xml

5ms RF-SWITCH 0.4ms after SRC
etc.

#### TEST
same as above


## Check global beam alignment:
A simple way to do this is to SRC then max delay ~3ms then apply 2PI MIX. 
Asymmetry in BEC transfer is then a sign of mis-alignment at rotation

do with WORST CASE SCENARIO: SRC @ 3ms, MIX @ +3ms.

### RESULT
* the bottom BEC (mf=1) hasn't done 2PI; top BEC does a beautiful 2PI
* beam is too high!


* pointed down by 2 mini-turns (~15deg?): I would guess this corr to ~0.3 mm
	* much more symmetric transfer!
	


* '''READY TO DO GLOBAL ROTATION'''


## experiment sequence design
 c:\remote\settings201811Apr122656.xml

### improvements/modifications
* SRC @ 2.5ms after trap SO
* Dynamics in B-field must be as short as possible after MIX
	* Nuller SO early: 1.5 ms after SRC
	* PUSH Z early (coincide with Nuller SO)

* SRC-MIX delay can be designed precisely as part of Keysight pulse
	* delay between MIX-(nuller SO/PUSH) can be controlled precisely
	* minimise: 0.1 ms --> SRC-MIX delay ~1.40(5) ms will do (> 1ms halo creation time)
	
c:\remote\settings201811Apr131410.xml

* INTERESTINGLY we need ~100 ms of Z-push...


#### TESTS
[x] get an estimate of number scattered
[x] set to ~20
[x] test for src no mix (i.e. theta=0)


## ideal_global_rotation_v1
c:\remote\settings201811Apr132816.xml

* counts in halo ~ 20/mf --> g2 ~ 20
* scan angles: pi*[0, 1/2, 1]

* PUSH drifting again...
	* duration decreased to 60 ms. (from ?? 85 ms)

* number drift down too	
	
* need ~1500 shots each --> 4500 tot --> 2 days

* LASER RELOCKED @ #3921


c:\remote\settings201813Apr080541.xml  (END)



## rotation characterisation on poorish Bragg
* T_MIX  in [0, 30] us. every 1us.

* couldn't get a efficient Bragg process (ideally +1 order)

* will run as long as possible until switching over to remaining *global rotation* experiment
	* could be max ~1200 shots tot -- > 40 shots per theta --> number comparable with "ideal Bragg" in trap!

c:\remote\settings201813Apr174820.xml (END --- VVVLOW modeocc - so not usual)



## ideal_global_rotation_v2
c:\remote\settings201813Apr080541.xml (end of v1)

* evap profile shifted DOWN ~kHz

* scanned angles = [pi/8, pi/4, 3*pi/8]; assuming Tpi = 10us

* need 6000 shots
* LN2 dewar needs to be replaced on 2018.04.14 

* Raman laser keeps encountering backreflection faults
* only ~100 shots OK from ~3000... so just startniga gain...


* Lots of unfortunate scenarios...
	* push drift...
	* Raman laser faults

* c:\remote\settings201816Apr154046.xml (END)


	
2018.04.16
## ideal_global_rotation_v3
c:\remote\settings201816Apr154046.xml (END OF v2)

* '''low mode-occupancy''' run at Theta = pi/2 --> maximally(?) distinguish between class of pair-states
	* pi/2 --> ~5us (need to get a good estimate for 5us)


c:\remote\settings201816Apr155613.xml
* evap reduced for low n

c:\remote\settings201817Apr133223.xml (END)



2018.04.17
## ideal_global_rotation_v4
* Minimal SRC-MIX delay = 0.8ms
* halo should (?) be small enough to be coherent for 0.8 ms ~1000*T_larmor
** spin-spin relaxation: is that a strong effect in our sys?

### GOAL
* Tmix-Tsrc = 0.8 ms ~ Tcollision ~ Tbec_sep
* Tsrc = 3.1 ms after trap SO
* RF switch

c:\remote\settings201817Apr135303.xml
* SEQUENCE set up!

* tweaking Ncounts in halo
c:\remote\settings201817Apr135543.xml


c:\remote\settings201818Apr115110.xml (END)

* Exp went superbly


-----------------------------------------------------------------------------------------------------
2018.04.25
## triplet_global_rotation_v1
* following from preliminary results from experiments "ideal_global_rotation"
	* (SRC 3.1 ms after trap SO; 0.8 ms delay to ROT)  seems like a triplet under global rotation
	* FIRST proper param search at *moderate mode occupancy*

### setup
#### Labview
* load ideal_global_rotation_v4 (0.8 ms)
	* c:\remote\settings201818Apr115110.xml (ORIG to load)

	* c:\remote\settings201825Apr090723.xml (push duration increased)
	
#### waveform
* copy from 'global_rotation_constB\ideal_global_rotation_v4\wfmgen_v0.8_20180413'
* scanning theta from 0 to pi - 9 lin-spaced rot duration


#### laser
* 2.20 A
* +3 GHz from 2^3P_0


## triplet_global_rotation_v2
c:\remote\settings201826Apr080915.xml
* the experiment is generally unstable and scanning for so many (9) params at once, when each takes ~1200 shots at ~low mode occupancy, is too risky.
* this version I scan the 4 intermediate points at medium mode occupancy (~>50 counts in mF=0 halo)
	* [1, 2, 4, 5] * pi/6 pulses

STOP to diagnose experiment
* c:\remote\settings201826Apr115207.xml


* OOPS... we don't want to scan up to PI (just want 0-pi/2)


## triplet_global_rotation_v3
* c:\remote\settings201826Apr115207.xml
* [1, 2, 4, 5] * (pi/2)/6 pulses
	par_T_mix=(0.5*T_pi)*[1,2,4,5]*1/6;

	
* SOLVED Push coil drift --> broken Kepco power supply 
	* swapped out P supply with a Tenma 5A x2 --> 20V, 10A same as Kepco.

* now 60 ms push is enough! (as it should be!)

* new setting with working push power supply: c:\remote\settings201826Apr170315.xml

* good amount of data collected (having higher mode occupancy helps! g2 is small though)
* EXP STOPPED @ ~2000 shots
* PUSH was stable!
* Raman laser was happy for ~day
* END c:\remote\settings201827Apr062712.xml


* switch over to unsearched theta (for 0.8ms) in [pi/2, pi]



## triplet_global_rotation_v4
2018.04.27

* c:\remote\settings201827Apr062712.xml
* get moderate mode occ

* pi/2 + [1, 2, 3, 4, 5, 6] * (pi/2)/6 pulses
	* i.e. 7/12, 8/12 ... 12/12 * pi pulse
	



## YAXIS - rotation!
* 0.8 ms, etc. otherwise nominal
* let's check fluoro and optimise IT cooling
* END: c:\remote\settings201828Apr125830.xml

* fluoro good 165mV, and IT cooling already looks good (large BEC on phospor) --> improve Bell with evap profile.

* check evap profile for Bell
* evap like NORMAL (+6kHz): c:\remote\settings201828Apr133158.xml

* mode occupancy max (x3) -- hopefully okay...
* c:\remote\settings201828Apr135329.xml

* ran ok. mode occuapncy not as high as I woudl have liked though. will be bordering on violating the equality.
* discharge source still a problem.

*  (END) c:\remote\settings201829Apr065443.xml



-----------------------------------------------------------------------------------------------------

# the last experiments!
see david's logbook pg 136-

2018-04-29T06:25
## EXP 4: demonstrating time evolution of entangled source
* Given time constraint, reduced the testable evolution times from tdelay=[0.1 0.3 0.8 1.1 1.4, 1.7] ms to:
	* tdelay = [1.1, 1.7] ms, since we already have tdelay=[0.8, 1.4] ms.	(a quadratic rate should be seen)

* other misc configs are:
	* rotation axis: X
	
### wfmgen settings (Matlab - Raman pulse)
See each experiment's folder for configs.
Notice the BEC splitting SRC pi/2 pulse is tailored for each "B-field". tweaking only the pulse duration gives good result.

	
### Labview settings
These define the experiment sequence and timings.


#### 1.1ms 
*	c:\remote\settings201829Apr065804.xml

* need 450 shots 
	* prelim analysis --> OK
* completed
	

#### 1.7ms
*	c:\remote\settings201829Apr070810.xml

* need ~450 shots (3H)
* completed


#### 0.95 ms
*	c:\remote\settings201829Apr112501.xml

* need ~450 shots (3H)
* completed


#### 1.25 ms
* the last one!
*	c:\remote\settings201829Apr112519.xml

* need ~450 shots (3H)
* exp number drifting up :) good thing
* then it drifts down... :p

* c:\remote\settings201829Apr195918.xml (END)
* completed


-----------------------------------------------------------------------------------------------------


## EXP 1 (1.1): demonstrating the violation of Bell inequality
* The big one. the long one. but g2~15 so it's going to be ok - budgeting ~2 days (have additional 1500 shots):
	* Expected run schedule: 2018-04-29T22:00 until 2018-05-01T18:00
	* Expected data: ~7500 shots total

* directory: global_rotation_constB\exp1_attack_inequality
	
	
### CONFIG
* 0.8 ms src-rot delay
* Labview:
	* c:\remote\settings201829Apr065443.xml (latest 0.8ms sequence from Y-axis rot)
* wfmgen:
	* ORIG get from triplet_global_rotation/y_v1
	* set to pi/2 (X) pulse: T_mix = 5us; dphi_2tonedelay=0

* [x] check sequence: good
	
* mode occupancy
	* tune evap profile offset:	
		* NOTE: within 5 kHz to trap bottom
	* c:\remote\settings201829Apr200706.xml
	
* evap bottom to be tuned once 10s of shots available count in [0.416,0.425]
	* Num in above region should be ~15
	
* better measure: num in [0.42,0.425] = 5
	

* 2018-05-01T14:00 - sudden reduction in number. raised evap bottom by ~7 kHz... to get same number in halo

* optimising for fluoro - number is too low.
	* ~ #6300
	* okay tweaked up
	* taking data again freerun (2)
* it's quite unstable... tried tweaking fluoro up again and also trying pressure ~16
	* HV discharge source may drop out sooonish
	* start: 17:50
* number kept drifiting down so don't have as much shots as planned. ~7000 good shots (which is not bad)

	
-----------------------------------------------------------------------------------------------------


## EXP 5: raw E - theta-dependency for low mode-occupancy source
* Data saved under "exp5_low_n_theta"


### Goal
* same mode occupancy as the original points (g2 ~ 25)
	* requires 5 counts in [0.42,0.425]
	* number needs to be stable
	* ~1000 shots per param
	* use data from prev exps at similar mode occ: theta = [pi/4, pi/2]
	* scanned theta = pi/16*[1,2,3,5,6,7]
	* + pair source: theta = 0
		

### Configs
* [x] Labview
	* can use sequence setting for Tdelay = 0.8 ms
		* c:\remote\settings201829Apr200706.xml (the one that worked)
		
		* <2018-05-02T18:00>: c:\remote\settings201802May175851.xml
			* works
* [x] Matlab
	* [x] source
	* [x] scan	


### EXP 5.1: Source
* 1_source
### EXP 5.2: theta scan
* 2_theta_scan

* actually, scan everything!
	* source + scan --> 

* Paused Freerun for M2 installation: c:\remote\settings201803May093636.xml
* <2018-05-03T12:25>: ready to continue
	* evap profile bottom reduced by 1kHz - watch for scattered number

* <2018-05-04T06:50>: experiment remotely paused (only stopped Free run on Labview) since atom number down too low
* <2018-05-04T09:20>: LN2 had depleted last night and shots 4730..4836 have too low atom number. switched to filled larger dewar which should last 3 full-days.
	* ignore ~10's of shots following which is a test for counts
	* possibly due to sudden cooling of HV source, same exp config gave high number when loaded --> evap profile offset by -1.0 kHz - will monitor for ~0.5 hour, as it most likely will drift back down.
* <2018-05-04T09:50>: source pressure: 19.3
* <2018-05-04T10:15>: evap offset -0.7 kHz
	* won't be updating minor evap adjustment logs anymore

* <2018-05-05T12:50>: EXP COMPLETED: c:\remote\settings201805May125657.xml


-----------------------------------------------------------------------------------------------------

* updated IT cooling:
	ORIG weak trap: c:\remote\settings201805May125657.xml
		IT seq seemed to be quite off
	NEW weak trap: c:\remote\settings201805May141048.xml
	
gobal rotation sequence (0.8ms): c:\remote\settings201805May142516.xml
	* high number (above)
	* lowish -- c:\remote\settings201805May143609.xml

	
-----------------------------------------------------------------------------------------------------
	
## cont: EXP 1: continue exps for violation
* <2018-05-05T13:00>: tweak flux + fluoro and optimise IT cooling


* using the exact wfmconfig in 'global_rotation_constB\exp1_attack_inequality\wfmgen_v0.9_20180428' as is.

* aim for 5 counts in [0.42,0.425]

* damn. seems like discharge src dropped out.
* <2018-05-06T10:30>: no the seed laser had unlocked.

* <2018-05-07T08:50>: found that the intensity controller was misbehaving and the power reading at ~3.2V
	* not sure what happened to it, but can't trust all the data taken in this run:
	* data in directory: global_rotation_constB\exp1_attack_inequality\5
		* try to find out when laser power moved (hopefully "jumped?") by looking at sudden change in BEC fraction
	* had to tweak the controller for a while to get it working like OK again... but need to watch it now - may malfunction any time.
	
* <2018-05-07T09:15>: running extra day of violation
	* increasing mode occupancy --> realised the runs 2-5 mode occupancy was too low due to a bright background moving into the ROI

	
-----------------------------------------------------------------------------------------------------


## cont: EXP 5: raw E - theta-dependency for low mode-occupancy source
* missing thetas:
	* theta = pi/16*[9,10,11,12,13,14,15,16] (pi/2..pi)

	
### v2
* par_T_mix=T_pi/16*[10,12,14,16]
	* need to figure out mode occupancy given new background	

* <2018-05-08T08:45>: free-run stopped - lost counts but discharge src ok. was the seed laser.
	* END: c:\remote\settings201808May085951.xml
	* Raman laser intensity lock is OK!

	
### v3
* par_T_mix=T_pi/16*[9,11,13,15]
* need ~7-8 counts in [0.42,0.425]

* Done @ END: c:\remote\settings201809May093135.xml
* saved: C:\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\global_rotation_constB\exp5_low_n_theta\2_theta_scan_v3
* source Press: 14.2
* mode occupancy higher than desired, but that's ok
* Raman laser intensity + freq: OK

* up next: moving onto theta scan + more violation



-----------------------------------------------------------------------------------------------------


## cont: EXP 1: continue exps for violation
* <2018-05-09T09:40>: violation v7
	* same common config: Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\global_rotation_constB\exp1_attack_inequality\wfmgen_v0.9_20180428
	* setting from end of prev EXP5.3 with reduced mocc
		* need ~8 counts in [0.42,0.425]

* EXP stopped for diagnostic
	* END c:\remote\settings201809May171437.xml
	
---

* <2018-05-09T17:40>: violation v8
	* same wfm config
	* number MUCH more stable after tweaking MOT2 alignment on fluoro
	
-----------------------------------------------------------------------------------------------------


## EXP 3: Theta characterisation
* a little tricky to do Bragg, so do best with the "downward source" Raman in-trap
* very-large halo then at ~3.9 ms 
* run from end of EXP1 to wednesday morning when project finishes.
* scan theta~[0,2*pi] ; i.e. T ~ [0,20] us; ~15 points? - can try squeezing points near the "highest gradient" region like pi/2, but then misses small thetas.
* or can scan theta~[0,pi] at 15 points

<2018-05-11T15:00>
* c:\remote\settings201811May160208.xml
	* Raman timing
		* Bragg source OK: SRC ttl; shutter; 
		* AO switching OK: RF switch
		* Bragg-Raman delay OK
	
	
par_T_mix=linspace(0,22e-6,23);
par_value=par_T_mix';




2018-05-12 Bryce
c:\remote\settings201812May143148.xml



-----------------------------------------------------------------------------------------------------

IF EXP TIME/STABILITY permitting...

## EXP 6: High number halo (mf=1, Bragg) mode occupancy scan
* goal is to get as high n as possible and hopefully see glimpse of stimulation effects in g2
* global_rotation_constB\exp6_bragg_halo_high_number
* NO rotation pulse
	* + Ram3 shutter closed
* SRC is Bragg pulse in trap

c:\remote\settings201813May061157.xml

(start) c:\remote\settings201813May064244.xml
- evap raised for high number but not too much thermal frac

* Raman laser error ~19:00 so experiment is officially completed.
END: c:\remote\settings201813May194632.xml


-----------------------------------------------------------------------------------------------------

>>> ALL EXPERIMENTS DONE <<<

-----------------------------------------------------------------------------------------------------

## EXP X: Extra experiment: EPR steering

* low mode occupancy; pi/2 Y-rotation (Jx);
* mode occupancy must be very close to EXP 5 (X-rotation scan)
	* g2 ~ 30

	
### settings

#### LV
INIT: from end of X-rotation: c:\remote\settings201809May171437.xml
* tweaked IT cooling: c:\remote\settings201809May171437.xml

* tweaked number
c:\remote\settings201819May193451.xml



### WFM
INIT: from low-n rot <Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\global_rotation_constB\exp5_low_n_theta\2_theta_scan_v3\wfmgen_v0.9_20180428>

* dphi_2tonedelay=pi/2;     % rotation around +Y-axis
* T_mix=5e-6;     % ~pi/2 rotation     


### versions
#### 1
2018-05-19
* too low number in BEC + thermal

c:\remote\settings201820May144038.xml



#### 2
2018-05-20
* optimising N again - flux --> fluoro --> IT cooling

c:\remote\settings201820May160852.xml

c:\remote\settings201821May165133.xml (end)


#### 3
2018-05-21

* check fluoro - exp unstable
	* it's good (200mV). most likely evap profile in shunted trap is too harsh
	
* redesigned shunted evap profile.
c:\remote\settings201821May171254.xml
c:\remote\settings201821May172355.xml

* try N~8 in t in [0.42,0.425]

* done

-----------------------------------------------------------------------------------------------------


## Supplementary EXPERIMENTS
2018-06-02

* [x] IT cooling update
	* end of last: c:\remote\settings201821May172355.xml (EPR)
	* new (it's the same...!):
		* A: 7.12	6.745	5.1
		* F: 8.45	7.848	7.73
	* don't need to update IT cooling from end of last	
	* c:\remote\settings201821May172355.xml (EPR)
	
	
-----------------------------------------------------------------------------------------------------
### S1. Characterising rotation - longer experiment: (expS1_theta_characterisation)
* end of last:	c:\remote\settings201811May160208.xml
* IT updated:	c:\remote\settings201802Jun124846.xml

* use same wfmgen config as in: DATADIR\exp3_theta_characterisation
* low Nsc...


* c:\remote\settings201802Jun173120.xml

* inevitably, cannot reproduce the experimental config exactly

* DONE
	* cleanup
		* deleted data when Raman laser died ([498,551], [5759,5905])
* END: c:\remote\settings201804Jun115154.xml


-----------------------------------------------------------------------------------------------------
### S2. Trap frequency
* Probably hard to use Raman
* try simple RF PAL

* Trap has not been stabilised! shot-to-shot variation in BEC oscillation too.
* USE same DLD-trig point as original so as not to confuse timing


### A. transient trap from Nuller osc trajectory (A_cont_osc_traj)
* this is the early t behaviour that's relevant for BEC
* t in [0.26,0.41] (from trig) (150ms)

* PAL at 1 kHZ - probing continuous trajectory. only approx freq will be necessary at this point.

* c:\remote\settings201802Jun142424.xml


* PAL config:
	* Wfm: Sqr wave; 2.50 MHz; 1.3 Vpp
	* Burst: 5 cycles
	* pulse period: 1kHz
	
* end: c:\remote\settings201802Jun152457.xml


### B. in-trap osc trajectory (longer monitor) (B_long_osc_traj)
* t in [0.25,1.25]
* trap seems to be evolving and BEC oscillation gets larger

* c:\remote\settings201802Jun153228.xml

* PAL config:
	* Wfm: Sqr wave; 2.50 MHz; 1.3 Vpp
	* Burst: 5 cycles
	* pulse period: 200 Hz

	
### C. No nuller ramp in-trap osc (C_stat_osc_traj) 
* trap osc should not grow like B-dynamics
* trap freq should not be affected strongly by a uniform bias field

* config
	* nuller Vset all 0 + INT
	* c:\remote\settings201802Jun164829.xml
	
* PAL config:
	* Wfm: Sqr wave; 3.50 MHz; 0.8 Vpp
	* Burst: 20 cycles
	* pulse period: 200 Hz


-----------------------------------------------------------------------------------------------------
## S3. Optical pulse
* MAINDIR\expS3_optical_pulse
* measure Bragg/Raman pulse power/intensity

* Bell: (raman L1,L2)
	* wfmgen from expX_epr

* theta charact: (bragg L1,L2)
	* wfmgen from 


	
-----------------------------------------------------------------------------------------------------
	
* LOG DELETED FROM HERE... 2018-06-06
* LOG ammended for S4 as much as possible

-----------------------------------------------------------------------------------------------------

## S4. Ramsey interferometry
* EXP completed

* NOTE: log I wrote yesterday was deleted :'( NEVER let DISK get full
	* and deleted logs... :O retrieved a 2 day old log backed up on Windows

* PRO: At least I still have wfm and LV!
	* wfmgen in `MAINDIR/expS4_ramsey`
	* c:\remote\settings201806Jun091303.xml

	
### Experimental settings (Recall from logbook)
#### sequence
1. Bragg pi/2 pulse in-trap
2. T_delay_mix for stab B (3.88 ms) 
3. 1st pi/2 pulse: dphi_1=0
4. Ramsey TOF: 10.18 us delay (~15*T_L; T_L ~ 0.7 us)
5. 2nd pi/2 pulse: dphi_2 is varied

* dphi_2 was scanned in [0,2*pi] 
	* lin-spaced 13 points
* Ramsey TOF was chosen to ±0.02 us s.t. when dphi_2=d_phi_1=0, halo completely (> 95%) flop to mf=0
	* NOTE: this is a new way to generate a close to pure mf=0 state
	

