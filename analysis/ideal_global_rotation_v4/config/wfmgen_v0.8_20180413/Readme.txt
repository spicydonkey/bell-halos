Version 0.8
* renamed variables
* tidied up code for Bell test


Version 0.7
*double_sin option for driving raman with single beam


Version 0.6
*allow interfacign with multiple waveform generators for up to 4 channel output


Version 0.4
*hope to allow for arbitrarily long pauses without using any memory by stringing waveforms together
*rewriten using using god-fearing PROCEDURAL PROGRAMING istead of this class bulshit

Version 0.3

*improved the keysight waveform code in particular the auto update features in dir_watch
**multi parameter eg amplitude and duration 
**log is stored with data
**data file name in log, date is of creation for data file
**parms in log
**email alert for source dropout
**wont update params if source has dropped out



To be improved
*change param file storage location
*do the params even need to be stored in a file as there is no mechanism to resume a param sweep if stopped half way ATM.
**should create a resume command
*gaussian pulse should have some kind of normalization so that the rabi phase is not changed (except for fourrier broadning)
*integration for pulse area should follow same form as rabi phase
*documention for user
*combine backend functions into single matlab file