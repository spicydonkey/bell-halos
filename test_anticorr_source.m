% Testing Bell correlation with simulated anti-correlated source
% 28/09/2017
% DKS

clear all; clc; close all;


%% Configurations
% TODO - currently needed since the analysis currently requires user to
% step through some parts of "run_zonecount.m" after running this script
configs.zone.nazim=100;
configs.zone.nelev=50;

configs.zone.binmethod=1;
configs.zone.binwidth=0.05;

configs.halo{1}.string='mF=0';
configs.halo{2}.string='mF=1';


%% Prepare source
% Generate ideal source
run('gen_ideal_corr_source');
% halo_k0 simulated

