% Raman beam intensity vs AOM waveform amplitude
%
% 2017/11/08
% DKS
%

% aa: K_R_mix, amplitude set on WaveformGenMain.m
% vv: max voltage on PD
% vv_sd: stdev of vv

%% data log
% 2.5 V on noise eater
aa_raw{1}=[0 0.1 0.2 0.3 0.4 0.5 0.55 0.05 0.15 0.25 0.35 0.45 0.525];
vv_raw{1}=[0 0.715 2.33 3.58 4.26 4.59 4.67 0.176 1.500 3.06 3.99 4.48 4.65];
vv_raw_sd{1}=1e-3*[0 4 10 20 20 30 20 1 5 10 20 20 30];
datastr{1}='2.50';

% 3.0 V on noise eater
aa_raw{2}=[0 0.1 0.2 0.3 0.4 0.5 0.55 0.05 0.15 0.25 0.35 0.45 0.525];
vv_raw{2}=[0 0.871 2.82 4.30 5.12 5.50 5.60 0.212 1.79 3.68 4.79 5.40 5.59];
vv_raw_sd{2}=1e-3*[0 3 10 30 30 30 30 1 20 20 30 30 20];
datastr{2}='3.00';

%% preproc
aa=cell(size(aa_raw));
vv=cell(size(vv_raw));
vv_std=cell(size(vv_raw_sd));

% sort in increasing amplitude
for ii=1:numel(aa_raw)
    [aa{ii},thisI]=sort(aa_raw{ii});
    vv{ii}=vv_raw{ii}(thisI);
    vv_sd{ii}=vv_raw_sd{ii}(thisI);
end


%% visualise data
pmarker='o^';
hfig=figure();
for ii=1:numel(aa)
    hold on;
    plot(aa{ii},vv{ii},pmarker(ii),...
        'DisplayName',datastr{ii});
    % TODO - ploterr for displaying error-bars, etc.
    
end
box on;
pleg=legend('Location','northwest');
pleg.Title.String='Laser intensity [V]';
xlabel('Raman amplitude [V]');
ylabel('Intensity [arb.]');