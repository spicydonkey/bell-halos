clear all;

% USER CONFIG
fileToWatch='Y:\d_txy_forc*.txt';
% ramanKminmax=[0.1,0.5];
nShotsAtSwitch=100;
% nSets=10;
tPause=2;
% ramanK=linspace(ramanKminmax(1),ramanKminmax(2),nSets);

nSets=4;
ramanK=[0.347,0.382,0.42,0.44];

% MAIN
oldSet=0;
while true
    % watch directory and report most recent SHOT number
    dfiles=dir(fileToWatch);
    dfile_ids=arrayfun(@(x) str2num(x.name(11:end-4)),dfiles,'UniformOutput',false);
    
    if isempty(dfile_ids)
        nShotNow=1;
    else 
        dfile_ids=vertcat(dfile_ids{:});
        nShotNow=max(dfile_ids)+1;
    end
    
    thisSet=floor((nShotNow-1)/nShotsAtSwitch)+1;
    if thisSet>nSets
        break;
    end
    
    if thisSet~=oldSet
        thisRamanK=ramanK(thisSet);
        fprintf('* switching Raman amplitude: %0.3g\n',thisRamanK);
        
        % update waveform
        run('WaveformGenMain.m');       
        % takes less than 10 sec - pulse comes on ~22 s into sequence
        
        oldSet=thisSet;
    end
        
    pause(tPause);
end

disp('Done! Please stop Experiment and save data');