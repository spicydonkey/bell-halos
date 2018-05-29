function plot_segments(dev_chan_wfm,if_zeros)
%   this plots the wavefroms defined by channels
%
%
% TODO
%   add a time accululator to handle variable sample rates (if the gen
%   supports this)
%   titles ect
%

figure(10);
for n=1:2
    chan_amp=[]; %point list
%     chan_time=[];
    for m=1:size(dev_chan_wfm{n},2)
        sub_wf=dev_chan_wfm{n}.waveform;
        if if_zeros
            chan_amp=[chan_amp,repmat(dev_chan_wfm{n}(m).waveform,1,dev_chan_wfm{n}(m).repeats)];
        else
            chan_amp=[chan_amp,dev_chan_wfm{n}(m).waveform];
        end
    end
    fprintf('total points:\t%d\n',size(chan_amp,2));
    subplot(2,1,n)
    plot(chan_amp)
    
end
clear chan_amp
pause(0.1)
end
