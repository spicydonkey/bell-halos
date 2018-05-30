function ch_processed=ch_to_waveforms(ch_params)
%
%   this function takes an array of waveform params and returns a (cell array/or structure) of wavefroms
%   this code deals with constants by looping the minimum number of points (32) in order to conserve
%   wavefrom memory and allow for longer sequences
%
%   ch_params: form like [type('sine','const'),freq(hz),phase(rad),duration,amplitude,sample rate,Gauss mod]
%
%
% TODO
%   if the number of repeats exceeds the max then make the sequence larger so it fits

global max_points  points_min repeats_max
    
n_points_tot=0;
seg_index=1; %need a semprate wf index because the constant parts will generate 1 or 2
zero_at_end=1;
ch_processed=[];

for n=1:size(ch_params,2)
    seg_raw=ch_params{n};
    seg_type=seg_raw{1};
    
    if strcmp(seg_type,'sine')
        seg_freq=seg_raw{2};
        seg_phase=seg_raw{3};
        seg_amp=seg_raw{4};
        seg_gmod=seg_raw{5};
        seg_sr=seg_raw{6}; %sample rate
        seg_duration=seg_raw{7};
        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        
        if n_points>max_points
            error('Requested segment length is too long\n')
        end
        
        t=linspace(0, n_points*dt, n_points);   
        envelope=seg_amp*gausswin(n_points,seg_gmod).';
        wf_out=sin(2*pi*seg_freq*t+seg_phase);    
        wf_out=wf_out.*envelope;
        PA=trapz(t, envelope);
        if zero_at_end & n==size(ch_params,2)
            wf_out=[wf_out,0];
        end
        
		fprintf('Pulse_area=%2.3e for duration %3.1f µs\n', PA, seg_duration*1e6)
        
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points;
    elseif strcmp(seg_type,'double_sine')
        seg_freq1=seg_raw{2};
        seg_freq2=seg_raw{3};
        seg_phase1=seg_raw{4};
        seg_phase2=seg_raw{5};
        seg_amp1=seg_raw{6};
        seg_amp2=seg_raw{7};
        seg_gmod1=seg_raw{8};
        seg_gmod2=seg_raw{9};
        seg_sr=seg_raw{10}; %sample rate
        seg_duration=seg_raw{11};
        
        dt=1/seg_sr;
        n_points=round(seg_duration/dt);

        
        if n_points>max_points
            error('Requested segment length is too long\n')
        end
        
        t=linspace(0, n_points*dt, n_points);   
        envelope1=seg_amp1*gausswin(n_points,seg_gmod1).';
        envelope2=seg_amp2*gausswin(n_points,seg_gmod2).';
        
        wf_out1=sin(2*pi*seg_freq1*t+seg_phase1);    
        wf_out1=wf_out1.*envelope1;
        
        wf_out2=sin(2*pi*seg_freq2*t+seg_phase2);    
        wf_out2=wf_out2.*envelope2;
        
        wf_out=wf_out1+wf_out2;
        PA=trapz(t, envelope1)+trapz(t, envelope2);
        if zero_at_end & n==size(ch_params,2)
            wf_out=[wf_out,0];
        end
        
		fprintf('Pulse_area=%2.3e for duration %3.1f µs\n', PA, seg_duration*1e6)
        
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points; 
        

    elseif strcmp(seg_type,'const') 
        seg_amp=seg_raw{2};
        seg_sr=seg_raw{3}; %sample rate
        seg_duration=seg_raw{4};
        

        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        %here we break the contant spacing into as many points_min blocks
        %as possible with one points_min+excess
        seg_repeats=floor(n_points/points_min-1); 
        points_rem=n_points-(seg_repeats-1)*points_min;
  
        if seg_repeats>repeats_max
            error('Requested repeats for constant exceeds %i \n',repeats_max)
        elseif seg_repeats>0
            ch_processed(seg_index).waveform=ones(1,points_min)*seg_amp;
            ch_processed(seg_index).sr=seg_sr;
            ch_processed(seg_index).repeats=seg_repeats;
            ch_processed(seg_index).vpp=0.001; %set range to min
            seg_index=seg_index+1;
            n_points_tot=n_points_tot+points_min;
        end
        
        if points_rem>0
            ch_processed(seg_index).waveform=ones(1,points_rem)*seg_amp;
            ch_processed(seg_index).sr=seg_sr;
            ch_processed(seg_index).repeats=1;
            ch_processed(seg_index).vpp=0.001;
            seg_index=seg_index+1;
            n_points_tot=n_points_tot+points_rem;
        end    
        
        
        if zero_at_end & n==size(ch_params,2) %if this is the last segment add a zero at the end
            ch_processed(seg_index-1).waveform=[ch_processed(seg_index-1).waveform,0]
        end
        
    else
        error('Unknown Waveform type\n');  
    end
       
end

%check that the total length is resonable
if n_points_tot>max_points
    error('Requested total waveform length is too long\n');
end

fprintf('done generating waveforms\n');


end
