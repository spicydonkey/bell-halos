function send_segments(dev_chan_wfm,id_dev)
%
%   this program sends the wavefrom segments to the keysight function generator and programs the
%   appropriate repeats as specified
%   ch1(1) is a data feild that has properties waveform(vector of doubles),sr(sample rate) and repeats
%
% DEBUG
%   see:
%       instrreset
%       instrfindall
%


%%START user data
if id_dev==1
    %usb_visa_id='USB0::0x0957::0x5707::MY53800526::0::INSTR';
%     usb_visa_id='USB0::2391::22279::MY53800526::0::INSTR';
    usb_visa_id = 'USB0::0x0957::0x5707::MY53800526::0::INSTR';
    defaultIP = '150.203.178.177';
    defaultPORT = 5025; 
elseif id_dev==2
%     usb_visa_id='USB0::0x0957::0x5707::MY53802248::0::INSTR';
%     usb_visa_id='USB0::2391::22279::MY530802248::0::INSTR';
   usb_visa_id='USB0::2391::22279::MY53800526::0::INSTR';
    defaultIP = '150.203.178.177';
end    
ARB_DIR='Int:\UserFiles\';
         

reset=1;        %resets the instrument, takes ~5sec 
use_usb=1;      %else use LAN
usb_visa_id = 'USB0::0x0957::0x5707::MY53800526::0::INSTR';
%%END user data

%initialize connection to the device
if use_usb
    instr=instrfind('Type', 'visa-usb', 'RsrcName',usb_visa_id, 'Tag', '');    
else
    instr=instrfind('Type', 'tcpip', 'RemoteHost', pars.IP, 'RemotePort', pars.PORT, 'Tag', '');   
end

%dont really see why this case strucute is here,perhaps so that if the
%instrument is not found just try to connect on the default
if isempty(instr) & ~use_usb
    % Create the tcpip object if it does not exist otherwise use the object that was found.
    instr = tcpip(defaultIP,defaultPORT);   
elseif isempty(instr) & use_usb                                                
   % Create the VISA-USB object if it does not exist, otherwise use the object that was found.                
    instr=visa('ni', usb_visa_id);                
else %if not empty take the first reference
    fclose(instr);
    instr= instr(1);
end            


% Connect to instrument object, obj.instr.
instr
fopen(instr)
fprintf(query(instr, '*IDN?'),'\n')
fprintf(instr, '*CLS');
if reset== 1
    fprintf(instr, '*RST');
elseif reset == 0 %Just switch off the outputs
    fprintf(instr, 'OUTP1 OFF');               
    fprintf(instr, 'OUTP2 OFF');
end

check_errors(instr);
fclose(instr);
          
send_sequence(instr,dev_chan_wfm,ARB_DIR)
switch_on(instr)

end


function switch_off(instr)
    fopen(instr);              
    fprintf(instr, 'OUTP1 OFF');
    fprintf(instr, 'OUTP2 OFF');               
    fclose(instr);
end


function check_errors(instr)
    error_string=query(instr,':SYST:ERR?') ;
    if isempty(strfind(error_string, '"No error"')) %&& isempty(strfind( error_string, 'value clipped to lower limit'))
        error('!!!==========INSTRUMENT ERROR==========!!!\n %s---------------------------', error_string);            
    end      
end


function switch_on(instr)
    fopen(instr);              
    fprintf(instr, 'OUTP1 ON');
    fprintf(instr, 'OUTP2 ON');               
    fclose(instr);
end


function delete(instr)
    fprintf('Disconnecting...\n')            
    fclose(instr);           
end


 function upload_arb_binary(instr,ARB_DIR, arb_fname, wf, SRate, Vpp)
    %Sends ASCII data and creates text .ARB file in the instruments
    %memory         
    wf_bin=int16(round(wf*(2^15-1)*2/Vpp));
    max(wf_bin-(2^15-1))
    min(wf_bin+(2^15-1))
    %             min(wf_bin)
    data_size=2*length(wf_bin); %one point takes up 2 int8 feilds
    %definite_length_header defines the number of points in the
    %wavefrom as #[feildwidth][data_size]
    scpi_header=sprintf('DATA:ARB:DAC myArb, %s', definite_length_header(data_size));
    %increase the buffer size to handle this large wavefrom
    instr.OutputBufferSize = 8*length(scpi_header)+data_size;
    %             obj.instr.InputBufferSize = 8*length(scpi_command);
    instr.Timeout=100; %in seconds
    fopen(instr);
    %             get(obj.instr,{'OutputBufferSize','BytesToOutput'})      
    %             fprintf(obj.instr, 'SOUR:VOLT %g mVpp', 0.0); 
    fprintf(instr, 'FORM:BORD SWAP'); %set block mode to use LSB
    check_errors(instr);
    
    fwrite(instr, [scpi_header,  typecast(wf_bin, 'uint8')], 'uint8');      
    fprintf(instr, 'FUNC:ARB:FILT STEP'); %filter types OFF NORM STEP which gives 3db n/a 0.27 0.13 sample rate
    fprintf(instr, 'FUNC:ARB myArb');    %select the wavefrom called myArb       
    %             fprintf(obj.instr, 'FUNC:ARB:PTP %g', max(abs(wf)));   
    fprintf(instr, 'FUNC:ARB:SRATE %g', SRate);
    check_errors(instr);
    fprintf(instr, 'SOUR:VOLT %g Vpp', Vpp); 
    check_errors(instr);
    fprintf(instr, 'SOUR:VOLT:OFFS 0.0');  
    fprintf(instr, 'MMEM:STOR:DATA "%s"', [ARB_DIR, arb_fname]);
    check_errors(instr);
    fprintf(instr, '*WAI');
    %             obj.query('FUNC:ARB:PTP?');
    fprintf(instr, ':DATA:VOL:CLEar');  
    %             obj.query(sprintf('MMEM:CAT? "%s"', obj.ARB_DIR));            
    %             obj.query(sprintf('MMEM:UPLoad? "%s"', arb_fname));            
    check_errors(instr); 
    flushinput(instr);
    fclose(instr);
 end             

 
 function send_sequence(instr,chanels,ARB_DIR)  
    fprintf('Uploading sequences: ');
    vpp_net=[0,0];
    for n=1:2
        name_string=['CH', int2str(n)];
        fprintf('%s:', char(name_string)); 
        barb_files={};
        repeat_count=[];
        
        for k=1:size(chanels{n},2)
            wf_data=chanels{n}(k);
            vpp_net(n)=max([vpp_net(n),wf_data.vpp]);
        end
        %vpp_net(n)
        for k=1:size(chanels{n},2)
            %fprintf('.%s..', mat_files{k});
            wf_data=chanels{n}(k);                    
            repeat_count(k)=wf_data.repeats;
            barb_files{k}=[int2str(k), char(name_string), '.barb'];                   
            upload_arb_binary(instr,ARB_DIR,barb_files{k}, wf_data.waveform, wf_data.sr, vpp_net(n));   
        end 
        fopen(instr);
        seq_name=sprintf('seq_%s',name_string);
        seq_setting=sprintf('"%s"',seq_name);
        for b=1:numel(barb_files)                     
            barbfile=[ARB_DIR, barb_files{b}];
            fprintf(instr, 'MMEM:LOAD:DATA "%s"', barbfile);                                        
            if b==numel(barb_files)
                play_control='onceWaitTrig';
            else
                play_control='repeat';
            end                                         
            arb_setting=sprintf('"%s",%i,%s,highAtStart,10', barbfile, repeat_count(b), char(play_control));
            seq_setting=[seq_setting,',', arb_setting];
        end      
    %                 seq_setting=[seq_setting,',', ['"',obj.ARB_DIR,'short_DC.barb"'],',1,"onceWaitTrig","highAtStart",10']                 
        fprintf(instr, 'SOUR1:BURS:STAT OFF');                 
        fprintf(instr, 'SOUR2:BURS:STAT OFF');
        fprintf(instr, 'SOUR:DATA:SEQ %s',  [definite_length_header(numel(seq_setting)), seq_setting]);                
        fprintf(instr, 'SOUR:FUNC:ARB "%s"', seq_name);                  
        fprintf(instr, 'MMEM:STORE:DATA "%s"', [ARB_DIR, seq_name, '.seq']);
        fprintf(instr, '*WAI');
        fprintf(instr, ':DATA:VOL:CLEar');
    %                 obj.check_errors() 
    
        check_errors(instr);
        fclose(instr);                                            
    end
    fprintf('...done\n')            
    set_run_settings_seq(instr,ARB_DIR, 'seq_CH1.seq', 'seq_CH2.seq', vpp_net(1), vpp_net(2))              
 end  

 
function set_run_settings_seq(instr,ARB_DIR, seq1_fname, seq2_fname, Vpp1, Vpp2)              
    %%SETTING SEQUENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    fopen(instr);             
    fprintf(instr, 'SOUR1:DATA:VOL:CLEar');              
    fprintf(instr, 'SOUR2:DATA:VOL:CLEar');                              

    fprintf(instr, 'MMEM:LOAD:DATA1 "%s"', [ARB_DIR, seq1_fname]);
    fprintf(instr, 'MMEM:LOAD:DATA2 "%s"', [ARB_DIR, seq2_fname]);                       

    fprintf(instr, 'SOUR1:FUNC:ARB "%s"', [ARB_DIR, seq1_fname]);                        
    fprintf(instr, 'SOUR2:FUNC:ARB "%s"', [ARB_DIR, seq2_fname]);          
    fprintf(instr, 'SOUR1:FUNC ARB');
    fprintf(instr, 'SOUR2:FUNC ARB');  
    fprintf(instr, 'SOUR1:FUNC:ARB:FILT STEP'); %OFF NORM STEP   
    fprintf(instr, 'SOUR2:FUNC:ARB:FILT STEP'); %OFF NORM STEP   
    
    fprintf(instr, 'SOUR1:FUNC:ARB:SRAT %g', 1e9);
    fprintf(instr, 'SOUR2:FUNC:ARB:SRAT %g', 1e9);
    
    fprintf(instr, 'FUNC:ARB:SYNC'); 
    fprintf(instr, 'SOUR1:VOLT:OFFset 0.0'); 
    fprintf(instr, 'SOUR2:VOLT:OFFset 0.0'); 
    fprintf(instr, 'SOUR1:VOLT %g Vpp', Vpp1);            
    fprintf(instr, 'SOUR2:VOLT %g Vpp', Vpp2); 
%             fprintf(obj.instr, 'TRIG1:SOUR TIMer');
%             fprintf(obj.instr, 'TRIG2:SOUR TIMer');      
    fprintf(instr, 'TRIG1:SOUR EXT');
    fprintf(instr, 'TRIG2:SOUR EXT');                

%             fprintf(obj.instr, 'OUTP:SYNC:MODE NORM');                

    check_errors(instr);
    fclose(instr);            
end


% UTILITIES:
%definite_length_header defines the number of points in the
%wavefrom as #[feildwidth][data_size]
%eg 523 long waveform would be #3523
%because generator is limited to 64e6 max never going to reach the 9 feild
%limit of this method
function h=definite_length_header(data_size)
    num_to_follow=length(sprintf('%i', data_size)); 
    h=sprintf('#%i%i', [num_to_follow, data_size]);
end
        
        
