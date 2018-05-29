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