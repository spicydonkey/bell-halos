% (c) Roman Khakimov, 2015
function file_end=update_action(new_files,FileDateDay,Repeat_Scan)
%%
%UPDATE_INSTRUMENT Summary of this function goes here
%   Detailed explanation goes here
    global pointer_fname
    global pars_fname
    global pars_log_fname
    global start_time %for ETA
    global repeats
    global progress
    file_end=0;
    fprintf('%s\n', 'Update action triggered');  
    %Read pointer -- number of parameter to use
    if exist(pointer_fname, 'file')
        fid1 = fopen(pointer_fname, 'r');        
        pointer=fscanf(fid1, '%i');
        fprintf('Param Pointer is %i\n',pointer)
        fclose(fid1);        
    else
        pointer=1;
    end
    %Read the parameter value according to pointer:
    pars=dlmread(pars_fname,',');


    if isempty(pars)
        error(sprintf('\nERROR %s is empty', pars_fname));
    end

    %Logging the used parameter:
    fid3 = fopen(pars_log_fname, 'a');
    fprintf(fid3,'%s, ', char(new_files));
    fprintf(fid3,'%s, ',datestr(FileDateDay,'YYYY-mm-dd_HH:MM:SS'));
    fprintf(fid3,'%20.10e, ', pars(pointer,:));
    fprintf(fid3,'\n');
    fclose(fid3); 
    
     if size(pars,1)<=pointer
        fprintf('End of Param list\n')
        if Repeat_Scan==1
            pointer=0;
            repeats=repeats+1;
        else
        file_end=1;
        end
    end
    
    if ~file_end==1 %not at the file end
        %Update pointer with next number:
        fid4 = fopen(pointer_fname, 'w');        
        fprintf(fid4,'%i', pointer+1);    fclose(fid4); 
        
        %want to check if the last set of params was identical
        if pointer~=0
            if pars(pointer+1,:)==pars(pointer,:)
               fprintf('Same parameters will not update\n')
            else
                %Instrument update and parsing the parameters
                instr_update(pars(pointer+1,:));
           end
        else
        %Instrument update and parsing the parameters
        instr_update(pars(pointer+1,:));
        end
        %display the percentage through the loop with an ETA
        %the ETA will break if the source drops out
        progress=pointer/size(pars,1);
        ETA_Days=(1-progress)*(datenum(clock-start_time)/(progress+repeats));
        fprintf('%3.1f %% complete. ',100*(progress+repeats))
        fprintf('Projected ETA %s HH:MM finishing at %s \n',...
            datestr(ETA_Days,'HH:MM'),datestr(datenum(clock)+ETA_Days,'dd mmm HH:MM'))
    
    end
    if file_end==1 && Repeat_Scan==1
        
    end
end



