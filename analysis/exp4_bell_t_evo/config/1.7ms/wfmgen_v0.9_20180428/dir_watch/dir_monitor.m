function dir_monitor(Low_Count_Size,Email_add,Repeat_Scan)
global watching_dir
global start_time
global pointer_fname
global pars_fname
global repeats
global progress
low_count_files=0;
high_count_files=0;
emailsent=0;
file_end=0;
%DIR_MONITOR executes action_fun upon a new file in the watching_dir
% is there any need to have action fun as variable?

%TO BE IMPROVED
%move the update code that i have used in the autoconvert and dld front
%pannel to this code with improvments such as
%wait till done writing
%check if first few lines pf file are valid not just size


%read the first param value and send it to the wavefrom gen
if exist(pointer_fname, 'file')
    fid1 = fopen(pointer_fname, 'r');
    pointer=fscanf(fid1, '%i');
    fclose(fid1);
else
    pointer=1;
end
%Read the parameter value according to pointer:
pars=dlmread(pars_fname,',');

%fid2 = fopen(pars_fname, 'r');
%pars=fscanf(fid2, '%f');fclose(fid2);
if isempty(pars)
    error(sprintf('\nERROR %s is empty', pars_fname));
end

fprintf('\nFirst params sent to waveform gen');
instr_update(pars(1,:)); %update the instrument with the first set of parameters

pause on
dir_init_content = dir(watching_dir);
initial_files = {dir_init_content.name};
fprintf('\nMonitoring %s', watching_dir);

while true
    dir_init_content = dir(watching_dir);
    filenames = {dir_init_content.name};
    new_files = setdiff(filenames,initial_files);
    %prevents updating on creation of _txy_forc files
    new_files=new_files(cellfun(@(x) isempty(findstr('_txy_forc',x)),new_files));
    new_files=new_files(cellfun(@(x) ~isempty(findstr('.txt',x)),new_files));
    
    if ~isempty(new_files)
        % deal with the new files
        %Bryce- added an error message if multiple new files found
        initial_files = filenames;
        if numel(new_files)~=1
            fprintf(2,'\nERROR multiple(%d) new files detected: ',numel(new_files));
            for k=1:numel(new_files)
                fprintf('\n %s\t', new_files{k});
            end
        end
        
        fprintf('%s\n', new_files{1});
        pause(6)%this must be set to the max timeout time on the TDC front
        %pannel, 5 seconds is rarely exceeded
        %find the creation time of this file to store in the params LOG
        FileInfo=dir(fullfile(watching_dir,new_files{1}));
        FileSize=FileInfo.bytes;
        FileSize=FileSize/(2^20);
        %I read the file size, if it is too small then i will not update
        if FileSize>Low_Count_Size
            FileDateDay=FileInfo.datenum;%the file mod date is returned in seconds
            file_end=update_action(new_files{1},FileDateDay,Repeat_Scan);
            if file_end==1
                break
            end
            %here i only pass the first entry from the
            %new_files list, this will rarely be a problem
            low_count_files=0;
            
            high_count_files=high_count_files+1;
            if high_count_files>5 && emailsent==1 && ~isempty(Email_add)
                emailsent=0; %hysterisis in the reset to prevent email spam
                fprintf(2,'Source Started Sending Email Alert\n')
%                 SendEmail(Email_add,'Source Started',...
%                     sprintf('The source has started up based on the file size for %d shots of the experiment. The time is %s. The scan is %4.1f %% complete.Projected ETA %s HH:MM finishing at %s \n',...
%                     high_count_files,datestr(now, 'dd mmm HH:MM'),100*(progress+repeats),datestr(ETA_Days,'HH:MM'),datestr(datenum(clock)+ETA_Days,'dd mmm HH:MM')));
            end
        else
            %if it is too low then a red text message is displayed
            %if this happens enough then i send an email to myself
            low_count_files=low_count_files+1;
            fprintf(2,'\n file is too small(%2.3f MiB) will not update. %d low files  \n',FileSize,low_count_files)
            if low_count_files>5 && emailsent==0 && ~isempty(Email_add)
                high_count_files=0;
                emailsent=1; %prevents spaming my inbox
                fprintf(2,'Sending Email Alert\n')
                ETA_Days=(1-progress)*(datenum(clock-start_time)/(progress+repeats));
%                 SendEmail(Email_add,'Source Dropout',...
%                     sprintf('The source has dropped out based on the file size for %d shots of the experiment. The time is %s. The scan is %4.1f %% complete.Projected ETA %s HH:MM finishing at %s \n',...
%                     low_count_files,datestr(now, 'dd mmm HH:MM'),100*(progress+repeats),datestr(ETA_Days,'HH:MM'),datestr(datenum(clock)+ETA_Days,'dd mmm HH:MM')));
            end
        end
        fprintf('Monitoring %s', watching_dir);
    else
        pause(.5) %little wait animation
        fprintf('.')
        pause(.5)
        fprintf('.')
        pause(.5)
        fprintf('.')
        pause(.5)
        fprintf('\b\b\b')
    end
end
fprintf('Scan Finished %s HH:MM',datestr(datenum(clock-start_time),'HH:MM'))
% SendEmail(Email_add,'Scan Finished',sprintf('The scan has finished. The time is %s. It took %s HH:MM to run.',...
%     datestr(now, 'YYYY-mm-dd_HH:MM:SS'),datestr(datenum(clock-start_time),'HH:MM')) );

end

