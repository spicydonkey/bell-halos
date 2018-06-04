% Raman/Bragg laser optical pulses
%
% DK Shin
% 2018-06-04


%% config
path_data='C:\Users\HE BEC\bell\raman-laser\expS3_optical_pulse\raw';


%% main
%%% parse raw data directory
dir_data=dir(path_data);
f_name_all={dir_data.name}';

% get csv data
b_file=cellfun(@(s) is_file(fullfile(path_data,s)),f_name_all);
f_name_csv=f_name_all(b_file);

b_csv=cellfun(@(s) strcmp(s(end-2:end),'csv'),f_name_csv);
f_name_csv=f_name_csv(b_csv);
ndata=numel(f_name_csv);


% get file id
f_id=cellfun(@(s) str2double(char(extractBetween(string(s),'NewFile','.csv'))),f_name_csv);
[f_id_sort,Is]=sort(f_id);
f_name_csv_sort=f_name_csv(Is);


%%% Read raw data
d=cell(ndata,1);

for ii=1:ndata
    d{ii}=read_osci_rigol(fullfile(path_data,f_name_csv_sort{ii}));
end


%% vis
figure();
ax=gca;

for ii=1:ndata
    hold on;
    td=d{ii};
    tt=td(:,1);     % time (s) since Keysight TRIG
    ttrig=td(:,2);  % trigg channel
    tv=td(:,3);     % PD output V
    
%     plot(1e3*tt,ttrig);
    plot(1e3*tt,tv);
    
end

box on;
ylabel('Power (arb. unit)');
xlabel('$t$ (ms)');
title('Raman/Bragg: Optical pulse');

y_lim=ax.YLim;
ylim([0,y_lim(2)]);
