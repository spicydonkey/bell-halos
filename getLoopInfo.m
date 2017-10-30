function [mf,ampRaman,nver] = getLoopInfo(pathdir)
% [MF, AMPRAMAN, NVER] = GETLOOPINFO(PATHDIR)
%
% Gets details about local operation data from name of database directory 
%

% initialise - this is the error state
mf=NaN;
ampRaman=NaN;
nver=NaN;

params=strsplit(pathdir,'_');
if length(params)~=3
    return;
end

mf=str2double(params{1});
ampRaman=str2double(params{2});
nver=str2double(params{3}(2:end));

end