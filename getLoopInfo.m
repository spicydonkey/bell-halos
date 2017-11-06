function [mf,ampRaman,nver,tdelay] = getLoopInfo(pathdir)
% [MF, AMPRAMAN, NVER] = GETLOOPINFO(PATHDIR)
%
% Gets details about local operation data from name of database directory 
%

% TODO
% [] output should be varargout - polymorphic?

% initialise - this is the error state
mf=NaN;
ampRaman=NaN;
nver=NaN;

params=strsplit(pathdir,'_');
if length(params)<3
    % file isn't a valid Loop directory. returns NaN
    return;
elseif length(params)==3
    % for MF/RAMANAMP/VERS type
    mf=str2double(params{1});
    ampRaman=str2double(params{2});
    nver=str2double(params{3}(2:end));
elseif length(params)==4
    % for MF/RAMANAMP/DELAY/VERS type
    mf=str2double(params{1});
    ampRaman=str2double(params{2});
    nver=str2double(params{3}(2:end));
    tdelay=str2double(params{4});
end

end