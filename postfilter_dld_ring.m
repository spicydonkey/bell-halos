function [txy_filt,bool_ring]=postfilter_dld_ring(txy,deadtime,deadxy)
% [txy_filt,bool_ring]=postfilter_dld_ring(txy,deadtime,deadxy)
%
% txy: Nx1 cell array of T,X,Y counts
% deadtime: dld dead-time in seconds for post-filtering (Default: 100ns)
% deadxy: spatial proximity [m] in which to rejects counts if two counts
% arrives within given deadtime. (Default: Inf)
%
% txy_filt: ringing filtered TXY cell-array
% bool_ring: cell array of boolean; marks the duplicates from detector ringing
%

% check inputs
if ~exist('deadtime','var')
    deadtime=100e-9;        % Default: 100 ns dead-time (Sean's thesis)
end
if ~exist('deadxy','var')
    deadxy=Inf;          % Default: reject all counts if dt in deadtime
end

% find duplicates from detector ringing
bool_ring=cellfun(@(x) findalias(x,deadtime,deadxy),txy,'UniformOutput',false);

% take care of shots with zero counts - can't index empty matrix
bool_empty=cellfun(@(x) isempty(x),txy);    % get empty shots
txy(bool_empty)={zeros(1,3)};   % fill empty shots with a NaN vector
txy_filt=cellfun(@(x,y) x(~y,:),txy,bool_ring,'UniformOutput',false);   % filter ringing
txy_filt(bool_empty)={zeros(0,3)};    % empty shots

end