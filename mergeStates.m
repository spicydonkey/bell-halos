function Xcomb=mergeStates(X)
% Merges speciated/distinguished atoms
%
% XCOMB = MERGESTATES(X)
%
% X: Nshots x Nstates cell-array of ncounts x 3 vector arrays
%
% NOTE: the function does not handle multiplicity! every atom should belong to ONE unique state
%

nShots=size(X,1);
Xcomb=cell(nShots,1);
for ii=1:nShots
    Xcomb{ii}=vertcat(X{ii,:});
end

end