function [configs_new,errout]=updateConfig(configs_old,verbose)
% [CONFIGS,ERROUT]=UPDATECONFIG(CONFIGS,VERBOSE)
% Updates configs
%

configs_new=configs_old;

if ~exist('verbose','var')
    verbose=0;
end

% 2017/11/06
if ~isfield(configs_new.flags,'savefigs')
    errout=1;
    if verbose>0
        warning('CONFIGS.flags.savefigs does not exist. setting to default 0.');
    end
    configs_new.flags.savefigs=0;
end


end