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

% 2017/11/09
% - major code update: halo_2bec.m
%   - elev_max feature (spherical polar) replaces zcap
for ii=1:numel(configs_new.halo)
    if ~isfield(configs_new.halo{ii},'elev_max')
        if ~isfield(configs_new.halo{ii},'zcap')
            errout=1;
            if verbose>0
                warning('Neither zcap nor elev_max exists. Setting elev_max to default pi/2: no elevation culling.');
            end
            configs_new.halo{ii}.elev_max=pi/2;
        else
            if verbose>0
                warning('elev_max is undefined. Will translate from zcap.');
            end
            % simple trigonometry
            configs_new.halo{ii}.elev_max=asin(configs_new.halo{ii}.zcap);
        end
    end
end
    
end