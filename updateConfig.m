function [configs_new,errout]=updateConfig(configs_old,verbose)
% [CONFIGS,ERROUT]=UPDATECONFIG(CONFIGS,VERBOSE)
% Updates configs
%

configs_new=configs_old;

if ~exist('verbose','var')
    verbose=0;
end

%%% 2017/11/06
if ~isfield(configs_new.flags,'savefigs')
    errout=1;
    if verbose>0
        warning('CONFIGS.flags.savefigs does not exist. setting to default 0.');
    end
    configs_new.flags.savefigs=0;
end

%%% 2017/11/09
% - major code update: halo_2bec.m
%   - elev_max feature (spherical polar): replaces zcap
%   - thermal fraction (r_th) culling: replaces dR_tail

% handle halo elevation filtering
for ii=1:numel(configs_new.halo)
    if ~isfield(configs_new.halo{ii},'elev_max')
        errout=1;
        if ~isfield(configs_new.halo{ii},'zcap')
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

% handle thermal fraction filtering
for ii=1:numel(configs_new.bec)
    if ~isfield(configs_new.bec{ii},'r_th')
        errout=1;
        if ~isfield(configs_new.bec{ii},'dR_tail')'
            if verbose>0
                warning('Size for thermal fraction is undefined. Setting r_th to default 0.01.');
            end
        else
            if verbose>0
                warning('r_th is undefined. Will evaluate from dR_tail.');
            end
            % evaluate thermal radius given the tail factor
            configs_new.bec{ii}.r_th=(1+configs_new.bec{ii}.dR_tail)*configs_new.bec{ii}.Rmax;
        end
    end
end

end