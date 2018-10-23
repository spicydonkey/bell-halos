%% Atom number mixing from spin,k processed halo - Ramsey interferometry
% DKS
% 20181022

%% load data
% get spin-momentum processed data
path_dir='C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\expS4_ramsey';
fname_data='expS4_20180823_1.mat';
path_data=fullfile(path_dir,fname_data);

load(path_data);        % load data

close all;  % tidy figures

%% CONFIGS
%%% warning some vars from loaded data clash and are overwritten by configs 
% general
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_ren='painters';

[c,cl,cd]=palette(4);
% c_gray=0.6*ones(1,3);
% line_sty={'-','--',':'};
% mark_typ={'o','s','^'};
mark_typ={'o','s','^','d'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$(1,1)$','$(0,0)$','$(1,0)$'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
str_ss={'$m_J=1$','$m_J=0$','$m_J=-1$','loss'};
mark_siz=7;
line_wid=1.5;
ax_lwid=1.2;
fontsize=12;


%% process data
K=k_par;
n_shot_par=cellfun(@(x) size(x,1),K);

N=cell(nparam,1);
for ii=1:nparam
    N{ii}=cellfun(@(x) size(x,1),K{ii});
end
Navg=cellfun(@(x) mean(x,1),N,'UniformOutput',false);
Navg=cat(1,Navg{:});
Nstd=cellfun(@(x) std(x,0,1),N,'UniformOutput',false);
Nstd=cat(1,Nstd{:});
Nse=Nstd./sqrt(n_shot_par)';

% NUMBER LOSS
Ntot=cellfun(@(x) sum(x,2),N,'UniformOutput',false);
Ntot_avg=cellfun(@(x) mean(x), Ntot);
Ntot_std=cellfun(@(x) std(x), Ntot);
Ntot_se=Ntot_std./sqrt(n_shot_par)';
Nloss=max(Ntot_avg)-Ntot_avg;     % decrease in counts from MAX

%%% collate loss term to mJ
Navg(:,4)=Nloss;
Nstd(:,4)=NaN(size(Nstd,1),1);
Nse(:,4)=NaN(size(Nse,1),1);


%%% smooth fit to data
ninterp=1e4;
phi_fit=linspace(min(par_dphi),max(par_dphi),ninterp);
Nfit=NaN(3,ninterp);

% % Model: polynomial
% polyfit_N=cell(3,1);
% for ii=1:4
%     polyfit_N{ii}=polyfit(par_dphi,Navg(:,ii),6);
%     Nfit(ii,:)=polyval(polyfit_N{ii},phi_fit);
% end

% Model 2: amp,offset,decayrate (i.e. freq, phase are fixed for simplicity)
modelfun=cell(1,4);
modelfun{1}='y~c+a*sin(x1-pi/2)';
modelfun{2}='y~c+a*sin(x1+pi/2)';
modelfun{3}=modelfun{2};
modelfun{4}=modelfun{2};
modelcoef={'a','c'};
par0={[300 150], [250 125], [1 1], [50 -25]};
fo = statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);

sinefit_N=cell(1,3);
for ii=1:4
    sinefit_N{ii}=fitnlm(par_dphi,Navg(:,ii),modelfun{ii},par0{ii},'CoefficientNames',modelcoef,'Options',fo);
    Nfit(ii,:)=feval(sinefit_N{ii},phi_fit);
end


%% vis: atom number Ramsey interferometry
hfig=figure('Name','N_ramsey','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;
pleg=NaN(1,4);      % data for legend
for ii=1:4
    tp=myploterr(par_dphi,Navg(:,ii),[],Nse(:,ii),mark_typ{ii},c(ii,:));  
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
    
    % fit
    tpfit=plot(phi_fit,Nfit(ii,:),'LineWidth',line_wid,'Color',0.5*(c(ii,:)+cl(ii,:)));
    uistack(tpfit,'bottom');
end

box on;
ax=gca;
axis tight;
ax.XTick=pi/2*[0:4];
ax.XTickLabel={'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};
lgd=legend(pleg);
lgd.FontSize=fontsize;
lgd.Location='East';
% lgd.Title.String='$m_J$';
xlabel('Relative phase $\phi$');
ylabel('Atom counts');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwid;

%% Misc number loss analysis
N0=Navg(:,2);

hfig=figure('Name','Nloss_ramsey','Units',f_units,'Position',f_pos,'Renderer',f_ren);

plot(par_dphi,Nloss./(N0+Nloss),'ko-','LineWidth',line_wid);

axis tight;
ylim([-0.1,0.3]);
ax=gca;
ax.XTick=pi/2*[0:4];
ax.XTickLabel={'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};
ax.FontSize=fontsize;
ax.LineWidth=ax_lwid;
xlabel('Relative phase $\phi$');
ylabel('$N_{\textrm{loss}}/(N_{0}+N_{\textrm{loss}})$');
grid on;
