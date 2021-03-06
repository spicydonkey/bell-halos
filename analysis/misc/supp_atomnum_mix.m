%% Atom number mixing from spin,k processed halo
% DKS
% 20181010

%% load data
% get spin-momentum processed data
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp3_yrot_char';
fname_data='exp3_s1_20180823_1.mat';
path_data=fullfile(path_dir,fname_data);

load(path_data);        % load data

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
Nloss=Ntot_avg(1)-Ntot_avg;     % decrease in counts from NO PULSE

%%% collate loss term to mJ
Navg(:,4)=Nloss;
Nstd(:,4)=NaN(size(Nstd,1),1);
Nse(:,4)=NaN(size(Nse,1),1);


%%% smooth fit to data
ninterp=1e4;
tau_fit=linspace(min(par_T),max(par_T),ninterp);
Nfit=NaN(3,ninterp);

% % Model: polynomial
% polyfit_N=cell(3,1);
% for ii=1:4
%     polyfit_N{ii}=polyfit(par_T,Navg(:,ii),6);
%     Nfit(ii,:)=polyval(polyfit_N{ii},tau_fit);
% end

% Model 2: amp,offset,decayrate (i.e. freq, phase are fixed for simplicity)
modelfun=cell(1,4);
modelfun{1}='y~c+a*exp(-x1/xc)*sin(2*pi*0.05*x1+pi/2)';
modelfun{2}='y~c+a*exp(-x1/xc)*sin(2*pi*0.05*x1-pi/2)';
modelfun{3}=modelfun{2};
modelfun{4}=modelfun{2};
modelcoef={'a','c','xc'};
par0={[120 130 50], [120 130 50], [1 1 50], [25 25 -40]};
fo = statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);

sinefit_N=cell(1,3);
for ii=1:4
    sinefit_N{ii}=fitnlm(1e6*par_T,Navg(:,ii),modelfun{ii},par0{ii},'CoefficientNames',modelcoef,'Options',fo);
    Nfit(ii,:)=feval(sinefit_N{ii},1e6*tau_fit);
end


%% vis
% grahics configs
% build color palette (blue,red,black,green)
[c,cl,cd]=palette(4);
% third color to black
c(4,:)=c(3,:);
cl(4,:)=cl(3,:);
c_gray=0.6*ones(1,3);
c(3,:)=[0,0,0]; 
cl(3,:)=c_gray;

mark_typ={'o','s','^','d'};
mark_siz=8;
line_wid=2;
str_ss={'$m_J=1$','$m_J=0$','$m_J=-1$','loss'};

hfig=figure('Units', 'normalized', 'Position', [0.2,0.2,0.3,0.3]);
hfig.Renderer='painters';
hold on;
pleg=NaN(1,4);      % data for legend
for ii=1:4
    tp=ploterr(1e6*par_T,Navg(:,ii),[],Nse(:,ii),'o','hhxy',0);
    set(tp(1),'color',c(ii,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'MarkerFaceColor',cl(ii,:),...
        'DisplayName',str_ss{ii});
    set(tp(2),'color',c(ii,:),'LineWidth',line_wid);
    pleg(ii)=tp(1);     % line data to show in legend
    
    % fit
    tpfit=plot(1e6*tau_fit,Nfit(ii,:),'LineWidth',line_wid,'Color',c(ii,:));        %0.5*(c(ii,:)+cl(ii,:))
    uistack(tpfit,'bottom');
end

ax=gca;
set(ax,'Layer','Top');     % graphics axes should be always on top
axis tight;
lgd=legend(pleg);
lgd.FontSize=11;
lgd.Location='East';
% lgd.Title.String='$m_J$';
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Atom counts');
ax.FontSize=11;
ax.LineWidth=1.2;

%% Misc number loss analysis
N0=Navg(:,2);

hfig=figure('Units', 'normalized', 'Position', [0.2,0.2,0.3,0.3]); 
hfig.Renderer='painters';

plot(1e6*par_T,Nloss./(N0+Nloss),'ko-','LineWidth',1.2);

axis tight;
ylim([-1,1]);
ax=gca;
ax.FontSize=11;
ylabel('$N_{\textrm{loss}}/(N_{0}+N_{\textrm{loss}})$');
xlabel('Pulse duration $\tau$ [$\mu$s]');
grid on;
