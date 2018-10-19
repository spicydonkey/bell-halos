%% Integrated g2 correlation
%   to be free of correlation length/volume
% DKS
% 20181019

%% CONFIG
dk_intg=0.1;     % half-width around dk=0 for g2 integration

flag_save_data=false;

path_dir='C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\exp5_bell_yrot';
% path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';

% get data (.mat) files to load
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

%% main
tau=cell(n_data,1);
G_intg=cell(n_data,1);
E=cell(n_data,1);

for ii=1:n_data
    %% load data
    tdataname=dataname{ii};
    tS=load(fullfile(path_dir,tdataname));
    tidx0=tS.idx_dk0;
    tau{ii}=tS.par_T;
    tg2=tS.g2;
    
    %% routine
    tdk_bin=tS.dk_cent_vec(2)-tS.dk_cent_vec(1);
    tIlim=round(dk_intg/tdk_bin);
    tidx_intg=tidx0+[-tIlim:tIlim];
    tV_intg=(tdk_bin*(2*tIlim+1))^3;    % integration vol (norm units)
    
    tG_int=cellfun(@(s) cellfun(@(g) sumall(g(tidx_intg,tidx_intg,tidx_intg)-1),s),tg2,'UniformOutput',false);
    tG_int=cat(1,tG_int{:});
    G_intg{ii}=tG_int;
    
    tE=Ecorr(tG_int(:,1),tG_int(:,2),tG_int(:,3),tG_int(:,3));
    E{ii}=tE;
end
clearvars tS;

tau=cat(1,tau{:});
G_intg=cat(1,G_intg{:});
E=cat(1,E{:});

%% vis
% configs
[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;

%% Integrated g2
h=figure('Name','integrated_g2',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');

hold on;

pleg=NaN(3,1);
for ii=1:3
    tp=myploterr(1e6*tau,G_intg(:,ii),[],[],mark_typ{ii},c(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
%     set(tp(2),'LineWidth',line_wid,'DisplayName','');
end
box on;
ax=gca;
lgd=legend(pleg);
lgd.FontSize=fontsize;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Integrated correlation $\mathcal{G}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
axis tight;

%% correlator
h=figure('Name','correlator',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');

hold on;

tp=myploterr(1e6*tau,E,[],[],mark_typ{1},c(1,:));
set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName','Data');
pleg=tp(1);     % line data to show in legend
%     set(tp(2),'LineWidth',line_wid,'DisplayName','');

box on;
ax=gca;
lgd=legend(pleg);
lgd.FontSize=fontsize;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Correlator $\mathcal{B}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
axis tight;
