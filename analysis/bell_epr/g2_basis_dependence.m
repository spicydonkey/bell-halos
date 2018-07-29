% Basis dependence of g2 correlations of interest 
%   - g2 types will be momentum-BB, spin-anticorr and corr (there are 11,00)
%   - theta: global rotation
%
% DKS
%

% demo data as a placeholder
theta = linspace(0,pi,1e3);

g2_amp = 30;
g2_00 = g2_amp * sin(theta).^2 + 1;
g2_11 = g2_amp * sin(theta).^2 + 1;
g2_01 = g2_amp * cos(theta).^2 + 1;        % is the same as 10

%% demo plot
%%% graphics config
font_siz_reg=12;
font_siz_sml=10;
font_siz_lrg=14;
% mark_siz=7;
lwid=1.5;
cols = palette(3);
lins = {'-','--','-.'};

% plot
figure;
hold on;
p00=plot(theta,g2_00,...
    'Color',cols(1,:),'LineStyle',lins{1},'LineWidth',lwid,...
    'DisplayName','$\uparrow\uparrow$');
p11=plot(theta,g2_11,...
    'Color',cols(2,:),'LineStyle',lins{2},'LineWidth',lwid,...
    'DisplayName','$\downarrow\downarrow$');
p01=plot(theta,g2_01,...
    'Color',cols(3,:),'LineStyle',lins{3},'LineWidth',lwid,...
    'DisplayName','$\uparrow\downarrow$');

xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});

xlabel('$\theta$');
ylabel('$g^{(2)}$');

xlim([0,pi]);
% Yll=ylim;
% axis tight;
% ylim([0,Yll(2)]);
% axis square;

lgd = legend([p00,p11,p01]);
set(lgd,'FontSize',font_siz_reg);
box on;

set(gca,'FontSize',font_siz_reg);
