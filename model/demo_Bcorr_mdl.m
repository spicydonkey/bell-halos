% Modelling experimental imperfections for Bell correlations
%   - state purity
%   - Rabi amplitude
%
% DKS
% 2019-07-10

%% DEMO
% CONFIG
theta = linspace(0,2*pi);
p = 0.9;
alpha = 0.8*pi/2;

%%% PLOT
H=figure('Name','demo');

B_theory = Bcorr_mdl([p,alpha],theta);


% figure;
plot(theta/pi,B_theory);

box on;
grid on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
titlestr = sprintf('$p =$ %0.2g; $\\alpha = %0.2g \\cdot (\\pi/2)$',p,alpha/(pi/2));
title(titlestr);


%% Purity dependency (perfect Rabi oscillation)
theta = linspace(0,2*pi);
p = linspace(0,1,5);
alpha = pi/2;

B_theory = arrayfun(@(x) Bcorr_mdl([x,alpha],theta),p,'uni',0);


%%% plot
% config
col = palette_cmap(length(p),@parula);      % without the max valued color

% main
H=figure('Name','psi+_purity');
hold on;

tp=[];
for ii=1:length(p)
    tp(ii) = plot(theta/pi,B_theory{ii});
    set(tp(ii),'DisplayName',num2str(p(ii),3));
    set(tp(ii),'Color',col(ii,:));
    set(tp(ii),'LineWidth',1);
end

grid on;
box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
title('purity vs. correlator')

lgd=legend(tp);
lgd.Title.String = '$\vert\Psi^+\rangle$ purity $p$';


%% Rabi amplitude (purity=1)
theta = linspace(0,2*pi);
p = 1;
alpha = linspace(0,1,10)*pi/2;

B_theory = arrayfun(@(x) Bcorr_mdl([p,x],theta),alpha,'uni',0);

%%% plot
% config
col = palette_cmap(length(alpha),@parula);      % without the max valued color


% main
H=figure('Name','rabi_amplitude');
hold on;

tp=[];
for ii=1:length(alpha)
    tp(ii) = plot(theta/pi,B_theory{ii});
    set(tp(ii),'DisplayName',num2str(alpha(ii)/(pi/2),3));
    set(tp(ii),'Color',col(ii,:));
        set(tp(ii),'LineWidth',1);
end

grid on;
box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
title('Rabi amplitude vs. correlator')


lgd=legend(tp);
lgd.Title.String = 'polar angle $\alpha/(\pi/2)$';