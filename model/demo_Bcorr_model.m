% Modelling experimental imperfections for Bell correlations
%   - state purity
%   - Rabi amplitude
%
% DKS
% 2019-07-10

%% Model
% see David's logbook for derivation
%
% theta: rotation angle (e.g. pi-rotation is theta=pi)
% p: triplet purity (e.g. pure triplet is p=1)
% alpha: polar angle of rotation axis (e.g. resonance is alpha=pi/2)

Bcorr_mdl = @(theta,p,alpha) -p * (cos(alpha)^2 + sin(alpha)^2*cos(theta)).^2 ...
    + p * ( (sin(alpha)*cos(alpha)*(1-cos(theta))).^2 + (sin(alpha)*sin(theta)).^2);


%% plot
H=figure('Name','demo');

theta = linspace(0,2*pi);
p = 0.9;
alpha = 0.85*pi/2;

B_theory = Bcorr_mdl(theta,p,alpha);


% figure;
plot(theta/pi,B_theory);

box on;
grid on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
titlestr = sprintf('$p =$ %0.2g; $\\alpha = %0.2g \\cdot (\\pi/2)$',p,alpha/(pi/2));
title(titlestr);

%% Purity 
theta = linspace(0,2*pi);
p = linspace(0,1,5);
alpha = pi/2;

B_theory = arrayfun(@(x) Bcorr_mdl(theta,x,alpha),p,'uni',0);


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


%% Rabi amplitude
theta = linspace(0,2*pi);
p = 1;
alpha = linspace(0,1,10)*pi/2;

B_theory = arrayfun(@(x) Bcorr_mdl(theta,p,x),alpha,'uni',0);

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