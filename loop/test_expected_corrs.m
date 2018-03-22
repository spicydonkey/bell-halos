%% Expected spin-correlations with a simple detuning model


%% Load data
% load('anal_rabi.mat');


%% config
tpulse=3.8e-6;    % pulse duration (s)


%% geometry
rabiGamma=acos(sqrt(rabiAmpConst));        % torque vector's tilt angle rel to XY-plane
rabiTheta=tpulse*rabiOmega;                 % Rabi phase

rabizTheta=2*asin(cos(rabiGamma).*sin(rabiTheta/2));    % end Bloch vector angle from pole


%% mixed state model

E_mixed=Emixed(rabizTheta,flip_bb(rabizTheta,1));

figure('Name','mixed state');
plotFlatMapWrappedRad(az,el,E_mixed);
cb=colorbar('southoutside');
cb.Label.String='E';
title(sprintf('tpulse = %0.2g',tpulse));


%% quantum entanglement model

S_bloch=NaN([size(rabiGamma),3]);       % bloch vectors for rotated spin-up

for ii=1:size(rabiGamma,1)
    for jj=1:size(rabiGamma,2)
        S_bloch(ii,jj,:)=rotBlochVector(rabiGamma(ii,jj),rabiTheta(ii,jj));
    end
end


% Plot Bloch vectors
figure('Name','Bloch vectors');
for ii=1:size(S_bloch,1)
    for jj=1:size(S_bloch,2)
        if rand()<0.1
            hold on;
            
            tS=S_bloch(ii,jj,:);
            plot3([0,tS(1)],[0,tS(2)],[0,tS(3)],'-');
            scatter3(tS(1),tS(2),tS(3),'o','filled');
        end
    end
end
% draw the sphere
hold on;
[xs,ys,zs]=sphere(30);
hSphSurf=surf(xs,ys,zs);
set(hSphSurf,'FaceColor',[0 0 1], ...
    'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none')
% axis config
axis equal;
axlim=1.1*[-1,1];
xlim(axlim);
ylim(axlim);
zlim(axlim);
xlabel('x');
ylabel('y');
zlabel('z');
box on;
view(3);
hold off;
camlight;


% expected correlation
dotSpair=sum( S_bloch.*(-flip_bb(S_bloch,1)) , 3);      % dot-prod entangled pairs
dth_bloch=acos(dotSpair);       % relative angle between pair

Eqm = @(dth) 2*cos(dth/2).^2 - 1;      % quantum mechanics - from relative angles between pairs

E_qm=Eqm(dth_bloch);

figure('Name','QM');
plotFlatMapWrappedRad(az,el,E_qm);
cb=colorbar('southoutside');
cb.Label.String='E';
title(sprintf('tpulse = %0.2g',tpulse));