
clear

addpath(genpath('utils'))
addpath(genpath('SeReM'))

% Here, I decided to separate the part of loading the data because I believe
% we will start to analyse the methodology parameter. And to compare the
% results with different parameterization we should use the same dataset.
% In this case, the following 2 lines should be commented. 
read_UNISIM
%close all

Phi = Phi(176:255,71:150);
Sw13 = Sw13(176:255,71:150);
Sw24 = Sw24(176:255,71:150);
Ip13 = Ip13(176:255,71:150);
Ip24 = Ip24(176:255,71:150);
VPVS13 = VPVS13(176:255,71:150);
VPVS24 = VPVS24(176:255,71:150);


cond_pos = [];
cond_value = [];

figure
imagesc(Sw24)
hold all
for well=1:length(WELLS)
    WELLS(well).i = WELLS(well).i - (176 - 1);
    WELLS(well).j = WELLS(well).j - (71 - 1);
    plot(WELLS(well).j,WELLS(well).i,'r+','LineWidth',2)
    text(WELLS(well).j,WELLS(well).i,WELLS(well).name)
    cond_pos = [cond_pos ; WELLS(well).i WELLS(well).j];
    cond_value = [  cond_value ; WELLS(well).Ip13 WELLS(well).VPVS13 WELLS(well).Ip24 WELLS(well).VPVS24 WELLS(well).Phi WELLS(well).Sw13 WELLS(well).Sw24  ];
end

%criticalporo = 0.3;

%% PRIOR SAMPLING
%Porous and no porous
n_sim = 100000;
Phi_train = 0.1 + 0.05*randn(n_sim/2,1);
Phi_train = [ Phi_train ; 0.25 + 0.05*randn(n_sim/2,1) ];
sw_train13 = rand(n_sim,1);
sw_train24 = rand(n_sim,1);

% Include Shale
Phi_train = [ Phi_train ; 0.02 + 0.005*randn(n_sim/2,1) ];
sw_train13 = [ sw_train13 ; 0.9 + 0.02*randn(n_sim/2,1) ];
sw_train24 = [ sw_train24 ; 0.9 + 0.02*randn(n_sim/2,1) ];
Phi_train(Phi_train<0) = 0.001;
Phi_train(Phi_train>0.4) = 0.4;
sw_train13(sw_train13>1) = 0.999;
sw_train24(sw_train24>1) = 0.999;

% Simulate observed data (elastic properties) with noise
[Vp, Vs, Rho] = RPM_unisim(Phi_train, sw_train13, criticalporo );
Vp = Vp + std_vp*randn(size(Vp));
Vs = Vs + std_vs*randn(size(Vs));
Rho = Rho + std_rho*randn(size(Rho));
Ip_train13 = Vp.*Rho;
VPVS_train13 = Vp./Vs;

% Simulate observed data (elastic properties) with noise
[Vp, Vs, Rho] = RPM_unisim(Phi_train, sw_train24, criticalporo );
Vp = Vp + std_vp*randn(size(Vp));
Vs = Vs + std_vs*randn(size(Vs));
Rho = Rho + std_rho*randn(size(Rho));
Ip_train24 = Vp.*Rho;
VPVS_train24 = Vp./Vs;

% Final training data:
mtrain = [Phi_train , sw_train13, sw_train24];
dtrain = [Ip_train13, VPVS_train13, Ip_train24, VPVS_train24];

%% FIGURES ABOUT PRIOR TRAINING

% figure
% subplot(211)
% histogram(POR,'Normalization','probability')
% set(gca,'YScale','log')
% xlabel('Reference Porosity')
% title('Training vs reference')
% subplot(212)
% histogram(Phi_train,'Normalization','probability')
% set(gca,'YScale','log')
% xlabel('Training Porosity')

% figure
% subplot(121)
% plot(Ip_train13,VPVS_train13,'k.')
% hold all
% plot(Ip13,VPVS13,'b.')
% title('Training vs reference')
% subplot(122)
% plot(Ip_train24,VPVS_train24,'k.')
% hold all
% plot(Ip24,VPVS24,'b.')

% figure
% subplot(121)
% plot(Phi_train,sw_train13,'k.')
% hold all
% plot(Phi,Sw13,'b.')
% title('Training vs reference')
% subplot(122)
% plot(Phi_train,sw_train24,'k.')
% hold all
% plot(Phi,Sw24,'b.')

%% INVERSION 

dcond = [Ip13(:), VPVS13(:), Ip24(:), VPVS24(:)];

% inverted_properties = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, 25);
% 
% Phi_inverted = reshape(inverted_properties(:,1), size(Phi));
% Sw13_inverted = reshape(inverted_properties(:,2), size(Phi));
% Sw24_inverted = reshape(inverted_properties(:,3), size(Phi));

%% Unconditioned DMS

ref_variables = [dtrain mtrain];
[I,J,K] = size(Ip13); 
[logs_simulated_all] = DMS(I, J, 15, 'sph', [0 0 0], 0.1, ref_variables, cond_pos, cond_value, 1, []);
simulation_DMS = logs_simulated_all{1};

generate_histograms(ref_variables)
generate_histograms(reshape(simulation_DMS,7,I*J)')

figure
imagesc(squeeze(simulation_DMS(7,:,:)))
hold all
for well=1:length(WELLS)
    scatter(WELLS(well).j,WELLS(well).i,100,WELLS(well).Sw24,'filled')
    hold all
    text(WELLS(well).j,WELLS(well).i,WELLS(well).name)
end

figure
imagesc(Sw24)
hold all
for well=1:length(WELLS)
    scatter(WELLS(well).j,WELLS(well).i,100,WELLS(well).Sw24,'filled')
    hold all
    text(WELLS(well).j,WELLS(well).i,WELLS(well).name)
end




%% FIGURES

% data
figure
ax1 = subplot(221);
imagesc(Ip13)
caxis([4000 16000])
title('IP13')
ax2 = subplot(222);
imagesc(VPVS13)
caxis([1.35 1.75])
title('VPVS13')
ax3 = subplot(223);
imagesc(Ip24)
caxis([4000 16000])
title('IP24')
ax4 = subplot(224);
imagesc(VPVS24)
caxis([1.35 1.75])
title('VPVS24')
linkaxes([ax1, ax2, ax3, ax4], 'xy');

% Estimates

figure
ax1 = subplot(231);
imagesc(Phi)
caxis([0 0.35])
title('Reference Porosity')
ax2 = subplot(232);
imagesc(Sw13)
caxis([0 1])
title('Reference Sw13')
ax3 = subplot(233);
imagesc(Sw24)
caxis([0 1])
title('Reference Sw24')
ax4 = subplot(234);
imagesc(Phi_inverted)
caxis([0 0.35])
title('Estimated Porosity')
ax5 = subplot(235);
imagesc(Sw13_inverted)
caxis([0 1])
title('Estimated Sw13')
ax6 = subplot(236);
imagesc(Sw24_inverted)
caxis([0 1])
linkaxes([ax1, ax2, ax3, ax4 ax5 ax6], 'xy');
title('Estimated Sw24')































