
close all
clear all

load('data_2TimeLapses.mat')
addpath('functions')

%% Classical petrophysical 4D inversion

% % % petrophysical domain discretization
ndiscr = 25;
phidomain = linspace(0, 0.4, ndiscr)';   
swdomain = linspace(0, 1, ndiscr)';   
mdomain = [phidomain swdomain];
% elastic domain discretization
AIdomain  = linspace(min(acousticimpedance(:)), max(acousticimpedance(:)),ndiscr)';
VPVSdomain = linspace(min(VPVS(:)), max(VPVS(:)),ndiscr)';
ddomain =[AIdomain VPVSdomain ];

% kernel bandwidths 
h = 5;
hm(1) = (max(phidomain)-min(phidomain))/h;
hm(2) = (max(swdomain)-min(swdomain))/h;
hd(1) = (max(AIdomain)-min(AIdomain))/h;
hd(2) = (max(VPVSdomain)-min(VPVSdomain))/h;


N_SubSamp = 10000;

indices_highSat = find(saturation10yrs>0.85);
indices_highSat = indices_highSat(randperm( numel(indices_highSat), N_SubSamp ));
indices_midSat = find(saturation10yrs<0.85 & saturation10yrs>0.4 );
indices_midSat = indices_midSat(randperm( numel(indices_midSat), N_SubSamp ));
indices_lowSat = find(saturation10yrs<0.4);
indices_lowSat = indices_lowSat(randperm(numel(indices_lowSat),N_SubSamp ));

indices = [indices_highSat' indices_midSat'  indices_lowSat' ];

mtrain = [porosity(indices); saturation10yrs(indices)]';
dtrain = [acousticimpedance10yrs(indices); VPVS10yrs(indices) ]';

phi_reference = porosity(:,:,5);
sw_reference = saturation10yrs(:,:,5);
data_AI = acousticimpedance10yrs(:,:,5);
data_VPVS = VPVS10yrs(:,:,5);
dcond = [data_AI(:) data_VPVS(:) ];

inverted_properties = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, 25);

Phimap = reshape(inverted_properties(:,1),150,200);
Swmap = reshape(inverted_properties(:,2),150,200);

figure
subplot(221)
imagesc(phi_reference)
caxis([0.05 0.3])
title('Reference Porosity')
subplot(222)
imagesc(sw_reference)
title('Reference Sw')
caxis([0 1])
subplot(223)
imagesc(Phimap)
title('Estimated Porosity')
caxis([0.05 0.3])
subplot(224)
imagesc(Swmap)
title('Estimated Sw')
caxis([0 1])


%% 4D petrophysical inversion

ndiscr = 25;
phidomain = linspace(0, 0.4, ndiscr)';   
swdomain = linspace(0, 1, ndiscr)';   
mdomain = [phidomain swdomain swdomain];
% elastic domain discretization
AIdomain  = linspace(min(acousticimpedance10yrs(:)), max(acousticimpedance(:)),ndiscr)';
VPVSdomain = linspace(min(VPVS(:)), max(VPVS(:)),ndiscr)';
ddomain =[AIdomain AIdomain VPVSdomain VPVSdomain ];

% kernel bandwidths 
h = 5;
hm(1) = (max(phidomain)-min(phidomain))/h;
hm(2) = (max(swdomain)-min(swdomain))/h;
hm(3) = (max(swdomain)-min(swdomain))/h;
hd(1) = (max(AIdomain)-min(AIdomain))/h;
hd(2) = (max(AIdomain)-min(AIdomain))/h;
hd(3) = (max(VPVSdomain)-min(VPVSdomain))/h;
hd(4) = (max(VPVSdomain)-min(VPVSdomain))/h;


N_SubSamp = 2000;%20000;
 
indices_highSat = find(saturation10yrs>0.85);
indices_highSat = indices_highSat(randperm( numel(indices_highSat), N_SubSamp ));
indices_midSat = find(saturation10yrs<0.85 & saturation10yrs>0.4 );
indices_midSat = indices_midSat(randperm( numel(indices_midSat), N_SubSamp ));
indices_lowSat = find(saturation10yrs<0.4);
indices_lowSat = indices_lowSat(randperm(numel(indices_lowSat),N_SubSamp ));

indices = [indices_highSat' indices_midSat'  indices_lowSat' ];

mtrain = [porosity(indices);  saturation5yrs(indices); saturation10yrs(indices)]';
dtrain = [acousticimpedance5yrs(indices); acousticimpedance10yrs(indices); VPVS5yrs(indices); VPVS10yrs(indices) ]';

slice = 5;
phi_reference = porosity(:,:,slice);
sw_reference5yrs = saturation5yrs(:,:,slice);
sw_reference10yrs = saturation10yrs(:,:,slice);
data_AI_5yrs = acousticimpedance5yrs(:,:,slice);
data_AI_10yrs = acousticimpedance10yrs(:,:,slice);
data_VPVS_5yrs = VPVS5yrs(:,:,slice);
data_VPVS_10yrs = VPVS10yrs(:,:,slice);
dcond = [data_AI_5yrs(:) data_AI_10yrs(:) data_VPVS_5yrs(:) data_VPVS_10yrs(:) ];
ns = size(dcond,1);

inverted_properties = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, 25);

Phimap = reshape(inverted_properties(:,1),150,200);
Sw1map = reshape(inverted_properties(:,2),150,200);
Sw2map = reshape(inverted_properties(:,3),150,200);

figure
subplot(231)
imagesc(phi_reference)
caxis([0.05 0.3])
title('Reference Porosity')
subplot(232)
imagesc(sw_reference5yrs)
caxis([0 1])
title('Reference Sw 5 years')
subplot(233)
imagesc(sw_reference10yrs)
caxis([0 1])
title('Reference Sw 10 years')
subplot(234)
imagesc(Phimap)
caxis([0.05 0.3])
title('Estimated Porosity')
subplot(235)
imagesc(Sw1map)
caxis([0 1])
title('Estimated Sw 5 years')
subplot(236)
imagesc(Sw2map)
caxis([0 1])
title('Estimated Sw 10 years')



























