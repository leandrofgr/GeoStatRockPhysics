
%close all
clear all

load('data/data_2TimeLapses.mat')
addpath(genpath('SeReM'))

%% Classical petrophysical inversion

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


N_SubSamp = 2000;

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
ns = size(dcond,1);


Ppost_classic = RockPhysicsKDEInversion(mtrain, dtrain, mdomain, ddomain, dcond, hm, hd);


for i=1:ns
    Ppostjoint=reshape(Ppost_classic(i,:),length(phidomain),length(swdomain));
    Ppostphi(i,:)=sum(squeeze(Ppostjoint),2);
    Ppostsw(i,:)=sum(squeeze(Ppostjoint),1)';
    Ppostphi(i,:)=Ppostphi(i,:)/sum(Ppostphi(i,:));
    Ppostsw(i,:)=Ppostsw(i,:)/sum(Ppostsw(i,:));
    [~,Phimapind]=max(Ppostphi(i,:));
    [~,Swmapind]=max(Ppostsw(i,:));
    Phimap(i)=phidomain(Phimapind);
    Swmap(i)=swdomain(Swmapind);
end

Phimap = reshape(Phimap,150,200);
Swmap = reshape(Swmap,150,200);

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


Ppost = RockPhysicsKDEInversion(single(mtrain), single(dtrain), single(mdomain), single(ddomain), single(dcond), hm, hd);


for i=1:ns
    Ppostjoint=reshape(Ppost(i,:),length(phidomain),length(swdomain),length(swdomain));
    Ppostphi(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),2);
    Ppostsw1(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),1);
    Ppostsw2(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),2)),1);
    Ppostphi(i,:)=Ppostphi(i,:)/sum(Ppostphi(i,:));
    Ppostsw1(i,:)=Ppostsw1(i,:)/sum(Ppostsw1(i,:));
    Ppostsw2(i,:)=Ppostsw2(i,:)/sum(Ppostsw2(i,:));
    [~,Phimapind]=max(Ppostphi(i,:));
    [~,Sw1apind]=max(Ppostsw1(i,:));
    [~,Sw2mapind]=max(Ppostsw2(i,:));
    Phimap(i)=phidomain(Phimapind);
    Sw1map(i)=swdomain(Sw1apind);
    Sw2map(i)=swdomain(Sw2mapind);
end

Phimap = reshape(Phimap,150,200);
Sw1map = reshape(Sw1map,150,200);
Sw2map = reshape(Sw2map,150,200);

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



























