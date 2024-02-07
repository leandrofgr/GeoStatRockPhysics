
%% Read models from petrel 
addpath(genpath('utils'))
addpath(genpath('SeReM'))

criticalporo = 0.3;

PermK = readNPY('data/UNISIM/PermK.npy');
NETGROSS = readNPY('data/UNISIM/NETGROSS.npy');
POR = readNPY('data/UNISIM/POR.npy');
POR(POR>0.4) = 0.399;
POR(POR<=0) = 0.02;
SW13 = readNPY('data/UNISIM/SW_2013.npy');
SW13(isnan(SW13)) = 0.9;
SW13(isnan(POR)) = nan;
SW24 = readNPY('data/UNISIM/SW_2024.npy');
SW24(isnan(SW24)) = 0.9;
SW24(isnan(POR)) = nan;

load('data/UNISIM/WellPositionInGrid.mat')

% 3D CASE:
%Phi = POR;
%Sw13 = SW13;
%Sw24 = SW24;

% 2D CASE:
Phi = POR(:,:,20);
Sw13 = SW13(:,:,20);
Sw24 = SW24(:,:,20);

%% Running RPM to simulate the observed data (elastic properties) for inversion 
std_vp = 50;
std_vs = 25;
std_rho = 0.05;
%correlation function to add an spatial correlated noise
if length(size(Sw24)) == 2
    correlation_function = construct_correlation_function_beta(20,10,Phi,2);
else
    correlation_function = construct_correlation_function_beta(20,10,Phi,2);
end

[Vp, Vs, Rho] = RPM_unisim(Phi, Sw13, criticalporo );
% Adding noise
Vp = Vp + std_vp*FFT_MA_3D(correlation_function, randn(size(Vp)));
Vs = Vs + std_vs*FFT_MA_3D(correlation_function, randn(size(Vp)));
Rho = Rho + std_rho*FFT_MA_3D(correlation_function, randn(size(Vp)));
Ip13 = Vp.*Rho;
VPVS13 = Vp./Vs;

[Vp, Vs, Rho] = RPM_unisim(Phi, Sw24, criticalporo );

% Adding noise
Vp = Vp + std_vp*FFT_MA_3D(correlation_function, randn(size(Vp)));
Vs = Vs + std_vs*FFT_MA_3D(correlation_function, randn(size(Vp)));
Rho = Rho + std_rho*FFT_MA_3D(correlation_function, randn(size(Vp)));
Ip24 = Vp.*Rho;
VPVS24 = Vp./Vs;


%% Extract well logs from the cubes

% Create wells strusct in main struct 
clear WELLS 
for well = 1:size(WellPositionInGrid,1)
   WELLS(well) = struct('name',WellPositionInGrid{well,1},'i',WellPositionInGrid{well,2},'j',WellPositionInGrid{well,3})
end
   
[WELLS] = add_welllog_from_cube(WELLS,Phi,'Phi');
[WELLS] = add_welllog_from_cube(WELLS,Sw13,'Sw13');
[WELLS] = add_welllog_from_cube(WELLS,Sw24,'Sw24');

[WELLS] = add_welllog_from_cube(WELLS,Ip13,'Ip13');
[WELLS] = add_welllog_from_cube(WELLS,Ip24,'Ip24');
[WELLS] = add_welllog_from_cube(WELLS,Ip24,'VPVS13');
[WELLS] = add_welllog_from_cube(WELLS,Ip24,'VPVS24');


    


%% FIGURES

figure 
imagesc(Sw24(:,:,1))
hold all
plot(WellPositionInGrid{1,3},WellPositionInGrid{1,2},'r+','LineWidth',2)


% data
figure
ax1 = subplot(221);
imagesc(Ip13(:,:,1))
caxis([4000 16000])
title('IP13')
ax2 = subplot(222);
imagesc(VPVS13(:,:,1))
caxis([1.35 1.75])
title('VPVS13')
ax3 = subplot(223);
imagesc(Ip24(:,:,1))
caxis([4000 16000])
title('IP24')
ax4 = subplot(224);
imagesc(VPVS24(:,:,1))
caxis([1.35 1.75])
title('VPVS24')
linkaxes([ax1, ax2, ax3, ax4], 'xy');

figure
subplot(131)
scatter(Ip24(:), VPVS24(:), 25, Sw24(:), 'filled');
grid
subplot(132)
scatter(Phi(:), Ip24(:), 25, Sw24(:), 'filled');
grid
subplot(133)
scatter(Sw24(:), VPVS24(:), 25, Phi(:), 'filled');
grid

figure
ax1 = subplot(121);
imagesc(Sw24(:,:,1))
ax2 = subplot(122);
imagesc(VPVS24(:,:,1) - VPVS13(:,:,1))