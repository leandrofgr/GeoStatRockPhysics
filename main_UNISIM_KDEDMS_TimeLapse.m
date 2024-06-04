
%clear
close all

addpath(genpath('utils'))
addpath(genpath('SeReM'))

%% read UNISIM data and apply optional crops for testing
%read_UNISIM

Phi_ = Phi;
Sw13_ = Sw13;
Sw24_ = Sw24;
Ip13_ = Ip13;
Ip24_ = Ip24;
VPVS13_ = VPVS13;
VPVS24_ = VPVS24;

% % Cropping
% Imin = 60;
% Jmin = 50;
% Size = 80;
% Phi_ = Phi(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% Sw13_ = Sw13(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% Sw24_ = Sw24(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% Ip13_ = Ip13(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% Ip24_ = Ip24(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% VPVS13_ = VPVS13(Imin:Imin + Size-1,Jmin:Jmin+Size-1);
% VPVS24_ = VPVS24(Imin:Imin + Size-1,Jmin:Jmin+Size-1);

cond_pos = [];
cond_value = [];

% Adjusting well position to the crop and constructing the conditional data for kriging
figure
imagesc(Sw24_)
hold all
for well=1:length(WELLS)
    % Cropping
    %WELLS(well).i = WELLS(well).i - (Imin - 1);
    %WELLS(well).j = WELLS(well).j - (Jmin - 1);
    plot(WELLS(well).j,WELLS(well).i,'r+','LineWidth',2)
    text(WELLS(well).j,WELLS(well).i,WELLS(well).name)
    cond_pos = [cond_pos ; WELLS(well).i WELLS(well).j];
    cond_value = [  cond_value ; WELLS(well).Ip13 WELLS(well).VPVS13 WELLS(well).Ip24 WELLS(well).VPVS24 WELLS(well).Phi WELLS(well).Sw13 WELLS(well).Sw24  ];
end


% Delete well out of the crop region
WELLS_ = WELLS;
list2delete = [];
for well=1:length(WELLS)
    if WELLS(well).i <0 || WELLS(well).i >80 || WELLS(well).j <0 || WELLS(well).j >80 
        list2delete = [list2delete well];
    end
end

% WELLS(list2delete) = [];

referece_variables = [Ip13(:) VPVS13(:) Ip24(:) VPVS24(:) Phi(:) Sw13(:) Sw24(:) ];
%generate_histograms(referece_variables)
generate_histograms([Phi(:) Sw13(:) Sw24(:) Ip13(:) VPVS13(:) Ip24(:) VPVS24(:) ])
sgtitle('Joint distribution of the reference models');

%%%% FIGURE - Reference petrophysical properties
figure
subplot(231)
imagesc(Phi_)
caxis([0 0.35])
hold all
plot_wells(WELLS)
title('Porosity - Reference')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(232)
imagesc(Sw13_)
caxis([0 0.8])
hold all
plot_wells(WELLS)
title('Sw 2013 - Reference')
xlabel('J')
grid
colorbar
subplot(233)
imagesc(Sw24_)
caxis([0 0.8])
hold all
plot_wells(WELLS)
title('Sw 2024 - Reference')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])
sgtitle('Reference models');

%%%% FIGURE - Elastic properties
figure
subplot(231)
imagesc(Ip13_)
%caxis([0 0.35])
hold all
plot_wells(WELLS)
title('P-impedance 2013')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(232)
imagesc(VPVS13_)
caxis([1.35 1.7])
hold all
plot_wells(WELLS)
title('Vp/Vs 2013')
xlabel('J')
grid
colorbar
subplot(234)
imagesc(Ip24_)
%caxis([0 0.8])
hold all
plot_wells(WELLS)
title('P-impedance 2013')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])
subplot(235)
imagesc(VPVS24_)
caxis([1.35 1.7])
hold all
plot_wells(WELLS)
title('Vp/Vs 2024')
xlabel('J')
grid
colorbar
sgtitle('Elastic properties');



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
mtrain = [Phi_train, sw_train13, sw_train24];
dtrain = [Ip_train13, VPVS_train13, Ip_train24, VPVS_train24];

prior_variables = [dtrain mtrain];
%generate_histograms(prior_variables)
names = {'\phi', 's_{w1}', 's_{w2}','I_{p1}','\alpha / \beta_1',  'I_{p2}','\alpha / \beta_2'};
generate_histograms([mtrain dtrain],[0 1500],names)
sgtitle('Joint distribution by Monte Carlo sampling');

dcond = [Ip13_(:), VPVS13_(:), Ip24_(:), VPVS24_(:)];

%% INVERSION - DMS for computing the mean value and the distributions. Kriging must be done as a prior. 

%inverted_properties = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, 20);
%Phi_inverted = reshape(inverted_properties(:,1), size(Phi_));
%Sw13_inverted = reshape(inverted_properties(:,2), size(Phi_));
%Sw24_inverted = reshape(inverted_properties(:,3), size(Phi_));

% When running with n_sim=1, it considers phi(g) = 0.5 to estimate the median model
[DME_result_median_uncond] = DMS(I, J, 35, 'sph', [0 0 0], 0.1, prior_variables, [], [], 1, cond_variables );

DME_result_median_uncond = cat(4, DME_result_median_uncond{:});
Phi_median_uncond = squeeze(DME_result_median_uncond(5,:,:,:));
Sw13_median_uncond = squeeze(DME_result_median_uncond(6,:,:,:));
Sw24_median_uncond = squeeze(DME_result_median_uncond(7,:,:,:));

% When running with n_sim=1, cond_pos and cond_value, it considers phi(g) as the kriging of the transformed hard data conditioning
[DME_result_median_cond] = DMS(I, J, 35, 'sph', [0 0 0], 0.1, prior_variables, cond_pos, cond_value, 1, cond_variables );

DME_result_median_cond = cat(4, DME_result_median_cond{:});
Phi_median_cond = squeeze(DME_result_median_cond(5,:,:,:));
Sw13_median_cond = squeeze(DME_result_median_cond(6,:,:,:));
Sw24_median_cond = squeeze(DME_result_median_cond(7,:,:,:));



%%%% Median by defining Phi(g) = 0.5 and kriging
figure
subplot(231)
imagesc(Phi_median_uncond)
caxis([0 0.35])
hold all
plot_wells(WELLS,'Phi')
title('Porosity - Unconditioned Median')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(232)
imagesc(Sw13_median_uncond)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw13')
title('Sw 2013 - Unconditioned Median ')
xlabel('J')
grid
colorbar
subplot(233)
imagesc(Sw24_median_uncond)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw24')
title('Sw 2024 - Unconditioned Median ')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])

subplot(234)
imagesc(Phi_median_cond)
caxis([0 0.35])
hold all
plot_wells(WELLS,'Phi')
title('Porosity - Conditioned Median')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(235)
imagesc(Sw13_median_cond)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw13')
title('Sw 2013 - Conditioned Median')
xlabel('J')
grid
colorbar
subplot(236)
imagesc(Sw24_median_cond)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw24')
title('Sw 2024 - Conditioned Median')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])
sgtitle('Median models');


%% INVERSION -  DMS for sampling spatially correlated simulations with conditioning of hard data
n_sims = 100;
cond_variables = dcond;
[I,J,K] = size(Ip13_); 
[logs_simulated_all] = DMS(I, J, 35, 'sph', [0 0 0], 0.1, prior_variables, cond_pos, cond_value, n_sims, cond_variables );

logs_simulated_all_ = cat(4, logs_simulated_all{:});

Phi_std = squeeze(std(logs_simulated_all_(5,:,:,:),[],4));
Phi_mean = squeeze(median(logs_simulated_all_(5,:,:,:),4));
[Phi_likely] = compute_most_likely((squeeze(logs_simulated_all_(5,:,:,:) )));

Sw13_std = squeeze(std(logs_simulated_all_(6,:,:,:),[],4));
Sw13_mean = squeeze(median(logs_simulated_all_(6,:,:,:),4));
[Sw13_likely] = compute_most_likely((squeeze(logs_simulated_all_(6,:,:,:) )));

Sw24_std = squeeze(std(logs_simulated_all_(7,:,:,:),[],4));
Sw24_mean = squeeze(median(logs_simulated_all_(7,:,:,:),4));
[Sw24_likely] = compute_most_likely((squeeze(logs_simulated_all_(7,:,:,:) )));


%%%% FIGURE - Histograms, reference, prior ans posterior
simulations2histogram = reshape(logs_simulated_all_, 7,I*J, n_sims);
for prop = 1:7
    for sim = 1:n_sims
        simulations2histogram(prop,isnan(reshape(Phi_,I*J,1)),sim) = nan;
    end
end
simulations2histogram = reshape(simulations2histogram,7,[])';

simulations2histogram = [simulations2histogram(:,5:7) simulations2histogram(:,1:4) ];
%simulations2histogram(simulations2histogram<1e-3) = nan;
simulations2histogram(isnan(reshape(Phi_,I*J,1)),:) = nan;
generate_histograms(simulations2histogram,[0 20000],names)


%% Median and std. dev. of simulations
figure
subplot(231)
imagesc(Phi_mean)
caxis([0 0.35])
hold all
plot_wells(WELLS,'Phi')
title('Porosity - Median of Simulations')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(232)
imagesc(Sw13_mean)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw13')
title('Sw 2013 - Median of Simulations')
xlabel('J')
grid
colorbar
subplot(233)
imagesc(Sw24_mean)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw24')
title('Sw 2024 - Median of Simulations')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])

subplot(234)
imagesc(Phi_std)
caxis([0 0.05])
hold all
plot_wells(WELLS)
title('Porosity - Std. Dev.')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(235)
imagesc(Sw13_std)
caxis([0 0.35])
hold all
plot_wells(WELLS)
title('Sw 2013 - Std. Dev.')
xlabel('J')
grid
colorbar
subplot(236)
imagesc(Sw24_std)
caxis([0 0.35])
hold all
plot_wells(WELLS)
title('Sw 2024 - Median Cond')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])
sgtitle('Median models and Std from simulations');

%%%% FIGURE - Multiple simulations

sim = 1;
Phi_sim = squeeze(logs_simulated_all_(5,:,:,sim));
Sw13_sim = squeeze(logs_simulated_all_(6,:,:,sim));
Sw24_sim = squeeze(logs_simulated_all_(7,:,:,sim));
Phi_sim(Phi_sim<0.01) = 0.01; Phi_sim(isnan(Phi_)) = nan;
Sw13_sim(Sw13_sim<0.01) = 0.01; Sw13_sim(isnan(Phi_)) = nan;
Sw24_sim(Sw24_sim<0.01) = 0.01; Sw24_sim(isnan(Phi_)) = nan;


figure
subplot(231)
imagesc(Phi_sim)
caxis([0 0.35])
hold all
plot_wells(WELLS,'Phi')
title('Porosity - Simulation 1')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(232)
imagesc(Sw13_sim)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw13')
title('Sw 2013 - Simulation 1')
xlabel('J')
grid
colorbar
subplot(233)
imagesc(Sw24_sim)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw24')
title('Sw 2024 - Simulation 1')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])

sim = 2;
Phi_sim = squeeze(logs_simulated_all_(5,:,:,sim));
Sw13_sim = squeeze(logs_simulated_all_(6,:,:,sim));
Sw24_sim = squeeze(logs_simulated_all_(7,:,:,sim));
Phi_sim(Phi_sim<0.01) = 0.01; Phi_sim(isnan(Phi_)) = nan;
Sw13_sim(Sw13_sim<0.01) = 0.01; Sw13_sim(isnan(Phi_)) = nan;
Sw24_sim(Sw24_sim<0.01) = 0.01; Sw24_sim(isnan(Phi_)) = nan;

subplot(234)
imagesc(Phi_sim)
caxis([0 0.35])
hold all
plot_wells(WELLS,'Phi')
title('Porosity - Simulation 2')
ylabel('I')
xlabel('J')
grid
colorbar
subplot(235)
imagesc(Sw13_sim)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw13')
title('Sw 2013 - Simulation 2')
xlabel('J')
grid
colorbar
subplot(236)
imagesc(Sw24_sim)
caxis([0 0.8])
hold all
plot_wells(WELLS,'Sw24')
title('Sw 2024 - Simulation 2')
xlabel('J')
grid
colorbar
colormap([1 1 1; parula])
sgtitle('Geostatistical Simulations');

%% RMSE calculation
index_non_nan = ~isnan(Phi);


Phi_rmse_inverted = sqrt( mean( (Phi(index_non_nan ) - Phi_median_uncond(index_non_nan )).^2 ) );
Sw13_rmse_inverted = sqrt( mean( (Sw13(index_non_nan ) - Sw13_median_uncond(index_non_nan )).^2 ) );
Sw24_rmse_inverted = sqrt( mean( (Sw24(index_non_nan ) - Sw24_median_uncond(index_non_nan )).^2 ) );

uncond_rmse = Phi_rmse_inverted + Sw13_rmse_inverted + Sw24_rmse_inverted 


Phi_rmse_mean = sqrt( mean( (Phi(index_non_nan ) - Phi_median_cond(index_non_nan )).^2 ) );
Sw13_rmse_mean = sqrt( mean( (Sw13(index_non_nan ) - Sw13_median_cond(index_non_nan )).^2 ) );
Sw24_rmse_mean = sqrt( mean( (Sw24(index_non_nan ) - Sw24_median_cond(index_non_nan )).^2 ) );

cond_rmse = Phi_rmse_mean + Sw13_rmse_mean + Sw24_rmse_mean 


Phi_rmse_mean_sim = sqrt( mean( (Phi(index_non_nan ) - Phi_mean(index_non_nan )).^2 ) );
Sw13_rmse_mean_sim = sqrt( mean( (Sw13(index_non_nan ) - Sw13_mean(index_non_nan )).^2 ) );
Sw24_rmse_mean_sim = sqrt( mean( (Sw24(index_non_nan ) - Sw24_mean(index_non_nan )).^2 ) );

cond_rmse_sim = Phi_rmse_mean_sim + Sw13_rmse_mean_sim + Sw24_rmse_mean_sim


















