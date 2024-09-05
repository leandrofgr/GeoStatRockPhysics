%% Rock physics inversion Driver %%
% TO DO  TO DO  TO DO  TO DO  TO DO  TO DO  TO DO  
% TO DO  TO DO  TO DO  TO DO  TO DO  TO DO  TO DO  

%% Load data
addpath(genpath('SeReM'));
load SeReM/Data/data4.mat
load data/cmap_3facies.mat
Facies3(Facies3==3) = 1.1;
Facies3(Facies3==1) = 3;

%% Monte carlo sampling of the joint distribution
% Prior of petrophysical properties p(r) with 2 facies:
n_sim = 10*[400 800 300];
petro_sim = [];
for facie = 1:length(unique(Facies))
    MU = mean([Phi(Facies==facie) Clay(Facies==facie) Sw(Facies==facie) ]) ;
    COV = cov([Phi(Facies==facie) Clay(Facies==facie) Sw(Facies==facie) ]);
    if facie == 1
        MU(3) = 0.2;
    end
    petro_sim = [ petro_sim; mvnrnd(MU, COV, n_sim(facie))];
end
% Additional facies for brine sand
MU = [0.24 0.07 0.95];
COV = diag([ 0.017 0.06 0.02 ]).^2;
petro_sim = [ petro_sim; mvnrnd(MU, COV, n_sim(3))];

Phi_sim = petro_sim(:,1);
Clay_sim = petro_sim(:,2);
Sw_sim = petro_sim(:,3);
Phi_sim(Phi_sim<0)=0; Phi_sim(Phi_sim>0.4)=0.4;
Clay_sim(Clay_sim<0)=0; Clay_sim(Clay_sim>0.8)=0.8;
Sw_sim(Sw_sim<0)=0.02; Sw_sim(Sw_sim>1)=0.98;

%% Rock-physics model - Linearization base on the logs
% Sampling elastic properties from p(m|r):
nd = size(Phi,1);
nv = 3;
R = zeros(nd,nv+1);
X = [Phi Clay Sw ones(size(Phi))];
R(1,:) = regress(Vp,X); 
R(2,:) = regress(Vs,X); 
R(3,:) = regress(Rho,X); 
[Vpsim, Vssim, Rhosim] = LinearizedRockPhysicsModel(Phi_sim, Clay_sim, Sw_sim, R);

noise_perc = 0.05;
Vpsim = Vpsim + noise_perc * std(Vpsim) * randn(size(Vpsim));
Vssim = Vssim + noise_perc * std(Vssim) * randn(size(Vssim));
Rhosim = Rhosim + noise_perc * std(Rhosim) * randn(size(Rhosim));

% training dataset
mtrain = [Phi_sim Clay_sim Sw_sim];
nv = size(mtrain,2);
%dtrain = [Vprpm Vsrpm Rhorpm];
dtrain = [Vpsim Vssim Rhosim];
nd = size(dtrain,2);
nf = max(unique(Facies));

% domain to evaluate the posterior PDF
phidomain = (0:0.005:0.4);   
cdomain = (0:0.01:0.8); 
swdomain = (0:0.01:1);
[P,V,S] = ndgrid(phidomain, cdomain, swdomain);
mdomain = [P(:) V(:) S(:)];

% measured data (elastic logs)
dcond = [Vp Vs Rho];
ns = size(dcond,1);

generate_histograms([mtrain dtrain])


%% Non-parametric case (Kernel density estimation) - Traditional approach
% % petrophysical domain discretization
% ndiscr = 25;
% phidomain = linspace(0, 0.4, ndiscr)';   
% cdomain = linspace(0, 0.8, ndiscr)';   
% swdomain = linspace(0, 1, ndiscr)';   
% mdomain = [phidomain cdomain swdomain];
% % elastic domain discretization
% vpdomain  = linspace(min(Vp), max(Vp),ndiscr)';
% vsdomain = linspace(min(Vs), max(Vs),ndiscr)';
% rhodomain = linspace(min(Rho), max(Rho),ndiscr)';
% ddomain =[vpdomain vsdomain rhodomain];
% % kernel bandwidths 
% h = 5;
% hm(1) = (max(phidomain)-min(phidomain))/h;
% hm(2) = (max(cdomain)-min(cdomain))/h;
% hm(3) = (max(swdomain)-min(swdomain))/h;
% hd(1) = (max(vpdomain)-min(vpdomain))/h;
% hd(2) = (max(vsdomain)-min(vsdomain))/h;
% hd(3) = (max(rhodomain)-min(rhodomain))/h;
% 
% % inversion
% tic
% Ppost = RockPhysicsKDEInversion(mtrain, dtrain, mdomain, ddomain, dcond, hm, hd);
% toc
% 
% % marginal posterior distributions
% Ppostphi = zeros(ns,length(phidomain));
% Ppostclay = zeros(ns,length(cdomain));
% Ppostsw = zeros(ns,length(swdomain));
% Phimap = zeros(ns,1);
% Cmap = zeros(ns,1);
% Swmap = zeros(ns,1);
% for i=1:ns
%     Ppostjoint=reshape(Ppost(i,:),length(phidomain),length(cdomain),length(swdomain));
%     Ppostphi(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),2);
%     Ppostclay(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),1);
%     Ppostsw(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),2)),1);
%     Ppostphi(i,:)=Ppostphi(i,:)/sum(Ppostphi(i,:));
%     Ppostclay(i,:)=Ppostclay(i,:)/sum(Ppostclay(i,:));
%     Ppostsw(i,:)=Ppostsw(i,:)/sum(Ppostsw(i,:));
%     [~,Phimapind]=max(Ppostphi(i,:));
%     [~,Cmapind]=max(Ppostclay(i,:));
%     [~,Swmapind]=max(Ppostsw(i,:));
%     Phimap(i)=phidomain(Phimapind);
%     Cmap(i)=cdomain(Cmapind);
%     Swmap(i)=swdomain(Swmapind);
% end
% 
% % plots
% figure
% subplot(131)
% %pcolor(phidomain, Depth, Ppostphi); 
% hold on; % shading interp; colorbar; 
% plot(Phi, Depth, 'k', 'LineWidth', 2);  
% ylabel('Depth (m)'); xlabel('Porosity (v/v)');  xlim([0 0.4]);
% plot(Phimap, Depth, 'r', 'LineWidth', 2);
% set(gca,'Ydir','reverse')
% grid
% subplot(132)
% %pcolor(cdomain, Depth, Ppostclay); 
% hold on; % shading interp; colorbar; 
% plot(Clay, Depth, 'k', 'LineWidth', 2); 
% xlabel('Clay volume (v/v)'); xlim([0 0.8]);
% plot(Cmap, Depth, 'r', 'LineWidth', 2);
% set(gca,'Ydir','reverse')
% grid
% subplot(133)
% %pcolor(swdomain, Depth, Ppostsw); 
% hold on; % shading interp; colorbar; 
% plot(Sw, Depth, 'k', 'LineWidth', 2); 
% plot(Swmap, Depth, 'r', 'LineWidth', 2);
% xlabel('Water saturation (v/v)');  xlim([0 1]);
% %hbc=colorbar; title(hbc, 'Probability');
% set(gca,'Ydir','reverse')
% grid
% sgtitle('Traditional Implementation');


%% Non-parametric case  - DMS approach - Simulating order Phi, vchaly, sw
I = length(Depth);
J = 1;

reference_variables = [dtrain mtrain];
[reference_variables] = extend_dateset_KDE(reference_variables,5,0.05);

% Median model:
tic
[DME_result_median] = DMS(I, J, 30, 'gau', [0 0 0], 0.06, reference_variables, [], [], 1, dcond );
toc
DME_result_median = cat(4, DME_result_median{:});
Phi_median_uncond = squeeze(DME_result_median(4,:,:,:));
Vc_median_uncond = squeeze(DME_result_median(5,:,:,:));
Sw_median_uncond = squeeze(DME_result_median(6,:,:,:));

% % Simulaitons:
n_sims = 50;
tic
[DME_result_sim] = DMS(I, J, 30, 'gau', [0 0 0], 0.06, reference_variables, [], [], n_sims, dcond );
toc
DME_result_sim = cat(4, DME_result_sim{:});
Phi_sim_uncond = squeeze(DME_result_sim(4,:,:,:));
Vc_sim_uncond = squeeze(DME_result_sim(5,:,:,:));
Sw_sim_uncond = squeeze(DME_result_sim(6,:,:,:));

% plots
figure
subplot(171)
imagesc(1,Depth,Facies3)
colormap(cmap_3facies)
ylabel('Depth (m)');

% ELASTIC PROPERTIES
subplot(172)
plot(Vp, Depth, 'k', 'LineWidth', 2);  
set(gca,'Ydir','reverse')
grid
xlabel('P-wave velocity (km/s)'); 
subplot(173)
plot(Vs, Depth, 'k', 'LineWidth', 2);  
set(gca,'Ydir','reverse')
grid
xlabel('S-wave velocity (km/s)'); 
subplot(174)
plot(Rho, Depth, 'k', 'LineWidth', 2);  
set(gca,'Ydir','reverse')
xlabel('Density (g/cm^3)'); 
grid

% INVERSION 
subplot(175)
plot(Phi_sim_uncond,Depth,'Color', [0.6, 0.6, 0.6])
hold on; %shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
plot(Phi_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Porosity (v/v)');  xlim([0 0.4]);
set(gca,'Ydir','reverse')
grid
subplot(176)
plot(Clay, Depth, 'k', 'LineWidth', 2); 
hold on; %shading interp; colorbar; 
plot(Vc_median_uncond, Depth, 'r', 'LineWidth', 2);
plot(Vc_sim_uncond,Depth,'Color', [0.7, 0.7, 0.7])
plot(Clay, Depth, 'k', 'LineWidth', 2); 
plot(Vc_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
set(gca,'Ydir','reverse')
grid
legend('Reference','Median','Simulations')
subplot(177)
plot(Sw_sim_uncond,Depth,'Color', [0.7, 0.7, 0.7])
hold on; %shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Sw_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
set(gca,'Ydir','reverse')
grid
% subplot(144)
% imagesc(1,Depth, Facies3)
% %pcolor(swdomain, Depth, Ppostsw); 
% set(gca,'Ydir','reverse')
% grid
%hbc=colorbar; title(hbc, 'Probability');
sgtitle('DMS Approach')


%% Non-parametric case  - DMS approach - Simulating order sw, vchaly, Phi
I = length(Depth);
J = 1;

mtrain = [Sw_sim Clay_sim Phi_sim];
reference_variables = [dtrain mtrain];
[reference_variables] = extend_dateset_KDE(reference_variables,5,0.05);

% Median model:
tic
[DME_result_median] = DMS(I, J, 30, 'gau', [0 0 0], 0.06, reference_variables, [], [], 1, dcond );
toc
DME_result_median = cat(4, DME_result_median{:});
Phi_median_uncond = squeeze(DME_result_median(6,:,:,:));
Vc_median_uncond = squeeze(DME_result_median(5,:,:,:));
Sw_median_uncond = squeeze(DME_result_median(4,:,:,:));

% % Simulaitons:
n_sims = 50;
tic
[DME_result_sim] = DMS(I, J, 30, 'gau', [0 0 0], 0.06, reference_variables, [], [], n_sims, dcond );
toc
DME_result_sim = cat(4, DME_result_sim{:});
Phi_sim_uncond = squeeze(DME_result_sim(6,:,:,:));
Vc_sim_uncond = squeeze(DME_result_sim(5,:,:,:));
Sw_sim_uncond = squeeze(DME_result_sim(4,:,:,:));

% plots
figure
% ELASTIC PROPERTIES
subplot(161)
plot(Vp, Depth, 'k', 'LineWidth', 2);  
ylabel('Depth (m)');
set(gca,'Ydir','reverse')
grid
xlabel('P-wave velocity (km/s)'); 
subplot(162)
plot(Vs, Depth, 'k', 'LineWidth', 2);  
set(gca,'Ydir','reverse')
grid
xlabel('S-wave velocity (km/s)'); 
subplot(163)
plot(Rho, Depth, 'k', 'LineWidth', 2);  
set(gca,'Ydir','reverse')
xlabel('Density (g/cm^3)'); 
grid

% INVERSION 
subplot(164)
plot(Phi_sim_uncond,Depth,'Color', [0.6, 0.6, 0.6])
hold on; %shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
plot(Phi_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Porosity (v/v)');  xlim([0 0.4]);
set(gca,'Ydir','reverse')
grid
subplot(165)
plot(Clay, Depth, 'k', 'LineWidth', 2); 
hold on; %shading interp; colorbar; 
plot(Vc_median_uncond, Depth, 'r', 'LineWidth', 2);
plot(Vc_sim_uncond,Depth,'Color', [0.7, 0.7, 0.7])
plot(Clay, Depth, 'k', 'LineWidth', 2); 
plot(Vc_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
set(gca,'Ydir','reverse')
grid
legend('Reference','Median','Simulations')
subplot(166)
plot(Sw_sim_uncond,Depth,'Color', [0.7, 0.7, 0.7])
hold on; %shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Sw_median_uncond, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
set(gca,'Ydir','reverse')
grid
% subplot(144)
% imagesc(1,Depth, Facies3)
% %pcolor(swdomain, Depth, Ppostsw); 
% set(gca,'Ydir','reverse')
% grid
%hbc=colorbar; title(hbc, 'Probability');
sgtitle('DMS Approach')








