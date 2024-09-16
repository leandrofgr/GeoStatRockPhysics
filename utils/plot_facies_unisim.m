
close all 
%% PLOTS FACIES AND WELL LOGS
load('.\data\UNISIM\well_logs_table.mat')
phi_well = well_logs_table(:,2);
sw_well = well_logs_table(:,5);
phi_well(sw_well==1) = phi_well(sw_well==1)*0.12;
sw_well(sw_well>=0.9) = 0.9;

sw_axis = linspace(0,1,100);
phi_axis = linspace(0,0.35,100);
[X,Y] = meshgrid(phi_axis, sw_axis );

% FACIES 1:
for i=1:size(Y,1)
    pdf_1(i,:) = normpdf(X(i,:),0.1,0.05);
end

% FACIES 2:
for i=1:size(Y,1)
    pdf_2(i,:) = normpdf(X(i,:),0.25,0.05);
end

% FACIES 3:
pdf_3 = mvnpdf([X(:) Y(:)],[0.02 0.9],[0.01 0; 0 0.02].^2);
pdf_3 = reshape(pdf_3,size(X));


figure
subplot(222)
scatter(phi_well, sw_well, 15, 'k','filled')
grid
hold all
contour(X, Y, pdf_1,4,'LineColor',[0.8 0.8 0],'LineWidth',1.5);
contour(X, Y, pdf_2,4,'LineColor',[0.0 0.7 0.4],'LineWidth',1.5);
contour(X, Y, pdf_3,4,'LineColor',[0.5 0.5 0.5],'LineWidth',1.5);
legend('Well data','Mid-porosiy sand','High-porosiy sand','Shale')

subplot(224)
histogram(phi_well,phi_axis,'Normalization','pdf')
hold all
plot(phi_axis,pdf_1(1,:),'Color',[0.8 0.8 0],'LineWidth',1.5)
plot(phi_axis,pdf_2(1,:),'Color',[0.0 0.7 0.4],'LineWidth',1.5)
plot(phi_axis,0.05*pdf_3(90,:),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
grid
xlim([0 0.35])
ylim([0 40])
xlabel('Porosity (v/v)')

uniform_pdf = 5*ones(size(sw_axis));
subplot(221)
histogram(sw_well,sw_axis,'Normalization','pdf')
hold all
plot(sw_axis,uniform_pdf ,'Color',[0.8 0.8 0],'LineWidth',1.5)
plot(sw_axis,uniform_pdf ,'Color',[0.0 0.7 0.4],'LineWidth',1.5)
plot(sw_axis,0.03*pdf_3(:,7),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
grid
xlim([0 1])
ylim([0 25])
ylabel('S_w (v/v)')