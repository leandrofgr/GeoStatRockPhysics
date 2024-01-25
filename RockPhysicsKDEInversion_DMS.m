function [ mean_marginal, conditioned_data ] = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, h, n_posterior_pts)

% ROCK PHYSICS KDE INVERSION computes the posterior distribution of
% petrophysical properties conditioned on elastic properties assuming a
% non-parametric distribution.
% The joint distribution of the Bayesian inversion approach is estimated
% from a training dataset using Kernel Density Estimation
% INPUT mtrain = training dataset of petrophysical properties (ntrain, nm)
%       dtrain = training dataset of elastic properties (ntrain, nd)
%       dcond = measured data (nsamples, nd)
%       h = Number of the "discretized grid" of the distribution, in
%       practce is used as a tolerance for the conditioning the list of
%       data.
% OUTUPT Ppost = joint posterior distribution

% Written by Leandro Passos de Figueiredo(June 2023)

grid_size = 1/h;

reference_variables = [dtrain mtrain];

number_conditional = size(dcond,2);
number_conditioned =  size(reference_variables,2) - size(dcond,2);

min2norm = min(reference_variables);

reference_variables = reference_variables - repmat(min2norm, size(reference_variables,1), 1);
max2norm = max(reference_variables);
reference_variables = reference_variables ./ repmat(max2norm, size(reference_variables,1), 1);

dcond = dcond - repmat(min2norm(1:number_conditional), size(dcond,1), 1);
dcond = dcond ./ repmat(max2norm(1:number_conditional), size(dcond,1), 1);



mean_marginal = nan*zeros( size(dcond,1), number_conditioned );
conditioned_data =[];
if nargin > 4
    conditioned_data = zeros( size(dcond,1), number_conditioned, n_posterior_pts);
end
num_point_without_statistic = 0;
for i = 1:1:size(dcond,1)
    
    if ~isnan(dcond(i,1))
        
        reference_variables_filtered = reference_variables;        
        
        for j =1:1:number_conditional
            
            index = abs(reference_variables_filtered(:,j) - dcond(i,j)) < grid_size; %filtra apenas valores em torno do valor amostrado
            reference_variables_filtered(~index,:) = [];
            
        end
        if size(reference_variables_filtered,1) > 0
            for j=1:number_conditioned
                mean_marginal(i,j) = mean(reference_variables_filtered(:,number_conditional+j));
                if nargin > 4
                    sorting = randperm( size(reference_variables_filtered,1), min([ size(reference_variables_filtered,1) n_posterior_pts]) ) ;
                    conditioned_data(i,j,1:length(sorting)) = reference_variables_filtered(sorting,number_conditional+j);
                end
            end
        else
            
            num_point_without_statistic = num_point_without_statistic + 1;
            disp('Not enough data for contitioning for ' + string(num_point_without_statistic) + ' The method will draw from the marginal. It might generate artifacts. Consider using KDE to increase the number of data points or increasing the grid_size parameter.')
            mean_marginal(i,1:number_conditioned) = mean(reference_variables(:,number_conditional+1:end));
            %conditioned_data(i,j,:)
            
        end
       
    end
    
    
end


mean_marginal = mean_marginal.* repmat(max2norm(number_conditional+1:end), size(mean_marginal,1), 1);
mean_marginal = mean_marginal + repmat(min2norm(number_conditional+1:end), size(mean_marginal,1), 1);

if nargin > 4
    indices_zeros = find(conditioned_data==0);
    conditioned_data = conditioned_data .* repmat(max2norm(number_conditional+1:end), size(conditioned_data,1), 1);
    conditioned_data = conditioned_data + repmat(min2norm(number_conditional+1:end), size(conditioned_data,1), 1);
    conditioned_data(indices_zeros) = 0;
end


end