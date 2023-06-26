function [ conditioned_data ] = RockPhysicsKDEInversion_DMS(mtrain, dtrain, dcond, h)                                                          

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




for i = 1:1:size(dcond,1) %para cada ponto repete;
    reference_variables_filtered = reference_variables;
    for j =1:1:number_conditional %para cada variavel repete
               
        index = abs(reference_variables_filtered(:,j) - dcond(i,j)) < grid_size; %filtra apenas valores em torno do valor amostrado                        
        reference_variables_filtered(~index,:) = [];                

    end
    
    for j=1:number_conditioned
        conditioned_data(i,j) = mean(reference_variables_filtered(:,number_conditional+j));
    end
    
end

conditioned_data = conditioned_data.* repmat(max2norm(number_conditional+1:end), size(conditioned_data,1), 1);
conditioned_data = conditioned_data + repmat(min2norm(number_conditional+1:end), size(conditioned_data,1), 1);


end