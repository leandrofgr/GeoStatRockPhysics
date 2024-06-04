function [most_likely] = compute_most_likely(simulations)

[I,J,K] = size(simulations);

most_likely = zeros(I,J);

parfor i=1:I
    for j=1:J
        if ~isnan(simulations(i,j))
            simulation = squeeze(simulations(i,j,:));
            [f, xi] = ksdensity(simulation); % Kernel density estimation
            [~, idx] = max(f);
            most_likely(i,j) = xi(idx); % Most likely value for each point
        end
    end
end












% for i=1:I
%     for j=1:J
%         if ~isnan(simulations(i,j))
%             simulation = squeeze(simulations(i,j,:));
%             min_value = min(simulation);
%             max_value = max(simulation);
% 
%             axis = linspace(min_value,max_value,bins);
%             epsilonToAvoidProblemInInterp = linspace(0,1,K)'*0.0001;
%             
%             Prob = linspace(0,1,K);
%             
%             inv_cdf = sort(simulation);
%             
%             cdf = interp1(inv_cdf+epsilonToAvoidProblemInInterp,Prob,axis) ;
%             
%             axis = axis(1:end-1);
%             axis = axis - mean(diff(axis))/2;            
%             
%             pdf = diff(cdf);
%             
%             [~,index_max] = max(pdf );
%             
%             most_likely(i,j) = axis(index_max);
%         end
%     end
% end