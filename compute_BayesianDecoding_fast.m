function [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap,numSpatialBins,spatialDim,windowSizeDecoding,return_fullPosterior)

    numTimeBins = size(numSpks,2);
    numClusters = size(numSpks,1);

    if return_fullPosterior == 1
        posterior = NaN(numTimeBins,prod(numSpatialBins));
    else
        posterior = NaN;
    end

    %emtpy storage matrices
        posteriorCOM = NaN(numTimeBins,spatialDim);
        posteriorSpread = NaN(numTimeBins,1);
        posteriorPeak = NaN(numTimeBins,1);
        
    %add baseline to tuning curves
        rateMap = rateMap + 1e-3;
        
    %decoding
        factor_exp = exp(-windowSizeDecoding*nansum(rateMap));
        
        %compute f(x)^k beforehand for each cell across all possible k's 
            Factor_prod_library = cell(1,numClusters);
            K = max(numSpks(:)); 
%             if K>50
%                 keyboard
%             end
            for j = 1:numClusters
                factor_prod_library = zeros(K,prod(numSpatialBins));
                for k = 1:K
                    factor_prod_library(k,:) = rateMap(j,:).^k;
                end
                Factor_prod_library{j} = factor_prod_library;
            end
            
        for k = 1:numTimeBins
            
%             if mod(k,round(numTimeBins/4))==0
%                 disp(strcat('elapsed portion: ',num2str(k/numTimeBins)))
%             end

        %new fast version
            ind_clusters_nonzero = find(numSpks(:,k)>0);
            factor_prod = ones(length(ind_clusters_nonzero),prod(numSpatialBins));
            for j = 1:length(ind_clusters_nonzero)
                factor_prod(j,:) = Factor_prod_library{ind_clusters_nonzero(j)}(numSpks(ind_clusters_nonzero(j),k),:);
            end  
            factor_prod = prod(factor_prod);
            
            likelihood_sub = bsxfun(@times,factor_prod,factor_exp)';
                likelihood_sub(isinf(likelihood_sub)) = 0;
            posterior_sub = likelihood_sub/nansum(likelihood_sub(:));
            posterior_sub_mat = reshape(posterior_sub,numSpatialBins);

        %old version
%             ind = find(numSpks(:,k)>0);
%             factor_exp = exp(-windowSizeDecoding*nansum(rateMap));
%             factor_prod = prod(bsxfun(@power,rateMap(ind,:),numSpks(ind,k)));
%             likelihood_sub = bsxfun(@times,factor_prod,factor_exp)';
%                 likelihood_sub(isinf(likelihood_sub)) = 0;
%             posterior_sub = likelihood_sub/nansum(likelihood_sub(:));
%             posterior_sub_mat = reshape(posterior_sub,numSpatialBins);
            

            if return_fullPosterior == 1
                posterior(k,:) = posterior_sub;
            end
            posteriorCOM(k,:) = compute_matLocToPlot(compute_centerOfMass(posterior_sub_mat));
            posteriorSpread(k) = compute_imageSpread(posterior_sub_mat,2);
            posteriorPeak(k) = max(posterior_sub(:));
                        
            %plot
%             imagesc(posterior_sub_mat,'alphadata',~isnan(posterior_sub_mat)), set(gca,'ydir','normal')
%             hold on, plot(posteriorCOM(k,1),posteriorCOM(k,2),'r*'), hold off
%             title(num2str(round(posteriorSpread(k)*100)/100))
%             drawnow
            
        end