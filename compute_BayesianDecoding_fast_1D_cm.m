function [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast_1D_cm(numSpks,spatialTuningCurve,numSpatialBins,windowSizeDecoding,return_fullPosterior,directionalDecoding)

    numTimeBins = size(numSpks,2);
    
    if directionalDecoding==1
        numSpatialBins = [numSpatialBins(1) 2*numSpatialBins(2)];
    end

    if return_fullPosterior == 1
        posterior = NaN(numTimeBins,prod(numSpatialBins));
    else
        posterior = NaN;
    end

    %emtpy storage matrices
    if directionalDecoding == 1
        posteriorCOM = NaN(numTimeBins,2);
        posteriorSpread = NaN(numTimeBins,2);
        posteriorPeak = NaN(numTimeBins,2);
    else
        posteriorCOM = NaN(numTimeBins,1);
        posteriorSpread = NaN(numTimeBins,1);
        posteriorPeak = NaN(numTimeBins,1);
    end

    %decoding
        factor_exp = exp(-windowSizeDecoding*sum(spatialTuningCurve));
        for k = 1:numTimeBins
            
            if mod(k,round(numTimeBins/10))==0
                disp(strcat('elapsed portion: ',num2str(k/numTimeBins)))
            end

            factor_prod = prod(bsxfun(@power,spatialTuningCurve,numSpks(:,k)));
                        
            likelihood_sub = bsxfun(@times,factor_prod,factor_exp)';
            likelihood_sub(isinf(likelihood_sub)) = 0;
            posterior_sub = likelihood_sub/nansum(likelihood_sub(:));
            posterior_sub_mat = reshape(posterior_sub,numSpatialBins);
            
            if return_fullPosterior == 1
                posterior(k,:) = posterior_sub;
            end
            
            if directionalDecoding == 1
                 
                posterior_sub_mat_left = posterior_sub_mat(1:numSpatialBins(2)/2);
                posterior_sub_mat_right = posterior_sub_mat(numSpatialBins(2)/2 + 1:end);

                
                posteriorCOM_sub = compute_matLocToPlot(compute_centerOfMass(posterior_sub_mat_left));
                posteriorCOM(k,1) = posteriorCOM_sub(1);
                posteriorCOM_sub = compute_matLocToPlot(compute_centerOfMass(posterior_sub_mat_right));
                posteriorCOM(k,2) = posteriorCOM_sub(1);


               posteriorSpread(k,1) = compute_imageSpread(posterior_sub_mat_left,1);
               posteriorSpread(k,2) = compute_imageSpread(posterior_sub_mat_right,1);
                
                posteriorPeak(k,1) = max(posterior_sub_mat_left(:));
                posteriorPeak(k,2) = max(posterior_sub_mat_right(:));
            else
                posteriorCOM_sub = compute_matLocToPlot(compute_centerOfMass(posterior_sub_mat));
                posteriorCOM(k) = posteriorCOM_sub(1);

                posteriorSpread(k) = compute_imageSpread(posterior_sub_mat);
                posteriorPeak(k) = max(posterior_sub(:));
            end
    
                                    
            %plot
%             plot(posterior_sub_mat(1:numSpatialBins(1)/2),'b')
%             hold on, plot(posterior_sub_mat(numSpatialBins(1)/2+1:end),'r'), hold off
%             vline(posteriorCOM(k,1),'b')
%             vline(posteriorCOM(k,2),'r')
%             ylim([0 0.5])
%             title(num2str(compute_round(posteriorSpread(k,1),100)))
%             drawnow
%             keyboard
        end