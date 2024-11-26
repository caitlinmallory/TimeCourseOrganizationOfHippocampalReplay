function [posterior_color,posterior_color_binarized,posterior_sum,posterior_sum_binarized] = compute_flattenedDecoding(posterior,spatialDim,colormap_num)   



[numDecodedTimePts,numBins] = size(posterior);
if colormap_num == 1
    colormap_choice = cbrewer2('Blues',numDecodedTimePts);
else
    colormap_choice = cbrewer2('Oranges',numDecodedTimePts);
end
% colormap_choice = cool(numDecodedTimePts);
posterior_sum = squeeze(nansum(posterior,1));

posterior_color = [];
posteriorThr = 0.01; % This is what I had this set ot originally.

posteriorThr = 20*1/numBins;

posterior_sum_binarized = heaviside(posterior_sum-posteriorThr);

if spatialDim==1
    %posterior_color
        %convert to RGB (avg RGBs from two directions)
        posterior_direction_1 = posterior(1:numDecodedTimePts/2,:);
        posterior_direction_2 = posterior(numDecodedTimePts/2+1:end,:);
        
        posterior_color(:,:,1) = (1-posterior_direction_1 + ones(size(posterior_direction_1)))/2;
        posterior_color(:,:,2) = (1-posterior_direction_1 + 1-posterior_direction_2)/2;
        posterior_color(:,:,3) = (ones(size(posterior_direction_1)) + 1-posterior_direction_2)/2;

        %change saturation
        HSV = rgb2hsv(posterior_color);
        HSV(:,:,2) = HSV(:,:,2)*4;
        posterior_color = hsv2rgb(HSV);
    
    %posterior_binarized
        posterior_color_binarized = heaviside(posterior-posteriorThr);
else
    C1 = 0;
    C2 = 0;
    a = 1-colormap_choice; 
    c1 = nanmean((posterior).*repmat(a(:,1),1,numBins),1);
    c2 = nanmean((posterior).*repmat(a(:,2),1,numBins),1);
    c3 = nanmean((posterior).*repmat(a(:,3),1,numBins),1);
    C1 = [c1; c2; c3];

    %convert black pixels to white
    C1 = 1-500*C1;
    posterior_color = C1;
    
%     %change saturation
%     posterior_color(posterior_color>1) = 1;
%     HSV = rgb2hsv(posterior_color');
%     HSV(:,2) = 2*HSV(:,2); HSV(HSV(:,2)>1,2) = 1;
%     posterior_color = hsv2rgb(HSV)';
    

    decodedData_binary = heaviside(posterior-posteriorThr);    
    C2 = zeros(3,numBins);
    a = colormap_choice; 
    for j = 1:numDecodedTimePts
        C_new = zeros(3,numBins);
        for k = 1:numBins
            C_new(:,k) = a(j,:)*decodedData_binary(j,k);
            if nansum(C2(:,k))==0
                C2(:,k) = C_new(:,k);
            elseif nansum(C2(:,k))~=0 && nansum(C_new(:,k))~=0
                C2(:,k) = (C2(:,k)+C_new(:,k))/2;
            end
        end
    end

    %convert black pixels to white
    C2(:,nansum(C2,1)==0) = 1;
    posterior_color_binarized = C2;
end






% [numDecodedTimePts,numBins] = size(posterior);
% 
% posterior_sum = squeeze(nansum(posterior,1));
% 
% posterior_color = [];
% posterior_color_binarized = [];
% posteriorThr = 0.01;
% posterior_sum_binarized = heaviside(posterior_sum-posteriorThr);
% if spatialDim==1
%     if strcmp(type,'both') || strcmp(type,'normal')
%     %posterior_color
%         %convert to RGB (avg RGBs from two directions)
%         posterior_direction_1 = posterior(1:numDecodedTimePts/2,:);
%         posterior_direction_2 = posterior(numDecodedTimePts/2+1:end,:);
%         posterior_color(:,:,1) = (1-posterior_direction_1 + ones(size(posterior_direction_1)))/2;
%         posterior_color(:,:,2) = (1-posterior_direction_1 + 1-posterior_direction_2)/2;
%         posterior_color(:,:,3) = (ones(size(posterior_direction_1)) + 1-posterior_direction_2)/2;
% 
%         %change saturation
%         HSV = rgb2hsv(posterior_color);
%         HSV(:,:,2) = HSV(:,:,2)*4;
%         posterior_color = hsv2rgb(HSV);
%     end
%     
%     if strcmp(type,'both') || strcmp(type,'binarized')
%     %posterior_binarized
%         posterior_color_binarized = heaviside(posterior-posteriorThr);
%     end
%     
% else
%     if strcmp(type,'both') || strcmp(type,'normal')
%         C1 = 0;
%         C2 = 0;
%         a = 1-colormap_choice(numDecodedTimePts); 
%         c1 = nanmean((posterior).*repmat(a(:,1),1,numBins),1);
%         c2 = nanmean((posterior).*repmat(a(:,2),1,numBins),1);
%         c3 = nanmean((posterior).*repmat(a(:,3),1,numBins),1);
%         C1 = [c1; c2; c3];
% 
%         %convert black pixels to white
%         C1 = 1-80*C1;
% 
%         posterior_color = C1;
%     end
%     %%%%
%     if strcmp(type,'both') || strcmp(type,'binarized')
%         decodedData_binary = heaviside(posterior-posteriorThr);    
%         C2 = zeros(3,numBins);
%         a = colormap_choice(numDecodedTimePts); 
%         for j = 1:numDecodedTimePts
%             C_new = zeros(3,numBins);
%             for k = 1:numBins
%                 C_new(:,k) = a(j,:)*decodedData_binary(j,k);
%                 if nansum(C2(:,k))==0
%                     C2(:,k) = C_new(:,k);
%                 elseif nansum(C2(:,k))~=0 && nansum(C_new(:,k))~=0
%                     C2(:,k) = (C2(:,k)+C_new(:,k))/2;
%                 end
%             end
%         end
% 
%         %convert black pixels to white
%         C2(:,nansum(C2,1)==0) = 1;
%         posterior_color_binarized = C2;
%         
% %         posterior_color_binarized(isnan(posterior_color_binarized(:))) = 1;
% 
%     end
% end
% 
