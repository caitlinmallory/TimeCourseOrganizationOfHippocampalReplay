function [spatialInfo] = compute_spatialInfo(X,P)

P_norm = P./sum(P(:));

%mean activity
X_mean = nansum(X(:).*P_norm(:));

%calculate spatial info sum p*r*log(r)
spatialInfo = nansum(P_norm(:).*(X(:)/X_mean).*log2(X(:)/X_mean));

% subplot(121)
% imagesc(X), title(spatialInfo)
% keyboard

% numShuffles = 1000;
% spatialInfo_shuffles = nan(numShuffles,1);
% T = load_timeBounds(Run_Times)-Run_Times(1);
% spkTime_shift = spkTime-Run_Times(1);
% for i = 1:numShuffles
%     
%     spkTime_rand = mod(spkTime_shift + rand(1)*diff(T),T(2))+Run_Times(1);
% %     spkTime_rand = mod(spkTime_shift,T(2))+Run_Times(1);
%     spkPos = compute_dataInterpolation(positions,spkTime_rand,[5 9]);
% %     spkPos(spkPos(:,4)<speedThr) = [];
% 
%     %histogram of spike position data
%     X = hist3([spkPos(:,3) spkPos(:,2)],'Edges',{y_edges x_edges}); X = X(1:end-1,1:end-1);
% 
%     %mean activity
%     X_mean = nanmean(X(:));
% 
%     %calculate spatial info sum p*r*log(r)
%     spatialInfo_shuffles(i) = nansum(P(:).*(X(:)/X_mean/30).*log2(X(:)/X_mean/30));
%     
% %     subplot(122)
% %     imagesc(X), title(spatialInfo_shuffles(i))
% %     keyboard
% end




