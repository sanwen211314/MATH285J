clear; close all;

m=500; n=500;
M = 2*randomSparseMat(m,n,10)-1;

minMSE = Inf;
numTrials = 20;
fprintf('Number of trials: %i\n',numTrials);
fprintf('Trials completed: ');
for l = 1:numTrials
    %X_ij is a normal r.v. with mean M_ij and variance sigma^2
    sigma = 0.5;
    X = sigma*randn(m,n) + M;
    
    %proportion of observed indices
    p = 0.2;
    
    %generate observation indices
    indices = randomObservationIndices(m,n,p);
    
    %set unobserved entries to zero
    X(indices == 0) = 0;
    
    %find SVD of observed part of X
    [U,S,V] = svd(X);
    
    %only keep singular values that are large enough
    q = p*sigma^2 + p*(1-p)*(1-sigma^2);
    threshold = 2*sqrt(q*n);
    S = S.*(S >= threshold);
    
    %construct estimate
    Mhat = (1/p)*U*S*V';
    
    %Mean square error
    MSE = norm(M - Mhat,'fro')^2/(m*n);
    
    if MSE < minMSE
        minMSE = MSE;
        BestMhat = Mhat;
        BestInd = indices;
        BestS = S;
    end
    
    fprintf('%i ',l);
end
fprintf('\n');

a = (1:m) - 0.5; b = (1:n)-0.5;
hFig = figure(110);clf;
set(hFig,'Position',[0 350 1600 400]);
subplot(1,3,1);
surf(a,b,M,'edgecolor','none'); view([0,90]);
axis([min(a) max(a) min(b) max(b)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
caxis([min(min(M)),max(max(M))]);
title('Original Matrix')
subplot(1,3,2);
surf(a,b,BestMhat,'edgecolor','none'); view([0,90]);
axis([min(a) max(a) min(b) max(b)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
str = sprintf('Reconstruction, MSE = %.2e', MSE);
title(str)
caxis([min(min(M)),max(max(M))]);
subplot(1,3,3);
surf(a,b,BestInd,'edgecolor','none');view([0,90]);
axis([min(a) max(a) min(b) max(b)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
title('Observation Indices (yellow)');
