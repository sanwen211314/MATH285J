clear;

%matrix size and rank
m = 20; n = 20; r = 3;
%M = randomSparseMat(m,n,r);
M = rand(m,n);

%observation rate
observationRate = 0.3;

% tStart = tic;
count = 0;
tol = 1e-7;
averageError = 0;
numTrials = 10;
MaxDiff = 0;
MinDiff = Inf;
for trial = 1:numTrials

    %   Generate a set of random observation indices
    indices = randomObservationIndices(m,n,observationRate);
    
    cvx_begin quiet
        variable X(m,n)
        minimize( norm_nuc(X) )
        subject to
            X(indices==1) == M(indices==1)
    cvx_end
    
    diff = norm(M - X,2);
    if(diff < tol)
       count = count + 1; 
    end
    averageError = averageError + diff;
    if(diff > MaxDiff)
        Xworst = X;
        MaxDiff = diff;
        IndWorst = indices;
    end
    if(diff < MinDiff)
        Xbest = X;
        MinDiff = diff;
        IndBest = indices;
    end
end
averageError = averageError/(numTrials);

fprintf('Number of solutions within %.0e: %i\n',tol,count);
fprintf('Average absolute error: %.4f\n',averageError);
%tElapsed = toc(tStart);
%fprintf('Time elapsed: %.2f sec\n',tElapsed);

a = (1:20) - 0.5; b= a;
hFig = figure(110);clf;
set(hFig,'Position',[50 50 1300 700]);
subplot(2,3,1);
surf(a,b,M); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
caxis([0,max(max(M))]);
title('Original Matrix')
subplot(2,3,2);
surf(a,b,Xbest); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
str = sprintf('Best Reconstruction, ||X - M|| = %.2e',MinDiff);
title(str)
caxis([0,max(max(M))]);
subplot(2,3,3);
surf(a,b,IndBest);view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
title('Observation Indices (yellow)');
subplot(2,3,4);
surf(a,b,M); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
caxis([0,max(max(M))]);
title('Original Matrix')
subplot(2,3,5);
surf(a,b,Xworst); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
str = sprintf('Worst Reconstruction, ||X - M|| = %.2e',MaxDiff);
title(str)
caxis([0,max(max(M))]);
subplot(2,3,6);
surf(a,b,IndWorst);view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
title('Observation Indices (yellow)');