clear; close all; tStart = tic;

%matrix size and correlation parameter
m = 20; n = 20;
t = round(n/2);

%create random matrix, make first t columns scalar multiples of a single
%vector
M = rand(m,n);
M(:,1:t) = (3/2)*rand(m,1)*rand(1,t);
Mbar = M(:,1:t) + (0.5*rand(m,t)-0.25);

% what proportion of entries do we observe?
observationRate = 0.2;
a = 0;
gamma = 1;

% run the matrix completion problem 'numTrials' number of times
% At the end, we will calculate the 'averageError' in the solution
% for each trial and we will 'count' the number of times we calculate
% a solution with 'tol' of the exact solution.
count = 0;
tol = 1e-5;
averageError = 0;
numTrials = 10;
minDiff = Inf;

for trial = 1:numTrials
    %   Generate a set of random observation indices
    indices1 = randomObservationIndices(m,t,(n/t)*observationRate*a);
    indices2 = randomObservationIndices(m,n-t,(n/(n-t))*observationRate*(1-a));
    indices = [indices1,indices2];
    
    cvx_begin quiet
        variable X(m,n)
        minimize( norm_nuc(X)+gamma*sum(sum((Mbar-X(:,1:10)).^2)))
        subject to
            %require the solution to match at the observed indices
            X(indices==1) == M(indices==1)
    cvx_end
    
    diff = norm(M - X,2);
    if(diff < tol)
       count = count + 1; 
    end
    if(diff < minDiff)
       minDiff = diff;
       Xbest = X;
       IndBest = indices;
    end
    averageError = averageError + diff;
end
averageError = averageError/numTrials;

fprintf('Proportion of observations from correlated columns: %.2f\n',a);
fprintf('Average relative error: %.4f\n',averageError);
tElapsed = toc(tStart);
fprintf('Time elapsed: %.2f sec\n',tElapsed);


a = (1:20) - 0.5; b = a;
hFig = figure(110);clf;
set(hFig,'Position',[0 350 1600 400]);
subplot(1,3,1);
surf(a,b,M); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
caxis([0,max(max(M))]);
title('Original Matrix')
subplot(1,3,2);
surf(a,b,Xbest); view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
str = sprintf('Best Reconstruction, ||X - M|| = %.2e',minDiff);
title(str)
caxis([0,max(max(M))]);
subplot(1,3,3);
surf(a,b,indices);view([0,90]);
axis([0.5 19.5 0.5 19.5]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
title('Observation Indices (yellow)');