clear; close all; tStart = tic;

%matrix size and correlation parameter
m = 200; n = 200;
t = 10;

%create random matrix, make first t columns scalar multiples of a single
%vector
M = rand(m,n);
coeffVector = rand(1,t);
M(:,1:t) = (3/2)*rand(m,1)*coeffVector;

% what proportion of entries do we observe?
observationRate = 0.2;
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
    if(diff < minDiff)
       minDiff = diff;
       Xbest = X;
       IndBest = indices;
    end
    averageError = averageError + diff;
end
averageError = averageError/numTrials;

fprintf('Average relative error: %.4f\n',averageError);
tElapsed = toc(tStart);
fprintf('Time elapsed: %.2f sec\n',tElapsed);


A = (1:m) - 0.5; B = (1:n)-0.5;
hFig = figure(110);clf;
set(hFig,'Position',[0 350 1600 400]);
subplot(1,3,1);
surf(A,B,M); view([0,90]);
axis([min(A) max(A) min(B) max(B)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
caxis([0,max(max(M))]);
title('Original Matrix')
subplot(1,3,2);
surf(A,B,Xbest); view([0,90]);
axis([min(A) max(A) min(B) max(B)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
str = sprintf('Best Reconstruction, ||X - M|| = %.2e',minDiff);
title(str)
caxis([0,max(max(M))]);
subplot(1,3,3);
surf(A,B,indices);view([0,90]);
axis([min(A) max(A) min(B) max(B)]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
title('Observation Indices (yellow)');