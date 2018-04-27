clear; close all; tStart = tic;

%matrix size and correlation parameter
m = 100; n = 100;
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
numTrials = 1;
minDiff = Inf;

for trial = 1:numTrials
    %   Generate a set of ordered, observation indices for the correlated
    % columns.
    indices1 = zeros(m,t);
    for numRows=1:m
        try
            indices1(1:numRows,:) = 1;
            indvar = [ones(size(M(1:numRows,1))) M(1:numRows,1)];
            b = mvregress(indvar, M(1:numRows,2:t));
            disp("The number of samples needed to obtain a positive-definite covariance:");
            numRows * t
            constraintBuilder = "";
            for row=2:t
                constraintBuilder = constraintBuilder + sprintf("X(:,%d) == X(:,1) .* b(2,%d) + b(1,%d)\n", row,row-1,row-1);
            end
            constraintBuilder
            break
        catch
        end
    end
    
    %   Generate a set of random observation indices
    disp("Proportion of observations from non-correlated columns")
    (observationRate-(numRows*t/(m*n)))/observationRate
    indices2 = randomObservationIndices(m,n-t,max(0, observationRate-(numRows*t/(m*n))));
    indices = [indices1,indices2];
    
    cvx_begin quiet
        variable X(m,n)
        minimize(norm_nuc(X))
        subject to
         
            %require the solution to match at the observed indices
            X(indices==1) == M(indices==1)
            
            % require the solution to satisfy the following constraints
            eval(constraintBuilder)
            
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

fprintf('Proportion of observations from correlated columns: %.2f\n',1 - ((observationRate-(numRows*t/(m*n)))/observationRate));
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