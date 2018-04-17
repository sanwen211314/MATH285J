clear; tStart = tic;

%matrix size and rank
m = 20; n = 20; r = 3;

% correlation parameters
% this will make the first 't' number of columns pairwise correlated 
% with expected correlation rho
t = round(n/2);
rho = 0.6;
corrMat = rho*ones(t,t) + (1-rho)*eye(t,t);
U = chol(corrMat);

% observation rate
observationRate = 0.5;

% run the matrix completion problem 'numTrials' number of times
% At the end, we will calculate the 'averageError' in the solution
% for each trial and we will 'count' the number of times we calculate
% a solution with 'tol' of the exact solution.
count = 0;
tol = 1e-5;
averageError = 0;
numTrials = 1000;

%proportion of observations from first 't' rows:
%   If a = 0.5 and t = n/2, we are just sampling uniformly
%   so we are back in the 'non-adaptive' case
a = 0.5;

for trial = 1:numTrials
    %Generate a random, relatively sparse matix of the desired rank
    % then correlate the first half of the rows
    A = randomSparseMat(m,n,r,0.3,0.3);
    A(:,1:t) = A(:,1:t)*U;
    
    %   Generate a set of random observation indices
    indices1 = randomObservationIndices(m,t,2*observationRate*a);
    indices2 = randomObservationIndices(m,n-t,2*observationRate*(1-a));
    indices = [indices1,indices2];
    
    cvx_begin quiet
        variable X(m,n)
        %minimize the nuclear norm
        minimize( norm_nuc(X) )
        subject to
            %require the solution to match at the observed indices
            X(indices==1) == A(indices==1)
    cvx_end
    
    diff = norm(A - X,2);
    if(diff < tol)
       count = count + 1; 
    end
    averageError = averageError + diff;
end
averageError = averageError/numTrials;

fprintf('Proportion of observations from first group: %.2f\n',a);
fprintf('Number of solutions within %.0e: %i\n',tol,count);
fprintf('Average error: %.4f\n',averageError);
tElapsed = toc(tStart);
fprintf('Time elapsed: %.2f sec\n',tElapsed);