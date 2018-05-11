clear; close all; tStart = tic;

%matrix size and correlation parameter
m = 200; n = 200;
t = 50;
NumChosen = [];
constraintBuilder = "";
for num=28:7:42
    %create random matrix, make first t columns scalar multiples of a single
    %vector
%     M = rand(m,n);
%    M(:,1:t) = (1/2)*rand(m,1)*rand(1,t);
%     load('M.mat');
%      save('Mrank1.mat', 'M')
    load('Mrank1.mat', 'M')
%     load('M30.mat');
    % load('M30Rank3.mat');
    k = rank(M(:,1:t));

    % what proportion of entries do we observe?
    observationRate = 0.3;

    %   Generate a set of ordered, observation indices for the correlated
    % columns.
   indices1 = zeros(m,t);
   i = [1:m];
   j = [1:t];
   possiblePairs = nchoosek(1:t, 2);
   usedPairs = [];
   
   for num2=1:num
       if(size(possiblePairs) == 0)
           break
       end
       randi = datasample(i,1);
       randjs = datasample(possiblePairs,1);
       usedPairs = vertcat(usedPairs,randjs);
       newPairs = [];
       for item=1:size(usedPairs)-1
           for index=1:2
               if (ismember(randjs(index), usedPairs(item,:)))
                  
                   if (randjs(index) == usedPairs(item,1))
                        otherElem = usedPairs(item,2);
                   else
                       otherElem = usedPairs(item,1);
                   end
                   if (index == 1)
                        newPairs = vertcat(newPairs,sort([randjs(2),otherElem]));
                   else
                       newPairs = vertcat(newPairs,sort([randjs(1),otherElem]));
                   end
               end
           end
       end
       usedPairs = vertcat(usedPairs, newPairs);
       usedPairs = unique(usedPairs, 'rows');
       possiblePairs = setdiff(possiblePairs,usedPairs,'rows');
       indices1(randi, randjs) = 1;
       B = M(randi,randjs(1)) \ M(randi,randjs(2));
       constraintBuilder = constraintBuilder + sprintf("X(:,%d) == X(:,%d) .* %d\n", randjs(2),randjs(1),B);
   end
    constraintBuilder
    %indices1(:,1) = 1;
    %indices1(numRows+1:end,1:t) = randomObservationIndices(m-numRows,t,1/t);
    % disp("Proportion of observations from non-correlated columns")
    % (observationRate-(USED/(m*n)))/observationRate


    % run the matrix completion problem 'numTrials' number of times
    % At the end, we will calculate the 'averageError' in the solution
    % for each trial and we will 'count' the number of times we calculate
    % a solution with 'tol' of the exact solution.
    Xavg = zeros(size(M));
    count = 0;
    tol = 1e-5;
    averageError = 0;
    numTrials = 1;
    minDiff = Inf;
    nnz(indices1)
    for trial = 1:numTrials

        %   Generate a set of random observation indices
        indices1 = correlatedObservationIndices(m,t,0.75*k/length(j),indices1);
        indices2 = randomObservationIndices(m,n,max(0, observationRate-(nnz(indices1)/(m*n))), t);
        indices = [indices1,indices2];
        cvx_begin quiet
            variable X(m,n)
            minimize(norm_nuc(X))
            subject to

                %require the solution to match at the observed indices
                X(indices==1) == M(indices==1);
                eval(constraintBuilder)

        cvx_end
        diff = norm(M - X,2)/norm(M,2);
        if(diff < tol)
           count = count + 1; 
        end
        if(diff < minDiff)
           minDiff = diff;
           Xbest = X;
           IndBest = indices;
        end
        Xavg = Xavg + X;
        averageError = averageError + diff;    
    %     if mod(trial-1,10) == 0 
    %         fprintf('Trial %i done\n',trial);
    %     end
    end
    averageError = averageError/numTrials
    NumChosen = [NumChosen averageError];
    Xavg = Xavg/numTrials;

    fprintf('------- ADAPTIVE SAMPLE -------\n');
    fprintf('Proportion of observations from correlated columns: %.2f\n',1 - ((observationRate-(nnz(indices1)/(m*n)))/observationRate));
    fprintf('Number of trials: %i\n',numTrials);
    fprintf('Average relative error: %.4f\n',averageError);
    fprintf('Best relative error: %.4f\n',minDiff);
    tElapsed = toc(tStart);
    fprintf('Time elapsed: %.2f sec\n',tElapsed);


    A = (1:m) - 0.5; B = (1:n)-0.5;
    hFig = figure(120);clf;
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
    str = sprintf('Adapt Reconstruction, ||X - M||/||M|| = %.2e',minDiff);
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
end
X = [7,14,21,28,35,42,49]
plot(X,NumChosen)