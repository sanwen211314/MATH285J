function [ indices ] = randomObservationIndices( m,n,observationRate,t)
% What goes in:
%   m = columns of matrix
%   n = rows of matrix
%   observationRate = the proportion of the matrix you want to observe
%
% What comes out:
%
%   indices = the list of observed indices in the form of a binary matrix
%                   1 means you observe this entry
%                   0 means you don't observe this entry
if nargin == 3
    if(observationRate >= 1)
        indices = ones(m,n);
        fprintf('Observation rate >= 1. You are observing all entries.\n');
    else
        indices = zeros(m,n);
        while(sum(sum(indices)) < m*n*observationRate)
            i = randi([1 m],1);
            j = randi([1 n],1);
            if indices(i,j) == 0
                indices(i,j) = 1;
            end
        end
    end
else
    if(observationRate >= 1)
        indices = ones(m,n-t);
        fprintf('Observation rate >= 1. You are observing all entries.\n');
    else
        indices = zeros(m,n-t);
        while(sum(sum(indices)) < m*n*observationRate)
            i = randi([1 m],1);
            j = randi([1 n-t],1);
            if indices(i,j) == 0
                indices(i,j) = 1;
            end
        end
    end
end

