function [ A, ML, MR ] = randomSparseMat( m,n,r,ML_sparsity, MR_sparsity)
% What goes in:
%
%   m = desired number of rows
%   n = desired number of colums
%   r = desired rank 
%   ML_sparsity  = optional parameter, 
%                   proportion of non-zero entries in ML
%   MR_sparsity  = optional parameter, 
%                   proportion of non-zero entries in MR
%   NOTE: If only ML_sparsity is specified, then it will be used for
%         both the ML and MR sparsity. If neither is entered,
%         default values ae 0.3 for both.
%         
%
% What comes out:
%
%   A  = sparse random m-by-n matrix of rank r (or < r if we're unlucky)
%   ML = left factor of A
%   MR = right factor of A
%

if nargin == 3
   ML_sparsity = 0.3; 
   MR_sparsity = 0.3;
end

if nargin == 4
   MR_sparsity = ML_sparsity;
end

assert(ML_sparsity <= 1, 'Sparsity proportion larger than 1 will lead to infinite loop');
assert(MR_sparsity <= 1, 'Sparsity proportion larger than 1 will lead to infinite loop');

ML = zeros(m,r); MR = zeros(r,n);

%generate ML and MR by choosing random indices 
% and filling them with random numbers
% drawn uniformly from [0,1]
while(sum(sum(ML~=0)) <= ML_sparsity*r*m)
   i = randi([1 m],1);
   j = randi([1 r],1);
   if(ML(i,j)==0)
      ML(i,j) = rand(1,1); 
   end
end
while(sum(sum(MR~=0)) <= MR_sparsity*r*n)
   i = randi([1 r],1);
   j = randi([1 n],1);
   if(MR(i,j)==0)
      MR(i,j) = rand(1,1); 
   end
end

A = ML*MR;

end

