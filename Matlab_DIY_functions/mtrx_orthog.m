function [X1_orth] = mtrx_orthog(X1_to_orth, X0_to_orth_to)
%[X1_orth] = mtrx_orthog(X1_to_orth, X0_to_orth_to)
% mtrx_orthog orthogonalizes the matrix X1_to_orth to the matrix
% X0_to_orth_to. The output X1_orth contains X1_to_orth after
% it being orthogonalized to X0_to_orth_to.
%
% If X0_to_orth_to is empty, then mtrx_orthog orthogonalizes each column
% of X1_to_orth to all the preceding columns.
%
% INPUTS
% X1_to_orth: matrix where each column is a vector to be orthogonalized
% to the vectors contained in X0_to_orth_to (or to the preceding columns
% of X1_to_orth if X0_to_orth_to is empty)
%
% X0_to_orth_to: matrix where each vector is a column (or all empty).
% You want to orthogonalize X1_to_orth to this matrix.
%
% OUTPUTS
% X1_orth: X1_to_orth orthogonalized to X0_to_orth_to (or
% each column orthogonalized compared to the previous ones if
% X0_to_orth_to is empty)
%
% See also VBA_orth in the VBA toolbox which makes a similar job
%
% Initially written by N.Clairis with the great help of J.Brochard - march
% 2018

%% extract size of X1_to_orthog
[nX1_lines, nX1_col] = size(X1_to_orth);

%% check that the regressors have been correctly entered
if (nX1_lines == 1 && nX1_col > 1) || (exist('X0_to_orth_to','var') && ~isempty(X0_to_orth_to) && size(X0_to_orth_to,1) == 1 && size(X0_to_orth_to,2) > 1)
    error('One of the two inputs has only 1 line but many columns. Please enter each regressor/vector as a different column and the values of your regressor as lines so that the script works properly');
end

%% case where each column orthogonalized to the previous ones
if ~exist('X0_to_orth_to','var') || isempty(X0_to_orth_to)
    % initiate resulting matrix
   X1_orth = NaN(size(X1_to_orth));
   X1_orth(:,1) = X1_to_orth(:,1); % 1st column is not orthogonalized to anything => just equal
    
   % loop through columns to orthogonalize each column to the previous ones
   for iVec = 2:nX1_col
       % each column should be orthogonalized to the previous ones (once
       % they have been orthogonalized as well) => 1st extract the previous
       % orthogonalized columns and the current column to be checked
       X1_bis_for_orthog    = X1_orth(:,1:(iVec - 1));
       curr_vector          = X1_to_orth(:,iVec);
       
       % perform the orthogonalization
       P1      = pinv(X1_bis_for_orthog'*X1_bis_for_orthog)*X1_bis_for_orthog';
       id_mtrx = eye(size(X1_bis_for_orthog,1)); % extract identity matrix
       R1      = id_mtrx - X1_bis_for_orthog*P1;
       X1_orth(:,iVec) = R1*curr_vector;
   end
    
    %% orthogonalize X1_to_orthog to X0_to_orthog_to
else
    P0      = pinv(X0_to_orth_to'*X0_to_orth_to)*X0_to_orth_to';
    id_mtrx = eye(size(X0_to_orth_to,1)); % extract identity matrix
    R0      = id_mtrx - X0_to_orth_to*P0;
    % perform the orthogonalization
    X1_orth = R0*X1_to_orth;
    
end

end

