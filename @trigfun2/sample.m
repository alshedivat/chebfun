function varargout = sample( f, varargin )
%SAMPLE      Values of f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on a tensor
%   product grid.
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = U * D * V'.
%
%   [U, D, V] = SAMPLE(F,M,N) returns the values of F on a M-by-N
%   tensor product grid.
%
% See also CHEBCOEFFS2, PLOTCOEFFS2. 

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% CDR factorization: 
[C, D, R] = cdr( f ); 

% Evaluate: 
if ( nargout <= 1 )
    varargout = { C.values * D * (R.values).'}; 
else
    varargout = {C.values , D, R.values}; 
end

end