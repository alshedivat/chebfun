function varargout = sample( f, varargin )
%SAMPLE      Values of f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on a tensor
%   product grid.
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = U * D * V'.
%
%   SAMPLE(F,M) returns the values of F on a M-by-M tensor product grid.
%
%   SAMPLE(F,M,N) returns the values of F on a M-by-N tensor product grid.
%
% See also CHEBCOEFFS2, PLOTCOEFFS2. 

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% Get grid size to evaluate on: 
[m, n] = length(f); 
if ( nargin > 1 ) 
    m = varargin{1}; 
    n = m; 
end
if ( nargin > 2 ) 
    n = varargin{2};
end

% Use CDR decomposition so we can keep it in low rank form: 
[C, D, R] = cdr( f ); 
% Cvals = sample(chebfun(C), m);
% Rvals = sample(chebfun(R), n);
% See https://github.com/chebfun/chebfun/issues/1808
Cvals = sample(C, m);
Rvals = sample(R, n);


% Evaluate: 
if ( nargout <= 1 )
    varargout = {Cvals * D * Rvals.'}; 
else
    varargout = {Cvals , D, Rvals}; 
end

end
