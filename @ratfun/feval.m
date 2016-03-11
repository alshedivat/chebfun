function out = feval(f, x)
%FEVAL   Evaluate a RATFUN.
%   FEVAL(F, X) evaluates a RATFUN F at the points in X.
%
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If f or x are empty, there's nothing to do.
if ( isempty(f) )
    out = [];
    return
elseif ( isempty(x) )
    % Return empty matrix with dimensions of the appropriate size.
    out = zeros(size(x));
    return
end


% output
out = feval(f, x);
    

end
