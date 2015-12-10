function f = compose( f, op, varargin )
% COMPOSE     compose command for TRIGFUN2 objects. 
% 
%  F = COMPOSE(F, OP )  returns the TRIGFUN2 that approximates OP(F).
% 
%  F = COMPOSE(F, OP, G )  returns the TRIGFUN2 that approximates OP(F,G).
%
% This command is a wrapper for the TRIGFUN2 constructor. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 && nargin(op) == 1)
    % OP has one input variable.
    
    % Call constructor: 
    f = trigfun2(@(x,y) op( feval(f, x, y) ), f.domain);
    
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     % promote
        g = trigfun2(g, f.domain); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = trigfun2(f, g.domain); 
    end
    
    % Call constructor: 
    f = trigfun2(@(x,y) op( feval(f, x, y), feval(g, x, y) ), f.domain);
    
else
    % Not sure what to do, error: 
    error('CHEBFUN:TRIGFUN2:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end 