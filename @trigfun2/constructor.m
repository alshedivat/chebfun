function g = constructor(varargin)
%CONSTRUCTOR   The main TRIGFUN2 constructor.

f = chebfun2( varargin{:} ); 

g.cols = f.cols; 
g.rows = f.cols; 

end