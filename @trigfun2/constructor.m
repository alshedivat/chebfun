function g = constructor(g, varargin)
%CONSTRUCTOR   The main TRIGFUN2 constructor.

% For now call the chebfun2 constructor: 
f = chebfun2( varargin{:}, 'periodic' ); 

% Translate the properties across to trigfun2: 
g.cols = f.cols; 
g.rows = f.rows; 
g.pivotValues = f.pivotValues;
g.pivotLocations = f.pivotLocations; 
g.domain = f.domain; 

end