function H = cross( F, G )
%CROSS   Vector cross product.
%   CROSS(F, G) returns a SPHEREFUNV representing the 3D cross product of
%   the SPHEREFUNV objects F and G.
%
%   See also DOT

%   CROSS(F, G, 'n') returns a SPHEREFUN representing the normal component
%   of the cross product of the SPHEREFUNV objects F and G.  
%   For two vector fields tangent to the sphere, the cross product is a
%   a vector that points in the normal direction.  CROSS returns a 
%   SPHEREFUN object representing the component in the normal (radial)
%   direction.  Mathematically, this is N dot (F x G), where N is the
%   normal to the sphere.

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    H = spherefun;
    return
end

% Get the components: 
Fc = F.components; 
Gc = G.components; 

% Do cross: 
H = [ Fc{2} .* Gc{3} - Fc{3} .* Gc{2} ; ...
      Fc{3} .* Gc{1} - Fc{1} .* Gc{3} ; ...
      Fc{1} .* Gc{2} - Fc{2} .* Gc{1} ];
  
end