classdef trigfun2 < separableApprox
%TRIGFUN2   TRIGFUN2 class for representing periodic 
%                    functions on [-pi,pi]x[-pi,pi].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = trigfun2(varargin)
            % The main TRIGFUN2 constructor!
            
            % Return an empty TRIGFUN2:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});
            
        end
        
    end        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % The main bulk of the TRIGFUN2 constructor:
        g = constructor(g, op, dom, varargin);
        
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % For example: 
        
        % Convert Chebyshev coefficients to values:
        % X = coeffs2vals(U); 
        
        % Convert values to Chebyshev coefficients:
        % X = vals2coeffs(U); 
        
        % Padua points to tensor grid:
        % [C, V, X, Y] = paduaVals2coeffs( F, dom ); 
        
        % Tensor product of Chebyshev points:
        % [xx, yy] = chebpts2(nx, ny, domain, kind);
        
        % Outer-product of two chebfuns:
        % F = outerProduct(f, g);   
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by TRIGFUN2 class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    
end
