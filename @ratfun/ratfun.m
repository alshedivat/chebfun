classdef ratfun
%RATFUN   RATFUN class for rational representation of functions on [a,b].
%
%   Class for approximating functions defined on finite, semi-infinite, or
%   doubly-infinite intervals [a,b]. Functions may be smooth, piecewise smooth,
%   weakly singular, or blow up on the interval.
%

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RATFUN Class Description:
%
% The RATFUN class is for representations of functions with poles on the 
% interval [a, b]
%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        baryWeights;
        points;
        values;
        funHandle;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = ratfun(varargin)
            % The main RATFUN constructor!
            
            % Return an empty RATFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Function handle for the function to be approximated:
            fh = varargin{1};
            if ( ~isa(fh, 'function_handle') )
                error( 'RATFUN:ratfun:ctor', 'First argument must be a function handle')                
            end
            
            n = 8;
            collocPoints = chebpts(n);            
            fk = fh(collocPoints);
            tol = 1e-12;
            [r_handle, min_deg, Xk, Fk, Wk] = ratfun.minratinterp(collocPoints, fk, tol);
            f.baryWeights = Wk;
            f.points = Xk;   
            f.values = Fk;
            f.funHandle = r_handle;
            
        end                       
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Hidden = true, Static = false )        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = false )
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        [r_handle, min_deg, Xk, Fk, Wk] = minratinterp(x, y, tol);
        [a, perm] = leja_ordering(a);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Hidden = true, Static = true )
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Class-related functions: private utilities for this m-file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = str2op(op)
    % Convert string inputs to either numeric format or function_handles.
    sop = str2num(op); %#ok<ST2NM> % STR2DOUBLE doesn't support str2double('pi')
    if ( ~isempty(sop) )
        op = sop;
    else
        depVar = symvar(op);
        if ( numel(depVar) ~= 1 )
            error('CHEBFUN:CHEBFUN:str2op:indepvars', ...
                'Incorrect number of independent variables in string input.');
        end
        op = eval(['@(' depVar{:} ')', op]);
    end
end

function [op, dom, data, pref, flags] = parseInputs(op, varargin)
    

end