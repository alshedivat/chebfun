function varargout = subsref(f, index)
%SUBSREF   RATFUN subsref.
% ( )
%   F(X) returns the values of the RATFUN F evaluated on the array X. 
%
% {}
%
%
% See also FEVAL

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document for array-valued CHEBFUN objects and quasimatrices.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'                
                  
        % Where to evaluate:
        x = idx{1};                         
            
        if ( length(idx) > 1 )
            error('RATFUN:RATFUN:subsref:dimensions', ...
                'Index exceeds RATFUN dimensions.')            
        end

        % Compute the output:
        if ( isnumeric(x) )
            % Call FEVAL():
            out = feval(f.funHandle, x);                               
        else
            error('RATFUN:RATFUN:subsref:nonnumeric',...
              'Cannot evaluate ratfun for non-numeric type.')          
        end                
        
        
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '{}'

        
end

% Convert to a cell:
varargout = {out};

end

