%POPULATE   Populate a FOURTECH class with values.
%   F = F.POPULATE(OP) returns a FOURTECH representation populated with values
%   F.VALUES of the function OP evaluated on an equally spaced grid. The fields
%   F.ISHAPPY and F.EPSLEVEL indicate whether the representation is deemed
%   'happy' and to what accuracy (see HAPPINESSCHECK.m). Essentially this means
%   that such an interpolant is a sufficiently accurate (i.e., to a relative
%   accuracy of F.EPSLEVEL) approximation to OP. If F.ISHAPPY is FALSE, then
%   POPULATE was not able to obtain a happy result.
%
%   OP should be vectorized (i.e., accept a vector input), and output a vector
%   of the same length. Furthermore, OP may be an array-valued function, in
%   which case it should accept a vector of length N and return a matrix of size
%   NxM.
%
%   F.POPULATE(OP, VSCALE, HSCALE) enforces that the happiness check is relative
%   to the initial vertical scale VSCALE and horizontal scale HSCALE. These
%   values default to 0 and 1 respectively. During refinement, VSCALE updates
%   itself to be the largest magnitude values to which (each of the columns in)
%   OP evaluated to.
%
%   F.POPULATE(OP, VSCALE, HSCALE, PREF) enforces any additional preferences
%   specified in the preference structure PREF (see FOURTECH.TECHPREF).
%
%   F.POPULATE(VALUES, ...) (or F.POPULATE({VALUES, COEFFS}, ...)) populates F
%   non-adaptively with the VALUES (and COEFFS) passed. These values are still
%   tested for happiness in the same way as described above, but the length of
%   the representation is not altered.
%
% See also FOURTECH, TECHPREF, HAPPINESSCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
function f = populate(f, op, vscale, hscale, pref)

if ( (nargin < 3) || isempty(vscale) )
    vscale = 0;
end
if ( (nargin < 4) || isempty(hscale) )
    f.hscale = 1;
else
    f.hscale = hscale;
end
if ( nargin < 5 )
    pref = chebtech.techPref();
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Non-adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%
% Values (and possibly coefficients) have been given.
if ( isnumeric(op) || iscell(op) )
    if ( isnumeric(op) )
        % OP is just the values.
        f.values = op;
        f.isReal = isreal(op);
        f.coeffs = f.vals2coeffs(op);
    else                 
        % OP is a cell {values, coeffs}
        f.values = op{1};
        f.isReal = isreal(op);
        f.coeffs = op{2};
        if ( isempty(f.values) )
            f.values = f.coeffs2vals(f.coeffs);
            if max(abs(imag(f.values))) > 1e3*eps
                f.isReal = false;
            else
                f.isReal = true;
            end
        end
    end
    
    % Update vscale:
    f.vscale = max(abs(f.values), [], 1);
    
    % We're always happy if given discrete data:
    f.ishappy = true;
    
    % Scale the epslevel relative to the largest column:
    vscale = f.vscale;
    f.epslevel = 10*eps(max(f.vscale));
    vscale(vscale <= f.epslevel) = 1;
    f.epslevel = f.epslevel./vscale;

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty values to pass to refine:
f.values = [];

% Check a random value of op in (-pi,pi) to see if the result is complex.
f.isReal = isreal(feval(op,pi*(2*rand-1)));

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine: (in PREF.REFINEMENTFUNCTION)
    [f.values, giveUp] = f.refine(op, f.values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale: (Only include sampled finite values)
    valuesTemp = f.values;
    valuesTemp(~isfinite(f.values)) = 0;
    vscale = max(vscale, max(abs(valuesTemp(:))));
    
    % Compute the Fourier coefficients:
    coeffs = f.vals2coeffs(f.values);
    
    % Check for happiness:
    f.coeffs = coeffs;
    f.vscale = vscale;
    [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref); 
        
    % We're happy! :)
    if ( ishappy ) 
        coeffs = f.alias(coeffs, cutoff);   % Alias the discarded coefficients.
        f.values = f.coeffs2vals(coeffs);   % Compute values on this grid.
        break
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update the vscale. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the 'true' vscale (as defined in FOURTECH classdef):
vscaleOut = max(abs(f.values), [], 1);
% Update vertical scale one last time:
vscaleGlobal = max(vscale, vscaleOut);

% Output the 'true' vscale (i.e., the max of the stored values):
vscale = vscaleOut;

% Adjust the epslevel appropriately:
% if ( any(vscaleOut > 0) )
%     epslevel = epslevel*vscaleGlobal./vscaleOut;
% else 
%     % Deal with zero vscale:
%     epslevel = epslevel./(1+vscaleOut);
% end
vscaleOut(vscaleOut < epslevel) = 1;
epslevel = epslevel*vscaleGlobal./vscaleOut;

% This may not be the best place to do this.
if f.isReal
    f.values = real(f.values);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Assign to FOURTECH object. %%%%%%%%%%%%%%%%%%%%%%%%%%
f.coeffs = coeffs;
f.vscale = vscale;
f.ishappy = ishappy;
f.epslevel = epslevel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouput. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ishappy )
    % We're done, and can return.
    return
end

end