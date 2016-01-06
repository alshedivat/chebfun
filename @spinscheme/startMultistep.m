function [uSol, dt] = startMultistep(K, adapTime, dt, L, LR, Nc, Nv, pref, ...
    S, uSol)
%STARTMULTISTEP  Get enough initial data when using a multistep scheme.
%    [USOL, DT] = STARTMULTISTEP(K, ADAPTIME, DT, L, LR, NC, PREF, S, USOL) does
%    a few steps of EXPRK5S8 with timestep DT to get enough initial data C to 
%    start the multistep SPINSCHEME K, using the linear part L, the linear part 
%    for complex means LR, the nonlinear part of the operator in coefficient
%    space NC, the SPINPREFERENCE object PREF, and the SPINOPERATOR S. ADAPTIME
%    is 1 if adpative in time, 0 otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
errTol = pref.errTol;
M = pref.M;
q = K.steps;
nVars = S.numVars;

% Number of points to discretize the domain:
N = size(uSol{1}, 1)/nVars;

% Create a cell-array to store the coefficients at the Q steps:
coeffs = cell(q, 1);

% Store the initial conidition in the last column:
coeffs{q} = uSol{1};

% Set-up EXPRK5S8 (5th-order one-step scheme):
K = spinscheme('exprk5s8');
schemeCoeffs = computeCoeffs(K, dt, L, LR, S);
if ( adapTime == 1 )
    LR2 = computeLR(S, dt/2, L, M, N);
    schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
end

% Do Q-1 steps of EXPRK5S8:
iter = 1;
while ( iter <= q-1 ) 
    
    cnew = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uSol);
     
    % If adpative in time, two steps in time with DT/2 and N points:
    if ( adapTime == 1 )
        cnew2 = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, uSol);
        cnew2 = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, cnew2);
        err = max(max(max(abs(cnew{1} - cnew2{1}))));
        % If successive step, store it:
        if ( err < errTol ) 
            coeffs{q-iter} = cnew2{1};
            uSol = cnew2;
            iter = iter + 1;
        % Otherwise, redo the step with DT/2 and N points:
        else
            dt = dt/2;
            schemeCoeffs = schemeCoeffs2;
            LR2 = computeLR(S, dt/2, L, M, N);
            schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
        end
        
    % If not adaptive in time, keep CNEW:
    else
        coeffs{q-iter} = cnew{1};
        uSol = cnew;
        iter = iter + 1;
    end

end
uSol = coeffs;

end