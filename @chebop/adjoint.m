function [N, adjcoeffs] = adjoint(N)
%ADJOINT   Compute the adjoint of a linear CHEBOP.
%   N = ADJOINT(N), where N is a CHEBOP, returns the adjoint CHEBOP of N.
%
%   [N, adjcoeffs] = ADJOINT(N) also returns a CHEBMATRIX ADJCOEFFS which stores 
%   the (variables) coefficients of the adjoint. The indexation is as follows:
%      N = adjcoeffs{1}*u^(n) + adjcoeffs{2}*u^(n-1) + ... + adjcoeffs{n+1}*u
%
% See also LINOP/ADJOINT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Tell CHEBOP/LINEARIZE to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, ~] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:adjoint:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'ADJOINT supports only linear CHEBOP instances.']);

% Get the value of the highest derivative:
n = L.diffOrder;

nbcs = length(L.constraint.values);

% Trivial case:
if ( n == 0 )
    return
end

if ( isempty(N.bc) )
    if ( isempty(N.lbc) && isempty(N.rbc) )
        bcType = 'nobc';
    elseif ( nbcs ~= n )
        bcType = 'bvp';
    elseif ( isempty(N.rbc) )
        bcType = 'ivp';
    elseif ( isempty(N.lbc) )
        bcType = 'fvp';
    else
        bcType = 'bvp';
    end
elseif ( ~strcmpi(N.bc, 'periodic') )
    bcType = 'periodic';
else
    error('CHEBFUN:CHEBOP:adjoint:boundaryconds', ...
        'ADJOINT does not supprot this case');
end

%% FORMAL ADJOINT

% Extract the blocks of the LINOP:
blocks = L.blocks;

% ADJOINT doesn't support system of equations for the moment.
% [TODO]: Support system of equations.
if ( max(size(blocks,2)) > 1 )
    error('CHEBFUN:LINOP:adjoint:system', ...
        'ADJOINT doesn''t support system of equations for the moment.');
end

% Get the first block and the domain:
block = blocks{1};
dom = L.domain;

% Get the coefficients:
if ( n == 0 )
    % This is a multiplication operator, nothing to do here:
    if ( nargout > 1 )
        adjcoeffs = conj(block);
    end
    return
else
    coeffs = toCoeff(block);
end

% Compute the coefficients of the adjoint:
adjcoeffs = cell(n+1,1);
for k = 0:n
    adjcoeffs{k+1} = chebfun('0', dom);
end
for k = 0:n
    for l = 0:k
        adjcoeffs{n+1-l} = adjcoeffs{n+1-l} + ...
            (-1)^k*nchoosek(k,l)*conj(diff(coeffs{n+1-k}, k-l));
    end
end

% Construct a CHEBOP from these new coefficients:
op = @(x, u) coeffs2func(u, adjcoeffs);
N = chebop(op, dom);

%% BOUNDARY CONDITIONS
str = '@(u) [u';
for k = 1:(n-1)
    str = [str, ' ; diff(u, ' int2str(k) ')']; %#ok<AGROW>
end
str = [str '];'];
f = str2func(str);

switch ( bcType )
    case 'periodic'
        N.bc = 'periodic';
    case 'nobc'
        N.lbc = f;
        N.rbc = f;
        
    case 'ivp'
        N.rbc = f;
        
    case 'fvp'       
        N.lbc = f;
        
    case 'bvp'
        % Do stuff. (Can of worms.)
        error('you suck')
        
end
   
if ( nargout > 1 )
    % Store the coefficients in a CHEBMATRIX:
    adjcoeffs = chebmatrix(adjcoeffs);
end

end

function out = coeffs2func(u, coeffs)
%COEFFS2FUNC   Output a function representing the operator field of the adjoint.

out = 0;
n = length(coeffs);
for k = 1:n
    out = out + coeffs{k}.*diff(u, n-k);
end

end
