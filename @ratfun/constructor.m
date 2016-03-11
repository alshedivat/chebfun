function [r, deg, xk, fk, wk] = ratfun(f, tol)
%RATFUN   Compute rational interpolant for a function on [-1, 1].
%   [R, DEG] = RATFUN(F) computes a rational interpolant to the function handle
%   F, and returns also the McMillan degree DEG of R. R is a function handle.
% 
%   [R, DEG] = RATFUN(F, TOL) also specifies the tolerance of approximation.
%   Default is TOL = 1e-14.
% 
%   [R, DEG, XK, FK, WK] = RATFUN(F) also returns the parameters of the
%   barycentric interpolation formula used to represent R, i.e., the
%   interpolation points XK, the function values at these points FK, and the
%   weights WK.

%% Basic checks

% Add (?): Assert that f is a function handle.


%% Tolerance
if ( ( nargin == 1 ) || isempty(tol) )
    tol = 1e-14;
end


%% main loop
maxdegree = 500; % Maximal degree of interpolant.
Nmax = 2*maxdegree+1; % Maximal number of points in which to interpolate.
for N = 3:2:Nmax
    xx = chebpts(N); % interpolation points
%     xx = linspace(-1, 1, N).';
    fx = f(xx); % interpolation values
    T = 1:2:N; % column indices
    S = 2:2:N; % row indices
    L = loewner(xx, fx, S, T); % build Loewner matrix
    
    % Compute parameters for barycentric interpolation formula
    xk = xx(T); % interpolation points
    fk = fx(T); % interpolation values
    C = null(L);
    wk = C(:,1); % weights
    r = @(x) bary(x, fk, xk, wk); % Interpolant in xk.
    
    if ( r_interpolates(r, f, tol) )
        deg = length(xk)-1;
        return % We found a sufficiently accurate approximation of f.
    end
    % The above construction looks for the unique minimal degree interpolant (of
    % degree rank(L)). We could add the construction of non-unique interpolants
    % (which have degree N-rank(L)). This could be a short-cut to an accurate
    % interpolant.
end

% If interpolant not yet accurate enough, return last one with a warning and the
% current approximation error.
deg = length(xk)-1;
E = linspace(-1, 1, 1000);
warning(['Interpolant not accurate enough. Final approximation error on [-1,1] is ', ...
    num2str(norm(r(E)-f(E), Inf))])
end % end of RATFUN()


%% Accuracy test for interpolant
function happy = r_interpolates(r, f, tol)
%R_INTERPOLATES tests if |r(x)-f(x)| < tol for x in [-1, 1].

E = linspace(-1, 1, 1000); % grid for [-1, 1]
err = norm(r(E) - f(E), Inf);
if ( err < tol )
    happy = 1; % absolute accuracy good enough
% elseif ( err < tol*norm(f(E), Inf) )
%     happy = 1; % relative accuracy good enough
else
    happy = 0;
end
end % end of R_INTERPOLATES()


%% Construction of Loewner matrix
function L = loewner(x, y, S, T)
% Build Loewner matrix of size lenght(S) by length(T) from data in x and y.
% S = row indices, T = column indices.

X = zeros(length(S), length(T));
Y = zeros(length(S), length(T));
% Compute the rowwise numerators and denominators in L:
for iRow = 1:length(S)
    Y(iRow,:) = y(S(iRow)) - y(T);  % y(s_i) - y(t_j)
    X(iRow,:) = x(S(iRow)) - x(T);  % x(s_i) - x(t_j)
end
L = Y./X ;  % L_{i,j} = (y(s_i)-y(t_j)) / (x(s_i)-x(t_j))
end % end of LOEWNER()
