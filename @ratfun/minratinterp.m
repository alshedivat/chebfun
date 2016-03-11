function [r_handle, min_deg, Xk, Fk, Wk] = minratinterp(x, y, tol)
%MINRATINTERP   Minimal degree rational interpolation.
%   [R_HANDLE, MIN_DEG] = MINRATINTERP(X, Y) computes a rational function
%   R_HANDLE interpolating the points in X and Y and of minimal McMillan degree
%   MIN_DEG. (The McMillan degree of P/Q is the maximum of deg(P) and deg(Q).)
%   The points in X must be distinct.
%
%   [R_HANDLE, MIN_DEG] = MINRATINTERP(X, Y, TOL) specifies the tolerance. By
%   default, TOL = 1e-12.
% 
%   [R_HANDLE, MIN_DEG, XK, FK, WK] = MINRATINTERP(X, Y, TOL) also returns the
%   parameters of the barycentric interpolation formula for R_HANDLE, i.e., the
%   interpolation points XK, the function values FK, and the weights WK.
%
%   Reference:
%   [1] A. Antoulas, B. Anderson, "On the scalar rational approximation
%       problem", IMA J. Math. Control Inf., 3:61-88, 1986.

%% basic checks
if ( ~isvector(x) || ~isvector(y) )
   error( 'minratinter:notvector', 'make sure x and y are both vectors' );
else
    x = x(:); % Ensure x is a column vector.
    y = y(:); % Ensure y is a column vector.
end

if ( length(x) ~= length(y) )
    error( 'minratinerp:length', 'Inconsistent number of interpolation points and values.');
else
    N = length(x);
end
   
if ( N == 1 )
    r_handle = @(z) y;
    min_deg = 0;
    Xk = 0;
    Fk = y;
    Wk = 1;
    return
end


    
%% Tolerance
if ( ( nargin == 2 ) || isempty(tol) )
    tol = 1e-12;
end

%% Get row and column indices (S and T) for Loewner matrix
m = floor(N/2);
[~, perm] = ratfun.leja_ordering(x);
T = perm(1:(N-m));
S = perm((N-m+1):N);

%% Compute the interpolant
L = loewner(x, y, S, T); % get Loewner matrix of size m x (N-m)
q = rank(L);

if ( q < N-m )
    % We are possibly in case (a) of [1, Theorem 2.25] (unique interpolant).
    % Build interpolant r_handle for x(T) with barycentric interpolation
    % formula:
    C = null(L); % Nontrivial since q = rank(L) < N-m.
    Xk = x(T);   % interpolation points
    Fk = y(T);   % function values
    Wk = C(:,1); % weights, any vector ~= 0 is fine, so pick first column
    r_handle = @(z) bary(z, Fk, Xk, Wk); % The interpolant for x(T).
    
    % Does r_handle interpolate all points in x?
    if ( r_interpolates(r_handle, Xk, Fk, x, y, tol) )
        min_deg = q;
        disp('Minimal degree interpolant is UNIQUE.')
        fprintf('The minimal degree is %d.\n', min_deg)
        return % We found the unique minimal degree interpolant.
    end
end

%% We are in case (b) of [1, Theorem 2.25] (non-unique interpolants).
min_deg = N-q;
disp('Minimal degree interpolant is NOT UNIQUE.')
fprintf('The minimal degree is %d.\n', min_deg)

% Build (q-1) x (N-q+1) Loewner matrix and attached rational function
% according to [1, Theorem 2.26].
T = perm(1:(N-q+1));
S = perm((N-q+2):N);

% Parameters for the barycentric interpolation formula:
Xk = x(T); % interpolation points
Fk = y(T); % function values

L = loewner(x, y, S, T);
C = null(L);
C = [sum(C, 2), C]; % add as first column the sum of all columns of C. 
% This should help against zero entries.
% Go through the columns of C (= possible weights)
for j = 1:size(C,2)
    Wk = C(:,j); % weights
    r_handle = @(z) bary(z, Fk, Xk, Wk);
    
    if ( r_interpolates(r_handle, Xk, Fk, x, y, tol) )
        disp(['Weight vector positive at index ' num2str(j)])
        return % We found a non-unique interpolant.
    end
end
warning('Naive handling of case 2.25(b) in [Antoulas & Anderson, 1986] failed.')
end % End of MINRATINTERP()




%%
function happy = r_interpolates(rk, Xk, Fk, x, y, tol)
% Check if abs(r(x) - y) is sufficiently small.
err = norm(rk(x) - y, Inf); % Error in interpolation points not in Xk:
% Bary returns r(Xk) = Fk, even if the rational function does not interpolate
% these points, which may happen if one of the weights is (almost) zero.
if ( ( err < tol ) || ( err < tol*norm(y, Inf) ) )
    happy = 1; % error (or relative error) in points other than xk is ok
else % not accurate enough
    happy = 0; 
    return
end

% Test accuracy in points xk:
err_Xk = abs(rk(Xk + 1e-13) - Fk);
if ( ( max(err_Xk) < 100*tol ) || all(err_Xk < 100*tol*Fk) )
    return % Absolute or all relative errors are sufficiently small.
else
    happy = 0;
end
end % end of R_INTERPOLATES()




%%
function L = loewner(x, y, S, T)
% Build Loewner matrix of size lenght(S) by length(T) from data in x and y.
% S = row indices, T = column indices.

% Compute the rowwise numerators and denominators in L:
for iRow = 1:length(S)
    Y(iRow,:) = y(S(iRow)) - y(T);  % y(s_i) - y(t_j)
    X(iRow,:) = x(S(iRow)) - x(T);  % x(s_i) - x(t_j)
end
L = Y./X ;  % L_{i,j} = (y(s_i)-y(t_j)) / (x(s_i)-x(t_j))
end % end of LOEWNER()