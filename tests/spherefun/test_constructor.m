function pass = test_constructor( ) 
% Test the spherefun constructor 

% Get tolerance: 
tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

f = @(x,y,z) x.^2 + y.^2 + z.^2;
f = redefine_function_handle(f);
g = spherefun(f);
pass(1) = ( SampleError(f, g) < tol );

f = @(x,y,z) exp(-cos(pi*(x+y+z)));
f = redefine_function_handle(f);
g = spherefun(f);
pass(2) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) 1-exp(x);
g = spherefun(f);
f = redefine_function_handle(f);
pass(3) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) exp(y);
g = spherefun(f);
f = redefine_function_handle(f);
pass(4) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) exp(z);
g = spherefun(f);
f = redefine_function_handle(f);
pass(5) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) cos(x.*y);
g = spherefun(f);
f = redefine_function_handle(f);
pass(6) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) sin(x.*y.*z);
g = spherefun(f);
f = redefine_function_handle(f);
pass(7) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) sin(x+ y.*z);
f = redefine_function_handle(f);
g = spherefun(f);
pass(8) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) sin(x+ y.*z) + 1;
f = redefine_function_handle(f);
g = spherefun(f);
pass(9) = ( SampleError(f, g) < tol ); 

f = @(x,y,z) 0*x;
g = spherefun(f);
pass(10) = ( norm(g, inf) == 0 ); 

% Test the vectorize flag is working: 
f = spherefun(@(x,y,z) cos(z));
g = spherefun(@(x,y,z) cos(z), 'vectorize');
pass(11) = ( norm(f - g) < tol );

f = spherefun(@(x,y,z) x.*y.*z);
g = spherefun(@(x,y,z) x*y*z, 'vectorize');
pass(12) = ( norm(f - g) < tol );

% Test construction from samples
f = spherefun(@(x,y,z) 1 + x.*sin(x.*y));
[m,n] = length(f);
F = sample(f, m + mod(m, 2), n);
g = spherefun(F);
pass(13) = ( norm(f - g) < tol );

f = spherefun(@(x,y,z) 1 + 0*x);
F = ones(2, 2);
g = spherefun(F);
pass(14) = ( norm(f - g) < tol );

F = ones(1, 2);
try
    g = spherefun(F);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier,'CHEBFUN:SPHEREFUN:constructor:poleSamples');
end

end

function f = redefine_function_handle(f)
% nargin(f) = 2, then we are already on the sphere, if nargin(f) = 3,
% then do change of variables:

if ( nargin(f) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f, lam, th, 0);
end

end

function sample_error = SampleError(h, g)
m = 6; 
n = m;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = linspace(0, pi, m).';

end