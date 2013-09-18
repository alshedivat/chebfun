% Test file for @chebfun/max.m.

function pass = test_max(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% Check empty case.
[y, x] = max(chebfun());
pass(1) = isempty(y) && isempty(x);

% Check operation without breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), pref);
[y, x] = max(f);
y_exact = 1.884217141925336;
pass(2) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation with breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), ...
    linspace(-1, 1, 10), pref);
[y, x] = max(f);
y_exact = 1.884217141925336;
pass(3) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for complex-valued chebfuns.
f = chebfun({@(x) exp((1 + 1i)*x), @(x) sec(1i*(x - 0.5))}, [-1 0 1], ...
    pref);
[y, x] = max(f);
y_exact = 1;
pass(4) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for impulses.
f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
[y, ignored] = max(f);
pass(5) = y == 2;

f.impulses(1,1,1) = 10;
f.impulses(3,1,1) = -10;
[y, x] = max(f);
pass(6) = (y == 10) && (x == -1);

f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
f.impulses(1,1,2) = -1;
f.impulses(2,1,3) = 1;
[y, x] = max(f);
pass(7) = (y == inf) && (x == 0);

f.impulses(2,1,2) = -1;
[y, x] = max(f);
pass(8) = (y == 2) && (x == 1);

% Check computation of local extrema.
f = chebfun(@(x) sin(x).^2 + sin(x.^2), [0, 4]);
[y, x] = max(f, 'local');
y_exact = [1.923771282655145
           1.117294907913736
           1.343997479566445];
x_exact = [1.323339426259694
           2.781195946808315
           3.776766383330969];
pass(9) = numel(y == 3) && norm(y - y_exact, inf) < 10*vscale(f)*epslevel(f);

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(10*x) cos(10*x) exp(x)], [-1 -0.5 0.5 1]);
[y, x] = max(f);
y_exact = [1 1 exp(1)];
fx = feval(f, x(:));
fx = [fx(1,1) fx(2,2) fx(3,3)];
pass(10) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun(@(x) [exp((1 + 1i)*x) sec(1i*(x - 0.5))], [-1 0 1], ...
    pref);
[y, x] = max(f);
y_exact = [exp(1 + 1i) 1];
fx = feval(f, x(:));
fx = [fx(1,1) fx(2,2)];
pass(11) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun({[-1 1 2], [1 -1 3]}, [-1 0 1]);
f.impulses(3, 1, 2) = 1;
f.impulses(2, 3, 3) = -1;
[y, x] = max(f);
pass(12) = isequal(y, [inf 1 3]) && isequal(x, [1 -1 0]);

op = @(x) sin(x).^2 + sin(x.^2);
f = chebfun(@(x) [op(x) op(x/2)], [0, 4]);
[y, x] = max(f, 'local');
y_exact = [1.923771282655145  1.923771282655145
           1.117294907913736  NaN
           1.343997479566445  NaN];
x_exact = [1.323339426259694  2.646678852519388
           2.781195946808315  NaN
           3.776766383330969  NaN];
fx1 = feval(f, x_exact(:,1));
fx2 = feval(f, x_exact(1,2));
pass(13) = isequal(size(y), [3 2]) && ...
    all(isnan(y(2:end,2))) && all(isnan(x(2:end,2)));
pass(14) = norm(y(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(y(1,2) - y_exact(1,2), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx1(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx2(:,2) - y_exact(1,2), inf) < 10*vscale(f)*epslevel(f);

% [TODO]:  Test the max(f, g) syntax, once it is implemented.

end
