fh = @(x) 1./(1+5*x.^2);
%%x = -1 + 2*rand(3, 1);
x = [-1, -.25, 0, 1];
y = fh(x);

%%
[r, deg, xk, fk, wk] = minratinterp(x, y);

%%
% Compute pole locations of r
a = wk/sum(wk);
b = xk;
P = diag(xk) - a*b.';
poles = eig(P);
poles(abs(poles) < 1e-13) = [];
poles(abs(poles) > 1 ) = []
%%
xx = -1:.01:1;
poleCheck = 1e-10*xx + poles;
[~, idx] = max(abs(r(poleCheck)));
poleCheck(idx)
poles
%%
xx = -1:.01:1;
plot( xx, fh(xx));
hold on
plot(xx, r(xx))
hold off
norm(r(xx)-fh(xx), inf)

%%
