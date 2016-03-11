function [a,perm] = leja_ordering(a)
% [a,perm] = leja_ordering(a)
%
% Compute a Leja ordering for a given set of points a.
%
% See "Accuracy and Stability of Numerical Algorithms"
% by Nicholas J. Higham

n = max(size(a));
perm= (1:n);
[~,i] = max(abs(a));
if i ~= 1
    a([1 i]) = a([i 1]) ;
    perm([1 i]) = perm([i 1]) ;
end
p = ones(n,1);
for k=2:n-1
    for i=k:n
        p(i) = p(i) *(a(i)-a(k-1));
    end
    [~, i]=max(abs(p(k:n)));
    i = i+k-1;
    if i ~=k
        a([k i]) = a([i k]);
        p([k i]) = p([i k]);
        perm([k i]) = perm([i k]);
    end
end