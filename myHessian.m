function [H] = myHessian(nu, x ,Lo, bw)
[d, n] = size(x); 
x0 = Lo*(x - repmat(nu, [1,n]))./repmat(bw, [1,n]);
K = exp(-0.5*sum(x0.^2, 1));
grad = nan(d,1); 
for j = 1:d
    grad(j) = mean(K.*sum(x0.*repmat(Lo(:,j)./bw, [1,n]), 1));
end
grad = grad/mean(K); 
H = nan(d); 
for j = 1:d
    for k = 1:j
        H(j,k) = mean(   K.*( sum(x0.*repmat(Lo(:,j)./bw, [1,n]), 1).*sum(x0.*repmat(Lo(:,k)./bw, [1,n]), 1) - ...
            repmat(sum(Lo(:,j).*Lo(:,k)./bw.^2), [1,n]) )   );
        H(k,j) = H(j,k) ;
    end
end
H = H/mean(K) - grad*grad';
end