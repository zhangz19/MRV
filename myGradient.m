function [grad] = myGradient(nu, x ,Lo, bw)
[d, n] = size(x); 
x0 = Lo*(x - repmat(nu, [1,n]))./repmat(bw, [1,n]);
K = exp(-0.5*sum(x0.^2, 1));
grad = nan(d,1); 
for j = 1:d
    grad(j) = mean(K.*sum(x0.*repmat(Lo(:,j)./bw, [1,n]), 1));
end
grad = grad/mean(K); 

% a = 0; c = zeros(d,1);
% for i = 1:n
%     b = Lo*(x(:,i)-nu)./bw;
%     a = a + exp(-0.5*(b'*b));
%     for j = 1:d
%         c(j) = c(j) + exp(-0.5*(b'*b))*(b'*(Lo(:,j)./bw));
%     end
% end
% norm(a/n - mean(K))
% norm(c/a - grad)

end