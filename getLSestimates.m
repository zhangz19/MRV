function [betas, gammas] = getLSestimates(lTstar, Xt, Z)
[~, L, p] = size(Xt);
Xall = [];
for j=1:p
    Xall = [Xall, Xt(:,:,j)];
end
if size(Z,1) > 1
    Xall = [Xall, Z];
end
betas = (Xall'*Xall)\(Xall'*lTstar);
gammas = betas((L*p+1):end);
if size(Z,1)==1
    gammas = 0;
end
betas = reshape(betas(1:(L*p)), [L,p]);

% % visualize
% load('datMRV.mat')
% Xnams = {'Intercept', 'Grade 1','Grade 2','Stage I','Stage II','Stage III','Black race','ER','PR'};
% % Xnams = {'Intercept', 'Grade 1','Grade 2','Stage I','Stage III','Black race'};
% for j = 1:p
%     fu = betas(:,j)'*B'; 
%     [~,I] = sort(U); 
%     subplot(3,3,j); plot(U(I),fu(I)', 'k.-')
%     title(Xnams(j))
%     axis tight
% end

% save('real_inits','betas','gammas')
end