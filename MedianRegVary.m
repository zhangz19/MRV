function [matpara, boot_weights] = MedianRegVary(V, Delta, Z, X, B)
% global simu  Xt Trues filenames  B U pilotRun X indmat performTest  %V Delta Z
% fprintf('case = %d, chain = %d:\n', [nvar, ch])

q = size(Z, 2);  p = size(X,2);  [n, L] = size(B);
Xt = zeros(n, L, p);  for j = 1:p; Xt(:,:,j) = B.*repmat(X(:,j), [1,L]);  end
cind = find(Delta==0); ncind = find(Delta==1);  %censored and non-censored
L_cind = length(cind);

%=============== basic running setup
simu = 1;
verbose = 0;
fixcluster = 0;
updateAll = 1;
reportrate = 0;
FDDP = 1; % use finite-dimensional DP approximation?  Recommended.
pilotRun = 1; %less run for testing code
performTest = 1;
% tot = 6e3; burn = 5e3; thin = 1;
tot = 5e4; burn = 4e4; thin =10;
if pilotRun == 1; tot = 50; burn = 0; thin = 1; end  % for testing if code can be executed
nsample = (tot-burn)/thin;

%----------------------- set prior
gamma_pre = 1e-4*eye(q);
gamma_mu = gamma_pre*zeros(q,1);
mean_IG = 0.01; var_IG = 10^4;
alpha_IG = 2+mean_IG^2/var_IG;
invbeta_IG = mean_IG*(alpha_IG-1);
A_v = alpha_IG*ones(1,p) + 0.5*L;
B_v = invbeta_IG*ones(1,p);
A_tau = alpha_IG*ones(1,2);
B_tau = invbeta_IG*ones(1,2);
A_a = 2; B_a = 0.01;   % gamma prior when FDDP=0
if FDDP == 1
    N = 150; %fix number of clusters to be large
    %ps = 1/N*ones(1,N);
    ps = [.9, (1-.9)/(N-1)*ones(1, N-1)];
    %labs = randsample(N, n, true, ps);
else
    N = 1; % start with 1 cluster
end
% use informative prior for base measure
mean_IG = 1./[0.3807, 0.2710].^2; % 0.3807, 0.2710: estimated nu(2,4): sqrt precesion
if simu==1;  mean_IG = 1./[1.3173 , 3.0886].^2; end  %for simulated data
var_IG = .01; %.01; % note this is prior for the variance
alpha_IG = 2+mean_IG.^2./var_IG;
A_nu = alpha_IG;
invbeta_IG = mean_IG.*(alpha_IG-1);
B_nu = invbeta_IG; % note B_nu is the scale parameter! In the Denominator
% mean( 1./gamrnd(A_nu(1), 1/B_nu(1), [1,1e4]) )
% disp(mean_IG)
% B_nu./(A_nu-1) % this is mean
nu_mean = [log(0.1411), log(0.1059)]; %zeros(1,2);
useMALA4alpha = 1;  %use Metropolis-adjusted Langevin algorithm (MALA)

%-----------------------  initial calculation
Pj = zeros(L);
Pj(1,1) = 1/1e4; Pj(2,2) = 1/1e4; % 1 over c^2_{j1} and c^2_{j2}
% Pj(1,1) = 1/1e1; Pj(2,2) = 1/1e1;
for l = 3:L
    Lambda_l = zeros(L,1);
    Lambda_l([l-2, l]) = 1; Lambda_l(l-1) = -2;
    Pj = Pj + Lambda_l*Lambda_l';
end
lTstar = log( V + unifrnd(0,1,[n,1]).*(Delta==0) );
[betas, gammas] = getLSestimates(lTstar, Xt, Z);
inv_vartheta = ones(1, p);
for j = 1:p
    inv_vartheta(j) = gamrnd( A_v(j),   1/(B_v(j)+0.5*sum(betas(:,j).*( Pj*betas(:,j) ))) );
end
mu = zeros(n,p); %mu1 = zeros(n,p);
for j = 1:p;  mu(:,j) = Xt(:,:,j)*betas(:,j);  end
Zpart = Z*gammas;
err = log(V) - sum(mu,2) - Zpart;
alpha =1;
labs = ones(n, 1); %start with single cluster
ms = histc(labs, 1:N);
nu = zeros(N,4);
invtau = 5*ones(1,2);
nu(:,2) = 1*sqrt(gamrnd( A_nu(1), 1/B_nu(1), [N,1] )); %sqrt precision
nu(:,4) = 1*sqrt(gamrnd( A_nu(2), 1/B_nu(2), [N,1] ));
nu(:,1) =  - exp(randn([N,1])./(nu(:,2)*invtau(1)) + nu_mean(1));
nu(:,3) =  exp(randn([N,1])./(nu(:,4)*invtau(2)) + nu_mean(2));

% estimates nu from err for singleton cluster
err_med = median(err);
indg0 = find(err < err_med);
nu(1,1) = min(-0.001,mean(err(indg0))); nu(1,2) = 1/std(err(indg0));
indg1 = find(err >= err_med);
nu(1,3) = max(0.001,mean(err(indg1))); nu(1,4) = 1/std(err(indg1));

nclabs = labs(ncind);
as = (0.5 - normcdf(-nu(labs,3).*nu(labs,4)))./( normcdf(-nu(labs,1).*nu(labs,2)) - normcdf(-nu(labs,3).*nu(labs,4)) ); % q(\theta)
bs = 1-as;
tmp = (err - nu(labs,1)).*nu(labs,2); %standardized residual
as(ncind) = as(ncind).*nu(nclabs,2).*exp( -0.5*tmp(ncind).^2);
as(cind) = as(cind).*normcdf( tmp(cind) );
tmp = (err - nu(labs,3)).*nu(labs,4);
bs(ncind) = bs(ncind).*nu(nclabs,4).*exp( -0.5*tmp(ncind).^2 );
bs(cind) = bs(cind).*normcdf( tmp(cind) );
gs = binornd( 1, as./(as+bs) );

% % initial becomes true
% if simu == 1
%     lTstar =   sum(mu,2) + Zpart + err;
% end

% Adaptive MCMC set-up
objrate = 0.44;
batchLen = 50; %min(50, tot);
batchNum = 0;
batchTot = tot/batchLen;
nlen = 4;
nu_accepts = zeros(N,nlen);
nu_rates = zeros(batchTot, nlen);
nu_tunings = repmat(log(0.1*abs([10000000*log(-nu(1,1)), -2*log(nu(1,2)), 1000000000*log(nu(1,3)), -2*log(nu(1,4))])), [N,1]); %ones(N, nlen);
if simu == 1
    nu_tunings = repmat(log(1e4*abs([-1e2*log(-nu(1,1)), -4*log(nu(1,2)), log(nu(1,3)), -2*log(nu(1,4))])), [N,1]); %ones(N, nlen);
end

% Adaptive MALA set-up
h_alpha = .01; alpha_M = 10;  alpha_accepts = zeros(1, tot);
alpha_objrate = 0.574;

% store results
matpara = nan(nsample, numel(betas)+numel(gammas)+numel(alpha));
Ms = zeros(1, nsample); %cell
% Ps = zeros(nsample, N);
% Theta = zeros(1, 4, nsample);
% BF = zeros(1, tot-burn);

% for bootstrap
if simu == 0
    load('indmat.mat')
else
    nboot = 500;
    %load('indmat_simu.mat')
    indmat = zeros(n,nboot);
    for b = 1:nboot
        inds = randsample(n,n,true);
        tab =tabulate(inds);
        indmat(tab(:,1), b) = tab(:,2);
    end
end
boot_weights = [];
for mytest =  1:performTest
    nboot = size(indmat,2);
    boot_weights = nan(nsample, nboot);
end

loop_save = 1; loop_counter = 1;

checkpoint = 0;
iter0 = 0; t0 = 0; completed = 0;
if checkpoint == 1
    if exist(nam, 'file')
        load(nam);  if completed == 1;  error('completed already'); end
    end
end

tic
for iter = (iter0+1):tot
    
    %tabulate(labs)
    
    if verbose == 1
        fprintf('%6d[%3.2f]', [sum(ms>0), alpha])
        %fprintf('%6d', iter)
        if(~mod(iter,10))
            fprintf('\n')
        end
    end
    
    % if sum(ms>0) == 1;  disp('singleton cluster'); end
    
    % step 1 and 3: update beta and inv_vartheta
    bs = reshape(nu(sub2ind(size(nu), labs, 3 - gs*2)), [n,1]); %mean
    err = lTstar - Zpart - sum(mu, 2) - bs;
    tmp = nu(:,[2,4]).^2; %sqrt(precision) becomes precision
    tmp = reshape(tmp(sub2ind(size(tmp), labs, 2 - gs)), [n,1]);
    as = repmat(tmp, [1, L]);
    
    %as = .1^2;
    for j = 1:p
        err = err + mu(:, j); %extract jth contribution
        %indnoj = 1:p; indnoj(j) = [];
        %err = lTstar - Zpart - sum(Trues.mu(:,indnoj), 2); % - bs;
        %plot(err, Trues.mu(:,j)+Trues.eps,'k.')
        Sigma = Xt(:,:,j).*as; % PX
        %         plot(Xt(:,:,j)*betas(:,j), Trues.mu(:,j),'k.')
        Sigma = Sigma'; %X'P
        Mu = Sigma* err;
        Sigma = Sigma*Xt(:,:,j) + Pj*inv_vartheta(j);
        Sigma = chol(Sigma, 'lower');
        Mu = Sigma\Mu;
        %         plot(Trues.mu(:,1), Xt(:,:,j)*betas(:,j), 'k.')
        %         plot(err, Xt(:,:,j)*betas(:,j)+Trues.eps,'k.')
        %         err = Xt(:,:,j)*betas(:,j) + 1e-4*randn([n,1]);
        %         disp((Xt(:,:,j)'*Xt(:,:,j))\(Xt(:,:,j)'*err))
        betas(:,j) = Sigma'\( randn(size(Mu)) + Mu );
        mu(:,j) = Xt(:,:,j)*betas(:,j);
        err = err - mu(:,j); %insert jth new contribution back
        inv_vartheta(j) = gamrnd( A_v(j),   1/(B_v(j)+0.5*sum(betas(:,j).*( Pj*betas(:,j) ))) );
    end
    
    %     % visualize Beta
    %     h=figure(1);
    %     for j = 1:p
    %         fu = betas(:,j)'*B';
    %         [~,I] = sort(U);
    %         subplot(2,4,j); plot(U(I),fu(I)', 'k.-')
    %         title(Xnams(j)) %Xnams = {'Intercept', 'Grade 1','Grade 2','Stage I','Stage II','Stage III','Black race'};
    %     end
    
    % step 2: update gamma
    if size(Z,1) > 1 %if there is no constant beta, code Z as a scalar
        err = err + Zpart; %extract Z contribution
        Sigma = Z.*repmat(tmp, [1,q]); % GZ ; Note we need q < p
        Sigma = Sigma'; %Z'G
        Mu = Sigma* err + gamma_mu;
        Sigma = Sigma*Z + gamma_pre;
        Sigma = chol(Sigma, 'lower');
        Mu = Sigma\Mu;
        gammas = Sigma'\( randn(size(Mu)) + Mu );
        Zpart = Z*gammas;
        err = err - Zpart; %insert Z's new contribution
    end
    
    % nclabs = labs(ncind);
    tmp =reshape(nu(sub2ind(size(nu), labs, 4 - gs*2)), [n,1]); %sqrt-precision
    
    err = - err + lTstar; % this is "Yhat"
    
    % step 4: sample censored Ttar
    R = rand([L_cind,1]);
    lTstar(cind) = err(cind) + 1./tmp(cind).*norminv( (1-R).*...
        min(0.9999, normcdf((log(V(cind)) - err(cind)).*tmp(cind))) + R );
    
    %     %scrsz = get(0,'ScreenSize');
    %     %width = .44; height = .34;
    %     h = figure(1);
    %     %set(h,'outerPosition',[1 scrsz(4)/5 scrsz(3)*width scrsz(4)*height]);
    %     %set(groot,'CurrentFigure',h);
    %     plot(log(V) , err, 'k.')
    %     getframe(h);
    
    %err(cind) = lTstar(cind) - err(cind); % becomes residual again
    % note for step 5, we now use residual from lTstar rather than log(V)
    % err = log(V) - err + bs; % exclude theta_i effect in the residual, using log(V)
    err = lTstar - err + bs; % exclude theta_i effect in the residual, using lTstar
    
    % step 5: sample gs
    qs = (0.5 - normcdf(-nu(:,3).*nu(:,4)))./( normcdf(-nu(:,1).*nu(:,2)) - normcdf(-nu(:,3).*nu(:,4)) ); % q(\upsilon): N*1
    as = qs(labs); % q(\theta)
    bs = 1-as;
    tmp = (err - nu(labs,1)).*nu(labs,2); %now tmp becomes the standardized residual
    %         as(ncind) = as(ncind).*nu(nclabs,2).*exp( -0.5*tmp(ncind).^2);
    %         as(cind) = as(cind).*normcdf( tmp(cind) );
    as = as.*nu(labs,2).*exp( -0.5*tmp.^2); % based on lTstar, no need to separate
    tmp = (err - nu(labs,3)).*nu(labs,4);
    %         bs(ncind) = bs(ncind).*nu(nclabs,4).*exp( -0.5*tmp(ncind).^2 );
    %         bs(cind) = bs(cind).*normcdf( tmp(cind) );
    bs = bs.*nu(labs,4).*exp(-0.5*tmp.^2); % based on lTstar, no need to separate s
    gs = binornd( 1, as./(as+bs) );
    %     if mean(gs)<=.01
    %         disp('here')
    %     end
    indg1 = find(gs==1); n1 = length(indg1);
    indg0 = find(gs==0); n0 = length(indg0);
    
    for myupdateAll = 1:updateAll
        
        for useFDDP = 1:FDDP % execute when FDDP == 1
            % step 6 b: sample ps
            ps = gamrnd(alpha/N + ms, 1);
            ps = ps'./sum(ps);
            
            % step 6 c: block sample  -- 2015-12-9
            ind0 = find(ms==0); d0 = numel(ind0);
            ind1 = find(ms>0); d1 = numel(ind1);
            nu(ind0, 2) = sqrt(gamrnd( A_nu(1), 1/B_nu(1), [d0, 1] )); %sqrt precision;
            nu(ind0, 4) = sqrt(gamrnd( A_nu(2), 1/B_nu(2), [d0, 1] )); %sqrt precision;
            nu(ind0, [1,3]) = exp(randn([d0,2])./(nu(ind0, [2,4]).*repmat(invtau,[d0,1])) + ...
                repmat(nu_mean, [d0,1])) .* repmat([-1, 1], [d0,1]);
            nu_accepts(ind0,:) = nu_accepts(ind0,:)+1; %always accept
            
            as = qs(labs); % q(\theta)
            bs = as;
            if n1 > 0
                labg1 = labs(indg1);
                tmp = (err(indg1) - nu(labg1,1)).*nu(labg1,2);
                bs(indg1) = log(as(indg1)) + log(nu(labg1,2)) -0.5*tmp.^2; % based on lTstar, no need to separate censored/noncensored
            end
            if n0 > 0
                labg0 = labs(indg0);
                tmp = (err(indg0) - nu(labg0,3)).*nu(labg0,4);
                bs(indg0) = log(1-bs(indg0)) + log(nu(labg0,4)) - 0.5*tmp.^2; % bs now becomes the log-likelihood
            end
            
            lik =  accumarray(labs, bs); %prod(bs(labs==j));
            for k = 1:4
                nu_new = nu; labg1 = labs(indg1); labg0 = labs(indg0);
                lik1 = zeros(N,1);
                switch k
                    case 1
                        nu_new(ind1,k) = -exp( exp(nu_tunings(ind1,k)).*randn([d1,1]) + log(-nu(ind1,k)) ); % propose and transform
                        ratio0 = - 0.5*((nu_new(ind1,k)-nu_mean(1)).*nu_new(ind1,2)*invtau(1)).^2 + 0.5*((nu(ind1,k)-nu_mean(1)).*nu(ind1,2)*invtau(1)).^2;
                    case 2 % sqrt(precision)
                        nu_new(ind1,k) = 1./sqrt(exp( exp(nu_tunings(ind1,k)).*randn([d1,1]) - 2*log(nu(ind1,k)) ));
                        ratio0 = - 0.5*((nu_new(ind1,1) - nu_mean(1)).*nu_new(ind1,k)*invtau(1)).^2 + ...
                            (2*(A_nu(1)+1)+1)*log(nu_new(ind1,k)) -nu_new(ind1,k).^2*B_nu(1) - 2*log(nu_new(ind1,k)) - ...
                            (     - 0.5*((nu(ind1,1) - nu_mean(1)).*nu(ind1,k)*invtau(1)).^2 + ...
                            (2*(A_nu(1)+1)+1)*log(nu(ind1,k)) -nu(ind1,k).^2*B_nu(1) - 2*log(nu(ind1,k))    );
                    case 3
                        nu_new(ind1,k) = exp( exp(nu_tunings(ind1,k)).*randn([d1,1]) + log(nu(ind1,k)) );
                        %disp([ nu_new(ind1,k) ,  nu(ind1,k) ])
                        ratio0 = - 0.5*((nu_new(ind1,k)-nu_mean(2)).*nu_new(ind1,4)*invtau(2)).^2 + 0.5*((nu(ind1,k)-nu_mean(2)).*nu(ind1,4)*invtau(2)).^2;
                    case 4 % sqrt(precision)
                        nu_new(ind1,k) = 1./sqrt(exp( exp(nu_tunings(ind1,k)).*randn([d1,1]) - 2*log(nu(ind1,k)) ));
                        ratio0 = - 0.5*((nu_new(ind1,3)-nu_mean(2)).*nu_new(ind1,k)*invtau(2)).^2 + ...
                            (2*(A_nu(2)+1)+1)*log(nu_new(ind1,k)) - nu_new(ind1,k).^2*B_nu(2) - 2*log(nu_new(ind1,k)) - ...
                            (    - 0.5*((nu(ind1,3)-nu_mean(2)).*nu(ind1,k)*invtau(2)).^2 + ...
                            (2*(A_nu(2)+1)+1)*log(nu(ind1,k)) - nu(ind1,k).^2*B_nu(2) - 2*log(nu(ind1,k))   );
                end
                % if isnan(nu_new(ind1,k)); disp('here'); end
                as = zeros(N,1);
                as(ind1) = (0.5 - normcdf(-nu_new(ind1,3).*nu_new(ind1,4)))./( normcdf(-nu_new(ind1,1).*nu_new(ind1,2)) - normcdf(-nu_new(ind1,3).*nu_new(ind1,4)) ); % q(\theta)'s
                % note ms(i) = sum(labs==j) must hold for the following
                if n1 > 0
                    tmp = accumarray(labg1, log(as(labg1).*nu_new(labg1,2)) -0.5*((err(indg1) - nu_new(labg1,1)).*nu_new(labg1,2)).^2 ); %num2str([tabulate(labg1), tmp],3)
                    tab = tabulate(labg1);
                    lik1(tab(:,2) ~= 0) = lik1(tab(:,2) ~= 0) + tmp(tmp~=0);
                end
                if n0 > 0
                    tmp = accumarray(labg0, log((1-as(labg0)).*nu_new(labg0,4)) -0.5*((err(indg0) - nu_new(labg0,3)).*nu_new(labg0,4)).^2 ); %num2str([tabulate(labg1), tmp],3)
                    tab = tabulate(labg0);
                    lik1(tab(:,2) ~= 0) = lik1(tab(:,2) ~= 0) + tmp(tmp~=0);
                end
                MH = lik1(ind1,:) - lik(ind1,:) + ratio0;
                inds = ind1(log(rand([d1,1])) <= MH);
                %disp([nu(inds,k), nu_new(inds,k)])
                nu(inds,k) = nu_new(inds,k);
                nu_accepts(inds,k) = nu_accepts(inds,k)+1;
                lik(inds,:) = lik1(inds,:);
            end
            
            % step 6 d & e: sample invtau
            invtau(1) = sqrt(gamrnd( A_tau(1)+d1, 1/(B_tau(1)+0.5*sum( ((nu(ind1,1) - nu_mean(1)).*nu(ind1,2)).^2)) ));
            invtau(2) = sqrt(gamrnd( A_tau(2)+d1, 1/(B_tau(2)+0.5*sum( ((nu(ind1,3) - nu_mean(2)).*nu(ind1,4)).^2)) ));
            
            % step 7: sample alpha
            sps = sum(log(ps));
            if useMALA4alpha == 1 % proposal using gradient information
                grad = alpha*(psi(alpha) - psi(alpha/N) + sps/N  + (1/alpha- 1) ); %1/alpha- 1 from d(log-prior) with Exp(1)
                alpha_t = log(alpha);
                alpha1_t = alpha_t + 0.5*h_alpha*grad + sqrt(h_alpha)*randn(1);
                alpha1 = exp(alpha1_t);
                grad1 = alpha1*(psi(alpha1) - psi(alpha1/N) + sps/N + (1/alpha1- 1) );
            else
                alpha1 = exprnd(1); % proposal density \pi(\alpha)
            end
            loglik =gammaln(alpha) - N*gammaln(alpha/N) +(alpha/N-1)*sps;
            loglik1 =gammaln(alpha1) - N*gammaln(alpha1/N) +(alpha1/N-1)*sps;
            MH = loglik1 - loglik;
            if useMALA4alpha == 1
                MH = MH - (alpha1-alpha) - ...
                    0.25*(grad1+grad)*(2*(alpha_t - alpha1_t) + 0.5*h_alpha*(grad-grad1));
            end
            if log(rand(1)) <= MH;  alpha = alpha1;  alpha_accepts(iter) = 1;  end %accept
            
            
            % tune MALA as suggested by Marshall & Robterts (2012),
            if iter > alpha_M
                eap = mean(alpha_accepts((iter-alpha_M+1):iter)); %empirical acceptance probability
                if reportrate == 1; disp(num2str(eap)); end
                h_alpha = h_alpha + sign((eap>=alpha_objrate)-0.5).*min(0.001*h_alpha, 1/iter); %define c(n) = 1/n
            end
            
            % step 6 a: sample labs
            % % for fast bootstrap sampling, we need to put it as last step
            as = repmat(err, [1,N]);
            bs = repmat(nu(:,2)', [n1,1]);
            %bs_log = repmat(log(nu(:,2))', [n1,1]);
            if n1 > 0
                %as(indg1,:) = repmat(log(qs)', [n1,1]) + bs_log - 0.5*((as(indg1,:) - repmat(nu(:,1)', [n1,1])).*bs).^2;
                as(indg1,:) = repmat(qs', [n1,1]).*bs.*exp(-0.5*((as(indg1,:) - repmat(nu(:,1)', [n1,1])).*bs).^2);
            end
            bs = repmat(nu(:,4)', [n0,1]);
            %bs_log = repmat(log(nu(:,4))', [n0,1]);
            if n0 > 0
                %as(indg0,:) = repmat(log(1-qs'), [n0,1]) + bs_log - 0.5*((as(indg0,:) - repmat(nu(:,3)', [n0,1])).*bs).^2;
                as(indg0,:) = repmat(1-qs', [n0,1]).*bs.*exp(-0.5*((as(indg0,:) - repmat(nu(:,3)', [n0,1])).*bs).^2);
            end
            %ps0 = repmat(log(ps), [n,1]) + as;
            
            %             % didn't improve it
            %             tic
            %             tmp = [1-qs'; qs'; nu(:,3)'; nu(:,1)';nu(:,4)';nu(:,2)']; %6*N matrix
            %             as1 = tmp(gs+1,:).*tmp(gs+5,:).*exp(-0.5*((repmat(err, [1,N]) - tmp(gs+3,:)).*tmp(gs+5,:)).^2);
            %             toc
            
            ps0 = repmat(ps, [n,1]).*as; %needed for calculating the importance sampling
            %maxp = max(ps0,[],2);
            %ps0 = exp(ps0-repmat(maxp,[1,N]))./repmat(sum(exp(ps0-repmat(maxp,[1,N])),2), [1,N]);
            %ps0 = cumsum(ps0./repmat(sum(ps0,2), [1,N]), 2);
            
            if fixcluster ~= 1
                labs = sum(cumsum(ps0./repmat(sum(ps0,2), [1,N]), 2) <= repmat(rand([n,1]), [1,N]), 2) + 1;
            end
            ms = histc(labs, 1:N);
            
        end
        
    end
    
    %     vec = [mean(betas,1), gammas', alpha];
    %     nvec = numel(vec);
    %     fprintf(strcat(repmat('%3.4f ',[1,nvec]),'\n'), vec)
    
    if iter > burn && ~mod(loop_counter, thin)
        % store results
        matpara(loop_save, :) = [reshape(betas, [numel(betas),1]); gammas; alpha];
        %         Ms{1, iter-burn} = ms;
        Ms(loop_save) = sum(ms>0);
        %Theta(:,:,loop_save) = nu(1,:);
        %         Ps(loop_save,:) = ps;
        
        % calculate the bootstrap estimates
        %ps0 = prod((repmat(sum(ps0, 2), [1,nboot])).^(indmat-1), 1);
        %boot_weights = boot_weights + [ps0; repmat(betas,[1,nboot]).*repmat(ps0, [p,1])];
        
        for mytest =  1:performTest
            boot_weights(loop_save, :) = sum((indmat-1).*repmat(log(sum(ps0, 2)), [1, nboot]), 1);
        end
        
        loop_save = loop_save + 1;
        loop_counter = 0;
        
        %         % calculate Bayes factor
        %         snu = sqrt(nu(cind));
        %         tmp =  Z(cind,:)*beta' + snu.*fulltheta(cind);
        %         loglik = sum(  1 - normcdf((log(V(cind)) - tmp)./snu)  );
        %         % note after step 5, err(ncind) is unchanged.
        %         if fixgamma0 == 0
        %             loglik = loglik + sum(...
        %                 - 0.5*(err(ncind) - fulltheta(ncind)).^2 - 0.5*lognu(ncind) - 0.5*log(2*pi) ...
        %                 );
        %             loglik = loglik + sum(...
        %                 - (eta-mu).^2./(2*tau2_eta) - 0.5*log(2*pi*tau2_eta) ...
        %                 );
        %         else
        %             loglik = loglik + sum(...
        %                 - 0.5*(err(ncind) - fulltheta(ncind)).^2 - lognu(ncind)/2 - 0.5*log(2*pi) ...
        %                 );
        %         end
        %         BF(iter-burn) = loglik;
    end
    
    if ~mod(iter,batchLen)
        batchNum = batchNum+1;
        nu_accepts = nu_accepts./batchLen;
        %         alpha_accepts = alpha_accepts/batchLen;
        nu_rates(batchNum,:) = mean(nu_accepts, 1);
        if reportrate == 1 && iter > burn
            if sum(ms>0) > 1
                disp(num2str(mean(nu_accepts(ms>0,:)), 2))  %alpha_accepts
            else
                disp(num2str(nu_accepts(ms>0,:), 2))
            end
        end
        %         if nu_accepts(1)==1
        %             disp('here')
        %         end
        nu_tunings = nu_tunings + sign((nu_accepts>objrate)-0.5).*min(0.01, 1/sqrt(batchNum));
        nu_accepts = zeros(N,nlen);
        %         alpha_accepts = 0;
    end
    
    loop_counter = loop_counter + 1;
    
    if toc > 14000 && checkpoint == 1
        iter0 = iter; CPUtime = toc; t0 = t0 + CPUtime/60;
        disp(iter0)
        save(nam,'matpara','Ms','iter0','t0','batchNum','nu_tunings','nu_accepts','alpha_accepts',...
            'betas', 'gammas', 'alpha','ps','ms','nu_rates','lTstar','nu','gs','labs','completed',...
            'boot_weights','loop_save','loop_counter')
        checkpoint = 0;
    end
    
end

runtime = t0 + toc/60;
completed = 1;
fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', tot, runtime)
disp(completed)
if simu == 0
    save(nam,'matpara','completed','runtime','boot_weights','Ms') %,'Theta'
end
end

% not run

