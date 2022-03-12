function [pvals, Ts] = getPval(mat, L, p)

bw_fac = 1; %[0.9, 0.95, 1, 1.05, 1.1]; 
bw_len = numel(bw_fac);
pvals = nan(2,p); Ts = nan(2,p); 
for nvar = 1:p
    %disp(nvar)
    x = mat(:,(nvar-1)*L+(1:L))';
    [d, n] = size(x);
    
    % sample covariance
    Sigma = x - mean(x,2)*ones(1,n);
    Sigma = Sigma*Sigma'/(n-1);
    [U,S,V] = svd(Sigma);
    
    Lo = U*(sqrt(S\eye(size(S,1))))*V';
    h0 = (4/((2*d+1)*n))^(1/(d+4));  %Hwang et al. (1994)

    %     norm(U-V)
    %     norm(Sigma - U*S*V') %U = V = B in the paper
    
    % % use kde, results not reasonable
    % x = U*x; %transformed space
    
    %     % using kde toolbox
    %     p1 = kde(x, 'rot');
    %     % p1 = kde( x, 'hall' ); % Plug-in type estimator (estimates each dim. separately)
    %     bw = getBW(p1); bw = bw(:,1);
    %     disp(std(x(1,:))*n^(-1/(d+4))) %check, = bw(1)
    %     [modes, ~, ~] = mymodes(p1); %open('kde/modes.m')
    %     % (-Hs{1}/n)\eye(L) %unreasonable
    %     x0 = modes(:,1);
    %     bw = std(x, [], 2)*n^(-1/(d+4));
    
    bw = h0*ones(d,1); %*1.181;
    
    % initial value
    % x0 = mean(x, 2); %posterior mean
    % use midpoint of the bin that corresponds to the peak of the histogram
    x_midPeak = zeros(d,1);
    for  j = 1:d
        [counts,centers] = hist(x(j,:), 50);
        ind = find(counts==max(counts));
        x_midPeak(j) = centers(ind(1));
    end
    
    auxfunc = @(nu) myGradient(nu, x, Lo, bw);
    options = optimoptions('fsolve','Display','off','MaxFunEvals',10000,'MaxIter',10000,'Algorithm','levenberg-marquardt'); %
    [x_postMode, ~,exitflag,~] = fsolve(auxfunc, mean(x, 2), options);
    
    if exitflag < 1
        x_postMode = x_midPeak; % use midpoint of optimal bin from histogram
        x_postMean = mean(x, 2); % use posterior mean
    end
    
    % compute Hessian for the mode x1
    H = 0; 
    for bw_loop = 1:bw_len
        H = H + myHessian(x_postMode, x, Lo, h0*bw_fac(bw_loop)*ones(d,1) );
    end
    H = H/bw_len; 
    %all(eig(H)<0)
    Sigma1 = (-H)\eye(L) ;
    if any(eig(Sigma1) < 0)
        fprintf('%d: Covariance is not p.d.\n', nvar)
    end
    F = [ones(L-1,1), -eye(L-1)];
    T = (F*x_postMode)' * ((F*Sigma1*F')\(F*x_postMode)); %chisq-(L-1)
    pvals(1, nvar) = 1- chi2cdf(T, L-1);
    Ts(1, nvar) = T;
    
    if exitflag < 1  %if not converged, use posterior mean
        H = 0;
        for bw_loop = 1:bw_len
            H = H + myHessian(x_postMean, x, Lo, h0*bw_fac(bw_loop)*ones(d,1) );
        end
        H = H/bw_len; 
        Sigma1 = (-H)\eye(L) ;
        if any(eig(Sigma1) < 0)
            fprintf('%d: Covariance not p.d.\n', nvar)
        end
        F = [ones(L-1,1), -eye(L-1)];
        T = (F*x_postMean)' * ((F*Sigma1*F')\(F*x_postMean)); %chisq-(L-1)
        pvals(2, nvar) = 1- chi2cdf(T, L-1);
        Ts(2, nvar) = T;
    end
    
end

% disp(pvals)

end
