function [] = runMRV(jobID)
% global simu V Delta Z Xt X Trues filenames B U pilotRun indmat performTest %#ok<*NUSED>

simu = 0;  % simulation or real data?
performTest = 1; % perform boostrap testing?
runHPC = 0;  if simu == 0; runHPC = 1; end  %runHPC=1 for real data
olderversionMATLAB = 0;
arrayID = str2double(num2str(jobID));  %simulation id
filenames = strcat('out',num2str(arrayID),'.mat');  % name of output file

% for simulation
for myloop = 1:(1-runHPC) %equivalant to "if runHPC == 0", easier to fold the code
    chs = 1:3; nch = length(chs);
    if simu == 1; load(strcat('sdat',num2str(arrayID),'.mat')); else   load('datMRV.mat'); end
    if olderversionMATLAB == 0; rng('default'); rng(arrayID*10); end
    q = size(Z, 2);  p = size(X,2);  [~, L] = size(B);
    % else rand('state', nvar*3);  randn('state', nvar*10);   end
    % construct bootstrap samples
    % boot_weightsAll = nan(tot, nboot);
    for ch = 1:nch
        [matpara, boot_weights] = MedianRegVary(V, Delta, Z, X, B);
        if ch==1  %initialize
            nsample = size(matpara,1);
            tot = nch*nsample;
            nboot = size(boot_weights,2);
            MCSamples = nan(nch, size(matpara,2),  nsample);  %used for psrf to check convergence
            matparaAll = nan(tot,  size(matpara,2));
            if performTest==1;  boot_weightsAll = nan(nch*nsample, nboot);  end
        end
        MCSamples(ch, :, :) = matpara';
        matparaAll((ch-1)*nsample+(1:nsample),:) = matpara;
        if performTest==1;   boot_weightsAll((ch-1)*nsample+(1:nsample),:) = boot_weights;  end
        %         MsAll((ch-1)*nsample+(1:nsample),:) = Ms;
        %         PsAll((ch-1)*nsample+(1:nsample),:) = Ps;
        %         ThetaAll((ch-1)*nsample+(1:nsample),:) = Theta;
    end
    R = psrf(MCSamples); %disp(R)  %check convergence
    %  disp(mean(matparaAll,1))
    %mat = matparaAll(:,1:(L*p));
    %save('mat_10.mat','mat','L','p')
    
    pvals1 = []; Ts1 = [];
    tot = size(boot_weightsAll,1);
    for mytest =  1:performTest
        [pvals1, Ts1] = getPval( matparaAll(:,1:(L*p)) , L, p);   disp(pvals1)%num2str(pvals1, 4)
        boot_weightsAll = exp(boot_weightsAll - repmat(max(boot_weightsAll,[],1), [tot,1]));
    end
    
    mBeta = reshape(mean(matparaAll(:,1:(L*p)), 1), [L, p]); % posterior mean
    % mBeta = reshape(quantile(matparaAll(:,1:(L*p)), .5, 1), [L, p]);
    lBeta = reshape(quantile(matparaAll(:,1:(L*p)), .025, 1), [L, p]);
    uBeta = reshape(quantile(matparaAll(:,1:(L*p)), .975, 1), [L, p]);
    
    mGamma = quantile(matparaAll(:,(L*p+1):end), [.5, .025,.975], 1);
    mGamma2 = mean(matparaAll(:,(L*p+1):end), 1);
    %disp(mGamma')
    
    pvals = []; Ts = zeros(1,p);
    for mytest =  1:performTest
        denom = sum(boot_weightsAll, 1);
        Bmat = nan(L*p, nboot);
        for j = 1:(L*p)
            Bmat(j,:) = sum(repmat(matparaAll(:,j), [1,nboot]) .* boot_weightsAll, 1)./ denom;
        end
        mBeta_b1 = reshape(mean(Bmat, 2), [L, p]);
        se = reshape(std(Bmat,[], 2), [L, p]);
        lBeta_b1 = mBeta_b1 - se*norminv(.975);
        uBeta_b1 = mBeta_b1 + se*norminv(.975);
        lBeta_b2 = reshape(quantile(Bmat, .025, 2), [L, p]);
        uBeta_b2 = reshape(quantile(Bmat, .975, 2), [L, p]);
        
        pvals = nan(1,p);
        for j = 1:p
            Sigma = zeros(L); b_bar = mean(Bmat((j-1)*L+(1:L),:),2);
            for i = 1:nboot
                b_jr = Bmat((j-1)*L+(1:L),i) - b_bar;
                Sigma = Sigma + b_jr*b_jr';
            end
            Sigma = Sigma/(nboot-1);
            b_hat = mBeta(:,j);
            F = [ones(L-1,1), -eye(L-1)];
            T = (F*b_hat)' * ((F*Sigma*F')\(F*b_hat)); %chisq-(L-1)
            pvals(j) =  1- chi2cdf(T, L-1) ;
            Ts(j) = T;
        end
        %num2str(pvals, 4)
        fprintf('case=%d, pvals = %3.4f, %3.4f, %3.4f.\n', [arrayID, pvals])
    end
    
    fu = cell(1,p);
    plotit = 1;
    for j = 1:p
        fu{j} = [U'; mBeta(:,j)'*B'; lBeta(:,j)'*B'; uBeta(:,j)'*B']; %; mBeta_b1(:,j)'*B'; lBeta_b1(:,j)'*B'; uBeta_b1(:,j)'*B'; lBeta_b2(:,j)'*B'; uBeta_b2(:,j)'*B'
        [~,I] = sort(U);
        fu{j} = fu{j}(:,I)';
        fu{j} = unique(fu{j}, 'rows');
        for myplot = 1:plotit
            subplot(1,p,j); h = plot(repmat(fu{j}(:,1),[1,size(fu{j},2)-1]), fu{j}(:,2:end), '.-','MarkerSize',12); %:end
            %set(h,'Color',.8*ones(1,3))
            %title(Xnams(j))
            hold on
            switch j
                case 1
                    h = ezplot(@(x) x.^3, [0,1]);
                case 2
                    h = ezplot(@(x) sin(2*pi*x), [0,1]);
                case 3
                    h = ezplot(@(x) 1, [0,1]);
            end
            title(''); set(h,'Color','k','LineWidth',1)
            %set(gca,'XTick',0.1:0.15:1,'XTickLabel',10:15:100); set(gca,'FontSize',15)
            hold off
            axis tight
        end
    end
    
    W0 = [.1,.3,.5,.7, .9]; lw = numel(W0);
    mse = zeros(p, lw); bias = zeros(p, lw); cp = zeros(p,lw); width = zeros(p,lw);
    for j = 1:p
        tmp = [mBeta(:,j)'*Trues.B0'; lBeta(:,j)'*Trues.B0'; uBeta(:,j)'*Trues.B0'];
        switch j
            case 1
                tmp1 = W0.^3;
            case 2
                tmp1 = sin(2*pi*W0);
            case 3
                tmp1 = 1;
        end
        bias(j,:) = tmp(1,:)-tmp1;
        mse(j,:) = (tmp(1,:)-tmp1).^2;
        cp(j,:) = (tmp(2,:)<=tmp1).*(tmp(3,:)>=tmp1);
        width(j,:) = tmp(3,:) - tmp(2,:);
    end
    
    disp(mse)
    
    save(filenames,'R','mse','bias','cp','width','pvals','pvals1','mBeta','mGamma','mGamma2','Bmat','Ts')
end

% for real data
for myloop = 1:runHPC % run HPC for real data, nvar becomes chain number
    rng('default'); rng(arrayID*23)
    load('datMRV.mat')
    %     induse = find(V>0);
    %     V = V(induse,:); Z = Z(induse,:); Delta = Delta(induse,:);
    % default B: knots = c(.1,.3,.5,.7,.9)
    % try 20 knots
    load('B20.mat');  B = B20;
    MedianRegVary(V, Delta, Z, X, B); 
end

end

% not run




