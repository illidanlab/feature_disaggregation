function [Wl,Wr,G,Xlr,Xlrtst,fval] =  MTML_L21i(y,Xl,Xltst,Xr_bar,lambda1,lambda2,lambda3,lambda4,coords,G_init,trnidx,tstidx,maxiter,toler)
    %solving MTMLa with regional variable imputation
    %input:
    %   y: cell of n*1 vector; Xl: cell of n*dl matrix, Xr_bar: r*dr matrix, 
    %   Xlr: cell of n*dr matrix, lambda1:hyperparameter, lambda2:hyperparameter
    %   lambda3:hyperparameter coords: the cell of n*2 matrix, which
    %   contains the spatial coordinates
    %formulation: 
    %minimize: sum_{1=i}^{n}||y_{i}-Xl_{i}*Wl_{i}-
    %Xlr_{i}*Wr{i}-diag(Xl{i}*G*Xr{i}')|| + lambda1*||[Wl,Wr]||_{1}+ lambda2*||G||_{1} +
    %lambda3*sum_{i}sum_{j}d_{ij}||Xlr{i}-Xlr{j}||_{2}^{2} +
    %lambda4*||Xlr-Xr_bar}||_{F}^{2}
    
    
    %args
    r = length(y);%number of regions
    ytmp = cell(1,r);
    for i = 1:r
        ytmp{i} = zeros(size(Xltst{i},1),1);
    end
    [~, sampletrn, ~, ~] = diagonalize(Xl, y);
    [~,sampletst,~,~] = diagonalize(Xltst,ytmp);
    
    
    dl = size(Xl{1},2);
    dr = size(Xr_bar,2);
    Xr_bar = (Xr_bar - repmat(mean(Xr_bar),[size(Xr_bar,1),1]))./repmat(std(Xr_bar),[size(Xr_bar,1),1]);
    %initialize the Xlr using the Xr_bar
    Xlrtrn_init = cell(1,r);
    Xlrtst_init = cell(1,r);
    Xlrall = cell(1,r);
    rng('shuffle')
    for i = 1:r
        ntrn = sampletrn(i);
        Xlrtrn_init{i} = repmat(Xr_bar(i,:),[ntrn,1]);
%         + randn(ntrn,dr);
        ntst = sampletst(i);
        Xlrtst_init{i} = repmat(Xr_bar(i,:),[ntst,1]);
%         + randn(ntst,dr);
        Xlrall{i} = repmat(Xr_bar(i,:),[ntrn+ntst,1]);
    end
    rng(2);
%     G_init = randn(dl,dr);
    coordG = cat(1,coords{:});
    D = squareform(pdist(coordG)); %distance matrix
    m = 1;
    for diagi = 1:r
        D(m:m+sampletrn(diagi)-1,m:m+sampletrn(diagi)-1) = D(m:m+sampletrn(diagi)-1,m:m+sampletrn(diagi)-1).*100;
        m = m + sampletrn(diagi);
    end
%     Xr_bar_cell = Xlr_init;
    Xr_bar_total = cat(1,Xlrall{:});%vertical concatenate regional mean value
%     R = sparse(dr,size(Xr_bar_total,1)*(size(Xr_bar_total,1)-1)/2);
    indR = 1;
    N = size(Xr_bar_total,1);
    ivec = zeros(N*(N-1)/2,1);
    jvec = zeros(N*(N-1)/2,1);
    vvec = zeros(N*(N-1)/2,1);
    DR = exp(-D);
    for Rrow = 1:N-1
        for Rcol = Rrow+1:N
            ivec(indR,1) = Rrow;
            jvec(indR,1) = indR;
            vvec(indR,1) = sqrt(DR(Rrow,Rcol));
            ivec(indR+1,1) = Rcol;
            jvec(indR+1,1) = indR;
            vvec(indR+1,1) = -sqrt(DR(Rrow,Rcol));
            indR = indR + 1;  
            
        end 
    end
    
    R = sparse(ivec,jvec,vvec,size(Xr_bar_total,1),size(Xr_bar_total,1)*(size(Xr_bar_total,1)-1)/2);
%     R = sparse(double(R));
    Xlr_oldtrn = Xlrtrn_init;
    Xlrall_old = Xlrall;
    G_old = G_init;
    fval = zeros(maxiter,1);
    fval(1) = primal_fval(y,Xl,Xlrall_old,Xr_bar_total,rand(1,r),rand(dl,r),rand(dr,r),G_old,lambda1,lambda2,lambda3,lambda4,D);
    %BCD starts here
    for i = 2:maxiter
        %fix Xlr, G, solving Wl,Wr;
        [W0,Wl,Wr] = solve_WlWr(y,Xl,Xlr_oldtrn,G_old,lambda1);
        %fix Wl,Wr,Xlr, solving G;
        W0_old = W0;
        Wl_old = Wl;
        Wr_old = Wr;
        G = solve_G(y,Xl,Xlr_oldtrn,W0_old,Wl_old,Wr_old,lambda2);
        [Xlr_total,~] = Xlrsolver_FISTA(y,Xl,G,R,Wl_old,Wr_old,Xr_bar,Xlrall_old,lambda3,lambda4,trnidx,500,0.01);
        Xlrall_old = cell(1,r);
        ind = 1;
        for j = 1:r
            Xlrall_old{j} = Xlr_total(ind:ind+sampletrn(j)+sampletst(j) -1,:);
            ind = ind + sampletrn(j)+sampletst(j);
        end
        G_old = G;
        fval(i,1) = primal_fval(y,Xl,Xlrall_old,Xr_bar_total,Wl,Wr,G,lambda1,lambda2,lambda3,lambda4,D);
        Xlr_oldtrn = cell(1,r);
        for regionid = 1:r
            Xlr_oldtrn{regionid} = Xlrall_old{regionid}(trnidx{regionid},:);
        end
        if abs(fval(i)- fval(i-1))/abs(fval(i-1)) < toler
            break;
        end        
    end
    Xlr = Xlrall_old;
    Xlrtst = cell(1,r);
    ind = 1;
    for j = 1:r
        Xlrtst{j} = Xlrall_old{j}(tstidx{j},:);
    end
    
    
% nested function
    %% 
    function [fval] = primal_fval(y,Xl,Xlrall,Xr_bar_total,Wl,Wr,G,lambda1,lambda2,lambda3,lambda4,~)
        fval = 0;
        Xlr = cell(1,r);
        for regionidx = 1:r
            Xlr{regionidx} = Xlrall{regionidx}(trnidx{regionidx},:);
        end
        Xlr_total = cat(1,Xlrall{:});
        %clustering loss
%         conttable = zeros(50*49/2,1)
        contRegu = norm((Xlr_total)'*R,'fro')^2;
        discRegu = norm(Xr_bar_total-Xlr_total,'fro')^2;%loss of variance within one region
        for region = 1:r
            rmse = norm(y{region}-Xl{region}*Wl(:,region) - Xlr{region}*Wr(:,region)-diag(Xl{region}*G*Xlr{region}'))^2;
            fval = fval + rmse  + lambda1 * sum(sqrt(sum(vertcat(Wl,Wr).^2, 2))) + lambda2 * l1_mat(G) + lambda3 * contRegu ...,
                    + lambda4 * discRegu; 
        end
      
    end

    function [Wl,Wr] = solve_WlWr(y,Xl,Xlr_old,G_old,lambda1)
        %minimize sum_||yi-Xl*G*Xlr'-[Xl,Xlr]*[Wl(:,i),Wr(:,i)]|| +
        %lambda1|[Xl,Xlr]|_{1}
        yhat = cell(1,r);
        X_all = cell(1,r);
        for region = 1:r
            yhat{region} = y{region} - diag(Xl{region}*G_old*Xlr_old{region}');
            X_all{region} = horzcat(Xl{region},Xlr_old{region});
        end
        W = Least_L21(X_all,yhat,lambda1);
        Wl = W(1:dl,:);
        Wr = W(dl+1:end,:);
    end
    
    function G = solve_G(y,Xl,Xlr_old,W0_old,Wl_old,Wr_old,lambda2)
        %solve global regression for G
        y_all = cell(1,r);
        for region = 1:r
            y_all{region} = y{region} - W0_old(region) - Xl{region}*Wl_old(:,region) - Xlr_old{region}*Wr_old(:,region);
        end
        y_all = cat(1,y{:});
        X_cross = createX_csi_local(Xl,Xlr_old);
        X_all = cat(1,X_cross{:});
        G_flat = lasso(X_all,y_all,'Lambda',lambda2);
%         G_flat = lasso(X_all,y_all,'Lambda',lambda2);
        G = reshape(G_flat,[dr,dl])';
    end
end