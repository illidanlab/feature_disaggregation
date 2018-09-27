addpath(genpath('./MALSAR-master/MALSAR/functions/')); % load function
addpath('./MALSAR-master/MALSAR/utils/'); % load utilities
addpath('./functions/'); % load utilities
addpath('./cvx/'); % load utilities
rng(2);
N = 20000;
%rcombo=[10];
rcombo=[2,4,8,10,20,40,100,200,500,800,1000,2000,5000];
rmsecombo = zeros(length(rcombo),1);
rmsesign = zeros(length(rcombo),1);
dl = 13; % 5 local feature
dr = 2; % 1 region feature for visualization
for iter = 1:length(rcombo)
    r = rcombo(iter);
%    r = 1600;% 5 region

    n = N/r;% 10 samples in each region


    coord = cell(1,r);%spatial coordinates
    center = randn(r,2)*10;
    std = 1;
    for i = 1:r
        coord{i}(:,1) = normrnd(center(i,1),std, [n,1]);
        coord{i}(:,2) = normrnd(center(i,2),std, [n,1]);
    end
    %generate local variable

    X_L = cell(1,r);
    for i = 1:r
        X_L{i} = rand(n,dl); 
    end
    %generate regional variable
    X_R_local = cell(1,r);

    %generate regional variable using coordinates information
    rng('shuffle');
    X_R = randn(r,dr);
    %normalize the data
    clear std ind;
    tmp1 = cat(1,X_L{:});
    tmp1 = (tmp1 - repmat(mean(tmp1),[n*r,1]))./repmat(std(tmp1),[n*r,1]);
    %generate cross term
    X_csi = createX_csi_local(X_L,X_R_local);

    X_csi_R = createX_csi(X_L,X_R);
    %local coefficient,region coefficient and vectorized CSI
    rng(3);
    alphatrue = randn(dl,r);
    betatrue = randn(dr,r);
    gammatrue = randn(dl*dr,1);
    % gammatrue = (abs(gammatrue) > 0.5).*gammatrue;

    % generate response variable

    Y = cell(1,r);
    rng('shuffle') 
    for i = 1:r
        %Y{i} = X_L{i}*alphatrue(:,i) + (X_R_local{i}*betatrue(:,i)) + X_csi{i}*gammatrue +  randn(n,1)./100;
        Y{i} = X_L{i}*alphatrue(:,i) + (repmat(X_R(i,:),[size(X_L{i},1),1])*betatrue(:,i)) + X_csi_R{i}*gammatrue +  randn(n,1);
    end
    Gtrue = reshape(gammatrue,[dr,dl])';
    Xtrn_csi = cell(1,r);
    for i = 1:r
    %     Xtrn_csi{i} = horzcat(ones(size(X_L{i},1),1),X_L{i},repmat(X_R(i,:),[size(X_L{i},1),1]),X_csi_R{i});
    end
    XG = cat(1,Xtrn_csi{:});
    YG = cat(1, Y{:});
    W = Least_Lasso(X_L, Y, 0);
    X_R = zscore(X_R);
    X_z = vertcat(X_R',ones(1,r));
    cvx_begin quiet
    cvx_precision high
    variable G(dl, dr+1)

    minimize(norm(W - G*X_z,'fro'))
    cvx_end
    G = G(:,1:dr)
    Gtrue
    rmsesign(iter) = CalcSign(Gtrue,G);
    rmsecombo(iter) = (norm(G-Gtrue,'fro'))/(dl*dr);
end

plot(rmsecombo)
hold on
plot(rmsesign)
