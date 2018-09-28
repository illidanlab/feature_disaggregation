Xreal = cat(1,X_R_local{:});
Xpred = cat(1,Xlr{:});
Xmean = zeros(300,1);
ind = 1;
for i = 1:r
    Xmean(ind:ind+29,1)=X_R(i);
    ind = ind+30;
end

% set(gca,'FontSize',50)
x = linspace(-5,5);
y = linspace(-5,5);
plot(x,y);
hold on;
xlabel('localized regional varible recovery','FontSize', 20)
ylabel('true localized regional varible','FontSize', 20)
scatter(Xpred,Xreal,100);
set(gca,'FontSize',20)

figure;
x = linspace(-5,5);
y = linspace(-5,5);
plot(x,y);
hold on;
scatter(Xmean,Xreal,100);
xlabel('low resolution regional data','FontSize', 20)
ylabel('true localized regional varible','FontSize', 20)
set(gca,'FontSize',20)