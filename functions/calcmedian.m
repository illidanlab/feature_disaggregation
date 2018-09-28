imp = result{10,7};
imp = cat(1,imp{:});
l = size(result{3,7},1);
chlaresult = zeros(l,5);
chlaresult(:,1) = imp;
chlaresult(:,2) = result{2,7};
chlaresult(:,3) = result{4,7};
chlaresult(:,4) = result{7,7};
chlaresult(:,5) = Ytstall;

absresult = zeros(l,4);

% chlaresult = exp(chlaresult)
% absresult(:,1) = abs(chlaresult(:,1)-chlaresult(:,5)./chlaresult(:,5));
% absresult(:,2) = abs(chlaresult(:,2)-chlaresult(:,5)./chlaresult(:,5));
% absresult(:,3) = abs(chlaresult(:,3)-chlaresult(:,5)./chlaresult(:,5));
% absresult(:,4) = abs(chlaresult(:,4)-chlaresult(:,5)./chlaresult(:,5));
absresultall = [absresultall;chlaresult];
% boxplot(absresult)
% 
% a = median(absresult,1)
% 
% b = mean(absresult,1)
% 
% scatter(chlasresult(:,1),chlaresult(:,5))
% [~,rmse] = rsquare(absresultall(:,1),absresultall(:,5))
% [~,rmse] = rsquare(absresultall(:,2),absresultall(:,5))
% [~,rmse] = rsquare(absresultall(:,3),absresultall(:,5))
% [~,rmse] = rsquare(absresultall(:,4),absresultall(:,5))
% 
expabsresultall = (absresultall)
median(abs(expabsresultall(:,1)-expabsresultall(:,5)./expabsresultall(:,5)))
median(abs(expabsresultall(:,2)-expabsresultall(:,5)./expabsresultall(:,5)))
median(abs(expabsresultall(:,3)-expabsresultall(:,5)./expabsresultall(:,5)))
median(abs(expabsresultall(:,4)-expabsresultall(:,5)./expabsresultall(:,5)))
error = expabsresultall(:,1:4) - repmat(expabsresultall(:,5),[1,4])
error = abs(error)
median(error)
median(abs(error./repmat(expabsresultall(:,5),[1,4])))

boxplot(error)
% ylim([0,20])

