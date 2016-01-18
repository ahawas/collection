filename = 'C:\Users\hawas1\Documents\MFE\Emperical\ass2\datahwk2_problem111.txt'
fid = fopen(filename);
out = textscan(fid,'%s %f %f %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fid);
hedge_returns = zeros(length(out{2}) , 8);
for i = 1:length(out{2})
    for j = 2:9
        hedge_returns(i,j-1) = out{j}(i);
    end
end
hedge_fund_dates = datetime(out{1}, 'inputformat', 'MM/dd/yyyy');
filename = 'C:\Users\hawas1\Documents\MFE\Emperical\ass2\datahwk2_problem2.csv'
fid = fopen(filename);
out_macro = textscan(fid,'%s %f %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fid);

snp = out_macro{7};
bond = out_macro{8};
usd = out_macro{3};
credit = out_macro{4};
cmdty = out_macro{5};
macro_dates = dateshift(datetime(out_macro{1}, 'inputformat', 'MM/dd/yyyy'), 'start',  'month', 'next');
indices_to_work_with = zeros();
p = 1;
x= [snp, bond, usd, credit, cmdty];
x=x(85:end-1,:);
rts = zeros(3,8);
betas = zeros(5,8);
lev_acfs = zeros(13, 8);
hedge_acfs = zeros(13,8);
ljung_decisions = zeros(1,8);
levys = zeros(202, 8);
for i =  1:8
   warning('off', 'optimlib:lsqlin:LinConstraints')
   y = hedge_returns(:,i)/100;
   options = optimoptions('lsqlin','Algorithm','active-set');
   [b,resnorm,residual,exitflag,output,lambda] = lsqlin(x,y, [],[],ones(1,5),1, [], [], [], options);
   %disp(strcat('Hedge #',int2str(i)))
   betas(:, i) = b;
   yhat = x*b;
   lev = sqrt( (1/(length(y)-1)) * sum( (y - mean(y)) .^ 2 ) );
   lev = lev / sqrt( (1/(length(y)-1)) * sum( (yhat - mean(yhat)) .^ 2 ) );
   levy = lev*yhat;
   rts(1,i) = mean(yhat);
   rts(2, i ) = mean(levy);
   rts(3, i )=  mean(y);
   %figure; autocorr(levy,12);
   levys(:,i) = levy;
   lev_acfs(:,i) = autocorr(levy, 12);
   hedge_acfs(:,i) = autocorr(y, 12);
   [h,pvalue] = lbqtest(levy, 'lags', [12]);
   ljung_decisions(1,i) = h;
  % pvalue,h
end
%%
% *Betas of factors for each hedge fund*
disp(betas)
%%
% *Average returns of Rstar, Rhat, and actual R for each hedge fund*
disp(rts)
%%
% *ACF values for Rhat, up to lag 12*
disp(lev_acfs)
%%
% *ACF values for actual returns, up to lag 12*
disp('ACF values for actual returns, up to lag 12')
disp(hedge_acfs);
%%
% *ACF Plot for leverage clones*
for i = 1:8
    figure;autocorr(levys(:,i), 12);
end
%%
% *Ljung decisions ( 1 per hedge fund, lag 12 )*
disp('Ljung decisions (1 per hedge fund, lag 12)');
disp(ljung_decisions);

%%
% *5)
 % It is clear that hedge funds do face an issue due to liquidity risk.
 % Through this example, we can see that it is possible to hedge some of
 % this risk by using more liquid assets such as the S&P500 index, Bonds,
  % USD, Commodities and Credit. We can see the liquidity risk is less since
 % most of the ACF values are decreased in our leveraged clone.
 % Nonetheless, we must also keep in mind that we can also notice a
 % substantial decrease in average returns than in the actual fund.* 


