clear all;clc

nino34 = dlmread('nino34.long.anom.data.txt');    %1870-2017
prec = xlsread('pr_1991_2015.xls');   %(mm,year1991-2015,month)
prec2 = dlmread('pr_1901_2015.txt');   %(mm,year1901-2015,month)

fire = xlsread('Fire_Alerts.xlsx');     %(#,year2012-2017,month)

%%
%use 2012-2015 to calculate corr btwn prec and fire
fire_sub1 = fire(1:48,1,1);
prec_sub = prec(253:end,1,1);
cor_fire_prec = corr(fire_sub1,prec_sub)
r2_fire_prec = cor_fire_prec^2
[R,P] = corrcoef(fire_sub1,prec_sub)
fitobj = fit(fire_sub1,prec_sub,'poly2')
%%
figure(11)
scatter(prec_sub,fire_sub1,'LineWidth',1.2)
title('Precipitation vs. Fires')
ylabel('Fire number')
xlabel('rainfall (mm)')
set(gca,'FontSize',14);
saveas(gcf,'11.png')

fire_sub1(find(fire_sub1>50000))=-999;
prec_sub(find(fire_sub1==-999))=[];
fire_sub1(find(fire_sub1==-999))=[];

cor_fire_prec = corr(fire_sub1,prec_sub)
r2_fire_prec = cor_fire_prec^2
[R,P] = corrcoef(fire_sub1,prec_sub)


figure(10)
scatter(prec_sub,fire_sub1,'LineWidth',1.2)
title('Precipitation vs. Fires')
ylabel('Fire number')
xlabel('rainfall (mm)')
set(gca,'FontSize',14);
saveas(gcf,'10.png')
%%

%use 2012-2017/9 to calculate corr btwn nino34 and fire
fire_sub2 = fire(1:end-1,1,1);
tmp = nino34(143:148,2:13);  %(2012-2017,1-12mth)
asd = reshape(tmp',[72,1]);
nino34_sub = asd(1:69);
cor_fire_nino = corr(fire_sub2,nino34_sub)
r2_fire_nino = cor_fire_nino^2
[R,P] = corrcoef(fire_sub2,nino34_sub)

figure(12)
scatter(nino34_sub,fire_sub2,'LineWidth',1.2)
title('Nino3.4 index vs. Fires')
ylabel('Fire number')
xlabel('Nino3.4')
set(gca,'FontSize',14);
saveas(gcf,'12.png')


fire_sub2(find(fire_sub2>50000))=-999;
nino34_sub(find(fire_sub2==-999))=[];
fire_sub2(find(fire_sub2==-999))=[];

cor_fire_nino = corr(fire_sub2,nino34_sub)
r2_fire_nino = cor_fire_nino^2
[R,P] = corrcoef(fire_sub2,nino34_sub)

figure(13)
scatter(nino34_sub,fire_sub2,'LineWidth',1.2)
title('Nino3.4 index vs. Fires')
ylabel('Fire number')
xlabel('Nino3.4')
set(gca,'FontSize',14);
saveas(gcf,'13.png')
%{
%lag (okay there is no lag)
for i=1:6
    cor_fire_nino_lag = corr(fire_sub2(i:end),nino34_sub(1:end-i+1))
end
%}
%%
%use 1901-2015 to calculate corr btwn nino34 and prec
prec_sub3 = prec2(:,1);
tmp = nino34(32:146,2:13);  %(1901-2015,1-12mth)
nino34_sub3 = reshape(tmp',[1380,1]);
cor_prec_nino = corr(prec_sub3,nino34_sub3)   %only -0.1999!!

%use 1991-2015 to calculate corr btwn nino34 and prec
prec_sub4 = prec(:,1);
tmp = nino34(122:146,2:13);  %(1901-2015,1-12mth)
nino34_sub4 = reshape(tmp',[300,1]);
cor_prec_nino = corr(prec_sub4,nino34_sub4)   %-0.31
r2_prec_nino = cor_prec_nino^2
[R,P] = corrcoef(prec_sub4,nino34_sub4)

%{
%lag:nino34 lags 2 months?doesn't make sense..
for i=1:6
    cor_prec_nino_lag = corr(prec_sub3(1:end-i+1),nino34_sub3(i:end))
end
%}
%%
%plot prec time series
tmp = reshape(prec2(:,1),[12,115]);
prec2_annual = sum(tmp,1);
asd = tmp(:,2);
year = 1901:2015;

figure(1)
plot(year,prec2_annual,'LineWidth',1.5)
title('Annual precipitation 1901-2015')
xlabel('year')
ylabel('precipitation (mm)')
set(gca,'FontSize',14);
saveas(gcf,'1.png')
%%
nino34_month = mean(nino34(32:end-2,2:13),1);   
nino34_month_ano = nino34_month-mean(nino34_month);
figure(20)
prec_month = mean(tmp,2);  %mean of 1901-2015
prec_month_ano = prec_month-mean(prec_month);
%plot(1:12,prec_month_ano)
%hold on
plot(1:12,nino34_month_ano)
%hold off
%%
nino_1d = reshape(nino34(32:end-2,2:13)',[1380 1]);
prec_1d = prec2(:,1);
figure(21)
plot(nino_1d(end-48:end))
figure(22)
plot(prec_1d(end-48:end),'LineWidth',1.5,'Color','blue')
%legend('rainfall')
set(gca,'FontSize',14);
saveas(gcf,'22.png')

%%
fire2 = reshape(fire(1:60,1),[12,5]);  %2012-2016
fire2_annual = sum(fire2,1);
year2 = 2012:2016;

figure(2)
plot(fire(1:end,1),'LineWidth',1.5,'Color','red')
%plot(year2,fire2_annual)
%plot(year2,fire(1:60,1))
%datetick('x','yyyy','keeplimits')
title('Monthly forest fires & precipitation 2012-2016')
xlabel('year')
ylabel('number')
%legend('fires')
set(gca,'FontSize',14);
saveas(gcf,'2.png')
%%
%climatology rainfall distribution 1901-2015
prec_bn = prctile(prec2_annual,33.33)
prec_an = prctile(prec2_annual,66.66)
%bar(prec2_annual)

min_cnt = min(prec2_annual);
max_cnt = max(prec2_annual);
binwidth = 50;
bins = [min_cnt:binwidth:max_cnt];

figure(3)
%histogram(prec2_annual,bins,'Normalization','pdf')
histogram(prec2_annual,bins,'Normalization','probability')
title('Annual precipitation 1901-2015 normal distribution')
xlabel('rainfall (mm)')
ylabel('%')
%%
%climatology rainfall distribution 1901-1980, 80yrs
prec_bn_pre = prctile(prec2_annual(1:80),33.3)
prec_an_pre = prctile(prec2_annual(1:80),66.6)

figure(4)
histogram(prec2_annual(1:80),bins,'Normalization','probability')
title('Annual precipitation 1901-1980 normal distribution')
xlabel('rainfall (mm)')
ylabel('%')
saveas(gcf,'4.png')

%climatology rainfall distribution 1981-2015, 35yrs
prec_bn_aft = prctile(prec2_annual(81:end),33.3)
prec_an_aft = prctile(prec2_annual(81:end),66.6)

figure(5)
histogram(prec2_annual(81:end),bins,'Normalization','probability')
title('Annual precipitation 1981-2015 normal distribution')
xlabel('rainfall (mm)')
ylabel('%')
saveas(gcf,'5.png')

%plot 2gether
figure(7)
histogram(prec2_annual(1:80),bins,'Normalization','probability','FaceAlpha',0.4,'FaceColor','blue','edgecolor','none')
hold on
histogram(prec2_annual(81:end),bins,'Normalization','probability','FaceAlpha',0.5,'FaceColor','yellow','edgecolor','none')
%xlim([-1 17])
%ylim([0 20])
title('Annual precipitation distribution','Fontsize',14)
xlabel('rainfall (mm)','Fontsize',12)
ylabel('%','Fontsize',12)
legend('1901-1980','1981-2015')
set(gca,'FontSize',14);
saveas(gcf,'7.png')

figure(8)
h1 = histogram(prec2_annual(1:80),bins,'Normalization','probability');
h1V = h1.Values;
h2 = histogram(prec2_annual(81:end),bins,'Normalization','probability');
h2V = h2.Values;
bar(h2V-h1V)
title('Annual precipitation difference (1981-2015) - (1901-1980)','Fontsize',14)
ylabel('difference of %','Fontsize',14)
set(gca,'FontSize',14);
saveas(gcf,'8.png')



%%
%nino34 plot
nino34_annualmean = mean(nino34(1:end-1,2:13),2);  %1870-2016
figure(6)
plot(1901:2016,nino34_annualmean(32:end))
title('Annual nino34 1901-2016 time series')
ylabel('Nino34 index')
xlabel('year')
saveas(gcf,'6.png')
%%
figure(9)
scatter(nino34_annualmean(32:end-1),prec2_annual,'LineWidth',1.2)
title('Nino3.4 index vs. Precipitation')
xlabel('Nino3.4 index')
ylabel('rainfall (mm)')
set(gca,'FontSize',14);
saveas(gcf,'9.png')
%%
%==========use nino34 to predict prec==================
%use prec2_annual(1901-2015) and nino34
std_err = sqrt(1-cor_prec_nino^2);   %0.9507

%cost table
percentage_basic = 0.35;
percentage_advan = 0.1;

%cost_table = [68266240 0;
%              38920000+68266240*percentage_basic 38920000;
%              38920000*1.5+68266240*percentage_advan 38920000*1.5];

cost_table = [1754 0;
              514+1754*percentage_basic 514;
              514*1.5+1754*percentage_advan 514*1.5];

%z of nino34
Z_nino34 = (nino34_annualmean-mean(nino34_annualmean))/std(nino34_annualmean);  %1870-2016
prec_annualmean = mean(prec2_annual);
prec_std = std(prec2_annual);
z_bn = (prec_bn-prec_annualmean)/prec_std;
z_an = (prec_an-prec_annualmean)/prec_std;

%choose one year within 1870-2016
for i = 1:147
    z_year(i) = cor_prec_nino*Z_nino34(i);
    z_bn_new(i) = (z_bn-z_year(i))/std_err;
    z_an_new(i) = (z_an-z_year(i))/std_err;
    P_bn(i) = normcdf(z_bn_new(i));
    P_n(i) = normcdf(z_an_new(i)) - normcdf(z_bn_new(i));
    P_an(i) = 1-normcdf(z_an_new(i));
  
end

%%
%EMV
P_year = [P_bn;P_n+P_an];
EMV_year = cost_table*P_year;

[M,I] = min(EMV_year);
one = find(I==1);
two = find(I==2);
three = find(I==3);
%nino34_annualmean(find(I>1))
year_all = 1870:2016;
year_2 = year_all(two);
year_3 = year_all(three);

%%
%EU -- risk-adversed one (sqrt)
Pmin = min(min(cost_table));
Pmax = max(max(cost_table));

for i=1:3
    for j=1:2
        P = cost_table(i,j);
        U_all(i,j)=sqrt((Pmax-P)/(Pmax-Pmin))*100;
    end
end

EU_year = U_all*P_year;

[M2,I2] = max(EU_year);
one2 = find(I2==1);
two2 = find(I2==2);
three2 = find(I2==3);
%nino34_annualmean(find(I2>1))
year_all = 1870:2016;
year_12 = year_all(one2)
nino34_12 = nino34_annualmean(one2)
year_22 = year_all(two2);
nino34_22 = nino34_annualmean(two2);
year_32 = year_all(three2);
nino34_32 = nino34_annualmean(three2);

%%
%EU -- risk-adversed one (3rd root)
Pmin = min(min(cost_table));
Pmax = max(max(cost_table));

for i=1:3
    for j=1:2
        P = cost_table(i,j);
        U_all(i,j)=((Pmax-P)/(Pmax-Pmin))^(1/3)*100;
    end
end

EU_year = U_all*P_year;

[M2,I2] = max(EU_year);
one2 = find(I2==1);
two2 = find(I2==2);
three2 = find(I2==3);
%nino34_annualmean(find(I2>1))
year_all = 1870:2016;
year_12 = year_all(one2)
nino34_12 = nino34_annualmean(one2)
year_22 = year_all(two2);
nino34_22 = nino34_annualmean(two2);
year_32 = year_all(three2);
nino34_32 = nino34_annualmean(three2);

%%
%okay it's a failed test :(
rcp85 = dlmread('pr_2020-2039_rcp85.txt');  %(12month,16model)
rcp85_annual = sum(rcp85);
rcp85_ensemble = mean(rcp85_annual);   %3.1118e+03

%=======ztest=========
%standard error of the mean
SE = prec_std/sqrt(2039-2020+1);

%z-score
zscore = (rcp85_ensemble-prec_annualmean)/SE   %rcp85: 4.5760 -> significant!

%%



