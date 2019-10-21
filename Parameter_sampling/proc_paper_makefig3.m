% a script that calculates the fraction of samples producing a certain
% parameter of interest, like glycogen, amino acids etc. and compares it to
% the general fraction of samples. We look specifically on staple uptakes
% such as net carbon, ammonia, phosphate and light uptakes. The script then
% plots all the fractions on one figure.

clear all;
load('iSO595c1.mat');
run='ct0_';
% conctatanting the matrices
f=[];
load('flux_iSO595c1_run_ct0_test2_100000.mat');
f=fluxmtx;
% for i=1:10
%     filename=['flux_iSO595c1_run_',run,num2str(i),'_100000.mat'];
%     load(filename);
%     f=[f,fluxmtx];
% end

%ranges
%HCO3
cmax=-1.35;
cmin=-0.35;
%CO2
co2max=0.818;
co2min=0;
%Ammonia
nmax=-0.56;
nmin=0;
%Light
hvmax=-181;
hvmin=-13.4;
%Phosphate
pmax=-0.227;
pmin=0;
%rubisco
rmax=4.7;
rmin=0.27;


nbins = 20;

% Fluxes of interest:
ind_growth = 63;
ind_glycg = 746;
ind_HCO3 = 749;
ind_CO2 = 740;
ind_ammonia = 734;
ind_light = 752;
ind_P = 745;
ind_RUBISCO = 86;
%Amino acids
ind_AA = [765 766 767 768 769 770 771 772 773 774 775 776 777 778 779 780 781 ...
782 783 784];
%Fatty acids
ind_FA = [717 786];
%Organic acids
ind_OA = [733 744 801 803 805 807 809];

%Find data points (indeces) for which there is glycogen production

% glycg_prod_ind = find((f(ind_glycg,:))>0);
[sortf1,sortf2]=sort(f(ind_glycg,:));
glycg_prod_ind = sortf2((end-1000):end);

% Plot nonzero value of glycogen production
figure; hold on; set(gca,'fontsize',[16])
plot(sort(f(ind_glycg,glycg_prod_ind)));
xlabel('index (sorted)')
ylabel('Nonzero Glycogen production flux')
title(['Nonzero/Total  = ' num2str(100*length(glycg_prod_ind)/20000) '%'])


% Exudates for datapoints producing glycogen:
% Amino acids
AA_prod_ind = f(ind_AA,glycg_prod_ind);
AA_exud = f(ind_AA,:);
[sumf1,sumf2]=sort(sum(f(ind_AA,:)));
ind_AA_exude = sumf2((end-1000):end);

% Fatty acids
FA_prod_ind = f(ind_FA,glycg_prod_ind);
FA_exud = f(ind_FA,:);
[sumf1,sumf2]=sort(sum(f(ind_FA,:)));
ind_FA_exude = sumf2((end-1000):end);

% Organic acids
OA_prod_ind = f(ind_OA,glycg_prod_ind);
OA_exud = f(ind_OA,:);
[sumf1,sumf2]=sort(sum(f(ind_OA,:)));
ind_OA_exude = sumf2((end-1000):end);

%ind_exud_data =  find(sum(f(ind_exud,:))>0.25);
%ind_exud_data =  find(sum(f(ind_exud,:))>0.2);

% ----------------------- PLOT ---------------------------------

% GROWTH levels that give rise to Exudate Production
this_min = min(f(ind_growth,:));
this_max = max(f(ind_growth,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
% Plot general fraction
[q,w]=hist(f(ind_growth,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_growth,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_growth,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_growth,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_growth,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('Growth Rate')
ylabel('Fraction')
filename=['EX_growth_',run,'.fig'];
savefig(filename);

% Nitrogen levels that give rise to Exudate Production
this_min = min(f(ind_ammonia,:));
this_max = max(f(ind_ammonia,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_ammonia,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_ammonia,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_ammonia,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_ammonia,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_ammonia,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('Ammonia Uptake')
ylabel('Fraction')
filename=['EX_nitrogen_',run,'.fig'];
savefig(filename);

% Light levels that give rise to Exudate Production
this_min = min(f(ind_light,:));
this_max = max(f(ind_light,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_light,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_light,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_light,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_light,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_light,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('Light Uptake')
ylabel('Fraction')
filename=['EX_light_',run,'.fig'];
savefig(filename);

% HCO3 levels that give rise to Exudate Production
this_min = min(f(ind_HCO3,:));
this_max = max(f(ind_HCO3,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_HCO3,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_HCO3,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_HCO3,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_HCO3,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_HCO3,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('HCO3 Uptake')
ylabel('Fraction')
filename=['EX_hco3_',run,'.fig'];
savefig(filename);

% CO2 outflux that give rise to Exudate Production
this_min = min(f(ind_CO2,:));
this_max = max(f(ind_CO2,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_CO2,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_CO2,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_CO2,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_CO2,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_CO2,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('CO2 secretion')
ylabel('Fraction')
filename=['EX_co2_',run,'.fig'];
savefig(filename);

% HCO3 - CO2 outflux that give rise to Exudate Production
this_min = min(f(ind_HCO3,:)-f(ind_CO2,:));
this_max = max(f(ind_HCO3,:)-f(ind_CO2,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
netc=f(ind_HCO3,:)-f(ind_CO2,:);
[q,w]=hist(netc,[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_HCO3,glycg_prod_ind)-f(ind_CO2,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_HCO3,ind_AA_exude)-f(ind_CO2,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_HCO3,ind_FA_exude)-f(ind_CO2,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_HCO3,ind_OA_exude)-f(ind_CO2,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('HCO3 - CO2')
ylabel('Fraction')
filename=['EX_net_c_',run,'.fig'];
savefig(filename);

% RUBISCO flux that give rise to Exudate Production
this_min = min(f(ind_RUBISCO,:));
this_max = max(f(ind_RUBISCO,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_RUBISCO,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]); 
% Plot glycogen fraction
[q,w]=hist(f(ind_RUBISCO,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_RUBISCO,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_RUBISCO,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_RUBISCO,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('RUBISCO flux')
ylabel('Fraction')
filename=['EX_rubisco_',run,'.fig'];
savefig(filename);

% Phosphorus uptake that give rise to Exudate Production
this_min = min(f(ind_P,:));
this_max = max(f(ind_P,:));
this_bin = (this_max - this_min)/nbins;
figure; hold on; %set(gca,'fontsize',[16])
[q,w]=hist(f(ind_P,:),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','k','linewi',[4]);
% Plot glycogen fraction
[q,w]=hist(f(ind_P,glycg_prod_ind),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color','r','linewi',[2]);
% Plot AA fraction
[q,w]=hist(f(ind_P,ind_AA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0,0.447,0.741],'linewi',[2]);
% Plot FA fraction
[q,w]=hist(f(ind_P,ind_FA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.929,0.694,0.125],'linewi',[2]);
% Plot OA fraction
[q,w]=hist(f(ind_P,ind_OA_exude),[this_min : this_bin : this_max]);
plot(w,q/sum(q),'color',[0.466,0.674,0.188],'linewi',[2]);
xlabel('P flux')
ylabel('Fraction')
filename=['EX_p_',run,'.fig'];
savefig(filename);






