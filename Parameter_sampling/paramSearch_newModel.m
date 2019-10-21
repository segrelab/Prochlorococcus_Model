%A script for parameter searching using random sampling:
%1. Picks random values for 6 parameters out of defines ranges, values will
%be recorded in a result matrix.
%2. Changes Lower and Upper bounds of selected exchanges in the model.
%3. Optimizes metabolic model to find FBA solution. The objective solution
%and the selected sampled exchange flux will be recorded in the result
%matrix.


%initCobraToolbox(); %To be activated if run first time.
clear;
tic();
%load model file
%glycogen created- not allowed back into biomass
% load('iJC568dcg_pro99_mustsecgly.mat');
% name='iJC568dcg_pro99_mustsecgly';
%check on the fakes
% model.ub(714)=0;
% model.ub(719)=0;
% model.ub(724)=0;
% model.ub(727)=0;
% model.ub(728)=0;
% model.ub(729)=0;
% model.ub(730)=0;
% model.ub(767:786)=0;
%%glycogen created, allowed back into biomass
% load('iJC568dcg_pro99_R02011on.mat');
% name='iJC568dcg_pro99_R02011on';
% %glycogen not created!
load('iSO595c1.mat');
runNo='ct0_test3';
name=['iSO595c1_run_',num2str(runNo)];

% model.lb(86)=4.7;
% model.ub(86)=4.7;
%Number of samples
samples=100000;
%general variables
i=1;
%ranges
%HCO3
cmax=-13.5;
cmin=-0.35;
%CO2
co2max=0.818;
co2min=0;
%Ammonia
nmax=-0.56;
nmin=0;
%Light
hvmax=181;
hvmin=13.4;
%Phosphate
pmax=-0.227;
pmin=0;
%rubisco
rmax=4.7;
rmin=0.27;


%result matrix
resMat=zeros(8,1);
%infeasible solution matrix
badSol=zeros(6,1);
%flux matrix- keep all the fluxes
fluxmtx=zeros(992,1);
badcounter=0;    
for i=1:samples
    %select random values for parameters and store in the result matrix.
    %HCO3
    c=selectParam(cmin,cmax);    
    %CO2
    co2=selectParam(co2min,co2max);    
    %Nitrogen (Ammonia)
    n=selectParam(nmin,nmax);    
    %Light
    hv=selectParam(hvmin,hvmax);    
    %Phosphate
    p=selectParam(pmin,pmax);
    %rubisco
    r=selectParam(rmin,rmax);
    
    %change lower and upper bounds of selected exchanges in the model.
    %HCO3 - fixed
    model.lb(749)=cmax;
    model.ub(749)=c;
    %CO2- range
    model.lb(740)=0;
    model.ub(740)=co2;
    %ammonia - range
    model.lb(734)=n;
    model.ub(734)=0;
    %light - fixed
    model.lb(752)=-hv;
    model.ub(752)=0;
    %phosphate - range
    model.lb(745)=p;
    model.ub(745)=0;
    %rubisco
    model.lb(86)=0;
    model.ub(86)=r;
  
    %solve for current conditions
    %optimize model to get a solution 
    sol=optimizeCbModel(model);
    % check if new changes are blocking 
%     BlockedReaction = findBlockedReaction(model);
%     for k=1:length(BlockedReaction)
%         temp=BlockedReaction{1,k}
%         temp=char(temp)
%         place=stridx(temp, model.rxns)
%         sol.x(place)
%         printRxnFormula(model, temp);
%     end
    
     
    %if strcmp(sol.origStat, 'INFEASIBLE')
    if sol.stat==0
        %save 'bad' parameter choices
        badSol(1,i)=c;
        badSol(2,i)=co2;
        badSol(3,i)=n;
        badSol(4,i)=hv;
        badSol(5,i)=p;
        badSol(6,i)=r;
        badcounter=badcounter+1;
        
        continue;
    end
    %assuming a feasible solution- get the glycogen flux
    x=sol.x(746);
    
    %keep all the solutions- populate the matrices
    %save parameter range maximum in a matrix
    resMat(1,i)=c;
    resMat(2,i)=co2;
    resMat(3,i)=n;
    resMat(4,i)=hv;
    resMat(5,i)=p;
    resMat(6,i)=sol.f;
    resMat(8,i)=x;
    resMat(7,i)=r;
    
    %save all lower flux bounds for currnet solution
    fluxmtx(:,i)=sol.x(1:992);
    disp(i);
end

%clean out all the all-zeroes columns from the result matrix
resMat=resMat(:,any(resMat,1));
%plotting

fig=figure('Position', get(0, 'Screensize'));
subplot(3,3,1)
histogram(resMat(1,:));
title('HCO3');
ylabel('Frequency');
xlabel('Flux (mmol/gdw h)');

subplot(3,3,2)
histogram(resMat(2,:));
title('CO2');
ylabel('Frequency');
xlabel('Flux (mmol/gdw h)');

subplot(3,3,3)
histogram(resMat(3,:));
title('Ammonia');
ylabel('Frequency');
xlabel('Flux (mmol/gdw h)');

subplot(3,3,4)
histogram(resMat(4,:));
title('Light');
ylabel('Frequency');
xlabel('Flux (mmol/gdw h)');

subplot(3,3,5)
histogram(resMat(5,:));
title('Phosphate');
ylabel('Frequency');
xlabel('Flux (mmol/gdw h)');

subplot(3,3,6)
histogram(resMat(7,:));
title('RuBisCo');
xlabel('Flux (mmol/gdw h)');
ylabel('Frequency');

subplot(3,3,7)
glyc=resMat(8,:);
keep = (glyc(:) ~= 0);
g=glyc(keep);
biomass=resMat(6,:);
b=biomass(keep);
sc=scatter(g,b);
sc.Marker='.';
sc.MarkerFaceColor='k';
sc.MarkerEdgeColor='k';
ylim([0 inf]);
ylabel('biomass (1/h)');
xlabel('Glycogen flux (mmol/gdw h)');

subplot(3,3,8:9)
plot(fluxmtx(746,:));
ylabel('Glycogen(mmol/gdw h)');
xlabel('Sample');

filename=[name,'_',num2str(samples),'.fig'];
tfilename=[name,'_',num2str(samples),'.tif'];
savefig(filename);

%plotting secreted exchanges
[selExc,selUpt]=findExcRxns(model);
ex=find(selExc==1);
secreted=fluxmtx(ex,:);
% ex=733:791;
%collect only lines with positive values
s=zeros(1,samples);
for j=1:length(ex)
    line=find(secreted(j,:)>0);
    if ~isempty(line)
        pos(j)=ex(j);
        s(j,:)=secreted(j,:);
    end
end
%clean out the zeros
nonzero=find(~all(s==0,2));
for z=1:length(nonzero)
    %clean out biomass and CO2
    marker(z)=model.rxns(ex(nonzero(z)));
    nz(z,:)=s(nonzero(z),:);
end
%finally we can plot
% figure('Position', get(0, 'Screensize'));
% plotrow=fix(length(nonzero)/4)+1;
% for plotId=1:length(nonzero)
%     figure;
%     %subplot(plotrow,4,plotId);
%     histogram(nz(plotId,:));
%     title(marker(plotId));
%     set(gca, 'YScale', 'log')
%     ylabel('Frequency');
%     xlabel('Value');
%     figure;
%     sz = 25;
%     % removing zero glycogen from the plot
%     glyc=resMat(8,:);
%     keep = (glyc(:) ~= 0);
%     g=glyc(keep);
%     vari=nz(plotId,:);
%     v=vari(keep);
%     sc=scatter(g,v);
%     sc.Marker='o';
%     sc.MarkerFaceColor='k';
%     sc.MarkerEdgeColor='k';
%     sc.LineWidth=1;
%     ylim([0 inf]);
%     ylabel(marker(plotId));
%     xlabel('Glycogen flux (mmol/gdw h)');
%     % calculate the mean
%     m=num2str(mean(nz(plotId,:)));
%     legend(m);
%     
% end
secname=[name,'_',num2str(samples),'_secreted','.fig'];
savefig(secname);

% %correlation matrix
% mtx=resMat';
% C = corrcoef(mtx);
% [R,PValue]=corrplot(mtx,'varNames',{'HCO3','CO2','NH3','Light','Phosphate','Biomass','RuBisCo','Glycogen'},'testR','on');
% corrname=[name,'_',num2str(samples),'_corr','.fig'];
% savefig(corrname);
% corrnamet=['corr_',name,num2str(samples),'.txt'];
% dlmwrite(corrnamet, mtx, 'delimiter','\t','newline','pc')
% 
% %correlation matrix- infeasible
% bs=badSol';
% corrplot(bs,'varNames',{'HCO3','CO2','NH3','Light','Phosphate','RuBisCo'},'testR','on');
% bcorrname=[name,'_',num2str(samples),'_inf_corr','.fig'];
% savefig(bcorrname);
% bcorrnamet=['inf_corr_',name,num2str(samples),'.txt'];
% dlmwrite(bcorrnamet, mtx, 'delimiter','\t','newline','pc')
% 
mtx=resMat';
matname=['result_mtx_',name,'_',num2str(samples),'.mat'];
save(matname,'mtx');
flbname=['flux_',name,'_',num2str(samples),'.mat'];
save(flbname,'fluxmtx');
badname=['infeasible_',name,num2str(samples),'.mat'];
save(badname,'badSol');
% 
% 
% 
% %make corrplot with excreted things
% %stack the matrices together
% [i1,j1] = ndgrid(1:size(nz,1),1:size(nz,2));
% [i2,j2] = ndgrid(1:size(resMat,1),(1:size(resMat,2))+size(nz,2));
% seccorrfull= accumarray([i1(:),j1(:);i2(:),j2(:)],[resMat(:);nz(:)]);
% % seccorrfull=[resMat;nz];
% % seccorr=[resMat;nz(8:end,:)];
% 
% %get the means and standard deviation for each parameter and excreted flux
% stats(1,:)=mean(seccorrfull,2);
% stats(2,:)=std(seccorrfull,0,2);
% stats=stats';
% statname=['mean_std_',name,num2str(samples),'.mat'];
% statnamet=['mean_std_',name,num2str(samples),'.txt'];
% dlmwrite(statnamet, stats, 'delimiter','\t','newline','pc')
% ncorr=seccorrfull';
% save(statname,'stats');
% corrplot(ncorr,'testR','on');
% scorrname=['corrwithsec_',name,num2str(samples),'.fig'];
% savefig(scorrname);
% ncorrname=['corrwithsec_',name,num2str(samples),'.mat'];
% ncorrnamet=['corrwithsec_',name,num2str(samples),'.txt'];
% save(ncorrname,'ncorr');
% dlmwrite(ncorrnamet, ncorr, 'delimiter','\t','newline','pc')


toc();
