% Calculating secretion means then summing them up to calculate the biomass
% to glycogen to excretion ratio.
% This script is using the flux matrices.

clear all;
%Load the model
load('iSO595c1.mat');
run='ct0_';
% conctatanting the matrices
wm=[];

for i=1:10
    filename=['flux_iSO595_run_',run,num2str(i),'_100000.mat'];
    load(filename);
    wm=[wm,fluxmtx];
end

%create a sub matrix of all the nonzero glycogen
line=find(wm(746,:)>0);
glycmtx=wm(:,line);
% prepare the plot matrix
% put biomass in
z(1,:)=wm(735,:);
%put glycogen in
z(2,:)=wm(746,:);




% getting the exchanges , need to leave out biomass and glycogen and only
% calculate the mean of the positive values. Also can disregard protons and
% water.

[selExc,selUpt]=findExcRxns(model);
ex=find(selExc==1);

% n=1;
% for f=1:500000
%     disp(f);
%     zsum=0;
%     for e=1:length(ex)
%         % is the line anything we want to disregard?
%         if (ex(e,:)==735) || (ex(e,:)==746) || (ex(e,:)==747) || (ex(e,:)==750)||(ex(e,:)==762)
%             continue
%         end
%         % are all the values positive? (What about part positive/part
%         % negative?)
%          line=find(wm(ex(e,1),:)>0);
%          l=wm(ex(e,1),line);
%         if ~isempty(line)
%             zsum=zsum+wm(ex(e,1),f);
%             
%         end
%         
%     end
%     z(3,f)=zsum;
% end

% Look into the exudates, use only positive nonzero values
exTable={};
n=1;
for e=1:length(ex)
    % is the line anything we want to disregard?
    if (ex(e,:)==735) || (ex(e,:)==746) || (ex(e,:)==747) || (ex(e,:)==750)||(ex(e,:)==762)
        continue
    end
    % are all the values positive? (What about part positive/part
    % negative?)
    
    line=find(wm(ex(e,1),:)>0);
    % find all the occurneces in glycogen producing configs
    gline=find(glycmtx(ex(e,1),:)>0);
    l=wm(ex(e,1),line);
    g=glycmtx(ex(e,1),gline);
    len=length(l);
    if ~isempty(line)
        exTable(n,1)={ex(e,1)};
        exTable(n,2)=model.rxns(ex(e,1));
        
        means(n,1)=mean(l);
        exTable(n,3)={mean(l)};
        exTable(n,4)={length(line)};
        exTable(n,5)={length(gline)};
        exTable(n,6)={mean(g)};
        disp(mean(l));
        n=n+1;
    end
end

exSum=sum(means)
% figure;
% subplot(2,2,1)
% scatter(z(1,:),z(3,:));
% subplot(2,2,2)
% scatter(z(2,:),z(3,:));

