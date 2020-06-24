%plot COMETS results per media and flux

tic();
j=1;
cellWeight=66E-15;

%load the output total biomass file into a matrix
files = dir('*.m');
%mediaFiles=dir('C:\COMETS\proc\dfba_results\media\*.m');
for i=1:length(files)
    file=load(files(i).name);
    n=strsplit(files(i).name,'_');
    %looking for the correct media file from the media directory
    mediafname=strcat('C:\COMETS\proc\dfba_results\media\media_',n(3));
    fluxfname=strcat('C:\COMETS\proc\dfba_results\flux\flux_',n(3));
    %loading the media file and converting media information to a table structure, easier for plotting
    media=parseMediaLog(mediafname{1});
    toc();
    flux=parseFluxLog(fluxfname{1});
    toc();
    flux{:,1}=flux{:,1}*0.00416667;
    media{:,1}=media{:,1}.*0.00416667;
    %getting fluxes of interest out of the flux table
    fcsource=flux(flux.rxn==756,:); %hco3
    fnsource=flux(flux.rxn==741,:); %ammonia
    fpsource=flux(flux.rxn==752,:); %Phosphate
    fssource=flux(flux.rxn==768,:); %Sulfate
    fgsource=flux(flux.rxn==753,:); %glycogen
    fbsource=flux(flux.rxn==759,:); %lightEX
    fggsource=flux(flux.rxn==759,:); %glycogengranule
    %collect reaction flux data
    r1=flux(flux.rxn==587,:); %R02110
    r2=flux(flux.rxn==588,:); %R02111
    r3=flux(flux.rxn==58,:); %BCarbonstorage
    r4=flux(flux.rxn==586,:); %R02109
    %converting from biomass/volume to cells/volume
    output(:,1)=file(:,1).*0.00416667;
    %converting from time step to days 0.1h=0.00416667days
    j=j+1;
    output(:,j)=file(:,2)./cellWeight;
end
%experimental data
xdot9312 = [0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27];
cell_data = [1419052.488   886758.691	2510054.533	4342450.579	6618183.367	9836826.858	13807174.51	22447085.89	40878834.36	55403033.4	70336826.86	NaN	101285105.7	110418967.3	96088786.64	96225630.54	NaN         87538854.81	83861877.98	75219154.74	71500937.29	NaN         66534417.34	56434688.35	47757181.57	39622764.23	34944173.44	29145799.46;
             1118204.628   1716336.772	2274463.34	3708809.59	6079871.759	8884165.04	11600641.2	18485503.21	25570393.09	40518817.95	NaN         NaN	90782408.7	112800111.5	99607053.25	105934650.2	77660181.2	97932422.4	81077528.59	96109906.43	87369523.24 84302539.73	78051091.64	61545522.06	53957077.08	45512995.69	34744244.76	30976533.49;
             1383908.046   1882758.621	2556475.096	3233409.962	5221302.682	7640459.77	NaN         18259003.83	25006896.55	39793256.7	NaN         NaN	91902222.22	116875862.1	119497777.8	118448582.4	NaN         101173842	85132970.03	87182425.07	80228201.63	73440871.93	61688692.1	49030517.71	34757765.67	32867574.93	28680926.43	23077588.56]*1000;
cell_data=cell_data./1000;
ydot9312 = nanmean(cell_data);
SD9312 = nanstd(cell_data);

Nx = [0 2 4	6	8	9	12	13	14	15	16	23];
N_data = [209.2235294	132.3294118	119.8117647	116.2352941	103.7176471	41.12941176	27.65803922	0           0           68.66823529	0	0;
          194.3215686	80.47058824	114.4470588	125.1764706	78.68235294	91.2        43.87137255	73.91372549	34.33411765	82.61647059	0	0;
          156.172549	82.25882353	118.0235294	155.5764706	121.6       44.10980392	82.02039216	75.64235294	21.45882353	15.03905882	0	0];
N_data=N_data./1000000;
Ny = nanmean(N_data);
N_SD = nanstd(N_data);


figure
subplot(1,6,1:2)
semilogy(output(:,1),output(:,2:(i+1)),'LineWidth',2);
ylabel('Cells/ml');
xlabel('Days');

hold on
errorbar(xdot9312,ydot9312,SD9312,'o','MarkerFaceColor','k','Color','k')


subplot(3,6,3:4)
plot(output(:,1),output(:,2:(i+1)),'LineWidth',2);
ylabel('Cells/ml');
xlabel('Days');
% xlim([0 40]);
hold on
errorbar(xdot9312,ydot9312,SD9312,'o','MarkerFaceColor','k','Color','k')

subplot(3,6,5)
sp1=plotMediaTimecourse(media,{'Photon[e]'});
sp1.LineWidth=2;
sp1.Color = [230, 81, 0]./255
ylabel('Mol');
xlabel('Days');
% ylim([0 0.2]);
% xlim([0 40]);

sub=subplot(3,6,6)
spg1=plot(fbsource.t,fbsource.flux)
spg1.LineWidth=2;
spg1.Color=[230, 81, 0]./255;
ylabel('mmol/gdw h');
xlabel('Days');
title('Flux');
% yyaxis right;


% ax=gca;
% ax.YColor = 'k';
% xlim([0 40]);
% legend ('Light');


subplot(3,6,7)
sp3=plotMediaTimecourse(media,{'HCO3[e]'});
sp3.LineWidth=2;
sp3.Color = [183, 28, 28]./255
ylabel('Mol');
xlabel('Days');
% xlim([0 40]);

subplot(3,6,8)
sp4=plot(fcsource.t,fcsource.flux);
sp4.LineWidth=2;
sp4.Color=[183, 28, 28]./255;
ylabel('mmol/gdw h');
xlabel('Days');
% xlim([0 40]);
title('Flux');

subplot(3,6,9)
sp5=plotMediaTimecourse(media,{'Ammonia[e]'});
sp5.Color = [156, 39, 176]./255
sp5.LineWidth=2;
ylabel('Mol');
xlabel('Days');
hold on
errorbar(Nx,Ny,N_SD,'o','MarkerFaceColor','r')
% xlim([0 40]);

subplot(3,6,10)
sp6=plot(fnsource.t,fnsource.flux);
sp6.LineWidth=2;
sp6.Color=[156, 39, 176]./255;
ylabel('mmol/gdw h');
xlabel('Days');
title('Flux');
% xlim([0 40]);

subplot(3,6,11)
sp7=plotMediaTimecourse(media,{'Orthophosphate[e]'});
sp7.Color = [0, 150, 136]./255
sp7.LineWidth=2;
ylabel('Mol');
xlabel('Days');
% xlim([0 40]);

subplot(3,6,12)
sp8=plot(fpsource.t,fpsource.flux);
sp8.LineWidth=2;
sp8.Color=[0, 150, 136]./255;
ylabel('mmol/gdw h');
xlabel('Days');
title('Flux');
% xlim([0 40]);

subplot(3,6,13)
sp9=plotMediaTimecourse(media,{'Sulfate[e]'});
sp9.Color = [253, 216, 53]./255
sp9.LineWidth=2;
ylabel('Mol');
xlabel('Days');
% xlim([0 40]);

subplot(3,6,14)
sp10=plot(fssource.t,fssource.flux)
sp10.LineWidth=2;
sp10.Color=[253, 216, 53]./255;
ylabel('mmol/gdw h');
xlabel('Days');
title('Flux');


%plot reaction fluxes
subplot(3,6,15);
sp11=plot(r1.t,r1.flux);
sp11.LineWidth=2;
sp11.Color=[216, 27, 96]./255;
ylabel('mmol/gdw h');
xlabel('Days');
hold on;
sp14=plot(r4.t,r4.flux);
sp14.LineWidth=2;
sp14.Color=[255, 128, 171]./255;
legend ('R02110','R02109');
title('Flux');

subplot(3,6,16);
sp12=plot(r2.t,r2.flux);
sp12.LineWidth=2;
sp12.Color=[236, 64, 122]./255;
ylabel('mmol/gdw h');
xlabel('Days');
legend('R02111');
title('Flux');


subplot(3,6,17);
sp13=plot(r3.t,r3.flux);
sp13.LineWidth=2;
sp13.Color=[124, 77, 255]./255;
ylabel('mmol/gdw h');
xlabel('Days');
legend('BCarbonstorage');
title('Flux');
% xlim([0 40]);

% subplot(3,6,18);
% spg=plot(fgsource.t,fgsource.flux)
% spg.LineWidth=2;
% spg.Color=[48, 79, 254]./255;
% 
% hold on;
% spgg=plot(fggsource.t,fggsource.flux);
% spgg.LineWidth=2;
% spgg.Color=[0, 184, 212]./255;
% ylabel('mmol/gdw h');
% xlabel('Days');
% legend ('Glycogen','Glycogen granule');
% subplot(3,6,15)
% % plotMediaTimecourse(media,{'L_Alanine[e]','L_Arginine[e]','L_Asparagine[e]','L_Aspartate[e]','L_Cystine[e]','L_Glutamate[e]','L_Glutamine[e]','L_Histidine[e]','L_Isoleucine[e]','L_Leucine[e]','L_Lysine[e]','L_Methionine[e]','L_Phenylalanine[e]','L_Proline[e]','L_Serine[e]','L_Threonine[e]','L_Tryptophan[e]','L_Tyrosine[e]','L_Valine[e]'},true);
% xlim([0 40]);
% sp11=plotMediaTimecourse(media,{'Oxygen[e]','H[e]'});
% % sp5.Color = {[0.26, 0.65, 0.96],[0.26, 0.65, 0]};
% % sp5.LineWidth={2,2};
% ylabel('Mol');
% xlabel('Days');


% subplot(3,6,16)
% sp12=plotMediaTimecourse(media,{'CO2[e]'});
% sp12.Color = [255, 111, 0]./255;
% sp12.LineWidth=2;
% ylabel('Mol');
% xlabel('Days');
% xlim=[0 40];

% subplot(3,6,15:18)
% sp13=plotMediaTimecourse(media,{'Oxygen[e]','H[e]','CO2[e]','7_Aminomethyl_7_carbaguanine[e]','Adenosine_3_5_bisphosphate[e]','Glycolaldehyde[e]','Methanol[e]','Nicotinate[e]','S_Adenosyl_4_methylthio_2_oxobutanoate[e]','S_Methyl_5_thio_D_ribose[e]'},true);
% ylabel('Mol (logscale)');
% xlabel('Days');
% xlim([0 40]);
% text({'Oxygen[e]','H[e]','CO2[e]','7_Aminomethyl_7_carbaguanine[e]','Adenosine_3_5_bisphosphate[e]','Glycolaldehyde[e]','Methanol[e]','Nicotinate[e]','S_Adenosyl_4_methylthio_2_oxobutanoate[e]','S_Methyl_5_thio_D_ribose[e]'});

toc();