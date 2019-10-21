%in case you are starting with an xml file, uncomment the two rows below
%and follow.
% fileName = 'kb_g.794_model.sbml.xml' %put your model here
% model = readCbModel(fileName);

%In case it is not done, set the objective function of biomass to 1. First find the reaction
%number for biomass
% stridx('biomass',model.rxns);
% model.c(447)=1;

%conveting metabolite names

%eco=normalizeMetNames(model);

layout=createLayout();
layout=addModel(layout,eco);

layout=setDims(layout,11,11);
layout=setInitialPop(layout,'colonies',1e-7,false);
LWM={'Ammonia'};


for i=1:length(layout.mets)

%         layout.media_amt(i)=cell2mat(newdata(i));
        layout=setMedia(layout,layout.mets{i},10);
    
    
end

%set static media values to 0
layout.static_media(11,11,11,11)=0;

%set time step
layout.params.timeStep=0.1

%set the death rate
layout.params.deathRate=0;

%set max cycles
layout.params.maxCycles=40000;

%set to max

%set max biomass allowance
layout.params.maxSpaceBiomass=1;

%set the volume cube in cm^3
layout.params.spaceWidth=1;

%set the output files
%total biomass
layout.params.writeTotalBiomassLog=1;

%media log
layout.params.writeMediaLog=1;
layout.params.mediaLogRate=400;
%flux log
layout.params.writeFluxLog=0;

%put a timestep
layout.params.useLogNameTimeStamp=1;


folder='.';
file='comets_layout.txt';

createCometsFiles(layout,folder,file,true);


