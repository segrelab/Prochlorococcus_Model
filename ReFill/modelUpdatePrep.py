''' Patching holes in the prochlorcoccus metabolic network with new annotated reactions.
the script uses a recursive function to build chains of possible reactions to be added.
 a chain of reaction can only be added if all the metabolites are connected to the network.
 the goal is to reduce reaction for manual curation.
 Written by Shany Ofaim - Segre Lab, Boston University, 2019'''

#-------- imports -------
from bioservices import *
import re
global prime
global met
global key
import scipy.io
import cobra
import cobra.test
global metgo
metgo=False

#-------- Functions -------
def addNetEntry(reaction,net):
    global neighbors
    # writes network entries and collect a set of added neighbors
    for m in reacBank[reaction]['mets']:
        if m in net:
            for r in net[m]:
                netout.write(reaction + '\t' + r +'\n')
                neighbors.add(r)

def isValid(reaction):
    ''' Checks if all the metabolites in a reaction are connected'''
    inModel = []
    global conn
    global network
    global missing
    conn=[]

    if 'mets' in reacBank[reaction]:
        l = len(reacBank[reaction]['mets'])
        # investigating metabolites in the reaction
        for met in reacBank[reaction]['mets']:
            # add name fixes for some metabolites following curation processes
            if '(R)_2,3_Dihydroxy_3_methylbutanoate[c]' in met:
                met = 'R_2_3_Dihydroxy_3_methylbutanoate[c]'
            if 'trans_Dodec_2_enoyl_[acp][c]' in met:
                met = 'trans_Dodec_2_enoyl[c]'
            if "9,15,9'_tricis_zeta_Carotene[c]" in met:
                met = '9_15_9_tricis_zeta_Carotene[c]'
            if 'Deamino_NAD[c]' in met:
                met = 'Deamino[c]'
            if "5'_Phosphoribosyl_N_formylglycinamide[c]" in met:
                met = '5_Phosphoribosyl_N_formylglycinamide[c]'
            if "5'_Phosphoribosylglycinamide[c]" in met:
                met = '5_Phosphoribosylglycinamide[c]'
            if "D_4'_Phosphopantothenate[c]" in met:
                met = 'D_4_Phosphopantothenate[c]'
            if "7,8_Dihydroneopterin_3'_triphosphate[c]" in met:
                met='7_8_Dihydroneopterin_3_triphosphate[c]'
            if met in modelMet or met in connected or 'e_[c]' in met:
                if met in connected:
                    inModel.append(1)
                    conn.append(met)
                else:
                    # build a binary vector of all the metabolites in the reaction.
                    # if a metabolite is in the model put 1 in the vector
                    # if a metabolites has already been connected to the network it is also valid as
                    # 'in the model'
                    inModel.append(1)

                if not sum(inModel)==l:
                    continue
                else:
                    return True
            else:
                # metabolite not in the model or in the connected dictionary.
                # validate it's not connected
                if met in connected:
                    missing.remove(met)
                else:
                    missing.add(met)


def updateDict(query, dict, value):
    # updates a dictionary record to a list
    if query in dict:
        dict[query].append(value)
    else:
        dict[query]=[]
        dict[query].append(value)


def printChain(m,rea,chains,level, marker):
    ''' traversing up the chain of reactions and collecting the candidates.'''
    global reacBank
    global pat0
    global ccount
    global bcount
    broken=False
    valid=False
    tcand={}
    tcand[pat0]=set()
    bcand={}
    bcand[pat0]=set()
    while level!=0:
        # looking to see if the queried reaction is someone's kid in any chain.
        if 'kids' in chains[m] and rea in chains[m]['kids']:
            # Now that this reaction is identified as a kid, let's look at its parents
            # check if all the metabolites in the parents are connected to the network or to other chains
            for p in chains[m]['parents']:
                if isValid(rea):
                    if rea in mreac:
                        print('Reaction '+ rea + ' is already in the model')
                        continue
                    print('Added: ' + rea)
                    tcand[pat0].add(rea)
                    valid=True
                    netout.write(rea +'\t')
                    # check if parent is also valid
                    if isValid(p):
                        if p in mreac:
                            print('Reaction ' + p + ' is already in the model')
                            continue
                        print('Added: ' + p)
                        tcand[pat0].add(p)
                        valid=True
                        # add network entries for both rea and parent
                        netout.write(p+'\n')
                        addNetEntry(rea,network)
                        addNetEntry(p,network)
                    else:
                        print('Chain broken when reached: ' + p)
                        level=0
                        broken=True
                        bcand[pat0].add(p)
                        bcand[pat0]=bcand[pat0]|tcand[pat0]
                        tcand[pat0]=set()
                        break
                    print(str(chains[m]['level']))
                else:
                    print('Chain broken when reached: '+ rea)
                    level = 0
                    broken=True
                    valid=False
                    break
            # the desired direction is going up levels in the recursion so if we are going down we are in a pot. loop
            if chains[m]['level'] == level and level != 0:
                print('Self loop warning for: ' + rea)
                level = 0
            if chains[m]['level']<= level:
                # only go up
                level=chains[m]['level']
                # check for self loops

            elif not broken:
                # you are going the wrong way, or have a loop
                print('suspected loop with: '+ rea + ' '+ p)
                loopout.write('suspected loop with: '+ rea + ' '+ p+ '\n')
                level=0

        # checking recursion level
        if level!=0:
            for par in chains[m]['parents']:
                for c in chains:
                    if c in marker:
                        # if reached the end of the dictionary, go up one level
                        level = 0
                    else:
                        if 'kids' in chains[c] and par in chains[c]['kids']:
                            printChain(c,par,chains,level,marker)
        else:
            # level is 0 - add the chains to their proper place. full chains go to the candidate list, broken go to the
            # broken.
            if not broken and valid:
                # this gets only the end products. add chain number
                if pat0 in cand and cand[pat0][1]:
                    ccount+=1
                    cand[pat0][ccount]=set()
                    cand[pat0][ccount]=cand[pat0][ccount]|tcand[pat0]
                else:
                    cand[pat0][ccount]=cand[pat0][ccount]|tcand[pat0]
            elif not valid:
                if pat0 in bad and bad[pat0][1]:
                    bcount+=1
                    bad[pat0][bcount]=set()
                    bad[pat0][bcount]=bad[pat0][bcount]|bcand[pat0]
                else:
                    bad[pat0][bcount]=bad[pat0][bcount]|bcand[pat0]





def connectReac(reac,reactionBank, step):
    ''' creates the chains of reactions to connect to the network'''
    inModel = []
    rflags={}
    global parentReac
    global pat0
    global metgo
    global bcount
    global ccount
    global mmets
    parentReac=[]
    hit=False
    conn=[]
    if 'mets' in reacBank[reac]:
        l = len(reactionBank[reac]['mets'])
        # investigating metabolites in the reaction
        for met in reactionBank[reac]['mets']:
            if met in model.metabolites or met in connected:
                if met in connected:
                    inModel.append(1)
                    conn.append(met)
                    if met in missing:
                        missing.remove(met)
                else:
                    # build a binary vector of all the metabolites in the reaction.
                    # if a metabolite is in the model put 1 in the vector
                    # if a metabolites has already been connected to the network it is also valid as
                    # 'in the model'
                    inModel.append(1)
                    # add the relevant reactions to the network dictionary
                    for m in mmets:
                        tm=str(m)
                        if met==tm:
                            r = m._reaction
                            for t in r:
                                updateDict(met,network,t._id)

                if not sum(inModel)==l:
                    continue
                else:
                    for c in conn:
                        if 'valid' in chains[c]:
                            chains[c]['valid'].append(reac)
                        else:
                            chains[c]['valid']=[]
                            chains[c]['valid'].append(reac)
                        # collapse back
                        marker =c
                        # traverse the chain up
                        printChain(c, reac, chains, step, marker)
                        step = 0
                        metgo = True


            else:
                # add met and reac as parent to the chain
                chains[met] = {}
                chains[met]['level']=step
                if step==0:
                    ccount=1
                    bcount=1
                    cand[met]={}
                    cand[met][ccount]=set()
                    bad[met]={}
                    bad[met][bcount]=set()
                    pat0=met
                chains[met]['parents'] = []
                chains[met]['parents'].append(reac)


                # What about kids? is this metabolite in other reactions in the bank?
                for entry in reactionBank:
                    if entry == reac:
                        continue
                    if 'mets' in reactionBank[entry] and met in reactionBank[entry]['mets']:
                        hit = True
                        if 'kids' in chains[met]:
                            # This metabolite has kid reactions, it's connected
                            chains[met]['kids'].append(entry)
                            connected[met] = True
                            if met in missing:
                                missing.remove(met)

                        else:
                            chains[met]['kids'] = []
                            chains[met]['kids'].append(entry)
                            connected[met] = True
                            if met in missing:
                                missing.remove(met)

                if 'kids' in chains[met]:
                    step += 1
                    # are the kids connected?
                    for kid in chains[met]['kids']:
                        connectReac(kid, reacBank, step)
                if not hit:
                    dead[met]=True

# ------------ start --------
missing=set()
# read input files
reac2gene={}
# read reaction to gene mapping
with open('pmm.txt') as gin:
    for line in gin:
        if line.startswith('gene'):
            continue
        line=line.split('\t')
        pml=line[0].split(':')
        pmm=pml[-1].strip()
        reac=line[3].strip()
        updateDict(reac,reac2gene,pmm)

# Storing the direction for reactions with R numbers from the modelSEED db
SeedReac={}
print('Reading Reaction directions')
with open('reactions2.19.tsv') as seedin:
    for line in seedin:
        line=line.split('\t')
        rname=line[1].strip()
        direction=line[9].strip()
        if rname.startswith('R') and not rname.startswith('rxn') and not rname.startswith('R_'):
            SeedReac[rname]=direction
# build a network for the reactions connected to metabolites visited

global chains
chains={}
global connected
connected={}
global dead
dead={}
global bad
bad={}
ccount=1
bcount=1
# open network file
netout=open('outputnet.txt','w')
neighbors=set()
print('Reading model reactions')

# read model mat file
mat = scipy.io.loadmat('iSO595c2.mat')
model=cobra.io.read_sbml_model('iSO595c2.xml')
print('Reading model metabolites')
global mmets
mmets=model.metabolites
# get the list of reactions for the model
mreac=model.reactions
genes=model.genes
# make a list for comfortable iterator.
gl=[]
for g in genes:
    gl.append(g.id)
global modelMet
modelMet=[]
for m1 in mmets:
    m1=str(m1)
    m1=m1.replace('__91__','[')
    m1=m1.replace('__93__',']')
    modelMet.append(m1)

bank=[]
print('Reading KEGG reaction bank')
with open('kegg_reaction_bank.txt') as bankin:
    for line in bankin:
        line=line.strip()
        bank.append(line)

# kegg reaction bank
kegg_con = KEGG()
kegg_con.organism='pmm'
paths=kegg_con.pathwayIds

# open an output m file
fout=open('modelRefresh.m','w')
added=set()

# loop log
loopout=open('loopLog.txt','w')


# Build the reaction bank + first pass
global reacBank
reacBank={}
for reaction in bank:
    reacBank[reaction]={}
    if reaction in SeedReac:
        reacBank[reaction]['direction']=SeedReac[reaction]
    else:
        reacBank[reaction]['direction']='='
    kegg_entry = kegg_con.parse(kegg_con.get(reaction))
    if isinstance(kegg_entry, int):
        continue
    definition=kegg_entry['DEFINITION']
    if 'PATHWAY' in kegg_entry:
        pathways=kegg_entry['PATHWAY']
        reacBank[reaction]['pathways']=pathways
    # Checking if all the components for this reaction exist in the model
    # splitting to reactants and products
    sides=definition.split('<=>')
    subs=sides[0].split('+')
    reacBank[reaction]['substrates']={}
    reacBank[reaction]['products']={}
    reacBank[reaction]['mets']=[]
    prods=sides[1].split('+')
    # Checking if the metabolite is in the model
    sflag=True
    pflag=True
    for s in subs:
        s=s.strip()
        # check for coefficient
        if re.match(r"^\d\s\w+",s):
            l=s.split(' ')
            coeff=str(-1*int(l[0].strip()))
            l.pop(0)
            l=' '.join(l)
            s=l
        else:
            coeff=-1
        # convert to underscore
        if s:
            # Cleaning
            if s.startswith('n'):
                s=s.replace('n ','')
            s = s.replace(',', '_')
            s=s.replace('(','')
            s=s.replace(')','')
            s=s.replace("'",'')
            s=s.replace(' ','_')
            s=s.replace('-','_')
            s=s+'[c]'
            reacBank[reaction]['substrates'][s]=coeff
            reacBank[reaction]['mets'].append(s)
        #print(s)
        if s in modelMet:
            index=modelMet.index(s)

            #print(modelMet.index(s))
            #print(modelMet[index])
        else:
            sflag=False
    for p in prods:
        p=p.strip()
        if re.match(r"^\d\s\w+",p):
            l=p.split(' ')
            coeff=l[0].strip()
            l.pop(0)
            l=' '.join(l)
            p=l
        else:
            coeff=1
        # convert to underscore
        if p:
            # cleaning and formatting
            if p.startswith('n'):
                p=p.replace('n ','')
            p = p.replace(',', '_')
            p=p.replace('(','')
            p=p.replace(')','')
            p = p.replace("'", '')
            p=p.replace(' ','_')
            p=p.replace('-','_')
            p=p+'[c]'
            reacBank[reaction]['products'][p]=coeff
            reacBank[reaction]['mets'].append(p)
        #print(p)
        if p in modelMet:
            index=modelMet.index(p)
            #print(modelMet.index(p))
            #print(modelMet[index])
        else:
            pflag=False
    if sflag and pflag:
        if reaction in mreac:
            print('Reaction '+ reaction +' is already in the model')
            continue
        added.add(reaction)
        print('Add: ' + reaction)

# Second pass for finding reactions in the bank that can complement other reactions


metgo=False
cand={}
brakes={}
for i,reac in enumerate(bank):
    # Reactions may become availble for one step additions as more metabolites get connected rendering more reactions
    # as valid.
    if isValid(reac):
        if reac in mreac:
            print('Reaction '+ reac +' is already in the model')
            continue
        print('Add: '+ reac)
        added.add(reac)
    if reac in added:
        i-=1
        continue
    print('Checking reaction: '+reac)
    #prime=True
    rank=0
    connectReac(reac,reacBank, rank)

output=open('chains.txt','w')
pout=open('pathways.txt','w')
dpout=open('detailed_pathways.txt','w')
pways={}
pways['No pathways']=set()
# add candidates to the m file
for obj in cand:
    if cand[obj]:
        for item in cand[obj]:
            added=added|cand[obj][item]
            output.write('\t'.join(list(cand[obj][item])) + '\n')

            # collect pathway information for distribution
            for f in cand[obj][item]:
                if 'pathways' in reacBank[f]:
                    for item in reacBank[f]['pathways']:
                        # discard global and overview maps
                        if any(L in item for L in (
                        'rn01100', 'rn01110', 'rn01120', 'rn01130', 'rn01200', 'rn01210', 'rn01212', 'rn01230',
                        'rn01220')):
                            continue
                        Pname=reacBank[f]['pathways'][item]
                        if Pname in pways:
                            pways[Pname].add(f)
                        else:
                            pways[Pname]=set()
                            pways[Pname].add(f)
                else:
                    pways['No pathways'].add(f)
# collect pathways for networks
npout=open('npathways.txt','w')
npways={}
npways['No pathways']=set()
for tenant in neighbors:
    kegg_entry = kegg_con.parse(kegg_con.get(tenant))
    if isinstance(kegg_entry, int):
        continue
    if 'PATHWAY' in kegg_entry:
        pathways=kegg_entry['PATHWAY']
        for pn in pathways:
            # discard global and overview maps
            if any(L in pn for L in ('rn01100','rn01110','rn01120','rn01130','rn01200','rn01210','rn01212','rn01230','rn01220') ):
                continue
            name=pathways[pn]
            if name in npways:
                npways[name].add(tenant)
            else:
                npways[name] = set()
                npways[name].add(tenant)
    else:
        npways['No pathways'].add(tenant)
for q in pways:
    pout.write(q+'\t'+str(len(pways[q]))+'\n')
    for p in pways[q]:
        dpout.write(p + '\t'+ q+'\n')
for d in npways:
    npout.write(d+'\t'+str(len(npways[d]))+'\n')
    for p in npways[d]:
        dpout.write(p + '\t'+ d+'\n')
# create a list of added reactions
aout=open('added_reactions.txt','w')
gout=open('new_added_genes.txt','w')
rout=open('gprs.txt','w')
for item in added:
    hit= False
    if 'R00124' in item:
        continue
    # add reaction function into the output m file
    aout.write(item+'\n')
    g=reac2gene[item]
    for cg,gene in enumerate(g):
        if g.__len__()==1:
            rout.write('x(')
            if gene not in gl:
                gout.write(gene + '\n')
                gl.append(gene)
                print(str(gl.index(gene)+1)+ ' is ' + gene)
                rout.write(str(gl.index(gene)+1)+')\n')
            else:
                print(str(gl.index(gene)+1)+ ' is ' + gene)
                rout.write(str(gl.index(gene)+1) + ')\n')

        else:
            if not hit:
                rout.write('( x(')
                hit=True
            if gene not in gl:
                if cg==g.__len__()-1:
                    if gene not in gl:
                        gout.write(gene + '\n')
                        gl.append(gene)
                        print(str(gl.index(gene) + 1) + ' is ' + gene)
                        rout.write(str(gl.index(gene) + 1) + ') )\n')
                    else:
                        print(str(gl.index(gene) + 1) + ' is ' + gene)
                        rout.write(str(gl.index(gene) + 1) + ') )\n')
                else:
                    gout.write(gene +'\n')
                    gl.append(gene)
                    print(str(gl.index(gene)+1)+ ' is ' + gene)
                    rout.write(str(gl.index(gene)+1) + ') & x(')
            else:
                if cg==g.__len__()-1:
                    if gene not in gl:
                        gout.write(gene + '\n')
                        gl.append(gene)
                        print(str(gl.index(gene) + 1) + ' is ' + gene)
                        rout.write(str(gl.index(gene) + 1) + ') )\n')
                    else:
                        print(str(gl.index(gene) + 1) + ' is ' + gene)
                        rout.write(str(gl.index(gene) + 1) + ') )\n')
                else:
                    print(str(gl.index(gene)+1)+ ' is ' + gene)
                    rout.write(str(gl.index(gene)+1) + ') & x(')


    fout.write("model=addReaction(model,'" + item + "','metaboliteList', {")
    # add metabolites
    subs=reacBank[item]['substrates']
    end=subs.__len__()
    for counter,s in enumerate(subs):
        s=s.replace("'_","_")
        fout.write("'"+ s + "',")
    prods=reacBank[item]['products']
    end=prods.__len__()
    for counter,p in enumerate(prods):
        if counter==(end-1):
            p=p.replace("'_","_")
            fout.write("'" + p + "'},")
        else:
            p=p.replace("'_","_")
            fout.write("'"+p+"',")
    # Add the stoichcoeffList
    fout.write("'stoichCoeffList',[")
    subs=reacBank[item]['substrates']
    end=subs.__len__()
    for counter,s in enumerate(subs):
            fout.write(str(reacBank[item]['substrates'][s]) + ",")
    prods=reacBank[item]['products']
    end=prods.__len__()
    for counter,p in enumerate(prods):
        if counter==(end-1):
            fout.write(str(reacBank[item]['products'][p]) + "],")
        else:
            fout.write(str(reacBank[item]['products'][p])+",")
    # Set reversibility
    if '=' in reacBank[item]['direction']:
        fout.write("'reversible', true);\n")
    else:
        fout.write("'reversible', false);\n")


output.close()
loopout.close()
fout.close()
pout.close()
netout.close()
aout.close()
gout.close()
rout.close()
# print out missing mets
misout=open('breaking_bad_mets.txt','w')
for m in missing:
    misout.write(m+'\n')

misout.close()
print('Done')
