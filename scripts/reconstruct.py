# -*- coding: utf-8 -*-
"""
Author: Snorre Sulheim
Email: snorre.sulheim@sintef.no
Date: 11/9/2019

This script contains functions used to edit and improve the genome-scale model of *Prochlorococcus marinus med*.

"""
import cobra
import re
import pandas as pd

# Make pandas print full dataframe
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def identify_duplicate_reactions(model, reaction_list = None):
    """
    Iterate over all reactions in a model to find duplciates based on metabolites and reactions bounds.

    Parameters
    ----------
    model: cobrapy model object
        The genome-scale metabolic model

    Notes
    -----
    From the model iSO595c1.xml at 11/9/2019 the following duplicated reactions were found:
        - R00546 and R00488
        - R00762 and R04780
        - R00959 and R08639
        - R01070 and R01068
        - R02736 and R00835
        - SuccinateEX and ThymidineEX
        - DethiobiotinEX and GuanineEX
    
    """
    if not reaction_list:
        reaction_list = [r for r in model.reactions]
    duplicate_reactions_list = []
    for i in range(len(reaction_list)):
        print(i, reaction_list[i].id)
        for j in range(i+1, len(reaction_list)):
            r1 = reaction_list[i]
            r2 = reaction_list[j]
            equal = _compare_reactions(r1, r2)
            if equal:
                duplicate_reactions_list.append([r1.id, r2.id])
    return duplicate_reactions_list

def _compare_reactions(r1, r2):
    """
    Compare if Reaction 1 is equal to Reaction 2 with respect to metabolites and reaction bounds.
    
    Parameters
    ----------
    r1: Cobra.Reaction,
        Reaction 1
    r2: Cobra.Reaction.
        Reaction 2
    """
    # If the number of metabolites are different they are not equal
    if len(r1.metabolites) != len(r2.metabolites):
        return False

    # If the compartments are different they are not equal
    if r1.compartments != r2.compartments:
        return False

    # Since the reactions can be defined in opposite directions we need to test twice
    for i in range(2):
        if i == 1:
            # Revert reaction
            r1.add_metabolites({m:-2*i for m,i in r1.metabolites.items()})
            # Revert reaction bounds
            r1.bounds = (r1.bounds[1]*-1, r1.bounds[0]*-1)
    
        if r1.metabolites == r2.metabolites:
            if r1.bounds == r2.bounds:
                return True
    return False

def delete_duplicate_reactions(model, reaction_tuple_list, save_fn = False):
    for keep, delete, grr in reaction_tuple_list:
        r_keep = model.reactions.get_by_id(keep)
        r_delete = model.reactions.get_by_id(delete)

        # Concatenate annotations
        for key, value in r_delete.annotation.items():
            try:
                existing_annotation = r_keep.annotation[key]
            except KeyError:
                r_keep.annotation[key] = value
            else:
                r_keep.annotation[key] = [existing_annotation, value]

        
        # Delete reaction
        print("Deleted {0}, duplicate of {1}.".format(r_delete.id, r_keep.id))
        r_delete.remove_from_model()

        # Concatenate genes
        print("Changed gene reaction rule of {0} from {1} to {2}".format(r_keep.id, r_keep.gene_reaction_rule, grr))
        r_keep.gene_reaction_rule = grr

    if save_fn:
        _export(model, save_fn)


def add_KEGG_annotations_from_IDs(model, save_fn = None):
    """
    Many of the reactions have the KEGG ID as the reaction identifier. Add these as annotations as well.
    
    Parameters
    ----------
    model: cobrapy model object
        The genome-scale metabolic model

    """

    for r in model.reactions:
        match = re.fullmatch("R\\d{5}", r.id)
        if match:
            r.annotation["kegg.reaction"] = r.id

    if save_fn:
        cobra.io.write_sbml_model(model, save_fn)

def add_metabolite_annotations_from_spreadsheet(model, spreadsheet_fn, sheetname = "METS", save_fn = None):
    """
    Reaction and metabolite identifiers are given in the spreadsheet (SI File 1) for the preceeding model iJC568 (https://dx.doi.org/10.1128/mSystems.00065-16). Add these annotations to the model.

    Parameters
    ----------
    model: cobrapy model object
        The genome-scale metabolic model
    spreadsheet_fn: str
        Path to the spreadsheet with model annotations
    """
    df = pd.read_excel(spreadsheet_fn, sheetname)
    
    n_kegg = 0
    n_inchi = 0
    n_pubchem = 0
    n_charge = 0
    n_formula = 0

    for i, row in df.iterrows():
        # Skip metabolites in spreadsheet from the "Boundary-compartment" b. 
        # These are not in the current model
        if row["COMPARTMENT"] == "b":
            continue

        # Get id from spreadsheet
        # For some of the metabolites the compartment is in the ID string, but for most they are not
        if row["ID"].rsplit("_")[-1] == row["COMPARTMENT"]:
            met_id = row["ID"]
        else:
            met_id = "{0}_{1}".format(row["ID"], row["COMPARTMENT"])

        try:
            m = model.metabolites.get_by_id(met_id)
        except KeyError:
            print("Metabolite {0} is not in the model".format(met_id))
        else:
            # Add annotations
            ## KEGG
            if isinstance(row["KEGG CPD"], str):
                kegg_match = re.search("C\\d{5}", row["KEGG CPD"])
                if kegg_match:
                    m.annotation["kegg.compound"] = kegg_match.group(0)
                    n_kegg+=1

            ## INCHI
            # The inchi key is a unique 25-character string that defines the molecular structur. 
            # There is also two dashes in the key making the total length 27
            if isinstance(row["InChiKey"], str) and len(row["InChiKey"]) == 27:
                m.annotation["inchikey"] = row["InChiKey"]
                n_inchi += 1

            ## PubChem
            if not pd.isna(row["PubChem CID"]):
                m.annotation["pubchem.substance"] = str(int(row["PubChem CID"]))
                n_pubchem += 1

            ## Molecular Charge
            if not pd.isna(row["Charge"]):
                m.charge = int(row["Charge"])
                n_charge += 1

            ## Chemical formula
            if isinstance(row["COMPOSITION"], str) and not len(m.formula) and len(row["COMPOSITION"]):
                m.formula = row["COMPOSITION"]
                n_formula += 1


    print("Added KEGG ID to {0} metabolites".format(n_kegg))
    print("Added InChi Key to {0} metabolites".format(n_inchi))
    print("Added PubChem CID to {0} metabolites".format(n_pubchem))
    print("Added molecular charge to {0} metabolites".format(n_charge))
    print("Added chemical formula to {0} metabolites".format(n_formula))

    if save_fn:
        _export(model, save_fn)

def _export(model, save_fn, yml_copy = True):
    cobra.io.write_sbml_model(model, save_fn)
    if yml_copy:
        yml_fn = save_fn.rsplit(".xml")[0] + ".yml"
        cobra.io.save_yaml_model(model, yml_fn)

if __name__ == '__main__':
    from pathlib import Path

    REPO_PATH = Path.cwd().parent
    print(REPO_PATH)
    model_path = REPO_PATH / "Model_files" / "iSO595c1.xml"
    model = cobra.io.read_sbml_model(str(model_path))
    if 0:
        # Find duplicate reactions
        duplicate_reactions = identify_duplicate_reactions(model)
        print("Duplicate reactions:")
        for lst in duplicate_reactions:
            print(lst)

    if 0:
        # Give KEGG annotations
        add_KEGG_annotations_from_IDs(model, str(model_path))

    if 0:
        # Get metabolite annotations from model iJC568 (preceeding model) 
        spreadsheet_fn = "C:\\Users\\snorres\\OneDrive - SINTEF\\Boston University\\Marine strains project\\GEM\\template_GEM_iJC568\\SI_File_1_iJC568_Excel.xls"
        add_metabolite_annotations_from_spreadsheet(model, spreadsheet_fn, save_fn = str(model_path))

    if 1:
        #                            Keep       Delete  Gene rule
        duplicate_reactions_list = [("R00488", "R00546", "PMM0222"),
                                    ("R00762", "R04780", "PMM0781 or PMM0767"),
                                    ("R08639", "R00959", "PMM0076 or PMM0278"),
                                    ("R01068", "R01070", "PMM0767 or PMM0781"),
                                    ("R00835", "R02736", "PMM1074 or PMM0771")]
        delete_duplicate_reactions(model, duplicate_reactions_list, save_fn = str(model_path))

