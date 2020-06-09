#!/usr/bin/env python

"""
Sample environments

Author: Snorre Sulheim / Shany Ofaim
Date: June 9, 2020
Email: snorre.sulheim@sintef.no
Affiliation: SINTEF; NTNU; Boston University

# Description
This script samples different environment within the range defined in the global variable PARAMETER_VALUES. These ranges are within biologically relevant ranges. The values are sampled using the "push-FBA" framework for light and bicarbonate, replicating the physiologically relevant conditions where the organism has to cope with excess light and bicarbonate since light absorbtion cannot be readily shut off.
    - Fixed influx rates are sampled for light and bicarbonate
    - CO2 efflux is limited to the maximal measurable value (0.818)
    - Upper bounds are sampled for phosphate and nitrogen uptake

The sampling runs FVA, FBA and pFBA.

"""

import cobra
import pandas as pd
import numpy as np
import argparse
import time

PARAMETER_VALUES = [#"Name",     "Reaction ID",          "Lower bound", "UpperBound"
                   ["HCO3",      "HCO3EXcar",            -8,            0],
                   ["Nitrogen",  "AmmoniaEX",            -0.56,         0],
                   ["Phosphate", "FAKEOrthophosphateEX", -0.1,          0],
                   ["Light",     "LightEX",              -150,          0]]
CO2MAX = 0.82

# Block fake reactions
FAKE_TRANSPORT = ["AminosugarsOUT", "FAKEAAOUT", "FAKEABPOUT", "FAKEacpTRANS", "FAKEApoacpTRANS", "FAKEThioredoxinTRANS", "FreefattyacidsOUT", "7NMeth7carbOUT", "ArtificialproteinOUT", "FADOUT", "LipoylproteinTRANS", "MenaquinoneOUT", "NicotinateOUT", "THFpolyglutOUT", "Thiamin_dpOUT"]

def random_sampling(model_fn, iterations = 3, reaction_ids = None, remove_blocked = True):
    """
    Function to randomly sample different environment at assess it using one of the selected methods
    """
    # Log time
    t0 = time.time()

    # Read model file
    model = cobra.io.read_sbml_model(model_fn)
    model.solver = 'gurobi'
    # Block H2S
    model.reactions.H2SEX.lower_bound = 0
    
    # Block fake transports
    for rid in FAKE_TRANSPORT:
        model.reactions.get_by_id(rid).bounds = (0,0)

    # Remove blocked reactions
    if remove_blocked:
        _remove_blocked(model, open_exchanges = True)

    # Block maximum CO2 production
    model.reactions.CO2EX.bounds = (0, CO2MAX)

    if not reaction_ids:
        reaction_ids = [r.id for r in model.exchanges] + ["R00024"]

    Nr = len(reaction_ids)
    flux_arr_FBA  = np.zeros((Nr, iterations))
    flux_arr_pFBA = np.zeros((Nr, iterations))
    flux_arr_min  = np.zeros((Nr, iterations))
    flux_arr_max  = np.zeros((Nr, iterations))
    environment_arr = np.zeros((4, iterations))
    total_samples = 0
  
    for i in range(iterations):
        print("{0}%".format(100*i/iterations))
        state = "infeasible"
        while state == "infeasible":
            total_samples += 1
            environment = sample_environment(model)
            solution = model.optimize()
            state = solution.status
        # Store the sampled values
        environment_arr[:, i] = environment

        # run pFBA
        with model:
            cobra.flux_analysis.parsimonious.add_pfba(model)
            pFBA_solution = model.optimize()

        df_fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_ids, fraction_of_optimum = 1)
        flux_arr_FBA[:, i] = solution[reaction_ids]
        flux_arr_pFBA[:,i] = pFBA_solution[reaction_ids]
        flux_arr_min[:, i] = df_fva.min(axis = 1)
        flux_arr_max[:, i] = df_fva.max(axis = 1)



    # Save to file
    timestr = time.strftime("%Y%m%d-%H%M")
    fn = "random_sampling_{0}_{1}.csv".format(iterations, timestr)

    
    fn_fba = "FBA_"+fn
    result_df = pd.DataFrame(flux_arr_FBA)
    result_df.index = reaction_ids
    result_df.index.name = "Reactions"
    result_df.to_csv(fn_fba)

    fn_pfba = "pFBA_"+fn
    result_df = pd.DataFrame(flux_arr_pFBA)
    result_df.index = reaction_ids
    result_df.index.name = "Reactions"
    result_df.to_csv(fn_pfba)

    fn_min = "minFVA_" + fn
    result_df_min = pd.DataFrame(flux_arr_min)
    result_df_min.index = reaction_ids
    result_df_min.index.name = "Reactions"
    result_df_min.to_csv(fn_min)

    fn_max = "maxFVA_" + fn
    result_df_max = pd.DataFrame(flux_arr_max)
    result_df_max.index = reaction_ids
    result_df_max.index.name = "Reactions"
    result_df_max.to_csv(fn_max)

    fn_env = "env_" + fn
    env_df = pd.DataFrame(environment_arr)
    env_df.index = [row[0] for row in PARAMETER_VALUES]
    env_df.index.name = "Exchanges"
    env_df.to_csv(fn_env)

    print("Time elapsed: {0:0.1f} hours".format((time.time()-t0)/3600))
    print("{0} feasible samples, out of {1} total".format(iterations, total_samples))

def _remove_blocked(model, open_exchanges = True):
    blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges = open_exchanges)
    print(len(blocked), blocked)
    model.remove_reactions([model.reactions.get_by_id(r_id) for r_id in blocked])

def sample_environment(model):
    environment = np.zeros(len(PARAMETER_VALUES))
    for i, row in enumerate(PARAMETER_VALUES):
        # Row: Name, Reaction ID, lower bound, upper bound
        key = row[0]
        reaction_id = row[1]
        lower_bound = row[2]
        upper_bound = row[3]
        flux = np.random.uniform(lower_bound, upper_bound)
        r = model.reactions.get_by_id(reaction_id)
        if key in ["HCO3", "Light"]:
            # Fix flux
            r.bounds = (flux, flux)
        elif key in ["Nitrogen", "Phosphate"]:
            # Set lower_bound only (i.e. define maximum uptake rate)
            r.bounds = (flux, 0)
        else:
            raise KeyError

        environment[i] = flux
    return environment


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Random sampling of environments')
    parser.add_argument('-i', "--iterations", help = "Number of iterations", default=1000, type = int)
    parser.add_argument('-f', "--filename", help = "Model filename", default = "iSO595v5.xml")
    args = parser.parse_args()
    random_sampling(args.filename, reaction_ids = None, iterations = args.iterations)
