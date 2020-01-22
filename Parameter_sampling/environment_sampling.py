#!/usr/bin/env python

"""
Sample environments

Author: Snorre Sulheim
Date: January 15, 2020

# Description
This script replicates the sampling of fluxes by Shany, however with exploring the space to check the robustness of the results



"""

import cobra
import pandas as pd
import numpy as np
from cobra.sampling import OptGPSampler
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

PARAMETER_VALUES = [#"Name",     "Reaction ID",          "Lower bound", "UpperBound", "Glycogen production mean"
                   ["HCO3",      "HCO3EXcar",            -1.35,          -0.35, -1.3], 
                   ["Nitrogen",  "AmmoniaEX",            -0.56,           0,    -0.0263],
                   ["Phosphate", "FAKEOrthophosphateEX", -0.227,          0,    -0.0152],
                   ["Light",     "LightEX",              -181,          -13.4,  -118.6]]
#                   ["Rubisco",   "R00024",                0.27,           4.7],-1.248
#                   ["CO2",       "CO2EX",                 0,              0.818,-0.22]]


def random_sampling(model_fn, iterations = 3, kind = "pFBA", reaction_ids = None, 
                    remove_blocked = True, processes = 2):
    """
    Function to randomly sample different environment at assess it using one of the selected methods
    """
    # Read model file
    model = cobra.io.read_sbml_model(model_fn)

    if not reaction_ids:
        reaction_ids = [r.id for r in model.exchanges] + ["R00024"]

    Nr = len(reaction_ids)
    flux_arr = np.zeros((Nr, iterations))
    flux_arr_min = np.zeros((Nr, iterations))
    flux_arr_max = np.zeros((Nr, iterations))

    if remove_blocked:
        _remove_blocked(model, open_exchanges = True)


    if kind == "pFBA":
        cobra.flux_analysis.parsimonious.add_pfba(model)
        for i in range(iterations):
            state = "infeasible"
            while state == "infeasible":
                _sample_environment(model)
                solution = model.optimize()
                state = solution.status
            flux_arr[:, i] = solution[reaction_ids]

    if kind == "FBA":
        for i in range(iterations):
            state = "infeasible"
            while state == "infeasible":
                _sample_environment(model)
                solution = model.optimize()
                state = solution.status
            flux_arr[:, i] = solution[reaction_ids]

    elif kind == "FVA":
        for i in range(iterations):
            print("{0}%".format(100*i/iterations))
            state = "infeasible"
            while state == "infeasible":
                _sample_environment(model)
                solution = model.optimize()
                state = solution.status
            df_fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_ids, fraction_of_optimum = 1)
            flux_arr[:, i] = df_fva.mean(axis = 1)
            flux_arr_min[:, i] = df_fva.min(axis = 1)
            flux_arr_max[:, i] = df_fva.max(axis = 1)

    elif kind == "random_sampling":
        for i in range(iterations):
            print("{0}%".format(100*i/iterations))
            state = "infeasible"
            while state == "infeasible":
                _sample_environment(model)
                solution = model.optimize()
                state = solution.status

            # Run the actual sampling
            model.reactions.BIOMASS.lower_bound = 0.99*solution.objective_value
            optgp = OptGPSampler(model, processes=processes)
            rs = optgp.sample(1000)
            model.reactions.BIOMASS.lower_bound = 0
            flux_arr[:, i] = rs[reaction_ids].mean(axis = 0)
            flux_arr_min[:, i] = rs[reaction_ids].min(axis = 1)
            flux_arr_max[:, i] = rs[reaction_ids].max(axis = 1)

    result_df = pd.DataFrame(flux_arr)
    result_df.index = reaction_ids
    result_df.index.name = "Reactions"
    # Save to file
    fn = "random_sampling_{0}_{1}.csv".format(kind, iterations)
    result_df.to_csv(fn)

    if kind in ["FVA", "random_sampling"]:
        fn_min = "min_" + fn
        result_df_min = pd.DataFrame(flux_arr_min)
        result_df_min.index = reaction_ids
        result_df_min.index.name = "Reactions"
        result_df_min.to_csv(fn_min)

        fn_max = "max_" + fn
        result_df_max = pd.DataFrame(flux_arr_max)
        result_df_max.index = reaction_ids
        result_df_max.index.name = "Reactions"
        result_df_max.to_csv(fn_max)

def _remove_blocked(model, open_exchanges = True):
    #_set_to_mean(model, lower_only = True)
    blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges = open_exchanges)
    print(len(blocked), blocked)
    model.remove_reactions([model.reactions.get_by_id(r_id) for r_id in blocked])

def _sample_environment(model):
    for row in PARAMETER_VALUES:
        flux = np.random.uniform(row[2], row[3])
        r = model.reactions.get_by_id(row[1])
        r.bounds = (flux, flux)

def _set_to_mean(model, lower_only = False):
    for row in PARAMETER_VALUES:
        #flux = (row[2] + row[3])/ 2.0
        r = model.reactions.get_by_id(row[1])
        if lower_only:
            flux = row[2]
            r.bounds = (flux, 0)
        else:
            flux = row[4]
            r.bounds = (flux, flux)

def compare_FVA_FBA(csv_fn_FBA, csv_fn_FVA_min, csv_fn_FVA_max):
    df_FBA = pd.read_csv(csv_fn_FBA, index_col = 0).T
    df_FVA_min = pd.read_csv(csv_fn_FVA_min, index_col = 0).T
    df_FVA_max = pd.read_csv(csv_fn_FVA_max, index_col = 0).T

    n_FBA = sum(df_FBA["GlycogenEX"]> 1e-7)
    n_FVA_max = sum(df_FVA_max["GlycogenEX"]> 1e-7)
    n_FVA_min = sum(df_FVA_min["GlycogenEX"]> 1e-7)

    fig, ax = plt.subplots(1)
    sns.kdeplot(df_FBA.loc[df_FBA["GlycogenEX"]> 1e-7, "LightEX"], ax = ax, label = "FBA", bw=2, lw = 2, shade = True)
    sns.kdeplot(df_FVA_max.loc[df_FVA_max["GlycogenEX"]> 1e-7, "LightEX"], ax = ax, label = "Max FVA", bw=2, lw = 2, shade = True)
    sns.kdeplot(df_FVA_max.loc[df_FVA_min["GlycogenEX"]> 1e-7, "LightEX"], ax = ax, label = "Min FVA", bw=2, lw = 2, shade = True)
    ax.set_xlabel("LightEX")
    ax.set_ylabel("Kernel density estimate of non-zero glycogen samples")
    plt.tight_layout()
    plt.show()

    print(n_FBA, n_FVA_max, n_FVA_min)



    # print(df_FBA.columns)
    # sns.jointplot(x = "LightEX", y = "GlycogenEX", data = df_FVA_max, ax = ax[0], label = "Max FVA")
    # sns.jointplot(x = "LightEX", y = "GlycogenEX", data = df_FVA_min, ax = ax[1], label = "Min FVA")  
    # sns.jointplot(x = "LightEX", y = "GlycogenEX", data = df_FBA, ax = ax[2], label = "FBA")
    # plt.show()
    
def plot_sampling_results(csv_fn, reactions):
    df = pd.read_csv(csv_fn, index_col = 0).T
    print(df.columns)
    #sns.pairplot()

    # Explore glycogen-producing samples
    df_nz = df.loc[df["GlycogenEX"] > 1e-10, :]

    sns.jointplot(x=df_nz["GlycogenEX"], y=df_nz["LightEX"], kind="hex")
    plt.tight_layout()
    plt.show()

    sns.jointplot(x=df_nz["GlycogenEX"], y=df_nz["AmmoniaEX"], kind="hex")
    plt.tight_layout()
    plt.show()
    
    sns.jointplot(x=df_nz["GlycogenEX"], y=df_nz["HCO3EXcar"], kind="hex")
    plt.tight_layout()
    plt.show()
    
    sns.jointplot(x=df_nz["GlycogenEX"], y=df_nz["BIOMASS"], kind="hex")
    plt.tight_layout()
    plt.show()

    # # Explore non glycogen-producing samples
    # df_z = df.loc[df["GlycogenEX"] < 1e-10, :]
    # # sns.jointplot(x=df_nz["AmmoniaEX"], y=df_nz["LightEX"], kind="hex")
    # g = sns.PairGrid(df_z[["LightEX", "BIOMASS", "AmmoniaEX", "HCO3EXcar"]])
    # g.map_diag(sns.distplot)
    # g.map_offdiag(sns.scatterplot);
    # plt.tight_layout()
    # plt.show()

    # # Explore all samples
    # g = sns.PairGrid(df[["LightEX", "BIOMASS", "AmmoniaEX", "HCO3EXcar"]])
    # g.map_diag(sns.kdeplot)
    # g.map_offdiag(sns.kdeplot, n_levels=6);
    # plt.tight_layout()
    # plt.show()

def pairplots_glycogen(csv_fn):
    df = pd.read_csv(csv_fn, index_col = 0).T
    # Explore glycogen-producing samples
    # df_z = df.loc[df["GlycogenEX"] > 1e-10, :]
    df["Glycogen production"] = (df["GlycogenEX"] > 1e-7).astype(str)
    reactions = ["LightEX", "BIOMASS", "AmmoniaEX", "HCO3EXcar", "Glycogen production"]

    g = sns.PairGrid(df[reactions], hue = "Glycogen production")
    g.map_diag(plt.hist, bins = 20, alpha = 0.5)
    g.map_lower(sns.kdeplot, alpha = 0.5)
    g.map_upper(plt.scatter, alpha = 0.3, s = 0.2)
    g.fig.suptitle("Environmental samples")
    # g.add_legend()
    # plt.tight_layout()
    plt.show()

    # # Explore non glycogen-producing samples
    # df_nz = df.loc[df["GlycogenEX"] < 1e-10, :]
    # g = sns.PairGrid(df_nz[reactions])
    # g.map_diag(sns.distplot)
    # g.map_offdiag(sns.kdeplot)
    # g.fig.suptitle("Glycogen producing samples")
    # plt.tight_layout()
    # plt.show()

    # # Explore all samples
    # g = sns.PairGrid(df[reactions])
    # g.map_diag(sns.distplot)
    # g.map_offdiag(sns.kdeplot)
    # g.fig.suptitle("All samples")
    # plt.tight_layout()
    # plt.show()

if __name__ == '__main__':
    # model_fn = "../Models/Prochlorococcus/iSO595_MarinePro.xml"
    model_fn = "../Models/Prochlorococcus/iSO595/Model_files/iSO595c1.xml"

    reaction_ids = ["GlycogenEX", "FumarateEX", "SuccinateEX", "PhospholipidEX", "PyruvateEX", "CitrateEX", 
                    "AcetateEX", "FormateEX", "4_Methyl_2_oxopentanoateEX", "HCO3EXcar", "CO2EX",
                    "R00024", "BIOMASS", "AmmoniaEX", "FAKEOrthophosphateEX", "LightEX"]

    if 0:
        model = cobra.io.read_sbml_model(model_fn)
        _set_to_mean(model)
        solution = model.optimize()
        print(model.summary(threshold = 1e-6))
        df_fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_ids, fraction_of_optimum = 1)
        print(df_fva)
    if 0:
        random_sampling(model_fn, reaction_ids = reaction_ids, iterations = 1, kind = "FBA", remove_blocked = True)

    if 0:
        parser = argparse.ArgumentParser(description='Random environment sampling')
        parser.add_argument('-i', "--iterations", help = "Number of iterations", default=1000, type = int)
        parser.add_argument('-k', "--kind", help = "FVA, random_sampling or pFBA", default = "pFBA")
        parser.add_argument('-f', "--filename", help = "model_filename", default = "iSO595c1.xml")
        parser.add_argument('-p', "--processes", help = "Number of processes", default = 4, type = int)
        args = parser.parse_args()
        random_sampling(args.filename, reaction_ids = reaction_ids, iterations = args.iterations, kind = args.kind, processes = args.processes)
    if 0:
        results_fn  = "../Data/random_sampling/max_random_sampling_FVA_100000.csv"
        plot_sampling_results(results_fn)

    if 1:
        results_fn  = "../Data/random_sampling/max_random_sampling_FVA_100000.csv"
        pairplots_glycogen(results_fn)

    if 0:
        FBA_fn =  "../Data/random_sampling/random_sampling_FBA_100000.csv"
        FVA_fn_min =  "../Data/random_sampling/min_random_sampling_FVA_100000.csv"
        FVA_fn_max =  "../Data/random_sampling/max_random_sampling_FVA_100000.csv"
        compare_FVA_FBA(FBA_fn, FVA_fn_min, FVA_fn_max)