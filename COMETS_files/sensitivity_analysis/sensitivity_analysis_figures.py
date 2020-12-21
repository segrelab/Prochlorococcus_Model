#!/usr/bin/env python

"""
Create figures for the sensitivity analysis for the Prochlorococcus manuscript 

Author: Snorre Sulheim
Date: December 14, 2020
Email: snorre.sulheim@sintef.no

"""
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 16})

# my_cmap = sns.hls_palette(n_colors = 10, as_cmap=True)
time_delta = 0.1
NIGHT_FACECOLOR = "#A9A9A9"

glycogen_mol_weight = 162.1406
cell_gDW = 66*1e-15


# RXN dict
RXN_NAME_DICT = {"R00024": "CBB (RuBisCO)",
                 "R05605": "ED (KDPG aldolase)",
                 "R02736": "PPP (G6PDH)",
                 "R00658": "Glycolysis (Enolase)",
                 # "R00762", 
                 "R01641": "PPP / CBB\n(Transketolase)",
                 "PyruvateEX": "Pyruvate\nexudation",
                 "FormateEX": "Formate\nexudation",
                 "GlycogenEX": "Glycogen storage"
                 }
def create_death_figures_growth_glycogen():
    matplotlib.rcParams.update({'font.size': 12})
    folder = Path("death/simulations")
    fig, [ax1, ax2] = plt.subplots(1,2, sharex = True, figsize = (6.8, 3))
    for i, results_folder in enumerate(folder.glob("*")):
        death_rate = str(results_folder).rsplit("_")[-1]
        print(results_folder)
        results_fn = results_folder / "results_table.csv"
        df_results = pd.read_csv(results_fn)
        print(df_results)
        ax1.plot(df_results["Time"], df_results["Biomass"], label = death_rate, lw = 1)
        
        # Convert to mg glycogen / gDW
        glycogen_n_lim = df1["Glycogen"]*1e-3*glycogen_mol_weight*cell_gDW*1e15/df1["Biomass"]
        ax2.plot(df_results["Time"], glycogen_n_lim, label = death_rate, lw = 1)
        # ax2.plot(df_results["Time"], df_results["Glycogen"]/df_results["Biomass"], label = death_rate, lw = 1)

    ax1.set_xlabel("Time [h]")
    ax2.set_xlabel("Time [h]")
    ax1.set_ylabel("Biomass [gDW]")
    # ax2.set_ylabel("Glycogen per biomass \n[mmol gDW$^{-1}$]")
    ax2.set_ylabel("Glycogen [fg cell$^{-1}$]")
    time = list(df_results["Time"])
    
    for ax in [ax1, ax2]:
        ax.axvspan(12,24, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(36,48, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(60,72, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(84,96, facecolor = NIGHT_FACECOLOR, alpha = 0.5)

    ax1.set_xlim(time[0], time[-1])
    ax2.set_xlim(time[0], time[-1])

    plt.legend(title = r"Death rate [d$^{-1}$]", ncol = 2)
    plt.tight_layout()
    plt.savefig("death_rate_growth_glycogen.svg")
    plt.close(fig)

def calculate_death_rate_deviations():
    matplotlib.rcParams.update({'font.size': 11})
    folder = Path("death/simulations")
    reactions = ["R00024", "R05605", "R02736", "R00658","GlycogenEX",# "R00762", 
                 "R01641", "PyruvateEX", "FormateEX"]
    flux_list = []
    death_rates = []
    for i, results_folder in enumerate(folder.glob("*")):
        death_rate = str(results_folder).rsplit("_")[-1]
        death_rates.append(death_rate)
        if i == 0:
            # Find rxn indexes
            rxn_names_fn = results_folder / "rxnIds.csv"
            rxn_df = pd.read_csv(rxn_names_fn, header = None)
            rxn_df.index += 1
            rxn_df.columns = ["Reaction name"]
            reaction_idxs = []
            reaction_names = []
            for j in rxn_df.index:
                if rxn_df.loc[j, "Reaction name"] in reactions:
                    reaction_idxs.append(j)
                    reaction_names.append(rxn_df.loc[j, "Reaction name"])

        # Get flux data
        death_rate = str(results_folder).rsplit("_")[-1]
        fn = results_folder / "fluxTable.csv"
        df = pd.read_csv(fn)



        for j, r_idx in enumerate(reaction_idxs):
            df_j = df[df["rxn"]==r_idx]
            reaction_name = rxn_df.loc[r_idx, "Reaction name"]
            lst = [death_rate, reaction_name]+list(df_j["flux"])
            flux_list.append(lst)

    full_df = pd.DataFrame(flux_list)
    hours = df_j["t"].values[:480]*time_delta
    fig, axes = plt.subplots(3,3, sharex = True, figsize = (6.6, 5))
    axes = axes.flatten()
    for j, rxn in enumerate(reaction_names):
        rx_df = full_df[full_df[1]==rxn]
        rx_ref = rx_df.loc[rx_df[0].astype(float)==0.10, 2:].values.T[:480]
        ax = axes[j]
        ax.set_title(RXN_NAME_DICT[rxn])
        ax.set_xlim(0,48)
        # if j in [3]:
        #     ax.set_ylabel("Flux [mmol gDW$^{-1}$ h$^{-1}$]")
        if j in [5,6,7]:
            ax.set_xlabel("Time [h]")

        ax.axvspan(12,24, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(36,48, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        print(death_rates)
        for k, d in enumerate(death_rates):
            d = float(d)
            # if d!=0.1:
            rx_d_df = rx_df.loc[rx_df[0].astype(float)==d, 2:].values.T[:480]
            diff = rx_d_df-rx_ref
            ax.plot(hours, rx_d_df, label = d, lw = 0.8)
            print(rxn , d)
            print(np.sum(np.abs(diff))/np.sum(np.abs(rx_ref)))

    plt.tight_layout()
    handles, labels = axes[0].get_legend_handles_labels()
    plt.legend(handles, labels, ncol = 2)
    # plt.show()
    plt.savefig("death_rate_deviation_fluxes_3_3.svg")

def calculate_glycogen_deviations():
    matplotlib.rcParams.update({'font.size': 11})
    folder = Path("glycogen/simulations")
    reactions = ["R00024", "R05605", "R02736", "R00658","GlycogenEX",# "R00762", 
                 "R01641", "PyruvateEX", "FormateEX"]
    flux_list = []
    labels = []
    i = 0
    for results_folder in folder.glob("*"):
        info = str(results_folder).split("_")
        if not len(info)>3:
            continue
        label = "Km: {0}; Vmax: {1}".format(info[1], info[3])
        labels.append(label)
        if i == 0:
            # Find rxn indexes
            rxn_names_fn = results_folder / "rxnIds.csv"
            rxn_df = pd.read_csv(rxn_names_fn, header = None)
            rxn_df.index += 1
            rxn_df.columns = ["Reaction name"]
            reaction_idxs = []
            reaction_names = []
            for j in rxn_df.index:
                if rxn_df.loc[j, "Reaction name"] in reactions:
                    reaction_idxs.append(j)
                    reaction_names.append(rxn_df.loc[j, "Reaction name"])

        # Get flux data
        fn = results_folder / "fluxTable.csv"
        df = pd.read_csv(fn)

        for j, r_idx in enumerate(reaction_idxs):
            df_j = df[df["rxn"]==r_idx]
            reaction_name = rxn_df.loc[r_idx, "Reaction name"]
            lst = [info[1], info[3], label, reaction_name]+list(df_j["flux"])
            flux_list.append(lst)

    full_df = pd.DataFrame(flux_list)
    hours = df_j["t"].values[:480]*time_delta
    fig, axes = plt.subplots(3,3, sharex = True, figsize = (6.6, 5))
    axes = axes.flatten()
    for j, rxn in enumerate(reaction_names):
        rx_df = full_df[full_df[3]==rxn]
        ax = axes[j]
        ax.set_title(RXN_NAME_DICT[rxn])
        ax.set_xlim(0,48)
        # if j in [3]:
        #     ax.set_ylabel("Flux [mmol gDW$^{-1}$ h$^{-1}$]")
        if j in [5,6,7]:
            ax.set_xlabel("Time [h]")

        ax.axvspan(12,24, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(36,48, facecolor = NIGHT_FACECOLOR, alpha = 0.5)

        for k, l in enumerate(labels):
            rx_d_df = rx_df.loc[rx_df[2]==l, 4:].values.T[:480]
            ax.plot(hours, rx_d_df, label = l, lw = 0.8)

    plt.tight_layout()
    handles, labels = axes[0].get_legend_handles_labels()
    plt.legend(handles, labels, ncol = 2)
    # plt.show()
    plt.savefig("glycogen_deviation_fluxes_3_3.svg")


def create_glycogen_figures_growth_glycogen():
    matplotlib.rcParams.update({'font.size': 12})
    folder = Path("glycogen/simulations")
    fig, [ax1, ax2] = plt.subplots(1,2, sharex = True, figsize = (6.8, 3))
    for i, results_folder in enumerate(folder.glob("*")):
        info = str(results_folder).split("_")
        if not len(info)>3:
            continue
        label = "Km: {0}; Vmax: {1}".format(info[1], info[3])
        print(results_folder)
        results_fn = results_folder / "results_table.csv"
        df_results = pd.read_csv(results_fn)
        print(df_results)
        ax1.plot(df_results["Time"], df_results["Biomass"], label = label, lw = 1)
        # Convert to mg glycogen / gDW
        glycogen_n_lim = df1["Glycogen"]*1e-3*glycogen_mol_weight*cell_gDW*1e15/df1["Biomass"]
        ax2.plot(df_results["Time"], glycogen_n_lim, label = death_rate, lw = 1)
        # ax2.plot(df_results["Time"], df_results["Glycogen"]/df_results["Biomass"], label = label, lw = 1)
    ax1.set_xlabel("Time [h]")
    ax2.set_xlabel("Time [h]")
    ax1.set_ylabel("Biomass [gDW]")
    # ax2.set_ylabel("Glycogen per biomass \n[mmol gDW$^{-1}$]")
    ax2.set_ylabel("Glycogen [fg cell$^{-1}$]")
    time = list(df_results["Time"])
    
    for ax in [ax1, ax2]:
        ax.axvspan(12,24, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(36,48, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(60,72, facecolor = NIGHT_FACECOLOR, alpha = 0.5)
        ax.axvspan(84,96, facecolor = NIGHT_FACECOLOR, alpha = 0.5)

    ax1.set_xlim(time[0], time[-1])
    ax2.set_xlim(time[0], time[-1])

    plt.legend(title = "Glycogen\ncoefficients", ncol = 1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("glycogen_growth_glycogen.svg")
    plt.close(fig)


def MM(Vmax, Km, concentration):
        return Vmax*(concentration)/(Km + concentration)

def light(I, X, a_b, a_w, dz):
    return I*a_b*(1-np.exp(-dz*(a_w+(a_b*X))))/(dz*(a_w+(a_b*X)))*3600

def light_ammonium_figures():
    matplotlib.rcParams.update({'font.size': 14})
    # http://www.aari.ru/docs/pub/060922/ree06.pdf
    # ammonium in the atlantic ocean
    nh4_conc = 0.109# uMol +-0.15 uMol
    N_growth_yield = 0.1166 #gDW/mmol NH3
    light_growth_yield = 0.00106 # gDW / mmol photons
    N_concentration = 0.1 #mM
    a_b = 0.2846 #m^2/gDW
    a_w = 0.465 
    dz = 0.01 #m
    grid_size = 1 #cm
    grid_volume_cm3 = grid_size**3 #Cubic centimeters
    grid_volume = grid_volume_cm3*1e-3 # Convert from cubic centimeters to L
    cell_weight = 66*1e-15 #g per cell (Cermak et al., 2017)
    initial_pop = 1.3e9*cell_weight*grid_volume # gDW/L

    light_amplitude = 0.04
    max_light = light(light_amplitude, initial_pop, a_b, a_w, dz)
    max_growth_light = max_light*light_growth_yield

    # Scatter points
    p1 = (0.39, 0.9) # This work
    p2 = (0.35*1e-3, 0.069) #  Maranon et al., 2013
    p3 = (0.35*1e-3, 0.66) #  Grossowiczc et al., 2017


    Vmax = np.arange(0, 1.2, 0.1) # mmol 
    Km = np.logspace(-4,0.1, 41)
    growth_rate_Nlim = np.zeros((len(Vmax), len(Km)))
    growth_rate_Nabund = np.zeros((len(Vmax), len(Km)))
    growth_rate_Natlantic = np.zeros((len(Vmax), len(Km)))
    growth_rate_N01 = np.zeros((len(Vmax), len(Km)))
    growth_rate_N001 = np.zeros((len(Vmax), len(Km)))
    growth_rate_N0001 = np.zeros((len(Vmax), len(Km)))
    max_growth = 0.0997
    for i, vi in enumerate(Vmax):
        for j, kj in enumerate(Km):
            growth_rate_Nlim[i,j] = MM(vi, kj, 0.1)*N_growth_yield
            growth_rate_Nabund[i,j] = MM(vi, kj, 0.8)*N_growth_yield
            growth_rate_Natlantic[i,j] = MM(vi, kj, 1e-4)*N_growth_yield
            growth_rate_N001[i,j] = MM(vi, kj, 1e-2)*N_growth_yield
            growth_rate_N0001[i,j] = MM(vi, kj, 1e-3)*N_growth_yield
                
    light_Nlim = growth_rate_Nlim/light_growth_yield
    light_Nabund = growth_rate_Nabund/light_growth_yield
    light_Natlantic = growth_rate_Natlantic/light_growth_yield
    light_N001 = growth_rate_N001/light_growth_yield
    light_N0001 = growth_rate_N0001/light_growth_yield



    # First plot
    xx, yy = np.meshgrid(Km, Vmax)
    fig, axes = plt.subplots(1,2, sharex = True, sharey=True,  figsize = (6, 3))
    CS1 = axes[0].contourf(xx,yy, growth_rate_Nabund, levels = [0, 0.042,1], alpha = 0.3, cmap = "Set1")
    CS2 = axes[1].contourf(xx,yy, growth_rate_Nlim, levels = [0, 0.042,1], alpha = 0.3, cmap = "Set1")
    for ax in axes:
        for i, p in enumerate([p1, p2, p3]):
            ax.scatter(p[0], p[1], c = "C{0}".format(i), s = 40)
        ax.set_xlabel("Km [mM]")
    # Format axes
    axes[0].set_xscale("log")
    axes[0].set_ylabel("Vmax [mmol gDW$^{-1}$ h$^{-1}$]")
    axes[0].set_title("Nitrogen-abundant")
    axes[1].set_title("Nitrogen-poor")
    plt.tight_layout()
    plt.savefig("km_vmax_sensitivity1.svg")
    plt.close(fig)

    # xx, yy = np.meshgrid(Km, Vmax)
    # fig, [ax1,ax2, ax3] = plt.subplots(1,3, sharex = True, figsize = (12, 8))
    # CS1 = ax1.contour(xx,yy, light_Nlim, levels = [0.1,1, 10, 40, 80], alpha = 1, colors = "k", linewidths = 2)
    # ax1.set_xscale("log")
    # plt.show()
    # Second plot
    fig, axes = plt.subplots(2,2, sharex = True, sharey = True, figsize = (6, 6))
    axes = axes.flatten()

    ax_titles = [0.1, 0.01, 0.001, 0.0001]
    data = [light_Nlim, light_N001, light_N0001, light_Natlantic]
    axes[0].set_xscale("log")        
    for i, ax in enumerate(axes):
        zz = data[i]
        CS = ax.contour(xx,yy, zz, levels = [0.1,1, 10, 40, 80], colors = "k", linewidths = 2)
        ax.clabel(CS, fontsize = 11,  fmt='%1.1f')
        for j, p in enumerate([p1, p2, p3]):
            ax.scatter(p[0], p[1], c = "C{0}".format(j), zorder = 3, s = 40)

        ax.set_title("[NH4$^+$] = {0:.1f}$\mu$M".format(ax_titles[i]*1e3))
        if i>1:
            ax.set_xlabel("Km [mM]")
    # Format axes
    axes[2].set_ylabel("Vmax [mmol gDW$^{-1}$ h$^{-1}$]")
    axes[0].set_ylabel("Vmax [mmol gDW$^{-1}$ h$^{-1}$]")
    plt.savefig("km_vmax_sensitivity_2.svg")


if __name__ == '__main__':
    if 0:
        create_death_figures_growth_glycogen()

    if 1:
        calculate_death_rate_deviations()

    if 0:
        create_glycogen_figures_growth_glycogen()
        calculate_glycogen_deviations()

    if 0:
        light_ammonium_figures()
