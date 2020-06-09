import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AffinityPropagation, SpectralClustering
import hdbscan
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn import decomposition
from matplotlib import pyplot as plt
from sklearn import manifold
import numpy as np

names_fn = "../Data/reaction_id_mapping.csv"
folder = "../Data/random_sampling_environment/iSO595v6/"
timestamp = "20200402-0956"

FBA_fn = folder + "FBA_random_sampling_10000_{0}.csv".format(timestamp)
pFBA_fn = folder + "pFBA_random_sampling_10000_{0}.csv".format(timestamp)
env_fn = folder + "env_random_sampling_10000_{0}.csv".format(timestamp)
FVA_min_fn = folder + "minFVA_random_sampling_10000_{0}.csv".format(timestamp)
FVA_max_fn = folder + "maxFVA_random_sampling_10000_{0}.csv".format(timestamp)

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)

CLASS_NAMES = ['Biomass','Amino acid','Nucleobases and nucleosides','Organic acid','Other','Uptake']

def figure_correlation_matrix_pFBA():
    df_names = pd.read_csv(names_fn, index_col = 0)
    df_env = pd.read_csv(env_fn, index_col = 0).T
    df = pd.read_csv(pFBA_fn, index_col = 0).T
    df[["LightEX", "HCO3EXcar", "AmmoniaEX", "FAKEOrthophosphateEX"]] = df_env[["Light", "HCO3", "Nitrogen", 
                                                                                "Phosphate"]].abs()
    drop_columns = ["tRNAEX", 'CadaverineEX', 'CadmiumEX','CalciumEX', 'ChlorideEX', 'CalciumEX', 'ChlorideEX', 
                    'IronEX','MagnesiumEX', 'PotassiumEX', 'SelenateEX', 'SodiumEX', 'StrontiumEX', 'ZincEX', 'R00024',
                    'H2OEX', "HEX"]

    df["SulfateEX"] *= -1

    df.drop(columns = drop_columns, inplace = True)
    df.drop(columns = df.columns[(df.abs()< 1e-3).all()], inplace = True)
    #df.columns = pd.merge(df.T, df_names["Name"], left_index = True, right_index = True, how = "left")["Name"]

    corr_matrix = df.corr()
    linked = linkage(corr_matrix, 'single', optimal_ordering = True)
    d = dendrogram(linked)
    df_corr = pd.DataFrame(corr_matrix)
    df_corr = df_corr.iloc[d['leaves'], d['leaves']]

    # Row colors
    row_class = pd.merge(df_corr, df_names, left_index = True, right_index = True, how = "left")
    row_class.set_index("Name", inplace = True)
    row_class = row_class["Class"]
    lut = dict(zip(CLASS_NAMES, sns.hls_palette(len(CLASS_NAMES))))
    row_colors = row_class.map(lut)


    df_corr.columns = pd.merge(df_corr.T, df_names, left_index = True, right_index = True, how = "left")["Name"]
    df_corr.index = df_corr.columns  

    sns.clustermap(df_corr, mask = np.triu(df_corr,1), yticklabels = True, xticklabels=True, 
                   row_cluster = False, col_cluster = False, cmap = "coolwarm", row_colors = row_colors, 
                   col_colors = row_colors, figsize = (16,16), vmin = -1, vmax = 1)#, annot = True, fmt='.1g')
    plt.savefig("pFBA_correlation.svg")
    sns.palplot(sns.hls_palette(len(df_names["Class"].unique())))
    plt.savefig("pFBA_colorbar.svg")

def figure_correlation_matrix_FVA_max():
    df_names = pd.read_csv(names_fn, index_col = 0)
    df_env = pd.read_csv(env_fn, index_col = 0).T
    df = pd.read_csv(FVA_max_fn, index_col = 0).T
    df[["LightEX", "HCO3EXcar", "AmmoniaEX", "FAKEOrthophosphateEX"]] = df_env[["Light", "HCO3", "Nitrogen", 
                                                                                "Phosphate"]].abs()
    drop_columns = ["tRNAEX", 'CadaverineEX', 'CadmiumEX','CalciumEX', 'ChlorideEX', 'CalciumEX', 'ChlorideEX', 
                    'IronEX','MagnesiumEX', 'PotassiumEX', 'SelenateEX', 'SodiumEX', 'StrontiumEX', 'ZincEX', 'R00024',
                    'H2OEX', "HEX"]
    df["SulfateEX"] *= -1

    df.drop(columns = drop_columns, inplace = True)
    df.drop(columns = df.columns[(df.abs()< 1e-3).all()], inplace = True)
    #df.columns = pd.merge(df.T, df_names["Name"], left_index = True, right_index = True, how = "left")["Name"]

    corr_matrix = df.corr()
    linked = linkage(corr_matrix, 'single', optimal_ordering = True)
    d = dendrogram(linked)
    df_corr = pd.DataFrame(corr_matrix)
    df_corr = df_corr.iloc[d['leaves'], d['leaves']]

    row_class = pd.merge(df_corr, df_names, left_index = True, right_index = True, how = "left")
    row_class.set_index("Name", inplace = True)
    row_class = row_class["Class"]
    lut = dict(zip(CLASS_NAMES, sns.hls_palette(len(CLASS_NAMES))))
    row_colors = row_class.map(lut)


    df_corr.columns = pd.merge(df_corr.T, df_names, left_index = True, right_index = True, how = "left")["Name"]
    df_corr.index = df_corr.columns                          

    sns.clustermap(df_corr, mask = np.triu(df_corr,1), yticklabels = True, xticklabels=True, 
                   row_cluster = False, col_cluster = False, cmap = "coolwarm", row_colors = row_colors, 
                   col_colors = row_colors, figsize = (16,16), vmin = -1, vmax = 1)#, annot = True, fmt='.1g')
    plt.savefig("max_FVA_correlation.svg")
    sns.palplot(sns.hls_palette(len(df_names["Class"].unique())))
    plt.savefig("max_FVA_colorbar.svg")
    fig, ax = plt.subplots(1, figsize = (16,16))
    sns.heatmap(df_corr, mask = np.triu(df_corr,1), yticklabels = True, xticklabels=True, 
                   cmap = "coolwarm", vmin = -1, vmax = 1, ax = ax, square = True)
    plt.subplots_adjust(bottom = 0.3)
    plt.savefig("max_FVA_correlation_heatmap.svg")

def figure_PCA_clustering(normalization = "log", n_clusters = 4):
    df_env = pd.read_csv(env_fn, index_col = 0).T
    df = pd.read_csv(pFBA_fn, index_col = 0).T
    df_names = pd.read_csv("../Data/reaction_id_mapping.csv", index_col = 0)
    df[["LightEX", "HCO3EXcar", "AmmoniaEX", "FAKEOrthophosphateEX"]] = df_env[["Light", "HCO3", "Nitrogen", 
                                                                                "Phosphate"]].abs()
    drop_columns = ["tRNAEX", 'CadaverineEX', 'CadmiumEX','CalciumEX', 'ChlorideEX', 'CalciumEX', 'ChlorideEX', 
                    'IronEX','MagnesiumEX', 'PotassiumEX', 'SelenateEX', 'SodiumEX', 'StrontiumEX', 'ZincEX', 'R00024',
                    'H2OEX', "HEX"]
    df.drop(columns = drop_columns, inplace = True)
    df.drop(columns = df.columns[(df.abs()< 1e-3).all()], inplace = True)
    
    df["SulfateEX"] *= -1
    if normalization == "log":
        df_pca = (df+1).transform(np.log)
        label = "log"
    elif normalization == "Z-score":
        df_pca = df.subtract(df.mean(axis = 0), axis = 1).divide(df.std(axis = 0), axis = 1)
        df_pca = df.divide(df.abs().max(axis = 0), axis = 1)
        label = "z_score"
    else:
        df_pca = df.divide(df.abs().max(axis = 0), axis = 1)
        label = "scaled_to_1"

    df_pca.columns = pd.merge(df_pca.T, df_names["Name"], left_index = True, right_index = True, how = "left")["Name"]

    k_fit = KMeans(n_clusters = n_clusters).fit(df_pca)
    df_cluster = pd.DataFrame(k_fit.cluster_centers_.T, index = df.columns)
    df_cluster.columns = np.arange(n_clusters, dtype = int)+1
    df_cluster.index = df_cluster.merge(df_names["Name"], left_on = "Reactions", right_index = True, how = "left")["Name"]


    temp = df_cluster.multiply(np.array(df.std(axis = 0)), axis = 0)
    temp = temp.add(np.array(df.mean(axis = 0)), axis = 0)
    df_c_scaled_to_1 = temp.divide(np.array(df.abs().max(axis = 0)), axis = 0)

    snsplot = sns.clustermap(df_c_scaled_to_1, z_score = False, row_cluster = True, col_cluster = False, yticklabels=True, figsize = (2, 16), vmin = 0, 
                             vmax = 1, cmap = 'Blues')#"coolwarm"
    snsplot.savefig('cluster_kmeans_{0}.svg'.format(label))
    snsplot.savefig('cluster_kmeans_{0}.pdf'.format(label))

    pca = decomposition.PCA(2)
    pca.fit(df_pca)
    fig, ax = plt.subplots(1)
    Y = pca.transform(df_pca)
    Yk = pca.transform(k_fit.cluster_centers_)
    colors_k = k_fit.predict(df_pca)
    plt.scatter(Y[:,0], Y[:,1], s = 0.1, c = colors_k, cmap=plt.get_cmap("Set2", 4))
    for i in range(4):
        plt.scatter(Yk[i,0], Yk[i,1], s = 20, c = "k")
        plt.text(Yk[i,0]+0.3, Yk[i,1]+0.3, "Cluster {0}".format(i+1), c = "k")
    plt.xlabel("Scores PC1")
    plt.ylabel("Scores PC2")
    plt.savefig('scores_{0}.svg'.format(label))
    plt.savefig('scores_{0}.pdf'.format(label))

    fig, ax = plt.subplots(1, figsize = (12, 4))
    reaction_ids = df_pca.columns
    row_class = pd.merge(df_pca.T, df_names, left_index = True, right_on = "Name", how = "left")
    row_class.set_index("Name", inplace = True)
    row_class = row_class["Class"]
    lut = dict(zip(CLASS_NAMES, sns.hls_palette(len(CLASS_NAMES))))
    row_colors = row_class.map(lut)

    for i in range(len(reaction_ids)):
        ax.scatter(pca.components_.T[i,0], pca.components_.T[i,1], 
            label = "{0}. {1}".format(i+1, reaction_ids[i]), color = row_colors[i])
        ax.text(pca.components_.T[i,0]+0.01, pca.components_.T[i,1]+0.01, i+1)
    plt.legend(ncol = 2, loc='upper left', bbox_to_anchor=(1.01, 1.02))
    plt.xlabel("Loadings PC1")
    plt.ylabel("Loadings PC2")
    #plt.tight_layout()
    plt.subplots_adjust(right = 0.5)
    plt.savefig('loadings_{0}.svg'.format(label))
    plt.savefig('loadings_{0}.svg'.format(label))

def tsne_clustering(normalization = "log", method = "KMeans", perplexity = 50, iterations = 3000, cluster = "center"):
    df_env = pd.read_csv(env_fn, index_col = 0).T
    df = pd.read_csv(pFBA_fn, index_col = 0).T
    df_names = pd.read_csv("../Data/reaction_id_mapping.csv", index_col = 0)
    df[["LightEX", "HCO3EXcar", "AmmoniaEX", "FAKEOrthophosphateEX"]] = df_env[["Light", "HCO3", "Nitrogen", 
                                                                                "Phosphate"]]
    drop_columns = ["tRNAEX", 'CadaverineEX', 'CadmiumEX','CalciumEX', 'ChlorideEX', 'CalciumEX', 'ChlorideEX', 
                    'IronEX','MagnesiumEX', 'PotassiumEX', 'SelenateEX', 'SodiumEX', 'StrontiumEX', 'ZincEX', 'R00024',
                    'H2OEX', "HEX"]

    #df["SulfateEX"] *= -1
    df.drop(columns = drop_columns, inplace = True)
    df.drop(columns = df.columns[(df.abs()< 1e-3).all()], inplace = True)
    #df.columns = pd.merge(df.T, df_names["Name"], left_index = True, right_index = True, how = "left")["Name"]

    #print(df[["LightEX", "HCO3EXcar", "AmmoniaEX", "FAKEOrthophosphateEX"]])

    if normalization == "log":
        df_pca = (df+1).transform(np.log10)
        label = "log"
        vmin = 0
        vmax = 2
        cmap = "Blues"
    elif normalization == "Z-score":
        df_pca = df.subtract(df.mean(axis = 0), axis = 1).divide(df.std(axis = 0), axis = 1)
        label = "z_score"
        vmin = -2
        vmax = 2
        cmap = "RdBu"
    else:
        df_pca = df.divide(df.abs().max(axis = 0), axis = 1)
        label = "scaled_to_1"
        vmin = -1
        vmax = 1
        cmap = "coolwarm"#"Blues"

    tsne = manifold.TSNE(n_components=2, init='random', random_state=0, perplexity = perplexity, n_iter = iterations)
    Y = tsne.fit_transform(df_pca)

    # Clustering of tsne data
    if method.lower() == "kmeans":
        n_clusters = 4
        k_fit_tsne = KMeans(n_clusters = n_clusters).fit(Y)
        classifier = k_fit_tsne.fit_predict(Y)
        Yk = k_fit_tsne.cluster_centers_
    elif method.lower() == "hdbscan":
        clusterer = hdbscan.HDBSCAN(min_cluster_size=200)
        classifier = clusterer.fit_predict(Y)
        n_clusters = len(clusterer.cluster_persistence_)
        Yk = np.array([clusterer.weighted_cluster_centroid(i) for i in range(n_clusters)])
        clusterer.condensed_tree_.plot(select_clusters = True)
        plt.savefig('tSNE_condensed_tree_{0}.svg'.format(label))
        plt.close()
        #plt.show()

    elif method.lower() == "spectral":
        n_clusters = 6
        k_fit_tsne = SpectralClustering(n_clusters = n_clusters).fit(Y)
        classifier = k_fit_tsne.fit_predict(Y)
        Yk = k_fit_tsne.cluster_centers_
    print(n_clusters)
    print(classifier)
    colors_y = [x if x != -1 else 7 for x in classifier]
    fig, ax = plt.subplots(1)
    l = ax.scatter(Y[:, 0], Y[:, 1], c=colors_y, cmap=plt.get_cmap("Set2", n_clusters +2), s= 10, alpha  =0.5)
    for i in range(n_clusters):
        ax.scatter(Yk[i,0], Yk[i,1], c = "k", s = 50)
        ax.text(Yk[i,0]+5, Yk[i,1]+5, i+1)
    #ax.scatter(Y[green, 0], Y[green, 1], c="g")
    #plt.colorbar(l)
    plt.savefig('tSNE_{0}_{1}_{2}_{3}.svg'.format(label, method, perplexity, iterations))

    if cluster == "center":
        cluster_idx = np.zeros(n_clusters)
        for i in range(n_clusters):
            dist = ((Y-Yk[i,:])**2).sum(axis = 1)
            cluster_idx[i] = dist.argmin()
            print("Cluster ", i, argmin)

        df_tsne = df_pca.iloc[cluster_idx, :]
        df_tsne.index = np.arange(1,n_clusters+1)


        print(df.iloc[cluster_idx, :])
        print(df_tsne)

        row_class = pd.merge(df_tsne.T, df_names, left_index = True, right_index = True, how = "inner")
        row_class.set_index("Name", inplace = True)
        row_class = row_class["Class"]
        lut = dict(zip(CLASS_NAMES, sns.hls_palette(len(CLASS_NAMES))))
        row_colors = row_class.map(lut)

        col_colors = sns.color_palette("Set2", n_clusters+2)[:-2]
        # col_colors = pd.Series(dict(zip(np.arange(1, n_clusters), sns.color_palette("Set2", n_clusters))))
        
        df_tsne.columns = pd.merge(df_tsne.T, df_names, left_index = True, right_index = True, how = "left")["Name"]

        sns.clustermap(df_tsne.T, row_cluster = True, col_cluster = False, yticklabels=True, figsize = (4, 16), 
                        vmin = vmin, vmax = vmax, cmap = cmap, z_score = None, 
                        square = True, row_colors = row_colors, col_colors = col_colors)
        plt.subplots_adjust(right = 0.5)
        plt.savefig('cluster_tSNE_{0}_{1}_{2}_{3}.svg'.format(label, method, perplexity, iterations))
    else:
        # Use mean as col in heatmaps
        cluster_means = []
        print("Classes: ", np.unique(classifier))
        for i in range(n_clusters):
            index = classifier == i
            print(i, sum(index))
            mean_i = df_pca.loc[index, :].mean()
            std_i = df_pca.loc[index, :].std()
            print(mean_i)
            print(df.loc[index, :].mean())
            #print(std_i)
            cluster_means.append(mean_i)

        df_tsne = pd.DataFrame(cluster_means)
        df_tsne.index = np.arange(1,n_clusters+1)

        print(df_tsne)

        row_class = pd.merge(df_tsne.T, df_names, left_index = True, right_index = True, how = "inner")
        row_class.set_index("Name", inplace = True)
        row_class = row_class["Class"]
        lut = dict(zip(CLASS_NAMES, sns.hls_palette(len(CLASS_NAMES))))
        row_colors = row_class.map(lut)

        col_colors = sns.color_palette("Set2", n_clusters+2)[:-2]
        # col_colors = pd.Series(dict(zip(np.arange(1, n_clusters), sns.color_palette("Set2", n_clusters))))
        
        df_tsne.columns = pd.merge(df_tsne.T, df_names, left_index = True, right_index = True, how = "left")["Name"]

        sns.clustermap(df_tsne.T, row_cluster = True, col_cluster = False, yticklabels=True, figsize = (4, 16), 
                        vmin = vmin, vmax = vmax, cmap = cmap, z_score = None, 
                        square = True, row_colors = row_colors, col_colors = col_colors)
        plt.subplots_adjust(right = 0.5)
        plt.savefig('cluster_means_tSNE_{0}_{1}_{2}_{3}.svg'.format(label, method, perplexity, iterations))


def high_bicarbonate_PPP():
    import cobra, cameo
    model_fn = "../../iSO595/Model_files/iSO595v6.xml"
    model = cobra.io.read_sbml_model(model_fn)

    light = -50
    model.reactions.AmmoniaEX.lower_bound = -0.3
    model.reactions.HCO3EXcar.lower_bound = -10
    model.reactions.FAKEOrthophosphateEX.lower_bound = -0.05
    model.reactions.LightEX.lower_bound = light
    model.reactions.CO2EX.bounds = (0, 0.83)
    model.reactions.R00024.bounds = (0, 4.7)
    #model.reactions.R06282.bounds = (-1000,0)
        

    with model:
        #model.reactions.AmmoniaEX.lower_bound = -0.1
        ppp = cameo.flux_analysis.phenotypic_phase_plane(model, ["HCO3EXcar", "GlycogenEX"], points = 40, source = "HCO3EXcar")
        wide_base = ppp.data_frame.pivot(index = "HCO3EXcar", columns = "GlycogenEX", values = "objective_upper_bound")# Possible with c_yield_upper_bound

    with model:
        r_atp = model.add_boundary(model.metabolites.get_by_id("ATP[c]"), "sink")
        r_adp = model.add_boundary(model.metabolites.get_by_id("ADP[c]"), "sink")
        r_atp.bounds = (-0.5, 0)
        r_adp.bounds = (0, 0.5)
        #model.reactions.AmmoniaEX.lower_bound = -0.1
        ppp = cameo.flux_analysis.phenotypic_phase_plane(model, ["HCO3EXcar", "GlycogenEX"], points = 40, source = "HCO3EXcar")
        wide_add_atp = ppp.data_frame.pivot(index = "HCO3EXcar", columns = "GlycogenEX", values = "objective_upper_bound")# Possible with c_yield_upper_bound

    with model:
        model.reactions.LightEX.lower_bound = light*2
        #model.reactions.AmmoniaEX.lower_bound = -0.1
        ppp = cameo.flux_analysis.phenotypic_phase_plane(model, ["HCO3EXcar", "GlycogenEX"], points = 40, source = "HCO3EXcar")
        wide_low_light = ppp.data_frame.pivot(index = "HCO3EXcar", columns = "GlycogenEX", values = "objective_upper_bound")# Possible with c_yield_upper_bound

    with model:
        model.reactions.LightEX.lower_bound = light*0.5
        #model.reactions.AmmoniaEX.lower_bound = -0.1
        ppp = cameo.flux_analysis.phenotypic_phase_plane(model, ["HCO3EXcar", "GlycogenEX"], points = 40, source = "HCO3EXcar")
        wide_lowlow_light = ppp.data_frame.pivot(index = "HCO3EXcar", columns = "GlycogenEX", values = "objective_upper_bound")# Possible with c_yield_upper_bound

    fig, axes = plt.subplots(2,2, figsize = (16, 16), sharey = True, sharex = True)
    axes = axes.flatten()
    title = ["Base model", "Additional ATP", r"200% light", r"50% light"]
    for i, wide in enumerate([wide_base, wide_add_atp, wide_low_light, wide_lowlow_light]):
        ax = axes[i]
        x = np.array(wide.index)
        y = np.array(wide.columns)
        m = ax.contourf(y, x, wide, levels = 20, cmap = "Greens", vmin = 0, vmax = 0.04)
        ax.set_xlabel("GlycogenEX  [mmol/gDW/h]")
        ax.set_ylabel("HCO3EX  [mmol/gDW/h]")
        ax.set_title(title[i])
        if i == 1:
            m0 = m

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(m0, cax=cbar_ax)
    cbar.set_label("Biomass")
    plt.savefig("high_bicarbonate_PPP.svg")
    #plt.show()

if __name__ == '__main__':
    if 0:
        figure_correlation_matrix_FVA_max()
    if 0:
        figure_correlation_matrix_pFBA()
    if 0:
        figure_PCA_clustering(normalization = "Z-score")
    if 1:
        tsne_clustering(normalization = "scaled_to_1", method = "hdbscan", cluster = "center")
    if 0:
        plt.rcParams.update({'font.size': 20})
        high_bicarbonate_PPP()