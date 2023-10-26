import os
import pandas as pd
import numpy as np
import xlsxwriter
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from openpyxl import load_workbook

""" ATTENTION, les billes DQA1*05:05	DQB1*03:19 et DQA1*05:05 DQB1*03:01 sont considéré comme strictement identique.
Demander à identifier les sérums avec DQA1*05:05	DQB1*03:19 positif revient à chercher les sérums DQA1*05:05	
DQB1*03:01 positifs également """


def run(CLASSE, RUN_NAME, SEUIL_POS, SEUIL_NEG, SEUIL_INTERMEDIAIRE, FORMAT, OPEN, CLUSTERMAP, TYPAGE, KMEAN,
        K_VALUE):
    all_param_str = ["CLASSE", "RUN_NAME", "SEUIL_POS", "SEUIL_NEG", "SEUIL_INTERMEDIAIRE", "FORMAT",
                     "OPEN", "CLUSTERMAP", "TYPAGE", "KMEAN",
                     "K_VALUE"]
    all_param = [CLASSE, RUN_NAME, SEUIL_POS, SEUIL_NEG, SEUIL_INTERMEDIAIRE, FORMAT, OPEN, CLUSTERMAP,
                 TYPAGE, KMEAN,
                 K_VALUE]

    if TYPAGE:
        df_typages = pd.read_csv(
            r"O:\Laboratoire_Immunologie_Histocompatibilitee\Biologistes\Typage Glims\Lymhoe_Glims.csv", sep=";",
            low_memory=True)

        if CLASSE == 1:
            file_equivalence = [i.replace("\n", "") for i in open(
                r"C:\DonneesLocales\FastHLA\data\serologie_c1.txt")]
            df_typages = df_typages[["code_fip", "A1", "A2", "B1", "B2", "C1", "C2"]]
        elif CLASSE == 2:
            file_equivalence = [i.replace("\n", "") for i in open(
                r"C:\DonneesLocales\FastHLA\data\serologie_c2.txt")]
            df_typages = df_typages[
                ["code_fip", "DR1", "DR2", "DQ1", "DQ2", "DP1", "DP2", "DRw1ded", "DRw2ded", "DRw1ded2", "DRw2ded2",
                 "DQA1", "DQA2", "DPA1", "DPA2"]]

        equivalence = {}
        for element in file_equivalence:
            equivalence[element.split(" ")[1]] = element.split(" ")[0]

    dir_path = r"C:\DonneesLocales\FastHLA\Profil Extractor\result\{}".format(
        RUN_NAME)

    if RUN_NAME not in os.listdir(r"C:\DonneesLocales\FastHLA\Profil Extractor\result"):
        os.mkdir(dir_path)

    os.startfile(dir_path)
    df = load_df(CLASSE)
    template = load_template(CLASSE)
    res = {}

    for sheetname, tmp_template in template.items():
        KMEAN_NAME = "KMEAN_{}_cl_{}_P_{}_N_{}".format(sheetname, str(CLASSE),
                                                       str(SEUIL_POS), str(SEUIL_NEG))
        try:
            os.mkdir(dir_path + "/" + KMEAN_NAME)
        except:
            print("Répertoire Kmean déjà créé")

        if CLASSE == 1:
            tmp_res, beads_of_interests = extract_classe_I(df, tmp_template, SEUIL_POS, SEUIL_NEG, SEUIL_INTERMEDIAIRE)
        elif CLASSE == 2:
            tmp_res, beads_of_interests = extract_classe_II(df, tmp_template, SEUIL_POS, SEUIL_NEG, SEUIL_INTERMEDIAIRE)

        tmp_res = reorder_columns(tmp_res, tmp_template, CLASSE)

        if CLUSTERMAP:
            cluster_map(tmp_res,
                        r"C:\DonneesLocales\FastHLA\Profil Extractor\result\{}\clustermap_{}_cl_{}_P_{}_N_{}.jpeg".format(
                            RUN_NAME, sheetname, str(CLASSE),
                            str(SEUIL_POS), str(SEUIL_NEG)),
                        CLASSE, SEUIL_POS, SEUIL_NEG)

        if KMEAN:
            kmean(tmp_res,
                  "C:/DonneesLocales/FastHLA/Profil Extractor/result/{}/{}".format(
                      RUN_NAME, KMEAN_NAME),
                  CLASSE, beads_of_interests, K_VALUE)

        if TYPAGE:
            tmp_res = add_typage_col(tmp_res, equivalence, df_typages)

        if FORMAT == "vertical":
            tmp_res = tmp_res.transpose()
            tmp_res.index.names = ['FIP']
            res[sheetname] = tmp_res

        elif FORMAT == "horizontal":
            res[sheetname] = tmp_res
            pass
        else:
            print("""La variable format peut prendre ces deux valeurs : "vertical" ou "horizontal".""")

    for sheet, result in res.items():
        output_excel = "C:/DonneesLocales/FastHLA/Profil Extractor/result/{}/{}_cl_{}_P_{}_N_{}.xlsx".format(
            RUN_NAME, sheet, str(CLASSE), str(SEUIL_POS),
            str(SEUIL_NEG))
        FC_liste = "C:/DonneesLocales/FastHLA/Profil Extractor/result//{}/FC_{}_cl_{}_P_{}_N_{}.txt".format(
            RUN_NAME, sheet, str(CLASSE), str(SEUIL_POS),
            str(SEUIL_NEG))
        colorized_df(result, output_excel, template, all_param_str, all_param)

        if "critère" not in sheet:
            write_FC_liste(result, FC_liste, FORMAT)

    if OPEN:
        subprocess.Popen(["C:\Program Files (x86)\Microsoft Office\Office16\EXCEL.EXE", output_excel]).pid

    for sheet, result in res.items():
        if FORMAT == "vertical":
            print(sheet + ": " + "Sur un total de {} sérums, ".format(str(df.shape[0])) +
                  "le pattern demandé est présent dans très exactement {} % des cas.".format(
                      str((result.shape[1] / df.shape[0]) * 100)))
        else:
            print(sheet + ": " + "Sur un total de {} sérums, ".format(str(df.shape[0])) +
                  "le pattern demandé est présent dans très exactement {} % des cas.".format(
                      str((result.shape[0] / df.shape[0]) * 100)))


def add_typage_col(tmp_res, equivalence, df_typage):
    typage = {}
    tmp_res.index = tmp_res.index.astype("float")
    df_typage = df_typage.set_index("code_fip")

    for patient in tmp_res.fip:
        try:
            tmp_typage = df_typage.loc[int(patient)]
            tmp_allele = []
            for allele, serologie in equivalence.items():
                if serologie in list(tmp_typage):
                    tmp_allele.append(allele)
            typage[patient] = tmp_allele
        except Exception as e:
            print("Typage pour " + str(int(patient)) + " non trouvé")

    for ligne, rows in tmp_res.iterrows():
        try:
            new_row = []
            for colonne in rows.keys():
                if colonne in typage[rows["fip"]]:
                    new_row.append(equivalence[colonne])
                else:
                    new_row.append(" ")

            tmp_res.loc[ligne - 0.1] = new_row
        except Exception as e:
            pass

    tmp_res = tmp_res.sort_index()
    return tmp_res


def kmean(df, output, CLASSE, beads_of_interests, K_VALUE):
    if CLASSE == 1:
        row_to_keep = ["A*", "B*", "C*"]
    else:
        row_to_keep = ["DR", "DQ", "DP"]

    locus_to_keep = set()
    for bead in beads_of_interests:
        for locus in row_to_keep:
            if locus in bead:
                locus_to_keep.add(locus)

    tmp_df = df[[i for i in df.columns if "*" in i]]
    tmp_df = tmp_df.transpose()

    if CLASSE == 2:
        dq0319_serie = tmp_df.loc["DQA1*05:05DQB1*03:19"].combine_first(tmp_df.loc["DQA1*05:05DQB1*03:01"])
        tmp_df.drop(["DQA1*05:05DQB1*03:19", "DQA1*05:05DQB1*03:01"], inplace=True)
        tmp_df.loc["DQA1*05:05DQB1*03:19"] = dq0319_serie

    tmp_df = tmp_df.dropna(axis=1)
    for locus in locus_to_keep:
        tmp_tmp_df = tmp_df.loc[[i for i in tmp_df.index if locus in i]]
        tmp_tmp_df = normalize(tmp_tmp_df)
        tmp_tmp_df = make_KM_and_append_cluster_to_normalized_df(10, tmp_tmp_df,
                                                                 (output + (
                                                                     "/{}_silhouette.jpeg".format(locus))).replace("*",
                                                                                                                   ""),
                                                                 K_VALUE)
        make_plot_by_clusters(tmp_tmp_df, output + ("/{}.jpeg".format(locus)), K_VALUE, locus)


def make_plot_by_clusters(df, output, K_VALUE, locus):
    clusters = set(df["Clusters"])
    for cluster in clusters:
        tmp_df = df[df["Clusters"] == cluster]

        tmp_df.drop("Clusters", axis=1, inplace=True)
        sns.boxplot(data=tmp_df, showfliers=False).set(
            title="Locus: " + locus + ". Cluster N°" + str(cluster + 1) + ": " + str(
                len(tmp_df.index)) + " échantillons")
        plt.xticks(rotation=90)
        plt.tight_layout()
        figure_name = output.replace(".jpeg", "cluster_ " + str(cluster + 1) + ".jpeg")
        if K_VALUE != "":
            figure_name = figure_name.replace(".jpeg", "K_" + str(K_VALUE) + ".jpeg")
        plt.savefig(figure_name.replace("*", ""), dpi=400)
        plt.close()


def make_KM_and_append_cluster_to_normalized_df(clusters, df_normalize, output, K_VALUE):
    if K_VALUE == "":
        best_clusters = kmean_best_cluster(clusters, df_normalize, output)
    else:
        best_clusters = int(K_VALUE)
    df_normalize = perform_kmean(best_clusters, df_normalize)
    return df_normalize


def kmean_best_cluster(clusters, df_normalize, output):
    res = np.arange(clusters, dtype="double")
    df_normalize.dropna(inplace=True)
    for i in df_normalize.columns:
        for j in range(len(df_normalize.index)):
            if str(df_normalize[i].iloc[j]) == "nan":
                print(i, j)

    for k in np.arange(clusters):
        km = KMeans(n_clusters=k + 2)
        km.fit(df_normalize)
        res[k] = metrics.silhouette_score(df_normalize, km.labels_)

    plot_kmean_silhouette(res, clusters, output)

    return list(res).index(np.max(res)) + 2


def perform_kmean(cluster, df_normalize):
    km = KMeans(n_clusters=cluster).fit(df_normalize)
    labels = km.labels_

    df_normalize["Clusters"] = labels
    return df_normalize


def plot_kmean_silhouette(res, clusters, output):
    plt.title("Silhouette")
    plt.xlabel("# of clusters")
    plt.plot(np.arange(2, clusters + 2, 1), res)
    plt.savefig(output, dpi=400)
    plt.close()


def cluster_map(df, output, CLASSE, SEUIL_MAX, SEUIL_MIN):
    if CLASSE == 1:
        row_to_keep = ["A*", "B*", "C*"]
    else:
        row_to_keep = ["DR", "DQ", "DP"]

    tmp_df = df[[i for i in df.columns if "*" in i]]
    tmp_df = tmp_df.transpose()

    if CLASSE == 2:
        dq0319_serie = tmp_df.loc["DQA1*05:05DQB1*03:19"].combine_first(tmp_df.loc["DQA1*05:05DQB1*03:01"])
        tmp_df.drop(["DQA1*05:05DQB1*03:19", "DQA1*05:05DQB1*03:01"], inplace=True)
        tmp_df.loc["DQA1*05:05DQB1*03:19"] = dq0319_serie

    tmp_df.columns = [str(i) for i in range(len(tmp_df.columns))]

    tokeep = []
    for i, rows in tmp_df.iterrows():
        if max(rows) < 1000:
            tokeep.append(False)
        else:
            tokeep.append(True)

    tmp_df = tmp_df[tokeep]
    tmp_df = tmp_df.dropna(axis=1)
    sns.clustermap(tmp_df, figsize=(25, 25), cmap="rocket_r", col_cluster=False, xticklabels=False,
                   vmax=int(SEUIL_MAX * 2.5), vmin=int(SEUIL_MIN / 2.5))
    plt.savefig(output, dps=400)
    plt.close()
    for locus in row_to_keep:
        tmp_tmp_df = tmp_df.loc[[i for i in tmp_df.index if locus in i]]
        sns.clustermap(tmp_tmp_df, figsize=(25, 15), cmap="rocket_r", col_cluster=False, xticklabels=False,
                       vmax=int(SEUIL_MAX * 1.5), vmin=int(SEUIL_MIN / 2.5))
        plt.savefig(output.replace("clustermap", "clustermap_" + locus.replace("*", "")))
        plt.close()


def normalize(df):
    result = df.copy()
    result = result / result.max().astype(int)
    result = result.astype(float).transpose()

    return result


def load_df(classe):
    if classe == 1:
        df = pd.read_csv(r"C:\DonneesLocales\FastHLA\data\CI_patient_with_fip_date_prvt_MOINS_ML_SINGLE_ONLY_CLEAN.csv",
                         sep=";")

    elif classe == 2:
        df = pd.read_csv(
            r"C:\DonneesLocales\FastHLA\data\CII_patient_with_fip_date_prvt_MOINS_ML_SINGLE_ONLY_CLEAN.csv", sep=";")
    else:
        print("Redéfinir les paramètres")

    # df.set_index("index", inplace=True)
    return df


def load_template(classe):
    if classe == 1:
        df = pd.read_excel(
            r"C:\DonneesLocales\FastHLA\template\template_CI.xlsx",
            sheet_name=None)
    elif classe == 2:
        df = pd.read_excel(
            r"C:\DonneesLocales\FastHLA\template\template_CII.xlsx",
            sheet_name=None)

    return df


def extract_classe_I(df, tmp, seuil_pos, seuil_neg, seuil_int):
    pos_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "O" or str(j).upper() == "P"]
    neg_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "X" or str(j).upper() == "N"]
    int_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "I"]

    bead_of_interest = pos_bead + neg_bead + int_bead

    df_pos = df[pos_bead]
    df_neg = df[neg_bead]
    df_int = df[int_bead]
    to_keep = []
    np_pos = df_pos.to_numpy()
    np_neg = df_neg.to_numpy()
    np_int = df_int.to_numpy()

    if np_pos.size == 0:
        min_mfi_pos = 1000000
    if np_neg.size == 0:
        max_mfi_neg = -10000000
    if np_int.size == 0:
        min_mfi_int = 10000000
        max_mfi_int = -10000000

    for i in range(len(df)):

        if np_pos.size != 0:
            min_mfi_pos = np.nanmin(np_pos[i])
        if np_neg.size != 0:
            max_mfi_neg = np.nanmax(np_neg[i])
        if np_int.size != 0:
            min_mfi_int = np.nanmin(np_int[i])
            max_mfi_int = np.nanmax(np_int[i])

        if (min_mfi_pos > seuil_pos) & (max_mfi_neg < seuil_neg) & (min_mfi_int > seuil_int[0]) & (
                max_mfi_int < seuil_int[1]):
            to_keep.append(True)
        else:
            to_keep.append(False)

    df_to_keep = df[to_keep]
    return df_to_keep, bead_of_interest


def extract_classe_II(df, tmp, seuil_pos, seuil_neg, seuil_int):
    pos_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "O" or str(j).upper() == "P"]
    neg_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "X" or str(j).upper() == "N"]
    int_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "I"]

    bead_of_interest = pos_bead + neg_bead + int_bead

    if "DQA1*05:05DQB1*03:01" in pos_bead:
        pos_bead.append("DQA1*05:05DQB1*03:19")
    if "DQA1*05:05DQB1*03:19" in pos_bead:
        pos_bead.append("DQA1*05:05DQB1*03:01")
    if "DQA1*05:05DQB1*03:01" in neg_bead:
        neg_bead.append("DQA1*05:05DQB1*03:19")
    if "DQA1*05:05DQB1*03:19" in neg_bead:
        neg_bead.append("DQA1*05:05DQB1*03:01")

    df_pos = df[pos_bead]
    df_neg = df[neg_bead]
    df_int = df[int_bead]

    to_keep = []
    np_pos = df_pos.to_numpy()
    np_neg = df_neg.to_numpy()
    np_int = df_int.to_numpy()

    if np_pos.size == 0:
        min_mfi_pos = 1000000
    if np_neg.size == 0:
        max_mfi_neg = -10000000
    if np_int.size == 0:
        min_mfi_int = 10000000
        max_mfi_int = -10000000

    for i in range(len(df)):

        if np_pos.size != 0:
            min_mfi_pos = np.nanmin(np_pos[i])
        if np_neg.size != 0:
            max_mfi_neg = np.nanmax(np_neg[i])
        if np_int.size != 0:
            min_mfi_int = np.nanmin(np_int[i])
            max_mfi_int = np.nanmax(np_int[i])

        if (min_mfi_pos > seuil_pos) & (max_mfi_neg < seuil_neg) & (min_mfi_int > seuil_int[0]) & (
                max_mfi_int < seuil_int[1]):
            to_keep.append(True)
        else:
            to_keep.append(False)

    df_to_keep = df[to_keep]
    return df_to_keep, bead_of_interest


def write_FC_liste(df, output, FORMAT):
    file = open(output, "w")

    if FORMAT == "vertical":
        for i in df.loc["ech"]:
            if len(i) > 2:
                file.write(str(i.split("_")[1]) + "/n")
    elif FORMAT == "horizontal":
        for i in df["ech"]:
            if len(i) > 2:
                file.write(str(i.split("_")[1]) + "/n")

    file.close()


def colorized_df_minimized(df_final, output):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(output, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df_final.to_excel(writer, sheet_name='Sheet1')

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    yellow = workbook.add_format({'bg_color': '#ffff00',
                                  'font_color': '#000000'})
    yellow.set_align("center")

    orange = workbook.add_format({'bg_color': '#ffaa00',
                                  'font_color': '#000000'})

    red = workbook.add_format({'bg_color': '#ff1e00',
                               'font_color': '#000000'})

    brique = workbook.add_format({'bg_color': '#960707',
                                  'font_color': '#ffffff'})
    # Calculate the range to which the conditional format is applied.
    (max_row, max_col) = df_final.shape
    min_row = 13  # Skip header.
    min_col = 1  # Skip index.
    max_row = min_row + max_row - 1
    max_col = min_col + max_col - 1

    # Apply a conditional format to the cell range.
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 500,
                                  'maximum': 2000,
                                  'format': yellow})

    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 2000,
                                  'maximum': 3000,
                                  'format': orange})
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 3000,
                                  'maximum': 5000,
                                  'format': red})
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 5000,
                                  'maximum': 30000000,
                                  'format': brique})

    worksheet.set_column(0, 0, 25)

    my_format = workbook.add_format()
    my_format.set_align('center')

    for i in range(11, 1000):
        worksheet.set_row(i, None, my_format)

    # Close the Pandas Excel writer and output the Excel file.
    try:
        writer.save()
    except:
        writer.close()


def colorized_df(df_final, output, template, all_param_str, all_param):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(output, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df_final.to_excel(writer, sheet_name='Sheet1')

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    yellow = workbook.add_format({'bg_color': '#ffff00',
                                  'font_color': '#000000'})
    yellow.set_align("center")

    orange = workbook.add_format({'bg_color': '#ffaa00',
                                  'font_color': '#000000'})

    red = workbook.add_format({'bg_color': '#ff1e00',
                               'font_color': '#000000'})

    brique = workbook.add_format({'bg_color': '#960707',
                                  'font_color': '#ffffff'})
    # Calculate the range to which the conditional format is applied.
    (max_row, max_col) = df_final.shape
    min_row = 13  # Skip header.
    min_col = 1  # Skip index.
    max_row = min_row + max_row - 1
    max_col = min_col + max_col - 1

    # Apply a conditional format to the cell range.
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 500,
                                  'maximum': 2000,
                                  'format': yellow})

    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 2000,
                                  'maximum': 3000,
                                  'format': orange})
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 3000,
                                  'maximum': 5000,
                                  'format': red})
    worksheet.conditional_format(min_row, min_col, max_row, max_col,
                                 {'type': 'cell',
                                  'criteria': 'between',
                                  'minimum': 5000,
                                  'maximum': 30000000,
                                  'format': brique})

    worksheet.set_column(0, 0, 25)

    my_format = workbook.add_format()
    my_format.set_align('center')

    for i in range(11, 1000):
        worksheet.set_row(i, None, my_format)

    for i, df in template.items():
        new_col1 = [i for i in all_param_str] + [None for i in range(len(df) - len(all_param_str))]
        new_col2 = [i for i in all_param] + [None for i in range(len(df) - len(all_param_str))]
        df = df.set_index(df.columns[0])
        df["Paramètres"] = new_col1
        df["valeur"] = new_col2
        df.to_excel(writer, sheet_name="critère_" + i)

    # Close the Pandas Excel writer and output the Excel file.
    try:
        writer.save()
    except:
        writer.close()


def reorder_columns(df, template, classe):
    new_index = []

    for i in df["fip"]:
        new_index.append(i)

    df.index = new_index

    new_order = []

    for col in df.columns:
        if "*" not in col:
            new_order.append(col)

    if classe == 1:
        for i in template["Allele"]:
            new_order.append(i)

    if classe == 2:
        for i, j in zip(template["alpha"], template["beta"]):
            bead = (str(i) + str(j)).replace("nan", "")
            new_order.append(bead)
            if bead == "DQA1*05:05DQB1*03:19":
                new_order.append("DQA1*05:05DQB1*03:01")

    df = df[new_order]
    return df
