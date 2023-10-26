import pandas as pd
import numpy as np
import xlsxwriter


def load_df(classe, MOINS_NC):
    if classe == 1:
        if MOINS_NC:
            df = pd.read_csv("data/patient_MOINS_NC_compilation_tout_CI_SAG_standard.csv", sep=";")
        else:
            df = pd.read_csv("data/patient_compilation_tout_CI_SAG_standard.csv", sep=";")
    elif classe == 2:
        if MOINS_NC:
            df = pd.read_csv("data/patient_MOINS_NC_compilation_tout_CII_SAG_standard.csv", sep=";")
        else:
            df = pd.read_csv("data/patient_compilation_tout_CII_SAG_standard.csv", sep=";")
    else:
        print("Redéfinir les paramètres")

    df.set_index("Unnamed: 0", inplace=True)
    return df


def load_template(classe):
    if classe == 1:
        df = pd.read_excel("input/template_CI.xlsx", sheet_name=None)
    elif classe == 2:
        df = pd.read_excel("input/template_CII.xlsx", sheet_name=None)

    return df


def extract_classe_I(df, tmp, seuil_pos, seuil_neg, seuil_int):
    pos_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "O" or str(j).upper() == "P"]
    neg_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "X" or str(j).upper() == "N"]
    int_bead = [i for i, j in zip(tmp.Allele, tmp.select) if str(j).upper() == "I"]

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
    return df_to_keep


def extract_classe_II(df, tmp, seuil_pos, seuil_neg, seuil_int):
    pos_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "O" or str(j).upper() == "P"]
    neg_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "X" or str(j).upper() == "N"]
    int_bead = [i + str(k).replace("nan", "") for i, j, k in zip(tmp.alpha, tmp.select, tmp.beta) if
                str(j).upper() == "I"]

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

        if (min_mfi_pos > seuil_pos) & (max_mfi_neg < seuil_neg) & (min_mfi_int > seuil_int[0]) & (max_mfi_int < seuil_int[1]):
            to_keep.append(True)
        else:
            to_keep.append(False)

    df_to_keep = df[to_keep]
    return df_to_keep


def write_FC_liste(df, output, FORMAT):
    file = open(output, "w")

    if FORMAT == "vertical":
        for i in df.loc["ech"]:
            file.write(str(i.split("_")[1]) + "\n")
    elif FORMAT == "horizontal":
        for i in df["ech"]:
            file.write(str(i.split("_")[1]) + "\n")

    file.close()


def colorized_df(df_final, output):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(output, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df_final.to_excel(writer, sheet_name='Sheet1')

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    # Add a format. Light red fill with dark red text.
    yellow = workbook.add_format({'bg_color': '#ffff00',
                                  'font_color': '#000000'})

    # Add a format. Green fill with dark green text.
    orange = workbook.add_format({'bg_color': '#ffaa00',
                                  'font_color': '#000000'})

    red = workbook.add_format({'bg_color': '#ff1e00',
                               'font_color': '#000000'})

    brique = workbook.add_format({'bg_color': '#960707',
                                  'font_color': '#ffffff'})
    # Calculate the range to which the conditional format is applied.
    (max_row, max_col) = df_final.shape
    min_row = 1  # Skip header.
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

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


def reorder_columns(df, template, classe):
    new_order = []

    for col in df.columns:
        if "*" not in col:
            new_order.append(col)

    if classe == 1:
        for i in template["Allele"]:
            new_order.append(i)

    if classe == 2:
        for i, j in zip(template["alpha"], template["beta"]):
            bead = (str(i) + str(j)).replace("nan","")
            new_order.append(bead)
            if bead == "DQA1*05:05DQB1*03:19":
                new_order.append("DQA1*05:05DQB1*03:01")

    df = df[new_order]

    return df
