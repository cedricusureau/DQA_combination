import pandas as pd
import itertools
import src.filter as filter
import time
import src.profil_extractor as profil_extractor
import os

DQA_alleles = ["DQA1*03:01", "DQA1*02:01", "DQA1*03:02", "DQA1*03:03", "DQA1*04:01", "DQA1*05:01",
               "DQA1*05:03", "DQA1*05:05", "DQA1*06:01"]
data_file = "data/CII_patient_with_fip_date_prvt_MOINS_ML_SINGLE_ONLY_CLEAN.csv"

minimum_pos = 1000
limite_neg = 500
ignore = ["DR", "DP"]
ignore_beads = True

job_name = "DQA1_comb_without_DR_DP_bis"


def load_df(file_name):
    df1 = pd.read_csv(file_name, sep=";")
    df1 = df1[(df1['NC'] < 500) & (df1['PC'] > 3000)]
    # Drop row with nan value in date_prvt column
    df1 = df1.dropna(subset=["date_prvt"])
    # Turn date_prvt into datetime
    df1["date_prvt"] = pd.to_datetime(df1["date_prvt"], dayfirst=True)
    df1 = df1.sort_values(by=["fip", "date_prvt"])
    return df1


if __name__ == "__main__":
    # Make all the combinations of DQA1 alleles for every possible lenght of the combination
    # (from 1 to 9)

    # Make_results_folder
    if not os.path.exists(job_name):
        os.makedirs(job_name)

    df = load_df(data_file)

    # if ignore_beads:
    #     df = filter.remove_some_bead(df, ignore)

    DQA_alleles_combinations = []

    for i in range(1, len(DQA_alleles) + 1):
        tmp_combination = list(itertools.combinations(DQA_alleles, i))
        DQA_alleles_combinations += tmp_combination

    beads_column = [i for i in df.columns if "*" in i]

    # For every beads combination in DQA_alleles_combinations, make a true/false list depending if the bead is in or not
    all_bool_list = []
    for DQA in DQA_alleles_combinations:
        tmp_bool_list = []
        # if any of the DQA is in the bead, then it's true
        for bead in beads_column:
            if ignore:
                if "DQ" in bead:
                    tmp_bool_list.append(any(x in bead for x in DQA))
                else:
                    tmp_bool_list.append("NA")
            else:
                tmp_bool_list.append(any(x in bead for x in DQA))

        all_bool_list.append(tmp_bool_list)

    results = {}
    count = 0

    folders = ["Très fréquent", "Fréquent", "Peu fréquent", "Rare", "Très rare"]

    # Make a folder for each element in folders in the results folder
    for folder in folders:
        if not os.path.exists(f"{job_name}/{folder}"):
            os.makedirs(f"{job_name}/{folder}")

    for i, DQA_comb in zip(all_bool_list, DQA_alleles_combinations):
        # if count < 10:
            count += 1
            start_time = time.time()
            tmp_sub_df = filter.filter_dataframe(df, i, beads_column, minimum_pos, limite_neg)
            # if the dataframe is empty, don't save it
            if len(tmp_sub_df) > 0:
                end_time = time.time()  # Enregistrer le temps de fin
                  # Calculer le temps écoulé

                DQA_str = "_".join([i.replace("*", "").replace(":", "") for i in DQA_comb])
                DQA_str_id = "_".join([i for i in DQA_comb])
                results[DQA_str_id] = len(tmp_sub_df)

                # Save the dataframe as an excel file
                tmp_sub_df = filter.reorder_columns(tmp_sub_df)
                # Depending of the len of tmp_sub_df, save it in the right folder. For example, if len > 500, save it in the "Très fréquent" folder

                if len(tmp_sub_df) > 500:
                    output = f"{job_name}/Très fréquent/{DQA_str}.xlsx"

                elif len(tmp_sub_df) > 100:
                    output = f"{job_name}/Fréquent/{DQA_str}.xlsx"
                elif len(tmp_sub_df) > 50:
                    output = f"{job_name}/Peu fréquent/{DQA_str}.xlsx"

                elif len(tmp_sub_df) > 10:
                    output = f"{job_name}/Rare/{DQA_str}.xlsx"

                elif len(tmp_sub_df) >= 1:
                    output = f"{job_name}/Très rare/{DQA_str}.xlsx"

                profil_extractor.colorized_df_minimized(tmp_sub_df.T, output)

                elapsed_time = end_time - start_time
                print(f"{DQA_comb} took {elapsed_time:.4f} seconds")
                print(f"{count}/{len(all_bool_list)}")
            else:
                DQA_str_id = "_".join([i for i in DQA_comb])
                results[DQA_str_id] = len(tmp_sub_df)
        # Save results as csv file

    results_df = pd.DataFrame.from_dict(results, orient="index", columns=["count"])
    results_df.to_csv(f"{job_name}/results_{limite_neg}_{minimum_pos}.csv")
