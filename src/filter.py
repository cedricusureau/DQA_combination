import pandas as pd

def filter_dataframe(df, bool_list, beads_column, minimum_pos, limite_neg):
    masks = []

    for idx, bead in enumerate(beads_column):
        if bool_list[idx] == True:  # Si True, on crée un masque pour les billes positives
            masks.append(df[bead] > minimum_pos)
        elif bool_list[idx] == "NA":
            masks.append(pd.Series([True] * len(df)))  # Crée un masque qui ne filtre rien
        else:  # Si False, on crée un masque pour les billes négatives
            masks.append(df[bead] < limite_neg)

    # On combine tous les masques avec l'opérateur logique 'et'
    final_mask = pd.concat(masks, axis=1).all(axis=1)

    return df[final_mask]


def reorder_columns(df):
    new_order = []

    for col in df.columns:
        if "*" not in col:
            new_order.append(col)

    for col in df.columns:
        if "*" in col:
            new_order.append(col)

    df = df[new_order]
    return df


def remove_some_bead(df, beads_to_remove):
    for col in df.columns:
        if "*" in col:
            for bead in beads_to_remove:
                if bead in col:
                    df = df.drop(col, axis=1)
    return df