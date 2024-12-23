"""
Load the data
"""
import uproot
import numpy as np
import pandas as pd
import torch
import sys
from parallel_pandas import ParallelPandas
import awkward as ak
import awkward_pandas #required for converting trees with vector data to pandas df's

from transformer_ee.dataloader.pd_dataset import Normalized_pandas_Dataset_with_cache


def get_sample_indices(sample_size: int, config) -> tuple:
    """
    A function to get the indices of the samples

    sample_size:    the number of samples
    config:         the configuration dictionary
    return:         the indices of train, validation and test sets
    """

    seed = config["seed"]

    _indices = np.arange(sample_size)
    np.random.seed(seed)
    np.random.shuffle(_indices)

    test_size = config["test_size"]
    valid_size = config["valid_size"]

    # test_size and valid_size can be either int or float
    if isinstance(test_size, float):
        test_size = int(sample_size * test_size)
    if isinstance(valid_size, float):
        valid_size = int(sample_size * valid_size)

    train_indicies = _indices[: sample_size - valid_size - test_size]
    valid_indicies = _indices[
        sample_size - valid_size - test_size : sample_size - test_size
    ]
    test_indicies = _indices[sample_size - test_size :]

    print("train indicies size:\t", len(train_indicies))
    print("valid indicies size:\t", len(valid_indicies))
    print("test  indicies size:\t", len(test_indicies))

    return train_indicies, valid_indicies, test_indicies


def get_train_valid_test_dataloader(config: dict):
    """
    A function to get the train, validation and test datasets
    Use the statistic of the training set to normalize the validation and test sets
    """
    
    with uproot.open(config["data_path"])[config["tree"]] as tree:
        df = tree.arrays(config["branches"], library = 'pd')
    
    ###apply desired cuts to the dataset
    df.dropna(subset=['W_truth'], inplace=True)  #drops events in W_truth' but null values seem to only appear in pandora track or shower related branches
    
    #select events with non-zero numbers of tracks or showers:
    QUERY = '(nshowers_pandoraShower != 0 or ntracks_pandoraTrack != 0)' #csv version uses 'and' instead of 'or'
    df = df.query(QUERY)
    
    #create branches with theta, phi from existing cosx, cosy, cosz 
    def _get_theta(cosyangle):
        theta = np.arccos(cosyangle)*180./np.pi
        return theta
    def _get_phi(cosxangle, coszangle):
        phi = np.arctan2(cosxangle, coszangle)*180/np.pi
        return phi
    df.loc[:,'nu_theta'] = df.apply(lambda x: _get_theta(x['nu_dcosy_truth']), axis =1)
    df.loc[:,'nu_phi'] =df.apply(lambda x: _get_phi(x['nu_dcosy_truth'], x['nu_dcosz_truth']), axis =1)

    genie_branches_to_filter = ['genie_primaries_pdg','genie_Eng', 'genie_Px','genie_Py','genie_Pz','genie_P','genie_mass','genie_status_code']
    
    def filter_genie(row, genie_columns):
        # Get the status codes
        status_codes = row['genie_status_code']
        
        # Filter each genie column based on the status codes
        filtered_genies = {col: [x for x, status in zip(row[col], status_codes) if status in [1]] for col in genie_columns}
        
        genie_Px_filtered = np.array(filtered_genies['genie_Px'], dtype=float)
        genie_Py_filtered = np.array(filtered_genies['genie_Py'], dtype=float)
        genie_Pz_filtered = np.array(filtered_genies['genie_Pz'], dtype=float)
        genie_Eng_filtered = np.array(filtered_genies['genie_Eng'], dtype=float)
        genie_P_filtered = np.sqrt(np.square(genie_Px_filtered) + np.square(genie_Py_filtered) + np.square(genie_Pz_filtered))
        filtered_genies['genie_P'] = genie_P_filtered

        return pd.Series(filtered_genies)

    filtered_genies = df.apply(filter_genie, axis=1, genie_columns=genie_branches_to_filter)
    df[genie_branches_to_filter] = filtered_genies[genie_branches_to_filter]

    #convert ak arrays from awkward_pandas into np arrays:
    df['enu_truth'] = df.apply(lambda x: ak.to_numpy(x['enu_truth']).astype(np.float64).item(), axis =1) #type must match 'nu_theta'
    
    train_idx, valid_idx, test_idx = get_sample_indices(len(df), config)
    print('valid_idx', valid_idx)
    print('train_idx', train_idx)

    train_set = Normalized_pandas_Dataset_with_cache(
        config, df.iloc[train_idx].reset_index(drop=True, inplace=False)
    )
    valid_set = Normalized_pandas_Dataset_with_cache(
        config,
        df.iloc[valid_idx].reset_index(drop=True, inplace=False),
        weighter=train_set.weighter,
    )
    test_set = Normalized_pandas_Dataset_with_cache(
        config,
        df.iloc[test_idx].reset_index(drop=True, inplace=False),
        weighter=train_set.weighter,
    )

    train_set.statistic()

    train_set.normalize()
    valid_set.normalize(train_set.stat)
    test_set.normalize(train_set.stat)

    batch_size_train = config["batch_size_train"]
    batch_size_valid = config["batch_size_valid"]
    batch_size_test = config["batch_size_test"]

    trainloader = torch.utils.data.DataLoader(
        train_set,
        batch_size=batch_size_train,
        shuffle=True,
        num_workers=10,
    )

    validloader = torch.utils.data.DataLoader(
        valid_set,
        batch_size=batch_size_valid,
        shuffle=False,
        num_workers=10,
    )

    testloader = torch.utils.data.DataLoader(
        test_set,
        batch_size=batch_size_test,
        shuffle=False,
        num_workers=10,
    )

    return trainloader, validloader, testloader, train_set.stat
