#!/usr/bin/env python3

# Taken from https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/cli/arboreto_with_multiprocessing.py

import sys
import time
import loompy as lp
import pandas as pd
from multiprocessing import Pool, cpu_count
import argparse
import tqdm

from arboreto.utils import load_tf_names
from arboreto.algo import genie3, grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network

from pyscenic.cli.utils import load_exp_matrix, suffixes_to_separator
from pathlib import PurePath


################################################################################
################################################################################
def run_infer_partial_network(target_gene_index):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=None
    )
    return n

def coexpression_modules(method, expression_mtx_fname, tf_names, num_workers,
    transpose='yes', sparse=True, cell_id_attribute='CellID', gene_attribute='Gene', seed=None):
        
    global gene_names
    global ex_matrix
    global method_params
    global tf_matrix
    global tf_matrix_gene_names
    
    if method == 'grnboost2':
        method_params = ['GBM', SGBM_KWARGS]  # regressor_type  # regressor_kwargs
    elif method == 'genie3':
        method_params = ['RF', RF_KWARGS]  # regressor_type  # regressor_kwargs

    start_time = time.time()
    ex_matrix = load_exp_matrix(
        expression_mtx_fname, (transpose == 'yes'), sparse, cell_id_attribute, gene_attribute
    )

    #print(f'sparse sparse: {sparse}', file=sys.stdout)
    if sparse:
        gene_names = ex_matrix[1]
        ex_matrix = ex_matrix[0]
        #print(f'gene_names sparse: {gene_names}', file=sys.stdout)
    else:
        gene_names = ex_matrix.columns
        #print(f'gene_names: {gene_names}', file=sys.stdout)

    end_time = time.time()
    print(
        f'Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time} seconds...',
        file=sys.stdout,
    )

    #tf_names = load_tf_names(tfs_fname)
    #print(f'Loaded {len(tf_names)} TFs...', file=sys.stdout)

    ex_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_names)

    print(f'starting {method} using {num_workers} processes...', file=sys.stdout)
    start_time = time.time()

    with Pool(num_workers) as p:
        adjs = list(
            tqdm.tqdm(
                p.imap(run_infer_partial_network, target_gene_indices(gene_names, target_genes='all'), chunksize=1),
                total=len(gene_names),
            )
        )

    adj = pd.concat(adjs).sort_values(by='importance', ascending=False)

    end_time = time.time()
    print(f'Done in {end_time - start_time} seconds.', file=sys.stdout)

    
    return adj
    
    
