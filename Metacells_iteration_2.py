#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:31:53 2024

@author: aphuong
"""
print('Import Packages')

import anndata as ad             # For reading/writing AnnData files
import matplotlib.pyplot as plt  # For plotting
import scanpy as sc
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sb             # For plotting
from scipy import io
from scipy import sparse
import scipy.sparse as sp        # For sparse matrices
import shutil                    # for filesystem operations
from math import hypot           # For plotting
from typing import *             # For type annotations
from tqdm import tqdm

mc.pl.set_max_parallel_piles(173)

print('Read Files')

cells = sc.read_h5ad("/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-2/cells.h5ad")
metacells = sc.read_h5ad("/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-2/metacells.h5ad")

def convey_cell_annotations_to_metacells() -> None:
    
    # Assign a single value for each metacell based on the cells.
    mc.tl.convey_obs_to_group(
        adata=cells, gdata=metacells,
        property_name="cellType", to_property_name="cellType",
        method=mc.ut.most_frequent  # This is the default, for categorical data
    )
    
    mc.tl.convey_obs_to_group(
        adata=cells, gdata=metacells,
        property_name="BC_type", to_property_name="BC_type",
        method=mc.ut.most_frequent  # This is the default, for categorical data
    )
    
    mc.tl.convey_obs_fractions_to_group(adata=cells, gdata=metacells, property_name="patient_id")
    mc.tl.convey_obs_fractions_to_group(adata=cells, gdata=metacells, property_name="scDblFinder.class")


def compute_next_iteration(next_iteration_index: int) -> None:

    print("# DIVIDE AND CONQUER...")
    global metacells
    metacells = None # So can be gc-ed
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456,target_metacell_size = 40)

    print("# COLLECT METACELLS...")
    metacells = mc.pl.collect_metacells(
        cells, name=f"Bassez.iteration-{next_iteration_index}.metacells", random_seed=123456
    )
    print(f"Iteration {next_iteration_index}: {metacells.n_obs} metacells, {metacells.n_vars} genes")

    print("# CONVEY CELL ANNOTATIONS...")
    convey_cell_annotations_to_metacells()

def finalize_next_iteration(next_iteration_index: int, *, with_types: bool) -> None:
    print("# COMPUTE FOR MCVIEW...")
    mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=123456)


    print("# SAVE CELLS...")
    cells.write_h5ad(f"/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-{next_iteration_index}/Bassez-it2.cells.h5ad")

    print("# SAVE METACELLS...")
    metacells.write_h5ad(f"/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-{next_iteration_index}/Bassez-it2.metacells.h5ad")


def next_iteration_without_types(next_iteration_index: int) -> None:
    compute_next_iteration(next_iteration_index)
    finalize_next_iteration(next_iteration_index, with_types=False)

print('Compute Next Iteration:')

next_iteration_without_types(2)
