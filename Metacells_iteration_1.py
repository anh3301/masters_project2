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

sc.set_figure_params(scanpy=True, dpi=200, dpi_save=600, frameon=True, vector_friendly=True, fontsize=14, figsize=None, color_map=None, format='pdf', facecolor=None, transparent=False, ipython_format='png2x')

cells = sc.read_h5ad("/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-1/cells.h5ad")

print('Parallel Piles')
max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)

print('Divide and Conquer')
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells,random_seed=123456,target_metacell_size = 40)

print('Collect Metacells')
metacells = mc.pl.collect_metacells(cells, name="Bassez.iteration-1.metacells",random_seed=123456)
print(f"Iteration 1: {metacells.n_obs} metacells, {metacells.n_vars} genes")

print('Compute for MCView')
with mc.ut.progress_bar():
    mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=123456)

print('Write files')

cells.write_h5ad("/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-1/Bassez-it1.cells.h5ad")
metacells.write_h5ad("/data/rds/DBC/UBCN/CANCDYN/aphuong/Bassez/output/iterative/iteration-1/Bassez-it1.metacells.h5ad")


 
