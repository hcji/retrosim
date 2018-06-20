# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 15:01:47 2018

@author: Ji
"""

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)

from __future__ import print_function
from rdkit.Chem.Draw import IPythonConsole, ReactionToImage, MolToImage, MolsToGridImage
from IPython.display import SVG, display, clear_output
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit import DataStructs
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import sys
sys.path.append('../../')
from retrosim.utils.draw import ReactionStringToImage, TransformStringToImage
from retrosim.utils.generate_retro_templates import process_an_example
from retrosim.data.get_data import get_data_df, split_data_df
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants

data = get_data_df('retrosim/data/data_processed.csv')

getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=False)
templates = {}
rcts_ref_fps = {}
for i in data.index:
    templates[i] = '(' + process_an_example(data['rxn_smiles'][i], super_general=True).replace('>>', ')>>')
    rcts_ref_fps[i] = getfp(data['rxn_smiles'][i].split('>')[0])

uni_templates = list(set(templates.values()))
rxns = {}
for i, template in enumerate(uni_templates):
    rxns[i] = rdchiralReaction(template)


class_ = 1
similarity_metric = DataStructs.BulkTanimotoSimilarity # BulkDiceSimilarity or BulkTanimotoSimilarity
similarity_label = 'Tanimoto'
getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=False)
getfp_label = 'Morgan2noFeat'

try:
    if prev_FP != getfp:
        raise NameError
except NameError:
    all_fps = []
    for smi in tqdm(data['prod_smiles']):
        all_fps.append(getfp(smi))
    data['prod_fp'] = all_fps
    prev_FP = getfp

if class_ != 'all':
    datasub = data.loc[data['class'] == class_]
else:
    datasub = data
fps = list(datasub['prod_fp'])

def do_one(ix):
    rct = rdchiralReactants(datasub['prod_smiles'][ix])
    outcomes = {}
    for i, rxn in rxns.items():
        try:
            outcomes[i] = rdchiralRun(rxn, rct, combine_enantiomers=False)
        except Exception as e:
            outcomes[i] = []
    
    
    
    
    
    
    
    
    
    

