SmallMolEval
----------------
Using machine learning to score potential drug candidates may offer an advantage over traditional imprecise scoring functions because the parameters and model structure can be learned from the data. However, models may lack interpretability, are often overfit to the data, and are not generalizable to drug targets and chemotypes not in the training data. Benchmark datasets are prone to artificial enrichment and analogue bias due to the overrepresentation of certain scaffolds in experimentally determined active sets. Datasets can be evaluated using spatial statistics to quantify the dataset topology and better understand potential biases. Dataset clumping comprises a combination of self-similarity of actives and separation from decoys in chemical space and is associated with overoptimistic virtual screening results. This code explores methods of quantifying potential biases and examines some common benchmark datasets.

Documentation
----------------
File: 
remove_AVE_bias2.py
slight modification on atomwise script to split data

run_remove_AVE_bias.py
example of running remove_AVE_bias2.py on DUDE dataset

main.py,main_activeonly.py,main.old.py
scripts that run the MUV spatial statistics

DescriptorSets.py
mostly contains functions used by MUV statistics and called in main files

gf.plot
plots MUV statistics

makegraphs.py
uses gf.plots to make whole dataset plots

analyze_AVE_bias.py
no revisions from atomwise, computes the bias score and AUC of ligand based models

aveanalyze.py
runs analyze_AVE_bias.py for directory of multiple directories containing splits on different receptors

Authors
----------------

SmallMolEval was written by Dr. Sally Ellingson.

Release
----------------

SmallMolEval is released under an MIT license.  For more details see the
NOTICE and LICENSE files.

``LLNL-CODE-759342``