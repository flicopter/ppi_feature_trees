ppi_feature_trees
=================
Protein-protein interaction modeling.

Libraries needed
-----------------

```bash
# [Prody] (http://prody.csb.pitt.edu)
sudo pip install -U ProDy

# Protein structure data
# [Zlab Benchmark] (http://zlab.umassmed.edu/benchmark/)
cd Data/
wget http://zlab.umassmed.edu/benchmark/benchmark4.tgz
tar xf benchmark4.tgz
```
Post-docking scripts
--------------------

```bash
python contactmap.py <docking_output_file> <partner_pdb_file>
```

*partner pdb file*

explanation comes here..

Note that a docking output file has a set of pdb files docked.
You need to have these files in the correct path.

```bash

```

 

