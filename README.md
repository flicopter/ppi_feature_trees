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
tar xzf benchmark4.tgz
```

It may be convenient to place all version benchmark data in one place, like below: 
```bash
cd benchmark4
mkdir unbound
mkdir bound
cp structures/*/*_[rl]_b.pdb bound/.
cp structures/*/*_[rl]_u.pdb unbound/.
```

Post-docking scripts
--------------------

```bash
python contactmap.py <docking_output_file> <partner_pdb_file>
```

**partner_pdb_file**

explanation comes here..

**Note**
A docking output file has a set of pdb files docked.
You need to have these files in the correct path.

```bash

```

 

