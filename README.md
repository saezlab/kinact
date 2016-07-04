# kinact
## Toolbox for Kinase Activity Scoring

This repository collects different computational methods to infer kinase activities from phosphoproteomics datasets. Each method will be accompanied by a protocol for explanation and in order to enable transfer to other datasets.

The package has been implemented by Jakob Wirbel in the [Saez-Rodriguez Group at RWTH-Aachen](http://www.combine.rwth-aachen.de/).

**Please be aware that the package is under on-going development.**

## Installation

In order to use `kinact`, clone the repository and type

```
python setup.py install
```
in your console.

## Quick start

The package contains example data from a publication by [de Graaf et al.](http://europepmc.org/abstract/MED/24850871), in which the phosphoproteome of Jurkat T-cells after prostaglandin E2 stimulation was analysed. The data can be loaded with:
```python
import kinact
data_fc, data_p_value = kinact.get_example_data()
print data_fc.head()

		"5min"		"10min"		"20min"		"30min"		"60min"
ID
A0AVK6_S71	−0.319306	−0.484960	−0.798082	−0.856103	−0.928753 
A0FGR8_S743	−0.856661	−0.981951	−1.500412	−1.441868	−0.861470 
A0FGR8_S758	−1.445386	−2.397915	−2.692994	−2.794762	−1.553398
A0FGR8_S691	0.2714580	0.2645960	0.5016850	0.4619840	0.655501 
A0JLT2_S226	−0.080786	1.0697100	0.5197800	0.5208830	−0.296040
```

Prior-knowledge information about kinase-substrate interactions can be loaded from the `pypath` package (see also the [documentation](http://omnipathdb.org/)).
```python
adjacency_matrix = kinact.get_kinase_targets(sources=['all'])
```

Finally, estimation of the kinase activities can be performed as follows for the example of the KSEA protocol:
```python
scores, p_values = kinact.ksea.ksea_mean(data_fc=data_fc["5min"].dropna(),
					 interactions=adjacency_matrix,
					 mP=data_fc["5min"].values.mean(),
					 delta=data_fc["5min"].values.std())
```

## Methods included in kinact:

+ KSEA

   KSEA is a method for the inference of kinase activities from phosphoproteomics data based on kinase substrate sets, which are constructed from curated information about kinase-substrate interactions from public resources like PhosphoSitePlus. The values of the fold changes of the phospho-sites in the substrate set of a given kinase are used to compute a score for the activity of this kinase. The main scoring system for KSEA is the mean or the median of the fold changes in the substrate set. KSEA was first proposed in the publication from [Casado et al.](http://europepmc.org/abstract/MED/23532336).

+ to be extended
