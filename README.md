# kinact
## Toolbox for Kinase Activity Scoring

This repository collects different computational methods to infer kinase activities from phosphoproteomics datasets. Each method will be accompanied by a protocol for explanation and in order to enable transfer to other datasets.

Please be aware that the package is under on-going development.

In order to use `kinact`, clone the repository and type

```
python setup.py install
```
in your console.

Implemented by Jakob Wirbel in the Saez-Rodriguez Group at RWTH-Aachen.

## Methods included in kinact:

+ KSEA

   KSEA is a method for the inference of kinase activities from phosphoproteomics data based on kinase substrate sets, which are constructed from curated information about kinase-substrate interactions from public resources like PhosphoSitePlus. The values of the fold changes of the phospho-sites in the substrate set of a given kinase are used to compute a score for the activity of this kinase. The main scoring system for KSEA is the mean or the median of the fold changes in the substrate set. KSEA was first proposed in the publication from Casado et al.: Sci Signal. 2013 Mar 26; 6(268):rs6. doi: 10.1126/scisignal.2003573.

+ to be extended
