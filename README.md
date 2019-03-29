# ELG_analysis

## Table of contents
* [Description](#description)
* [Dependencies](#dependencies)
* [Examples](#examples)

## Description

This repository contains a number of routines written as `Python` scripts and `jupyter` notebooks. Overall, all routines are used to select objects in the *J-PLUS* database, classify them as *Emission-line galaxies* and perform some analysis over them. A number of tools are used to achieve this end, including routines to write and retreive data from `SQL` queries, the Machine learning library `scikit-learn` and the Deep-learning package `Keras`.

## Dependencies

This code makes use of public data from the [J-PLUS Data Release 1](http://j-plus.es/datareleases/data_release_dr1), accessible after registration in the J-PLUS portal. 
Most of the code interacting with the database itself was developed by [Raul Angulo](mailto:reangulo@gmail.com).

In addition, make sure you have the following packages installed:
- numpy
- matplotlib
- scipy
- astropy

The machine learning analysis requires the following packages:

- [scikit-learn](https://scikit-learn.org/stable/index.html)
- [Keras](https://keras.io/)
- [TensorFlow](https://www.tensorflow.org/)

## Examples

[Keras_class.ipynb](https://github.com/aaorsi/ELG_analysis/blob/master/Keras_class.ipynb). This notebook loads a catalogue of ELGs to train a deep-learning sequential neural network with *Keras*. To scan for the best set of parameters, the code performs a Cross-validation Grid search using `scikit-learn`.

[3D-HSTxJ-PLUS.ipynb](https://github.com/aaorsi/ELG_analysis/blob/master/3D-HSTxJ-PLUS.ipynb). This notebook performs a cross-match between J-PLUS data and that from the 3D-HST project over the AEGIS field. The result can be used to study the optical properties (given by J-PLUS) of Near-Infrared, grism selected objects by 3D-HST.

[ELG_stacked.ipynb](https://github.com/aaorsi/ELG_analysis/blob/master/ELG_stacked.ipynb). This notebook stacks J-PLUS objects with a given selection over the regions corresponding to *RedMAPPER* clusters within a given redshift range. The stack of objects around clusters in a narrow redshift range around *z~0.3* should include the contribution of Emission-line galaxies at that redshift.





