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

The machine learning analysis requires the following packages:

- [scikit-learn](https://scikit-learn.org/stable/index.html)
- [Keras](https://keras.io/)
- [TensorFlow](https://www.tensorflow.org/)

## Examples

`3D-HSTxJ-PLUS.ipynb`:

This notebook performs a cross-match between J-PLUS data and that from the 3D-HST project over the AEGIS field. The result can be used to study the optical properties (given by J-PLUS) of Near-Infrared, grism selected objects by 3D-HST.



