# ELG_analysis

## Table of contents
* [Description](#description)
* [Dependencies](#dependencies)
* [Varying parameters](#varying-parameters)

## Description

This repository contains a number of routines written as ´Python´ scripts and ´jupyter´ notebooks. Overall, all routines are used to select objects in the *J-PLUS* database, and classify them as *Emission-line galaxies*. A number of tools are used to achieve this end, including routines to write and retreive data from ´SQL´ queries, the Machine learning library ´scikit-learn´ and the Deep-learning package ´Keras´.

### Dependencies

This code makes use of public data from the [J-PLUS Data Release 1](http://j-plus.es/datareleases/data_release_dr1), accessible after registration in the J-PLUS portal. 
Most of the code interacting with the database itself was developed by [Raul Angulo](mailto:reangulo@gmail.com).

In addition, make sure you have the following packages installed:
- numpy
- matplotlib
- scipy

The machine learning analysis requires the following packages:

- scikit-learn
- Keras
- TensorFlow




