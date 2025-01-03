# Open KMAP-C Toolkit

## Overview of KMAP

The **Kinetic Modeling and Analysis Package (KMAP)** is an open-source software environment designed to implement and apply various tracer kinetic models for analyzing dynamic positron emission tomography (PET) data. It particularly focuses on addressing the challenges associated with total-body PET kinetic modeling, for example, considering time delay correction. The main goal of this open-source toolkit is to provide developers of tracer kinetic modeling with a foundational library to build upon, saving them from starting from scratch. The initial version of **KMAP** was developed at the University of California, Davis, and its open-source version was launched to support the [Open Kinetic Modeling Initiative](https://www.openkmi.org/).

## KMAP-C Toolkit

The **KMAP-C Toolkit** is a core C/C++ source code library to define and implement the input function, kinetic models, optimization algorithms, and utility functions. For parametric imaging (voxel-wise kinetic modeling), it adopts the OpenMP parallel programming model for acceleration. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/31573cb0-b1f5-4c50-8a51-8da9490eb214" width="1000" >
</div>

Currently, this C toolkit is not aimed at providing a comprehensive modeling solution but instead, a development backbone to support those time-consuming tracer kinetic models for high-level programming languages such as MATLAB and Python. Therefore the toolkit also includes wapper functions designed to integrate the kinetic modeling implementations for other programming languages, primarily MATLAB but may be extended for Python in the future. Please find the [**KMAP-M**](https://github.com/sharekm/KMAP-M) toolbox that offers more kinetic modeling approaches to be used with MATLAB.

## Contributors
See the full list of [Contributors](CONTRIBUTORS.md).

## Licensing

This repository is licensed under [MIT License](KMAP-C/LICENSE).
