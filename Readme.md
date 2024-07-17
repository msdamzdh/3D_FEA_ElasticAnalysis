# 3D Solid Finite Element Analysis

This repository contains MATLAB code for performing 3D solid finite element analysis (FEA) and visualizing the results.

## Table of Contents
1. [Overview](#overview)
2. [Files](#files)
3. [Getting Started](#getting-started)
4. [Usage](#usage)
5. [Dependencies](#dependencies)
6. [Contributing](#contributing)
7. [Results](#results)
   
## Overview

This project implements a 3D solid finite element analysis solver with visualization capabilities. It supports multiple rectangular domains, 8-node brick elements, and can calculate displacements, stresses, and Von Mises stress.

## Files

1. `SolidFEA.m`: The core FEA solver
2. `PlotI3DterationHistory.m`: Function for plotting 3D iteration history
3. `LBeam.m`: Example script for analyzing an L-shaped beam

## Getting Started

1. Clone this repository: git clone https://github.com/msdamzdh/3D_FEA_ElasticAnalysis.git
2. Open MATLAB and navigate to the repository directory.
3. Run the `LBeam.m` script to see an example analysis.

## Usage

### SolidFEA.m

This file contains the main FEA solver and helper functions. It's not typically run directly but is called by other scripts.


### PlotI3DVonmisesCountour.m

Use this function to visualize the 3D Von Mises stress contours:

```mattlab
PlotI3DVonmisesCountour(nor, NodeCoord, IC, nelx, nely, nelz, scale, Phi, Vonmises, order)
```

## Results


https://github.com/user-attachments/assets/2640d3aa-2430-418c-bedb-39954f7a02f7


