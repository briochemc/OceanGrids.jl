# OceanGrids

<p>
  <img src="https://img.shields.io/badge/stability-experimental-orange.svg">
</p>
<p>
  <a href="https://travis-ci.com/briochemc/OceanGrids.jl">
    <img alt="Build Status" src="https://travis-ci.com/briochemc/OceanGrids.jl.svg?branch=master">
  </a>
  <a href='https://coveralls.io/github/briochemc/OceanGrids.jl'>
    <img src='https://coveralls.io/repos/github/briochemc/OceanGrids.jl/badge.svg' alt='Coverage Status' />
  </a>
</p>
<p>
  <a href="https://ci.appveyor.com/project/briochemc/OceanGrids-jl">
    <img alt="Build Status" src="https://ci.appveyor.com/api/projects/status/briochemc/OceanGrids-jl.svg?branch=master">
  </a>
  <a href="https://codecov.io/gh/briochemc/OceanGrids.jl">
    <img src="https://codecov.io/gh/briochemc/OceanGrids.jl/branch/master/graph/badge.svg" />
  </a>
</p>

This package is a dependency of [AIBECS](https://github.com/briochemc/AIBECS.jl.git).
It defines two types for objects used by AIBECS:
- the grids
- the transport matrices
The goal of [OceanGrids](https://github.com/briochemc/AIBECS.jl.git) is to standardize the format of these objects in order for AIBECS to use without confusion, regardless of the ocean grid or circulation that it uses.

For example, the transport matrix from products of the Ocean Circulation Inverse Model (OCIM, see [*DeVries et al*., 2014](https://doi.org/10.1002/2013GB004739)) are traditionally built as flux convergences with unit yr<sup>-1</sup>.
This is different from the format adopted in AIBECS, which expects transport matrices as flux divergences and in the SI unit of s<sup>-1</sup>.
In a similar vein, units for the grid data in OCIM products are not documented, so that it is easy to get confused and carry dimensional inconsistencies in one's model.
[OceanGrids](https://github.com/briochemc/AIBECS.jl.git) attempts to fix these discrepancies by always using the same format and provide tests to ensure some level of consistency.

 