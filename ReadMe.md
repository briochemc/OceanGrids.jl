<a href="https://gist.github.com/briochemc/10e891bdb7da49fc4bf5467a5876434f">
  <img src="https://user-images.githubusercontent.com/4486578/59238897-0a004c80-8c43-11e9-861c-5fe00069af92.png", align="right", width="50%">
</a>

# OceanGrids

<p>
  <a href="https://github.com/briochemc/OceanGrids.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/briochemc/OceanGrids.jl/mac.yml?label=OSX&logo=Apple&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/OceanGrids.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/briochemc/OceanGrids.jl/linux.yml?label=Linux&logo=Linux&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/OceanGrids.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/briochemc/OceanGrids.jl/windows.yml?label=Windows&logo=Windows&logoColor=white&style=flat-square">
  </a>
  <a href="https://codecov.io/gh/briochemc/OceanGrids.jl">
    <img src="https://img.shields.io/codecov/c/github/briochemc/OceanGrids.jl/master?label=Codecov&logo=codecov&logoColor=white&style=flat-square">
  </a>
</p>

This package is a dependency of [AIBECS](https://github.com/briochemc/AIBECS.jl.git).
It defines types for grids used by AIBECS.

The goal of [OceanGrids](https://github.com/briochemc/OceanGrids.jl.git) is to standardize the format of grids in order for AIBECS to avoid confusion when swapping the circulation it uses for another.

For example, units for the grid data in the Ocean Circulation Inverse Model (OCIM, see [*DeVries et al*., 2014](https://doi.org/10.1002/2013GB004739)) products are not documented, so that it is easy to get confused and carry dimensional inconsistencies in one's model.
[OceanGrids](https://github.com/briochemc/OceanGrids.jl.git) attempts to fix these discrepancies by always using the same format and provide tests to ensure some level of consistency.
