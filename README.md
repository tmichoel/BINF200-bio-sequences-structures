# BINF200 Biological Sequences and Structures

## Installation

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> BINF200

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that [raw data](#raw-data) are not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "BINF200"
```
which auto-activate the project and enable local path handling from DrWatson.

## Notebooks

The `.jl` files in the `notebooks` folder are [Pluto](https://plutojl.org/) notebooks. To run these notebooks, assuming you completed the above [installation instructions](#installation), open a Julia console in the project directory and do:
```
julia> using Pluto
julia> Pluto.run()
```
Find more details at the bottom of the [Pluto homepage](https://plutojl.org/)  or on the [MIT Introduction to Computational Thinking page](https://computationalthinking.mit.edu/Fall23/installation/).

## Raw data

To run the scripts and notebooks in this repository you first need to create a folder `data` in the repository's root folder. Then download the following datasets and store them in the right location.

### Long read sequencing data of expressed antigens in Trypanosoma brucei infections

1. Create a subfolder `PacBio_VSG` of the `data` folder.

2. Download the files in [this dataset](https://doi.org/10.18710/FFANM0) and store them in the `PacBio_VSG` folder, keeping the folder structure of the OneDrive folder intact.

**Note:** If you only want to run the [notebooks](#notebooks), you *don't* need to download the contents of the `original` folder of [the dataset](https://doi.org/10.18710/FFANM0).

### Coronavirus data

Download the files available [at this OneDrive URL](https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/Coronavirus?csf=1&web=1).


Store the files in a subfolder `Coronavirus` of the `data` folder.

## Releases

[![DOI](https://zenodo.org/badge/673794071.svg)](https://zenodo.org/doi/10.5281/zenodo.10043221)