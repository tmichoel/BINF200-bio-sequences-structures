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

The `.jl` files in the `notebooks` folder are [Pluto](https://plutojl.org/) notebooks. Follow the instructions at the bottom of the [Pluto](https://plutojl.org/) homepage or on the [MIT Introduction to Computational Thinking page](https://computationalthinking.mit.edu/Fall23/installation/) to install [Pluto](https://plutojl.org/) and run these notebooks.

## Raw data

To run the scripts and notebooks in this repository you need to download was data and store it in the right location.

### Long read sequencing data of expressed antigens in Trypanosoma brucei infections

Download the files available [at this OneDrive URL](https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/PacBio_VSG?csf=1&web=1).

Store the files in a subfolder `PacBio_VSG` of the `data` folder, keeping the folder structure of the OneDrive folder intact.

### Coronavirus data

Download the files available [at this OneDrive URL](https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/Coronavirus?csf=1&web=1).


Store the files in a subfolder `Coronavirus` of the `data` folder.