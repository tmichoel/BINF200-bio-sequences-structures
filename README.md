# BINF200 Biological Sequences and Structures

## Notebooks

The `.jl` files in the `notebooks/pluto` folder are [Pluto](https://plutojl.org/) notebooks. To run these notebooks, first [install Pluto](https://plutojl.org/#install). Then open a Julia console in the project directory and do:
```
julia> using Pluto
julia> Pluto.run()
```

See the [Pluto documentation](https://plutojl.org/en/docs/) for more information.

## Additional software

### BLAST

You need to [download and install BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). Installers are available from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

On Linux, you can install BLAST+ with the command.

```
sudo apt install ncbi-blast+
```

### EMBOSS

You need to [download and install EMBOSS](http://emboss.open-bio.org/), the European Molecular Biology Open Software Suite. Instructions are available from [the EMBOSS documentation](http://emboss.open-bio.org/html/adm/ch01s01.html).

On Linux, you can install EMBOSS with the command.

```
sudo apt install emboss
```

## Raw data

To run the scripts and notebooks in this repository you first need to create a folder `data` in the repository's root folder. Then download the following datasets and store them in the right location.

### Long read sequencing data of expressed antigens in Trypanosoma brucei infections

1. Create a subfolder `PacBio_VSG` of the `data` folder.

2. Download the files in [this dataset](https://doi.org/10.18710/FFANM0) and store them in the `PacBio_VSG` folder, keeping the folder structure intact.

**Note:** If you only want to run the [notebooks](#notebooks), you *don't* need to download the contents of the `original` folder of [the dataset](https://doi.org/10.18710/FFANM0).

### Coronavirus data

Download the files available [at this OneDrive URL](https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/Coronavirus?csf=1&web=1).


Store the files in a subfolder `Coronavirus` of the `data` folder.

## Releases

[![DOI](https://zenodo.org/badge/673794071.svg)](https://zenodo.org/doi/10.5281/zenodo.10043221)