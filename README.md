# BrookingsPC2020

Replication files for
["What's Up with the Phillips Curve?"](https://www.brookings.edu/bpea-articles/whats-up-with-the-phillips-curve/)
by Marco Del Negro, Michele Lenza, Giorgio E. Primiceri, and Andrea Tambalotti, *Brookings Papers on Economic Activity*,
Sprgin 2020. Conference Draft

## Required software

- Matlab 18
- Julia v1.0 or v1.1
- brookings_pc branch of [DSGE.jl](https://github.com/FRBNY-DSGE/DSGE.jl) v1.0.0
- [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl) v0.1.6

**Download instructions**

1. Download Julia from `https://julialang.org/downloads/`.
2. Open the Julia REPL and type:

   a. Press `]` to enter Pkg mode.

   b. Enter `add https://github.com/FRBNY-DSGE/DSGE.jl#brookings_pc`

   c. Enter `add SMC`, followed by `pin SMC@v0.1.6`.

   d. If, after running this replication code, you would like to use the most
      current version of DSGE.jl, go into Pkg mode, enter `rm DSGE`, and enter `add DSGE`.
      To use the most current version of smc, go into Pkg mode and enter `free SMC`.

## Installing this repository

Git users are welcome to fork this repository or clone it for local
use. Non-Git users will probably find it easiest to download the zip
file by clicking on the green `Clone or download` button on the right
hand side of this screen, and then clicking "Download ZIP".


## Directory structure


## How to run the VAR code
See the readme.docx in src/ReplicationFilesVAR.

## How to run the DSGE code
