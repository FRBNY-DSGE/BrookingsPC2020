# BrookingsPC2020

Replication files for
["What's Up with the Phillips Curve?"](https://www.brookings.edu/bpea-articles/whats-up-with-the-phillips-curve/)
by Marco Del Negro, Michele Lenza, Giorgio E. Primiceri, and Andrea Tambalotti, *Brookings Papers on Economic Activity*,
Spring 2020. Conference Draft

## Required software

VAR results
- Matlab 18

DSGE results
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


## How to run the [VAR code](https://github.com/FRBNY-DSGE/BrookingsPC2020/tree/master/src/ReplicationFilesVAR)
See the [readme.docx](https://github.com/FRBNY-DSGE/BrookingsPC2020/blob/master/src/ReplicationFilesVAR/readme.docx) in src/ReplicationFilesVAR.

## How to run the [DSGE code](https://github.com/FRBNY-DSGE/BrookingsPC2020/tree/master/src/ReplicationFilesDSGE)
See the [readme.txt](https://github.com/FRBNY-DSGE/BrookingsPC2020/blob/master/src/ReplicationFilesDSGE/readme.txt) in src/ReplicationFilesDSGE.

Disclaimer
------
The VAR code is copyright of Giorgio Primiceri. The DSGE code is copyright of the authors and the Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.
