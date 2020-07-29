# Black stain root disease landscape management model: BLASTMAN
## About
The **BLack STain root disease landscape MANagement model (BLASTMAN)** is an agent-based model that simulates forest management disturbances (i.e., thinnings and clear-cut regeneration harvests) and the spread of a tree disease (black stain root disease, BSRD) in forested landscapes.

Adam J. Bouché developed this model under the guidance of Dr. Klaus J. Puettmann and Dr. David C. Shaw as a thesis for a Master of Science (MS) degree in Forest Ecosystems and Society at Oregon State University between 2017 and 2020. The thesis text can be found at [ScholarsArchive@OSU](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/c247f0268?locale=en).

This ABM wa implemented in NetLogo (v.6.0.4), R (v.3.5.3), and Go (v.1.14). The most recent version of the 'spread infection' programis bundled in this repo, but the repository for the program implemented by Mario Vega is found [here](https://github.com/mariowhowrites/spread-infection).

## Using the model
This model runs in [NetLogo (v.6.0.4)](http://ccl.northwestern.edu/netlogo/). R must be installed on the computer (v.3.5.3 and later should work) and the R extension for NetLogo must be [configured](http://ccl.northwestern.edu/netlogo/docs/r.html). 

### Directories
In addition, several directories must be located in the same directory as the model ('.nlogo'):

  **/supporting_files** - A directory containing all of the '.nls' files used by the main model file, including:

  * basic_functions_1.10.0.0.nls
  * insect-attraction_1.10.0.0.nls
  * probability_parameters_1.10.0.0.nls
  * ring_lists_1.10.0.0.nls
  * spread-infection_Go_1.10.0.0.nls
  * spread-infection
  * spread-infection.exe

  **/supporting_files/programs** - A directory called 'programs' used to store programs and files used during model runs.
  
  **/output** - A directory where all model outputs (results) are created.
  
  **/output/rasters** - A directory specifically for raster model outputs ('maps' of the landscape)

### Necessary modifications to model code
While NetLogo can work on relative paths (relative to the directory), R needs the absolute paths and the working directory for the R instance running in NetLogo cannot be changed.
* Change in main model
* Change in spread-infection_Go_1.10.0.0.nls
