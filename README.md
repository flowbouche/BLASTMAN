# Black stain root disease landscape management model: BLASTMAN
The BLack STain root disease landscape MANagement model (BLASTMAN) is an agent-based model implemented in NetLogo (v.6.0.4), R, and Go that simulates forest management disturbances (i.e., thinnings and clear-cut regeneration harvests) and the spread of a tree disease in forested landscapes.

Adam J. Bouché developed this model under the guidance of Dr. Klaus J. Puettmann and Dr. David C. Shaw as a thesis for a Master of Science (MS) degree in Forest Ecosystems and Society at Oregon State University between 2017 and 2020. The thesis text can be found at [ScholarsArchive@OSU](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/c247f0268?locale=en).

This model also relies on the ['spread infection' program](https://github.com/mariowhowrites/spread-infection) written in Go designed by Adam Bouché and developed by Mario Vega .

## USING THE MODEL ON YOUR COMPUTER
Several directories must be located in the same directory as the model ('.nlogo'):

  #### /supporting_files
  A directory containing all of the '.nls' files used by the main model file, including:

  * basic_functions_1.10.0.0.nls
  * insect-attraction_1.10.0.0.nls
  * probability_parameters_1.10.0.0.nls
  * ring_lists_1.10.0.0.nls
  * spread-infection_Go_1.10.0.0.nls

  #### /supporting_files/programs/program_og
  A directory called 'programs' containing a subdirectory 'program_og' that contains the binary executables of the Go program 'spread-infection' designed by Adam Bouché and written by Mario Vega.
  
  * spread-infection
  * spread-infection.exe

  #### /output
  A directory where all model outputs (results) are created.
  
  #### /output/rasters
  A directory specifically for raster model outputs ('maps' of the landscape)
