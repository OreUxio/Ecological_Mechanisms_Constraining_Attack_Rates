- Parameter_space_generator: 

	- Creates directories (within ../Data/) and subsequent configuration file with which to run the full assembly model.
	- Creates the "list_of_files_run.txt" file which is used by Full_assembly_run.sh to loop over paramter combinations.

- NewAsMini.cfg: A template configuration file. 
 
- Full_assembly_run.sh - Uses list_of_files_run.txt to loop over configuration files and run subsequent food-web assembly simulations.
 
- NewWeb* - Contains the C++ source code and compiled programs for simulating the full food-web assembly model.

- Invasion_Probability_Full.R - Simulates the invasion probaility and birth rate of mutants attempting to invade a food-web. 
