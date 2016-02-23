# LAMMPS

#Update 1.1.3 (2/23/2016 @ 12:00pm)
-Minor code cleanup

-File reader will terminate program if input files are not found, rather than running program with incorrect values.

#Update 1.1.255 (1/20/2016 @ 11:58pm)
-Fixed an issue where the uncertainty coefficients were not being read in correctly due to a spelling error.

#Update 1.1.25 (1/22/2016 @ 10:53am)
-The uncertainty can now be written to file via the dump command

-Fixed issue with force symmetrization in pair_agni

#update 1.1.2 (1/20/2016 @1:22pm)
-Updated the check for periodicity in pair_agni

-Updated compute_uncertainty (this is not ready for implementation though!)

#Update 1.1.15 (1/18/2016 @11:30am)
-Updated all user arrays from vectors to pointers.

-Added an uncertainty class, and implemented uncertainty calculations.

-Updated the test.agni file. 

-Multi-elemental support has been removed for future update.

#Update 1.1.1 (12/21/2015 @10:38am)
-Massive code clean-up.

#Update 1.1.05 (12/7/2015 @9:38am)
-The file reader now supports multi-elemental systems. 

#Update 1.1.0 (11/25/2015 @12:15PM)
-MPI support has been established in the pair_agni files.

#Update 1.0.025 (11/3/2015 @10:25AM)
-The user can now call their input file for Pair_Agni whatever they wish, as long as they follow the outline in the included in.eam file.

#Update 1.0.02 (11/2/2015 @9:18AM)
-The user can now place values in test.agni whever they wish, as long as they follow the outline for placeing endVar.

#Update 1.0.015 (10/12/2015 @2:56PM)
-The user no longer has to label nTrain values for each line.

#Update 1.0.01 (10/9/2015 @11:25PM)

-Fixed issue with file reader breaking after nTrain > 3.

#Update 1.0.005 (10/9/2015 @10:58PM)
-Changed protected variables seetings in pair_user_temp.h to non-static reference to avoid compiler compatibility issues. 


