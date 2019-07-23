REM Set this to the directory where the executable and input files are
set root=...\input_files

REM The name of the compiled program
set exe=%root%\XYZ.exe

REM File listing which other files to import
set parms=default_parms.txt
REM File specifying parameters for a particular location
set site=site_files\site_1.txt

REM The name of the folder where the output will go
set output=output

set run=%exe% %root% %parms% %site% %output%    

REM Run the model using the options and parameters in the input files - run0.txt will be the name of the output file                                 
%run% run0 

REM interventions.txt contains a list of interventions that can be switched on from time 0 onwards - they are off by default (value of 0)
REM With ITNs after time 0
%run% run1 itn 1

REM With ITNs after time 0, at 70% coverage
%run% run2 itn 1 itn_coverage 0.7
   
REM SMC starting after 2 years
%run% run3 smc 1 smc_start 2   

REM MDA at 90% coverage
%run% run4 mda 1 mda_coverage 0.9

REM MSAT
%run% run5 mda 1 mda_screen 1

REM Any of the other parameters in the input files can be changed in the same way
REM Just add the pair parameter_name parameter_value to the command line as in these examples
REM It's not a good idea to change individual parameters in model_parms.txt or the larval parameters in mosq_species_parms.txt


REM The above use an exponential age distribution with a mean of 21 years
REM To use the age distribution for Tanzania:
set demog=demog\demog_Tanzania.txt demog\pop_Tanzania.txt
set run_demog=%exe% %root% %parms% %site% %demog% %output% 

%run_demog% run6
%run_demog% run7 itn 1

pause



