Repository containing all programs that are useful for observing the
crab pulsar with Leuschner. Documentation for code is below:

If you're a new user, do this everytime you log into crab till you
find someone who can show you how to modify your own bashrc file.

Source /home/nch/.bashrc
source /home/vishal/.bashrc

**remove_RFI.c** :Combines filterbank files in a given directory and
then removes RFI from data by simple SD threshold test.

*Prereqs (header files):
	 *rand_num_generator.h for creating gaussian distributions to replace RFI with.

	 *filterbank.cpp for opening the original filterbank file and
	 adjusting header values and writing the new cleaned filterbank file.

	 *make_cmd : Make file to compile remove_RFI.c

*To run: 
Type ./remove_RFI.o in the directory where the filterbank files are.
The macros at the top can be set to TRUE and FALSE accordingly
depending on what needs to be done

**removal.c** :This is like remove_RFI.c except only removes RFI from a single
file who's name must be manually edited in the code itself.  Requires
same header files as remove_RFI.c

**create_fil-single.py** :This is the RFI removal code/filterbank
  creation code that takes the raw binary data from Deepthi's acq data
  program, and converts it to filterbank format.

Three arguments need to be passed in when the program is called:
1. The output filename
2. The start channel for the filterbank file
3. The end channel for the filterbank file

**adc_snap_MOPG2.py** :shows the output of the ADC

**spectrum_monitor.py** shows the output spectrum for the four Stokes parameters from the SNAP

