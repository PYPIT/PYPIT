The following is a list of items that are currently in development

* All items in the code with msgs.work("message") need to be implemented

* Write sorted files to fits table (let the user pick which method)
* Store master calibrations for each master frame (to save recalculating the solutions -- including slit trace, tilts, and wavelength)
* Quick reduction methods (for use in real time at the telescope)
* In function pypit.PYPIT, check for successful reduction and print reduction status (near line 125)
* Avoid re-making Bias, Arc, etc. frames as an option
* Organize lamp line lists better (1 file!!)
* Add code to confirm rejected arc lines are in NIST lists (e.g. when the NIST wavelengths change)
* Code to parse settings file into a Web Page, especially image type conditions
* Flat indexing in "setup", e.g. arspecobj.init_exp
* What should we do about multiple standard stars?
* Allow both bias and overscan to be used
* How much effort for backwards-compatability (e.g. old detectors)?  RJC -> 1 settings file per detector (e.g. detector upgrade -> new settings file)?
* Generate tests for adding new instrument
  extinction
  wavelengths
  settings file
