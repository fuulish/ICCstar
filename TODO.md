[ ] fix contrast calculation (i.e., do it only once in the beginning and then always exchange/communicate the data between processes)
[ ] fix errors in compute_efield (see valgrind memcheck log)
[ ] fix non-working pre_force calculation in setup_pre_force
[ ] find correct unit conversion factor (also double-check current numerical value)
[ ] detect wether kspace is used (and qsum_qsq needs to be called)
