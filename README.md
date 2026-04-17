
### Timing notes

Orig vs Brute force

Original 
1
locate time =    3.9999999999970615E-006
 lon corner index =          375         376         376         375
 lat corner index =           13          13          14          14


1000
locate time =    1.0590000000000044E-003

PE 0: initialize_mpi_utilities:  Running with            1  MPI processes.
 Nx =         5000
 Ny =         4500
 PE 0: init_irreg_interp  to determine (minimum) max_reg_list_num values for new grids ...
 PE 0: init_irreg_interp ...  interp_handle%ii%grid_num is           59
 locate time =    3.0199999999998006E-003

Brute force
1
 locate time =    2.5100000000000122E-004
 lon corner index =          375         376         376         375
 lat corner index =           13          13          14          14

1000
locate time =   3.1

Nx =         5000
 Ny =         4500
 PE 0: init_irreg_interp  to determine (minimum) max_reg_list_num values for new grids ...
 PE 0: init_irreg_interp ...  interp_handle%ii%grid_num is           59
 locate time =    318.55002400000001     
