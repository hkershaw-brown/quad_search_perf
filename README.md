
# Usage

`./test_quad_search num_repeats grid p (print last lon lat)`  

e.g. 
  `./test_quad_search 1 nx_4500_ny_5000_irregular_grid.nc p`  

## Different versions of quad\_utils\_mod.f90

```
brute_quad_utils_mod

cp brute_quad_utils_mod quad_utils_mod.f90
./quickdbuild.sh
```

Use original

```
rm quad_utils_mod.f90
./quickdbuild.sh
```

two pass version

```
cp quad_utils_twopass quad_utils_mod.f90
./quickbuild.sh

```
## Memory usage

On mac

```
/usr/bin/time -l ./test_quad_search 1 nx_4500_ny_5000_irregular_grid.nc
```

### Time vs Memory notes

```
Program test_quad_search v11.22.0					 Program test_quad_search v11.22.0
 --------------------------------------					 --------------------------------------
									
  set_nml_output Echo NML values to log file only			  set_nml_output Echo NML values to log file only
 PE 0: initialize_mpi_utilities:  Running with           		 PE 0: initialize_mpi_utilities:  Running with           
 num_reps  =            1						 num_reps  =            1
 grid_file = nx_4500_ny_5000_irregular_grid.nc				 grid_file = nx_4500_ny_5000_irregular_grid.nc
 Nx =         5000							 Nx =         5000
 Ny =         4500							 Ny =         4500
 locate time =   0.59411800000000003     			|        PE 0: init_irreg_interp  to determine (minimum) max_reg_
								>        PE 0: init_irreg_interp ...  interp_handle%ii%grid_num i
								>        locate time =    1.3999999999958490E-005
 lon corner index =         4450        4451        4451 		 lon corner index =         4450        4451        4451 
 lat corner index =         4184        4184        4185 		 lat corner index =         4184        4184        4185 
									
 --------------------------------------					 --------------------------------------
 Finished ... at YYYY MM DD HH MM SS = 					 Finished ... at YYYY MM DD HH MM SS = 
                 2026  4 17 16 20 26				|                        2026  4 17 16 19 23
 Program test_quad_search v11.22.0					 Program test_quad_search v11.22.0
 --------------------------------------					 --------------------------------------
									
        1.13 real         0.64 user         0.16 sys		|               1.70 real         1.14 user         0.47 sys
           747995136  maximum resident set size			|                 6265372672  maximum resident set size
                   0  average shared memory size			                   0  average shared memory size
                   0  average unshared data size			                   0  average unshared data size
                   0  average unshared stack size			                   0  average unshared stack size
               45839  page reclaims				|                     382633  page reclaims
                  59  page faults				|                         44  page faults
                   0  swaps						                   0  swaps
                   0  block input operations				                   0  block input operations
                   0  block output operations				                   0  block output operations
                   0  messages sent					                   0  messages sent
                   0  messages received					                   0  messages received
                   0  signals received					                   0  signals received
                   4  voluntary context switches		|                          1  voluntary context switches
                 192  involuntary context switches		|                        250  involuntary context switches
         15718200888  instructions retired			|                17365909574  instructions retired
          2624557598  cycles elapsed				|                 5436004588  cycles elapsed
           726862592  peak memory footprint			|                 6247118080  peak memory footprint
```

# Results CSR vs main

```
./test_quad_search.main  100 regional_4500_5000.nc  p > main.out  
./test_quad_search  100 regional_4500_5000.nc r p > csr.out
```

```diff
(py-env) [hkershaw:runs]() > diff csr.out main.out 
4,5c4,5
<                  2026  6  4  8 32 59
<  Program test_quad_search v11.24.0-3-g24d04cfdf-dirty
---
>                  2026  6  4  8 32 39
>  Program test_quad_search v11.24.0-dirty
17,24c17,26
<  init time =    1.2169999999999959E-003
<  PE 0: init_irreg_interp  two-pass: max candidates per coarse box =           25
<  PE 0: init_irreg_interp  two-pass: min candidates per coarse box =            9
<  PE 0: init_irreg_interp  two-pass: empty coarse boxes =            0  of      2812329
<  PE 0: init_irreg_interp  two-pass: total coarse-index entries =     45616723
<  PE 0: init_irreg_interp two-pass: mean candidates per coarse box =   16.2
<  PE 0: init_irreg_interp  two-pass: boxes with >2x mean =            0  of      2812329
<  random number creation time =    10.450315000000000     
---
>  init time =    6.1400000000000343E-004
>  PE 0: init_irreg_interp  to determine (minimum) max_reg_list_num values for new grids ...
>  PE 0: init_irreg_interp ...  interp_handle%ii%grid_num is           48
>  PE 0: init_irreg_interp init_irreg_interp: max candidates per coarse box = 48
>  PE 0: init_irreg_interp init_irreg_interp: min candidates per coarse box = 30
>  PE 0: init_irreg_interp init_irreg_interp: empty coarse boxes = 0  of 810000
>  PE 0: init_irreg_interp init_irreg_interp: total coarse-index entries = 32809050
>  PE 0: init_irreg_interp init_irreg_interp: mean candidates per coarse box =   40.5
>  PE 0: init_irreg_interp init_irreg_interp: boxes with >2x mean = 0  of 810000
>  random number creation time =    10.543972000000000     
225c227
<  locate time =    3.8900000000019475E-004
---
>  locate time =    6.0099999999962961E-004
229,230c231,232
<                  2026  6  4  8 33 12
<  Program test_quad_search v11.24.0-3-g24d04cfdf-dirty
---
>                  2026  6  4  8 32 51
>  Program test_quad_search v11.24.0-dirty
```

