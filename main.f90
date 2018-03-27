	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Simple Cloud Model (SCM): 
	!>Solves a simple cloud model for investigating the hallet mossop process
	!>\f$ F\left(t,z \right)
	!>   = initialisation,microphysics,etc \f$
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>
	!> known bug: the lower halo is too large (i.e. wasted memory, should go from
	!> -o_halo+1 to kp+o_halo, rather than -o_halo to kp+o_halo
	!> o_halo might need 1 adding in order to work with advection scheme
	!> for now this bug does not affect any results.



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use nrtype
        use variables
        use initialisation
        use drivers
        implicit none

        character (len=200) :: nmlfile = ' '

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /sounding_spec/ adiabatic_prof, adiabatic_frac, &
        					psurf, tsurf, t_cbase, t_ctop, &
        					n_levels_s, q_read, theta_read, rh_read, z_read
        ! define namelists for environment
        namelist /run_vars/ outputfile, runtime, kp, dz, dt, &
                    nq, nprec, &
                    ord, halo, monotone, microphysics_flag,hm_flag, theta_flag, &
        			ice_init, num_ice, mass_ice, &
        			updraft_type, t_thresh, &
        			t_thresh2,w_peak
        namelist /run_vars2/ q_type, q_init
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        call allocate_arrays(nq,nlevels_r,q_type,q_init,q_read) 
        read(8,nml=run_vars2)
        read(8,nml=sounding_spec)
        close(8)
        o_halo=ord+2 !ord+1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_profile_1d(nq,nprec,n_levels_s,psurf,tsurf,t_cbase,t_ctop, &
        					adiabatic_prof, adiabatic_frac, &
        					q_type,q_init, z_read,theta_read, &
                            q_read,kp,o_halo,dz,grid1%q, grid1%precip, &
                            grid1%theta, grid1%p, &
                            grid1%z,grid1%t,grid1%rho,grid1%u, num_ice, mass_ice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call model_driver_1d(nq,nprec,kp,ord,o_halo,runtime,dt,updraft_type,t_thresh,w_peak, &
        					grid1%q,grid1%precip,grid1%theta, &
                            grid1%p,dz,grid1%z,grid1%t,grid1%rho,grid1%u,io1%new_file, &
                            micro_init,monotone,microphysics_flag,hm_flag,theta_flag, &
                            mass_ice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



