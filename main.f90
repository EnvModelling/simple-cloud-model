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
        use micro_module, only : set_qnames
        use w_micro_module, only : read_in_wmm_bam_namelist
        use p_micro_module, only : read_in_pamm_bam_namelist, p_initialise_aerosol_1d
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
                    advection_scheme, &
                    ord, halo, monotone, microphysics_flag, &
                    bam_nmlfile, aero_nmlfile, aero_prof_flag, hm_flag, theta_flag, &
        			drop_num_init, num_drop, ice_init, num_ice, mass_ice, &
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
        ! read in bulk aerosol namelists                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        grid1%nq=nq
        grid1%nprec=nprec
        select case(microphysics_flag)
            case(0:1) ! standard
                call set_qnames(grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,&
                    grid1%nq,grid1%ncat,grid1%nprec, &
                    grid1%iqv, grid1%iqc, grid1%ini, grid1%iqi)
            case(2) ! wmm
                call read_in_wmm_bam_namelist(bam_nmlfile, &
                    grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,grid1%nq,&
                    grid1%ncat, &
                    grid1%nprec, &
                    grid1%iqv, grid1%iqc, grid1%inc)
            case(3) ! pamm
                call read_in_pamm_bam_namelist(bam_nmlfile,aero_nmlfile, &
                    aero_prof_flag, &
                    grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,grid1%nq,&
                    grid1%ncat, &
                    grid1%nprec, grid1%n_mode, &
                    grid1%iqv, grid1%iqc, grid1%inc, grid1%cat_c, grid1%cat_r)    
			case default
				print *, 'error'
				stop
		end select
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_profile_1d(grid1%nq,grid1%nprec,n_levels_s,psurf,tsurf,t_cbase,t_ctop, &
        					adiabatic_prof, adiabatic_frac, &
        					grid1%q_type,q_init, z_read,theta_read, &
                            q_read,kp,o_halo,dz,grid1%dz2,grid1%q, &
                            grid1%iqv,grid1%iqc,grid1%inc,grid1%iqi,grid1%ini,&
                            grid1%precip, &
                            grid1%theta, grid1%p, &
                            grid1%z,grid1%t,grid1%rho,grid1%u, &
                            drop_num_init, num_drop, ice_init, num_ice, mass_ice, &
                            microphysics_flag)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise aerosol for microphysics_flag==3                          !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(microphysics_flag .eq. 3) then
            call p_initialise_aerosol_1d(aero_prof_flag, &
                        grid1%nq,grid1%ncat,grid1%c_s,grid1%c_e, &
                        grid1%inc, &
                        kp,o_halo, &
                        grid1%z,grid1%rho,&
                        grid1%p, grid1%t, &
                        grid1%q)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call model_driver_1d(grid1%nq,grid1%nprec, grid1%ncat,grid1%n_mode, &
                            kp,ord,o_halo,runtime,dt,updraft_type,t_thresh,w_peak, &
                            grid1%c_s,grid1%c_e, &
                            grid1%inc,grid1%iqc, &
                            grid1%cat_c, grid1%cat_r, &
                            grid1%q_name, &
        					grid1%q,grid1%precip,grid1%theta, &
                            grid1%p,dz,grid1%dz2,grid1%z,grid1%t,grid1%rho,grid1%u,io1%new_file, &
                            micro_init,advection_scheme, &
                            monotone,microphysics_flag,hm_flag,theta_flag, &
                            mass_ice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *,'SCM has finished'

    end program main



