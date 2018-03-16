	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers and physics code for the simple cloud model
    module drivers
    use nrtype
    use variables
    private
    public :: model_driver_1d
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs model microphysics and advection for a number of steps
	!>writes to output file, solves for the microphysics over one time-step, 
	!> then advects particles
	!>@param[in] nq: number of q fields
	!>@param[in] kp: number of vertical levels
	!>@param[in] ord: order of advection scheme
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] runtime
	!>@param[in] dt - timestep
	!>@param[in] updraft_type - updraft_type
	!>@param[in] t_thresh - t_thresh
	!>@param[in] w_peak - w_peak
	!>@param[in] dz - grid spacing
	!>@param[inout] q, theta, pressure, z, temperature, rho, u
	!>@param[inout] precip
	!>@param[inout] new_file
	!>@param[inout] micro_init - flag to initialise microphysics
	!>@param[in] monotone - flag for monotonic advection
	!>@param[in] microphysics_flag - flag for calling microphysics
	!>@param[in] hm_flag - flag for switching on / off hm process
	!>@param[in] theta_flag - flag for advecting theta dry
	!>@param[in] mass_ice - mass of a single ice crystal (override)
    subroutine model_driver_1d(nq,kp,ord,o_halo,runtime, &
                               dt,updraft_type,t_thresh,w_peak, &
                               q,precip,theta,p,dz,z,t,rho,u,new_file,micro_init,monotone, &
                               microphysics_flag,hm_flag,theta_flag,mass_ice)

    use nrtype
    use advection
    use micro_module

    implicit none
    integer(i4b), intent(in) :: nq,kp, ord, o_halo, updraft_type
    real(sp), intent(in) :: runtime, dt, dz, t_thresh,w_peak
    real(sp), dimension(nq,-o_halo:kp+o_halo), intent(inout) :: q
    real(sp), dimension(4,1:kp), intent(inout) :: precip
    real(sp), dimension(-o_halo:kp+o_halo), intent(inout) :: theta, p, z, t,rho,u
    logical, intent(inout) :: new_file, micro_init
    logical, intent(in) :: monotone,microphysics_flag, hm_flag,theta_flag
    real(sp), intent(in) :: mass_ice

    ! local variables
    integer(i4b) :: nt, i, j, nsteps, iter
    real(sp) :: time


    nt=ceiling(runtime / real(dt,kind=sp) )
    do i=1,nt
    	time=real(i-1,sp)*dt
        print *,'time-step ',i,' of ',nt, ' time=',time

		
    	! kinematics:
		call kinematics_1d(nq,kp,o_halo,updraft_type, t_thresh, w_peak, &
							time,q,theta,p,z,t,u)    	
    
    
        ! output:
        call output_1d(time,nq,kp,q(:,1:kp),precip(:,1:kp),theta(1:kp),p(1:kp), &
                       z(1:kp),t(1:kp),u(1:kp),new_file)



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! advection                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		nsteps=ceiling(dt/(0.7_sp*dz)*(max(maxval(u),abs(minval(u)))))
		do iter=1,nsteps
			do j=1,nq
				! set halos
				call set_halos_1d(kp,ord,o_halo,q(j,:))            
				! advection
				call bott_scheme_1d(kp,ord,o_halo,dt/real(nsteps,sp), &
									dz,z,u,q(j,:),monotone)
			enddo
			if(theta_flag) then
				! set halos
				call set_halos_1d(kp,ord,o_halo,theta)        
				! advection
				call bott_scheme_1d(kp,ord,o_halo,dt/real(nsteps,sp),&
									dz,z,u,theta,monotone)
			endif
		enddo
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        ! solve microphysics. initialise constants that are commonly used, if needed
        if (microphysics_flag) then
			call microphysics_1d(nq,kp,o_halo,dt,dz,q,precip,theta,p, &
						   z,t,rho,u,micro_init,hm_flag,mass_ice)
			! calculate precipitation diagnostics
			
		endif       






    enddo



    end subroutine model_driver_1d


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>sets the vertical wind speed
	!>@param[in] nq: number of q fields
	!>@param[in] kp: number of vertical levels
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] updraft_type - updraft_type
	!>@param[in] t_thresh - t_thresh
	!>@param[in] w_peak - w_peak
	!>@param[in] time
	!>@param[inout] q, theta, pressure, z, temperature, u
    subroutine kinematics_1d(nq,kp,o_halo,updraft_type, t_thresh, w_peak, &
                               time,q,theta,p,z,t,u)

    use nrtype
    use advection

    implicit none
    integer(i4b), intent(in) :: nq,kp, o_halo, updraft_type
    real(sp), intent(in) :: time, t_thresh, w_peak
    real(sp), dimension(nq,-o_halo:kp+o_halo), intent(inout) :: q
    real(sp), dimension(-o_halo:kp+o_halo), intent(inout) :: theta, p, z, t,u

	select case (updraft_type)
		case(0)
			u(:)=w_peak
		case(1)
			if(time.lt.t_thresh) then
				u(:)=w_peak*sin(1._sp*pi*time/t_thresh)
			else
				u(:)=0._sp
			endif
		case(2)
			if(time.lt.t_thresh) then
				u(:)=0._sp
			endif
			if(time.ge.t_thresh) then
				u(:)=w_peak*sin(2._sp*pi*(time-t_thresh)/t_thresh2)
			endif
			if(time.ge.(t_thresh+t_thresh2)) then
				u(:)=0._sp
			endif
		case default
			print *,'select defined updraft_type ',updraft_type
	end select
	
    end subroutine kinematics_1d
    
    
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>set halos 1-d 
	!>@param[in] kp:  number of grid points
	!>@param[in] ord: order of interpolation scheme
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[inout] psi: field to set
    subroutine set_halos_1d(kp,ord,o_halo,psi)  

    use nrtype
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: kp, ord, o_halo
    real(sp), dimension(-o_halo:kp+o_halo), intent(inout) :: psi

	! top
	psi(kp+1:kp+o_halo)=psi(kp)
	! bottom
	psi(-o_halo:0)=psi(1)
	
	end subroutine set_halos_1d

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>output 1 time-step of model
	!>@param[in] time in seconds
	!>@param[in] nq number of q fields
	!>@param[in] kp number of vertical levels
	!>@param[in] q, precip, theta, pressure, z, temperature,u
	!>@param[inout] new_file
    subroutine output_1d(time,nq,kp,q,precip,theta,p,z,t,u,new_file)

    use nrtype
    use netcdf
    use variables, only : io1

    implicit none
    real(sp), intent(in) :: time
    integer(i4b), intent(in) :: nq,kp
    real(sp), dimension(nq,kp), intent(in) :: q
    real(sp), dimension(4,kp), intent(in) :: precip
    real(sp), dimension(kp), intent(in) :: theta, p, z, t, u
    logical, intent(inout) :: new_file

    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_create(outputfile, NF90_CLOBBER, io1%ncid) )

        ! define dimensions (netcdf hands back a handle)
        call check( nf90_def_dim(io1%ncid, "times", NF90_UNLIMITED, io1%x_dimid) )
        call check( nf90_def_dim(io1%ncid, "nq", nq, io1%nq_dimid) )
        call check( nf90_def_dim(io1%ncid, "nprec", 4, io1%nprec_dimid) )
        call check( nf90_def_dim(io1%ncid, "kp", kp, io1%k_dimid) )


        ! close the file, freeing up any internal netCDF resources
        ! associated with the file, and flush any buffers
        call check( nf90_close(io1%ncid) )


        ! now define some variables, units, etc
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
        ! define mode
        call check( nf90_redef(io1%ncid) )




        ! define variable: time
        call check( nf90_def_var(io1%ncid, "time", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "time", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "seconds") )

        ! define variable: q
        call check( nf90_def_var(io1%ncid, "q", NF90_DOUBLE, &
                    (/io1%nq_dimid, io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "q", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg or number per kg") )

        ! define variable: precip
        call check( nf90_def_var(io1%ncid, "precip", NF90_DOUBLE, &
                    (/io1%nprec_dimid, io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "precip", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "mm hr-1") )



        ! define variable: theta
        call check( nf90_def_var(io1%ncid, "theta", NF90_DOUBLE, &
                    (/io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "theta", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "K") )


        ! define variable: p
        call check( nf90_def_var(io1%ncid, "p", NF90_DOUBLE, &
                    (/io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "p", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "Pa") )


        ! define variable: z
        call check( nf90_def_var(io1%ncid, "z", NF90_DOUBLE, &
                    (/io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "z", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "m") )


        ! define variable: t
        call check( nf90_def_var(io1%ncid, "t", NF90_DOUBLE, &
                    (/io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "t", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "K") )

        ! define variable: u
        call check( nf90_def_var(io1%ncid, "u", NF90_DOUBLE, &
                    (/io1%k_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "u", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                    "units", "ms-1") )

        call check( nf90_enddef(io1%ncid) )
        call check( nf90_close(io1%ncid) )

        new_file=.false.

    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write data to file                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
    ! write variable: q
    call check( nf90_inq_varid(io1%ncid, "time", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, time, &
                start = (/io1%icur/)))

    ! write variable: q
    call check( nf90_inq_varid(io1%ncid, "q", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, q, &
                start = (/1,1,io1%icur/)))

    ! write variable: precip
    call check( nf90_inq_varid(io1%ncid, "precip", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, precip, &
                start = (/1,1,io1%icur/)))

    ! write variable: theta
    call check( nf90_inq_varid(io1%ncid, "theta", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, theta, &
                start = (/1,io1%icur/)))

    ! write variable: p
    call check( nf90_inq_varid(io1%ncid, "p", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, p, &
                start = (/1,io1%icur/)))

    ! write variable: z
    call check( nf90_inq_varid(io1%ncid, "z", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, z, &
                start = (/1,io1%icur/)))

    ! write variable: t
    call check( nf90_inq_varid(io1%ncid, "t", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, t, &
                start = (/1,io1%icur/)))

    ! write variable: u
    call check( nf90_inq_varid(io1%ncid, "u", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, u, &
                start = (/1,io1%icur/)))


    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine output_1d







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HELPER ROUTINE                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check(status)
    use netcdf
    use nrtype
    integer(I4B), intent ( in) :: status

    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
    end if
    end subroutine check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module drivers
