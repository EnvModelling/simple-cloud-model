	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the simple cloud model and solves the 
	!>hydrostatic equation for pressure:
	!>\f$ \frac{\partial P}{\partial z} = - \frac{P}{R_aT}g \f$
    module initialisation
    use numerics_type
!    use variables
    private
    public :: calc_profile_1d, allocate_arrays
    contains
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for q_type and q_init
	!>@param[in] nq number of q fields
	!>@param[in] n_levels number of levels for reading in sounding
	!>@param[inout] q_type: integer array
	!>@param[inout] q_init: logical array
	!>@param[inout] q_read: real array
    subroutine allocate_arrays(nq,n_levels,q_type,q_init,q_read)
    use numerics_type
    implicit none
    integer(i4b), intent(in) :: nq, n_levels
    integer(i4b), dimension(:), allocatable, intent(inout) :: q_type
    logical, dimension(:), allocatable, intent(inout) :: q_init
    real(wp), dimension(:,:), allocatable, intent(inout) :: q_read
    ! local variables:
    integer(i4b) :: AllocateStatus
    
    ! allocate arrays
    allocate( q_type(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_init(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_read(1:nq,1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    
    
    end subroutine allocate_arrays



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>interpolates the sounding to the grid
	!>@param[in] nq number of q fields
	!>@param[in] nprec number of precipitation arrays
	!>@param[in] n_levels number levels for sounding
	!>@param[in] psurf surface pressure
	!>@param[in] tsurf surface temperature
	!>@param[in] t_cbase cloud base temperature
	!>@param[in] t_ctop cloud top temperature
	!>@param[in] adiabatic_prof: flag if we want an adiabatic profile
	!>@param[in] adiabatic_frac: fraction of adiabatic liquid water content in cloud
	!>@param[in] q_type flag for type of q-field
	!>@param[in] q_init flag for whether q-field initialised
	!>@param[in] z_read vertical levels for sounding
	!>@param[in] theta_read theta on vertical levels for sounding
	!>@param[in] q_read q fields on vertical levels for sounding
	!>@param[in] kp number of vertical levels of grid
	!>@param[in] o_halo number of extra grid levels required for advection
	!>@param[in] dz vertical resolution of grid
	!>@param[inout] dz2 vertical resolution of grid
	!>@param[inout] q, 
	!>@param[in] iqv,iqc,inc,iqi,ini
	!>@param[inout] precip, theta, pressure, z, temperature, rho,u
	!>@param[in] drop_num_init: flag to initialise number of drops where liquid water>0
	!>@param[in] number conc of drops #/kg
	!>@param[in] ice_init: flag to initialise ice crystals in model
	!>@param[in] number conc of ice crystals #/kg
	!>@param[in] mass of a single ice crystal kg.
	!>@param[in] microphysics_flag: flag for kind of microphysics used
    subroutine calc_profile_1d(nq,nprec,n_levels,psurf,tsurf,t_cbase, &
    						t_ctop, adiabatic_prof, adiabatic_frac, q_type,q_init, &
                             z_read,theta_read,q_read, &
                             kp,o_halo,dz,dz2,q,&
                             iqv, iqc, inc, iqi, ini, &
                             precip,theta,p,z,t,rho,u, &
                             drop_num_init, num_drop, ice_init,num_ice, mass_ice, &
                             microphysics_flag)
    use numerics_type
    use numerics, only : find_pos, poly_int, vode_integrate, zeroin
    use constants
    use variables, only : theta_surf, theta_q_sat, w_cb, t1old, p111

    implicit none
    ! inputs
    integer(i4b), intent(in) :: n_levels, nq,nprec,o_halo
    real(wp), dimension(n_levels), intent(in) :: z_read, theta_read
    real(wp), dimension(nq,n_levels), intent(in) :: q_read
    integer(i4b), dimension(nq), intent(in) :: q_type
    logical, dimension(nq), intent(in) :: q_init
    integer(i4b), intent(in) :: kp, microphysics_flag
    real(wp), intent(in) :: dz, psurf, tsurf, t_cbase, t_ctop
    logical, intent(in) :: adiabatic_prof, ice_init, drop_num_init
    real(wp), intent(in) :: adiabatic_frac
    real(wp), intent(in) :: num_drop, num_ice, mass_ice
    ! inouts
    real(wp), dimension(:), allocatable, intent(inout) :: theta, p, z, t, rho,u,dz2
    real(wp), dimension(:,:), allocatable, intent(inout) :: q, precip
    integer(i4b), intent(in) :: iqv, iqc, inc, iqi, ini

    ! local variables:
    integer(i4b) :: i,j, iloc, AllocateStatus, istore,istore2
    real(wp) :: var, dummy
	! variables for odesolver:
	real(wp), dimension(1) :: z1, p11,p22
	real(wp) :: htry,hmin,eps2,p1,p2,p_ctop,z11,z22,theta1

    ! allocate arrays
    allocate( precip(1:kp,1:nprec), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q(-o_halo+1:kp+o_halo,1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( theta(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( p(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( z(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( t(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( rho(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( u(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( dz2(-o_halo+1:kp+o_halo), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"


	precip=0._wp

    ! set up vertical level array
    z=dz*(/(i,i=-o_halo,kp+o_halo-1)/)+0.5_wp*dz ! note this is actually zn
    dz2=dz
	q=0._wp
	if (adiabatic_prof) then
		! calculate the dry adiabat:
		theta_surf=tsurf*(1.e5_wp/psurf)**(ra/cp)
		p1=psurf
		z1=0._wp
		p2=1.e5_wp*(t_cbase/theta_surf)**(cp/ra)
		htry=p2-psurf
		eps2=1e-5_wp
		call vode_integrate(z1,p1,p2,eps2,htry,hmin,hydrostatic1)
		p1=p2
		
		! integrate going downwards - dry adiabatic layer
		p(1)=psurf
		do i=1,-o_halo+2,-1
			z11=z(i)
			z22=z(i-1)
			p11=p(i)
			htry=-dz
			hmin=-1.e-2_wp
			call vode_integrate(p11,z11,z22,eps2,htry,hmin,hydrostatic1a)
			p(i-1)=p11(1)
		enddo
		! integrate going upwards - dry adiabatic layer
		p(1)=psurf
		do i=1,kp+o_halo-1
			z11=z(i)
			z22=z(i+1)
			if(z22.gt.z1(1)) exit
			p11=p(i)
			htry=dz
			hmin=1.e-2_wp
			call vode_integrate(p11,z11,z22,eps2,htry,hmin,hydrostatic1a)
			p(i+1)=p11(1)
		enddo
		istore=i-1
		! adiabatic temperature
		t(-o_halo+1:istore)=theta_surf*(p(-o_halo+1:istore)/1.e5_wp)**(ra/cp)
		! adiabatic vapour mixing ratio
		q(-o_halo+1:istore,iqv)=eps1*svp_liq(t_cbase)/(p1-svp_liq(t_cbase))


		! now calculate the moist adiabat
		w_cb=eps1*svp_liq(t_cbase)/(p1-svp_liq(t_cbase))
		theta_q_sat=t_cbase*(1.e5_wp/p1)**(ra/cp)*exp(lv*w_cb/cp/t_cbase)	
		! theta_q_sat is conserved. Use it to calculate the new temperature
		t1old=t_ctop

		p_ctop=zeroin(p1,3000._wp,calc_theta_q2,1.e-5_wp)

!		stop
		! integrate going upwards - moist adiabatic layer
		t1old=t_ctop
		istore=istore+1
		! now find the temperature
		p111=p(istore)
		t(istore)=theta_surf*( &
			(p(istore)-dz*p(istore)/ra/t_ctop)/1.e5_wp)**(ra/cp) ! a temperature colder than
		                                                   ! next level

		t(istore)=zeroin(1.01_wp*t(istore),t_ctop,calc_theta_q,1.e-5_wp)

		t1old=t(istore)
		do i=istore,kp+o_halo-1
			
			z22=z(i+1)
			if (i.eq.istore) then
				z11=z1(1)
				p11=p1
			else
				z11=z(i)
				p11=p(i)
			endif
			htry=dz
			hmin=1.e-2_wp
			!print *,t1old,p11,z11,z22
			call vode_integrate(p11,z11,z22,eps2,htry,hmin,hydrostatic2a)
			p(i+1)=p11(1)
			t(i+1)=theta_surf*(p(i+1)/1.e5_wp)**(ra/cp)
			t1old=t(i)
			p111=p(i+1)
			t(i+1)=zeroin(t(i+1),t1old*1.01_wp,calc_theta_q,1.e-5_wp)
			if(t(i+1).lt.t_ctop) exit
		enddo
		istore2=i-1
		do i=istore,istore2
			q(i,iqv)=eps1*svp_liq(t(i))/ &
								(p(i)-svp_liq(t(i)))
			q(i,iqc)=adiabatic_frac* &
					max(q(1,iqv)-eps1*svp_liq(t(i))/(p(i)-svp_liq(t(i))) ,0._wp)
			if(drop_num_init .and. &
			    ((microphysics_flag .eq. 2) .or. (microphysics_flag .eq. 3))) then
			    q(i,inc) = num_drop
			endif
		enddo

		! integrate going upwards - dry adiabatic layer
		theta1=t_ctop*(1.e5_wp/p_ctop)**(ra/cp)
		istore2=istore2+1
		do i=istore2,kp+o_halo-1
			z11=z(i)
			z22=z(i+1)
			p11=p(i)
			htry=dz
			hmin=1.e-2_wp
			call vode_integrate(p11,z11,z22,eps2,htry,hmin,hydrostatic1a)
			p(i+1)=p11(1)
		enddo
		t(istore2:kp+o_halo)=theta1*(p(istore2:kp+o_halo)/1.e5_wp)**(ra/cp)

		! initialise ice crystals
		if(ice_init .and. (microphysics_flag .eq. 1)) then
            where(t(istore:istore2).lt.ttr)
                q(istore:istore2,iqi)=num_ice*mass_ice
                q(istore:istore2,ini)=num_ice
            end where
        endif
		
	else
		! use linear interpolation to put sounding on grid:
		do i=-o_halo+1,kp+o_halo
			iloc=find_pos(z_read(1:n_levels),z(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(z(i),z_read(n_levels)), var,dummy)
			theta(i)=var


			! q-fields
			do j=1,nq
				if(q_init(j)) then
					! linear interp q fields
					call poly_int(z_read(iloc:iloc+1), q_read(j,iloc:iloc+1), &
								min(z(i),z_read(n_levels)), var, dummy)
					q(i,j)=var
				else
					if(q_type(j).eq.2) then
						q(i,j) = .0_wp
					else
						q(i,j) = 0.0_wp
					endif
				endif
			enddo

		enddo
	
		p(1)=psurf
		do i=1,kp+o_halo
			! solve the hydrostatic equation to get pressure and temperature
			t(i)=theta(i)*(p(i)/psurf)**(ra/cp)
			if(i.lt.(kp+o_halo)) then
				p(i+1)=p(i)-dz*p(i)/(ra*t(i))*grav
			endif
		enddo

		do i=0,-o_halo+1,-1
			! solve the hydrostatic equation to get pressure and temperature
			p(i)=p(i+1)+dz*p(i+1)/(ra*t(i+1))*grav
			t(i)=theta(i)*(p(i)/psurf)**(ra/cp)
		enddo
	endif

    u(:)=0
    theta=t*(1.e5_wp/p)**(ra/cp)
    rho=p/(ra*t)
    !q(6,:)=0._wp
    end subroutine calc_profile_1d


	subroutine hydrostatic1(p,z,dzdp)
	use numerics_type
	use constants
	use variables, only : theta_surf
	implicit none
	real(wp), intent(in) :: p
	real(wp), dimension(:), intent(in) :: z
	real(wp), dimension(:), intent(out) :: dzdp
	real(wp) :: t
	
	t=theta_surf*(p/1.e5_wp)**(ra/cp)
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic1

	subroutine hydrostatic1a(z,p,dpdz)
	use numerics_type
	use constants
	use variables, only : theta_surf
	implicit none
	real(wp), intent(in) :: z
	real(wp), dimension(:), intent(in) :: p
	real(wp), dimension(:), intent(out) :: dpdz
	real(wp) :: t
	
	t=theta_surf*(p(1)/1.e5_wp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a

	subroutine hydrostatic2(p,z,dzdp)
	use numerics_type
	use numerics, only : zeroin
	use constants
	use variables, only : theta_surf,theta_q_sat, w_cb, t1old, p111
	implicit none
	real(wp), intent(in) :: p
	real(wp), dimension(:), intent(in) :: z
	real(wp), dimension(:), intent(out) :: dzdp
	real(wp) :: t
	
	p111=p
	t=theta_surf*(p111/1.e5_wp)**(ra/cp)
	t=zeroin(t,t1old*1.01_wp,calc_theta_q,1.e-5_wp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic2

	subroutine hydrostatic2a(z,p,dpdz)
	use numerics_type
	use numerics, only : zeroin
	use constants
	use variables, only : theta_surf,theta_q_sat, w_cb, t1old, p111
	implicit none
	real(wp), intent(in) :: z
	real(wp), dimension(:), intent(in) :: p
	real(wp), dimension(:), intent(out) :: dpdz
	real(wp) :: t
	
	p111=p(1)
	t=theta_surf*(p111/1.e5_wp)**(ra/cp)
	t=zeroin(t,t1old*1.01_wp,calc_theta_q,1.e-5_wp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2a
	

	function calc_theta_q(t111)
	use numerics_type
	use constants
	use variables, only : theta_q_sat, p111
	implicit none
	real(wp), intent(in) :: t111
	real(wp) :: calc_theta_q
	real(wp) :: ws
	ws=eps1*svp_liq(t111)/(p111-svp_liq(t111))
	calc_theta_q=t111*(1.e5_wp/p111)**(ra/cp)*exp(lv*ws/cp/t111)-theta_q_sat

	end function calc_theta_q     

	function calc_theta_q2(p)
	use numerics_type
	use constants
	use variables, only : theta_q_sat, t1old
	implicit none
	real(wp), intent(in) :: p
	real(wp) :: calc_theta_q2
	real(wp) :: ws
	ws=eps1*svp_liq(t1old)/(p-svp_liq(t1old))
	calc_theta_q2=t1old*(1e5_wp/p)**(ra/cp)*exp(lv*ws/cp/t1old)-theta_q_sat

	end function calc_theta_q2    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp_liq(t)
		use numerics_type
		use constants, only : ttr
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_liq
		svp_liq = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-ttr)/ 234.5_wp)* &
			  (t-ttr)/(257.14_wp + (t-ttr)))
	end function svp_liq


    end module initialisation

