!  no skimmer version


!----- bounds_mod -----
! This fortran module defines all of the boundary condition stuff
module bounds_mod
  use core_types_mod
  implicit none



  ! REQUIRED INTERFACE DATA:

  ! the list of ghost cells
  type(cell_type), allocatable, dimension(:) :: gc

  ! a list of integers indicating ghost cell type
  ! 1=left, 2=top, 3=right   (Note: there is no need for a bottom)
  integer, allocatable, dimension(:) :: gc_type

  ! the number of ghost cells
  integer num_gc
  integer jstart,jend
  integer istart,iend



  ! DATA USED INTERNALLY:

  integer start_gc_capac
  parameter(start_gc_capac=10)

  real(8),allocatable,dimension(:) :: vz_top
  real(8),allocatable,dimension(:) :: vr_top
  real(8),allocatable,dimension(:) :: vmag_top
  real(8),allocatable,dimension(:) :: T_top
  real(8),allocatable,dimension(:) :: P_top

  real(8),allocatable,dimension(:) :: vz_end
  real(8),allocatable,dimension(:) :: vr_end
  real(8),allocatable,dimension(:) :: vmag_end
  real(8),allocatable,dimension(:) :: T_end
  real(8),allocatable,dimension(:) :: P_end

  real(8),allocatable,dimension(:) :: ur
  real(8),allocatable,dimension(:) :: uz

  real(8) Ptarget,ufix
  parameter(Ptarget=131.6,ufix=184.d0)

contains


  !--- initialize_boundaries
  ! initializes the ghost cells in a default configuration
  subroutine initialize_boundaries
    use constants_mod
    use core_data_mod
    use inflow_mod
    use parallel_info_mod
    implicit none
    integer i
    real(8) margin

    ! "margin" specifies how far out the ghost cells extend normal
    ! to the computing region

    margin = sqrt(kb*simT/m)*30*tau

    ! initializing the ghost cells

        ! for each mpi_rank figure out the number of boundary edges for
        ! this processor

        ! here the the test keys for each kind of ghost cell

        ! Left edge:  mod(mpi_rank,num_proc_z)=0
        ! Right edge: mod(mpi_rank+1,num_proc_z)=0
        ! Top edge:   mpi_rank+1 >= num_proc_z*(num_proc_r-1)

    ! allocate gc and gc_type to the maximum possible of 3 (top, left, right)
    allocate(gc(3))
    allocate(gc_type(3))

    ! start it at zero
    num_gc=0

    ! check to see if this processor touches the left edge
    ! if so increment num_gc 1 and set gc_type(1) to 1
    if(mod(mpi_rank,num_proc_z).eq.0) then
        num_gc=num_gc+1
        gc_type(num_gc)=1
    end if

    ! check to see if this processor touches the top
    ! if so increment num_gc to 1 and set gc_type(1) to 2
    if(mpi_rank.ge.num_proc_z*(num_proc_r-1)) then
        num_gc=num_gc+1
        gc_type(num_gc)=2
    end if

    ! check to see if this processor touches the right edge
    ! if so increment num_gc 1 and set gc_type(1) to 3
    if(mod(mpi_rank+1,num_proc_z).eq.0) then
        num_gc=num_gc+1
        gc_type(num_gc)=3
    end if

    ! now initialize the cells
    do i=1,num_gc
        call initialize_gc(i, margin)

    end do

    ! Initialize the inflow_mod
    call init_inflow(mpi_rank)

    ! allocating the arrays for saved boundary condition information
    allocate(vz_top(num_s_cells_z+1))
    allocate(vr_top(num_s_cells_z+1))
    allocate(vmag_top(num_s_cells_z+1))
    allocate(T_top(num_s_cells_z+1))
    allocate(P_top(num_s_cells_z+1))

    allocate(vz_end(num_s_cells_r+1))
    allocate(vr_end(num_s_cells_r+1))
    allocate(vmag_end(num_s_cells_r+1))
    allocate(T_end(num_s_cells_r+1))
    allocate(P_end(num_s_cells_r+1))

    ! loading the saved boundary condition information
     open(unit=900,file='top_bound.txt',status='old',action='read')
     do i=1,num_s_cells_z
       read(900,*) vz_top(i),vr_top(i),T_top(i),P_top(i)
     end do
     close(unit=900)


     open(unit=900,file='end_bound.txt',status='old',action='read')
     do i=1,num_s_cells_r
       read(900,*) vz_end(i),vr_end(i),T_end(i),P_end(i)
     end do
     close(unit=900)


  end subroutine initialize_boundaries

  !--- initialize_gc ---
  ! initializes a single ghost cell
  subroutine initialize_gc(i, margin)
    use helpers_mod
    use simulation_mod
    use constants_mod
    implicit none
    integer i
    real(8) margin

    if(gc_type(i) .eq. 1)then !left gc
      gc(i)%min_z = zmin_local - margin
      gc(i)%max_z = zmin_local
      gc(i)%min_r = max(rmin_local - margin, 0.d0)
      gc(i)%max_r = rmax_local + margin

    else if(gc_type(i) .eq. 2)then !top gc
      if(zmin_local .eq. zmin_global)then
        gc(i)%min_z = zmin_local
      else
        gc(i)%min_z = zmin_local - margin
      end if
      if(zmax_local .eq. zmax_global)then
        gc(i)%max_z = zmax_local
      else
        gc(i)%max_z = zmax_local + margin
      end if
      gc(i)%min_r = rmax_local
      gc(i)%max_r = rmax_local + margin

    else if(gc_type(i) .eq. 3)then !right gc
      gc(i)%min_z = zmax_local
      gc(i)%max_z = zmax_local + margin
      gc(i)%min_r = max(rmin_local - margin, 0.d0)
      gc(i)%max_r = rmax_local + margin

    end if

    gc(i)%min_th = -pi
    gc(i)%max_th = pi
    gc(i)%volume = (gc(i)%max_z-gc(i)%min_z) * (gc(i)%max_th-gc(i)%min_th)/2.d0 * (gc(i)%max_r**2 - gc(i)%min_r**2)
    gc(i)%mask = 0
    gc(i)%num_parts = 0
    gc(i)%partition = 0
    gc(i)%num_left = 0
    gc(i)%capacity = start_gc_capac
    allocate(gc(i)%ps(start_gc_capac))

! initialize everything else to zero
   gc(i)%efield_r=0.d0
   gc(i)%efield_z=0.d0



  end subroutine initialize_gc


  !--- update_bc_data ---
  ! Looks at the sampling cells in the exit window at rmax and
  ! adaptively produces consistent boundary conditions in this
  ! window. After much toil and trouble we found a couple of methods
  ! that work, maybe.

  ! Ambitious approach:

  ! Specify an average temperature across the window. Let the code
  ! choose the profile T(z), but rescale it every time this subroutine
  ! is called so that it has the desired average T. This target profile
  ! is called T_top (T_top is stored in the file (top_bound.txt).

  ! Specify the average pressure and do the same thing as for T, creating
  ! P_top. This profiles is not stored but is created in each update
  ! cycle by using the simulation density profile and T_top to find
  ! the profile of P that the simulation is choosing. This P profile is
  ! then used with the T profile to determine the density of particles
  ! in the ghost cell outside of the exit window.

  ! Let the code choose the profiles and magnitudes of Vr and Vz.

  ! Less ambitious approach:

  ! Specify both the pressure and the temperature at the top and
  ! at the end. These are in Ptarget and Ttarget

  ! let the velocities choose themselves by watching what particles
  ! are doing at the edge and updating when update_bc_data is called,
  ! on the output skip. This means that when you choose the output skip,
  ! you are also choosing how often the boundary conditions are updated.
  ! Not ideal, and even dangerous. An output skip less that 1000 makes
  ! for jumpy boundary conditions in the velocity

! NOTE: This routine could be called by any processor, but the sampling
! cell data will only have been accumulated for the sampling cells belonging
! that processor. This will defeat the purpose of this routine taking care
! of updating boundary conditions all around the edge of the simulation.
! To prevent this, this routine will only be called by processor 0, which
! has had the sum of all sampling data over all processors loaded into its
! sampling cells by write_collected_data.




  subroutine update_bc_data
    use core_data_mod
    use constants_mod
    use core_types_mod
    use simulation_mod
    use parallel_info_mod
    implicit none
    integer i,j,itop,switch
    real(8) n, vz, r, theta
    real(8) temp,Tavg,Ttarget,Pavg,PTarget,Pcorrect

    real(8) twopi,under

    type(sampling_cell_type) c

    data twopi/6.283185307d0/

! if this isn't processor 0, die
if(mpi_rank.ne.0) then
   write(*,*) 'mpi_rank isn''t 0 in update_bc_data. Dying.'
   stop 456
end if


 ! set the starting and ending sampling cells along the top for computing
 ! the outward flux and for updating boundary condition quantities

 ! This has to be looked at everytime you try a new geometry

! use the polygon information read in

  ! next to the last point in the sampling polygon, sampling cone
  ! intersection with the top boundary
    jstart = polys(1)%pts(polys(1)%num_pts-1)%z/s_cells_dz   + .1

  ! first point in the skimmer polygon, gas-facing edge
  ! of the skimmer intersecting the top boundary
    jend   = num_s_cells_z   !   the last j is at the very end of the computing region in z for this
                   !   no-skimmer simulation


  ! top edge, gotten from the first point in the skimmer polygon
    itop = num_s_cells_r  ! top edge

! Start of Temperature control


   Ttarget = 2740.d0 ! the target magnitude for Tavg

! let the code choose the distribution of T(z) along the top, but
! control its magnitude

! compute the average temperature across the outflow segment at the top
! and along the right downstream edge

   Tavg=0.d0
   i=itop

     if(mpi_rank.eq.0) then
   write(*,*) 'jstart,jend: ',jstart,jend
   end if

   do j=jstart,jend
     c = s_cells(i,j)

     if(c%sum_n.gt.0.d0) then
        temp=m/3.d0/kb*( (c%sum_vx2+c%sum_vy2+c%sum_vz2)/c%sum_n - &
                      (c%sum_vx**2+c%sum_vy**2+c%sum_vz**2)/c%sum_n**2 )
!       T_top(j) = temp
        T_top(j) = Ttarget ! nail the temperature
     end if
        ! if there are no particles out here just use what was read in via T_top

     Tavg = Tavg + T_top(j)

!    if(mpi_rank.eq.0) then
!    write(*,"(i5,1pe12.4,' ',a5)")  j,temp,'$$'
!  end if

   end do

   j=jend
     if(mpi_rank.eq.0) then
   write(*,*) 'itop: ',itop
   end if
   do i=1,itop
     c = s_cells(i,j)

     if(c%sum_n.gt.0.d0) then
        temp=m/3.d0/kb*( (c%sum_vx2+c%sum_vy2+c%sum_vz2)/c%sum_n - &
                      (c%sum_vx**2+c%sum_vy**2+c%sum_vz**2)/c%sum_n**2 )
!       T_end(i) = temp
        T_end(i) = Ttarget ! nail the temperature
     end if
        ! if there are no particles out here just use what was read in via T_top

     Tavg = Tavg + T_end(i)
!     if(mpi_rank.eq.0) then
!     write(*,"(i5,1pe12.4,' ',a5)")  i,temp,'%%'
!     end if

   end do


   Tavg = Tavg/( real(jend-jstart+1) + real(itop) )
   if(mpi_rank.eq.0) then
      write(*,*) 'Tavg: ',Tavg
   end if


!  adjust T_top and T_end to hit this target
!  if(Tavg.gt.0.d0) then
!     do j=jstart,jend
!       T_top(j) = T_top(j)*Ttarget/(Tavg+1d-8) ! protect against division by zero
!     end do
!  end if

!  if(Tavg.gt.0.d0) then
!     do i=1,itop
!       T_end(i) = T_end(i)*Ttarget/(Tavg+1d-8) ! protect against division by zero
!     end do
!  end if

! extend beyond the window to pick up the ghost cell extension regions
! on either side
!   T_top(jstart-1)=2*T_top(jstart)-T_top(jstart+1)
!   T_top(jstart-2)=3*T_top(jstart)-2*T_top(jstart+1)
!   T_top(jend+1)=2*T_top(jend)-T_top(jend-1)
!   T_top(jend+2)=3*T_top(jend)-2*T_top(jend-1)

! End of Temperature control

! Start of Pressure control

! let the code choose the distribution of P(z) along the top, but
! control its magnitude

! to make sure that we don't add particles in places other than
! the exit window, set P_top to zero inside the sampler


!   do j=1,max(1,jstart-3)
!    P_top(j)=0.d0
!   end do



! compute the average pressure across the outflow region at the top
!won't converge   under=.2 ! underrelaxation parameter
!won't converge   Pavg=0.d0
!won't converge   i=itop
!won't converge   do j=jstart,jend
!won't converge     c = s_cells(i,j)
!won't converge
!won't converge     if(c%sum_n.gt.0.d0) then
!won't converge        ! density
!won't converge        temp=c%sum_n*Nef/c%volume/num_samples
!won't converge        P_top(j) = (1.d0-under)*P_top(j) + under*temp*kb*T_top(j) ! P = nkT
!won't converge     end if
!won't converge        ! just keep the P_top that was read in
!won't converge
!won't converge     Pavg = Pavg + P_top(j)
!won't converge
!won't converge   end do
!won't converge   Pavg = Pavg/real(jend-jstart+1)
!won't converge
!won't converge   PTarget = P       ! the target magnitude for Pavg
!won't converge   if(mpi_rank.eq.0) then
!won't converge      write(*,*) 'PTarget, Pavg',Ptarget,Pavg
!won't converge   end if
!won't converge
!won't converge
!won't converge   Pcorrect=Ptarget/(Pavg+1.d-9) ! protect against division by zero
!won't converge
!won't converge!  adjust P_top to hit this target
!won't converge   if(Pavg.gt.0.d0) then
!won't converge      do j=jstart,jend
!won't converge        P_top(j) = P_top(j)*Pcorrect
!won't converge      end do
!won't converge   end if

! punt and peg it. The code won't have a pressure that is constant like this, but
! it will be in the neighborhood of it and that seems to be the best I can do
!    do j=jstart,jend
!     P_top(j)=Ptarget
!    end do

!    do i=1,itop
!     P_end(i)=Ptarget
!    end do

! comment this loop if you just want to use the pressure profile read in

! extend beyond the window to pick up the ghost cell extension regions
! on either side along the top by the sampler
!   P_top(jstart-1)=2*P_top(jstart)-P_top(jstart+1)
!   P_top(jstart-2)=3*P_top(jstart)-2*P_top(jstart+1)


! End of Pressure control


! now load vz_top and vr_top
    switch=0
    i=itop


    do j=jstart,jend

         c = s_cells(i,j)

! Let vr and vz do what the code wants
! if there is not data just keep what we had before from top_bound.txt

! try underrelaxing
        under=0.8;

! kluge
   write(*,507) mpi_rank,c%sum_n,c%sum_vx,c%sum_vz
507 format(' rank, n, vx, vz: ',i4,3(1x,1pe12.4))


        if(c%sum_n.gt.0.d0) then
          vr_top(j) = c%sum_vx/c%sum_n*under + (1.d0-under)*vr_top(j)
          vz_top(j) = c%sum_vz/c%sum_n*under + (1.d0-under)*vz_top(j)
          switch=1
        end if
        ! if there are no particles, just keep what was read in


if(c%sum_n.ne.0) then
   write(*,401) mpi_rank,j,c%sum_n,c%sum_vx,c%sum_vz
401 format(' $top: ',i4,i4,3(1x,1pe12.4))
end if

   end do

    switch=0
    j=jend
    do i=1,itop

         c = s_cells(i,j)

! Let vr and vz do what the code wants
! if there is not data just keep what we had before from top_bound.txt

        if(c%sum_n.gt.0.d0) then
          vr_end(i) = c%sum_vx/c%sum_n*under + (1.d0-under)*vr_end(i)
          vz_end(i) = c%sum_vz/c%sum_n*under + (1.d0-under)*vz_end(i)
          switch=1
        end if
        ! if there are no particles, just keep what was read in

   end do



! extrapolate to the margins
!    vr_top(jstart-1)=2*vr_top(jstart)-vr_top(jstart+1)
!    vr_top(jstart-2)=3*vr_top(jstart)-2*vr_top(jstart+1)
!    vr_top(jend+1)=2*vr_top(jend)-vr_top(jend-1)
!    vr_top(jend+2)=3*vr_top(jend)-2*vr_top(jend-1)
!
!    vz_top(jstart-1)=2*vz_top(jstart)-vz_top(jstart+1)
!    vz_top(jstart-2)=3*vz_top(jstart)-2*vz_top(jstart+1)
!    vz_top(jend+1)=2*vz_top(jend)-vz_top(jend-1)
!    vz_top(jend+2)=3*vz_top(jend)-2*vz_top(jend-1)

! no velocity beyond the margins - connect to zero flow at the wall
! also zero the velocity in the cell next to the wall
    vr_top(jstart)=0.d0
    vz_top(jstart)=0.d0


    vr_top(jstart-1)=0.d0
    vr_top(jstart-2)=0.d0


    vz_top(jstart-1)=0.d0
    vz_top(jstart-2)=0.d0

    if(mpi_rank.eq.0) then
      write(*,*) ' T_top(167) :',T_top(167), ' T_end(100) :',T_end(100)
      write(*,*) ' P_top(167) :',P_top(167)
      write(*,*) ' dens_top(167) :',s_cells(itop,167)%sum_n/s_cells(itop,167)%volume*Nef/num_samples
117   format(' n(167): ',1pe12.4)
    end if


! write these quantities if there are particles at the outer edge to work with
! you can tell when there aren't any because Tavg will be zero

     if(switch.gt.0.d0) then

       open(unit=900,file="top_bound.txt",status="replace")
       do j=1,num_s_cells_z
         write(900,101) vz_top(j),vr_top(j),T_top(j),P_top(j)
101      format(4(1x,1pe12.4))
!kluge
         write(*,407) mpi_rank,j,vr_top(j)
407 format(' proc, j, vr_top ',2(1x,i4),1pe12.4)
       end do
       close(unit=900)

       open(unit=900,file="end_bound.txt",status="replace")
       do i=1,itop
         write(900,101) vz_end(i),vr_end(i),T_end(i),P_end(i)
       end do
       close(unit=900)

     end if

     end subroutine update_bc_data









  !--- do_bc ---
  ! actually flushes and repopulates the ghost cells
  subroutine do_bc
    implicit none
    integer i,j

    do i=1,num_gc
      !flush the ghost cell
      gc(i)%num_parts = 0
      gc(i)%partition = 0

      if(gc_type(i) .eq. 1)then
        call do_left_bc(gc(i))
       else if(gc_type(i) .eq. 2)then
         call do_top_bc(gc(i))
       else if(gc_type(i) .eq. 3)then
         call do_right_bc(gc(i))
      end if

    end do

  end subroutine do_bc









  !-- do_left_bc ---
  ! this subroutine repopulates a left-type ghost cell
  subroutine do_left_bc(c)
    use inflow_mod
    implicit none
    type(cell_type) c

    if(c%min_r .ge. 0.519e-3) return

    call fill_inflow_cell(c)
   call fill_inflow_cell_trace(c)

  end subroutine do_left_bc




  !-- do_top_bc ---
  ! this subroutine repopulates a top-type ghost cell
  subroutine do_top_bc(c)
    use core_data_mod
    use core_types_mod
    use helpers_mod
    use constants_mod
    use simulation_mod
    use parallel_info_mod

    implicit none
    type(cell_type) c
    integer Ncandidates, i, Nloaded
    real(8) ::  n, T, x, y, z, r, temp, nmax, delr2, fracandidates=0.d0, rN
    real(8) zstart,zend,margin,vztmp,vrtmp,Ttmp,ntmp,vth
    real(8) vperp,theta,twopi
    type(particle_type) p0
    integer idx

    data twopi/6.283185307d0/



! just do vacuum here
!
!    c%num_parts = 0
!
!    c%partition = 0
!    return


    ! This routine adds particles at the top edge (maximum r)
    ! and must be controlled so that particles don't get injected
    ! into metal at the edge. This you have to handle with your
    ! own logic. In this version particles are loaded between
    ! zstart-margin and zend+margin. This loads a few particles
    ! into the metal, but not many. In this routine margin
    ! is simply the z-width of one sampling cell

    ! note that jstart and jend are given values in update_bc_data

    zstart=c%min_z
    zend=c%max_z



! find the maximum density along this cell's portion of the top
    jstart = (zstart-zmin_global)/s_cells_dz + 1
    jend   = (zend-zmin_global)/s_cells_dz + 1
    nmax=0.d0
    do idx=jstart,jend
       nmax = max( nmax,P_top(idx)/kb/T_top(idx) )
    end do

    nmax=nmax*2.d0 ! just make sure it's big enough

    nmax=4.d22


! compute the number of candidates to use


    rN = c%volume * nmax / real(Nef)
    Ncandidates = rN + fracandidates
    fracandidates = fracandidates + rN-Ncandidates


    delr2=c%max_r**2-c%min_r**2

    Nloaded=0

    do i=1,Ncandidates

        ! choose the x-position assuming a probability proportional to r
        call rand_mt(temp)
        x=sqrt(c%min_r**2+temp*delr2)
        y=0.d0

        r = x
        call rand_mt(temp)
        z = zstart  + temp*(zend-zstart)

! find the desired local density at the candidate particle location
        idx = (z-zmin_global)/s_cells_dz + 1
        idx = max( idx, 1 )
        idx = min( idx, num_s_cells_z )

        n = P_top(idx)/kb/T_top(idx)

! now do acceptance-rejection density and on the load window
        call rand_mt(temp)
        if(nmax*temp .lt. n .and. idx.gt.jstart-1 .and. idx.lt.jend+1)then
            Nloaded = Nloaded+1

! create the particle

! calculate vth at this location
            vth=sqrt(kb*T_top(idx)/m)

            p0%prev_coll = -1 ! no collisions have happened yet

! create the velocity of the particle
! use a fixed speed and a spherical expansion velocity
            call rand_mt(temp)
            vperp = vth*sqrt(-2.d0*log(temp))
            call rand_mt(temp)
            theta = twopi*temp
! chosen    p0%vx =  vperp*cos(theta) + x/sqrt(x**2+z**2)*ufix
            p0%vx =  vperp*cos(theta) + vr_top(idx)
            p0%vy = vperp*sin(theta)

            call rand_mt(temp)
            vperp = vth*sqrt(-2.d0*log(temp))
            call rand_mt(temp)
            theta = twopi*temp
! chosen    p0%vz = vperp*sin(theta) + z/sqrt(x**2+z**2)*ufix
            p0%vz = vperp*sin(theta) + vz_top(idx)

! create its position
            p0%x = x
            p0%y = 0.d0
            p0%z = z
            p0%r = x

            p0%flag=99 ! flag it as a particle born as a ghost particle

            p0%element=0  ! argon particle

! put the particle where it should go


            call cell_insert_particle(c, p0)

! we set the input trace element density to the input argon density
! to get good statistics. So here in the boundary we can just insert
! a trace particle every time we insert an argon particle and we will
! keep the densities proportional. If a different density is needed,
! use acceptance rejection here

! turn off the outside load and just let the nozzle populate the simulation
!           p0%element=1 ! trace particle
!           call cell_insert_particle(c, p0)


!             write(311,101) mpi_rank,p0%x,p0%z,p0%vx,p0%vy,p0%vz
        end if
  101   format(i5,5(1x,1pe12.4))


    end do

! kluge
!    write(*,102) nmax, Ncandidates,Nloaded,mpi_rank,c%volume
102  format(' nmax, Ncandidates, Nloaded, mpi_rank: ',1(1x,1pe12.4),3(2x,i12),2x,1pe12.4)


  end subroutine do_top_bc






  !-- do_right_bc ---
  ! this subroutine repopulates a right-type ghost cell
  subroutine do_right_bc(c)

    use core_data_mod
    use core_types_mod
    use helpers_mod
    use constants_mod
    use simulation_mod
    use parallel_info_mod


    implicit none
    type(cell_type) c
    integer Ncandidates, i, Nloaded
    real(8) :: n, T, x, y, z, r, temp, nmax, delr2, delz, fracandidates=0.d0, rN
    real(8) rstart,rend,margin,vztmp,vrtmp,Ttmp,ntmp,vth
    real(8) vperp,theta,twopi,Vlocal,rmaxlocal
! note that rmin_local and rmax_local are the min and max r-values in the
! physical region (no ghost cell extensions) of the processor

    type(particle_type) p0
    integer idx

    data twopi/6.283185307d0/




!   write(*,107) c%min_z,c%max_z,c%min_r,c%max_r,mpi_rank,' %%'
107 format('zstart,zend, rmin, rend, mpi_rank',4(1pe12.4),i5,a5)

    ! just do vacuum
!    type(cell_type) c
!    c%num_parts = 0
!    c%partition = 0

!    return


    ! This routine adds particles at the right edge (maximum z)
    ! Note that this is called by every processor that has
    ! a right-edge ghost cell, so you have to use cell variables,
    ! i.e., c%xxx variables

    rstart=c%min_r
    rend=c%max_r


! find the maximum density along this cell's portion of the top
    istart = (rstart-rmin_global)/s_cells_dr + 1
    iend   = (rend-rmin_global)/s_cells_dr


    nmax=0.d0
    do idx=istart,iend
       nmax = max( nmax,P_end(idx)/kb/T_end(idx) )
    end do

    nmax=nmax*2.d0 ! just make sure it's big enough

    nmax=4.d22


! compute the number of candidates to use
! don't use the full cell volume because we don't want
! to load the upper right-hand corner twice.



!    if(c%max_r.gt.rmax_global) then
!      rmaxlocal=rmax_global
!      Vlocal=pi*(rmaxlocal**2-c%min_r**2)*delz
!    end if
!

     delz =c%max_z - c%min_z  ! the z-width of the ghost region
     delr2=rend**2-rstart**2
     Vlocal=pi*(rend**2-rstart**2)*delz

    rN = Vlocal * nmax / real(Nef)
    Ncandidates = rN + fracandidates
    fracandidates = fracandidates + rN-Ncandidates


    Nloaded=0

    do i=1,Ncandidates

        ! choose the z-position assuming a uniform density
        call rand_mt(temp)
        z=c%min_z + temp*delz

        y=0.d0

        ! choose the r-position assuming a uniform density
        call rand_mt(temp)
        x=sqrt(c%min_r**2+temp*delr2)

        r = x

! find the desired local density at the candidate particle location
        idx = r/s_cells_dr + 1
        idx = max( idx, 1 )
        idx = min( idx, num_s_cells_r )

        n = P_end(idx)/kb/T_end(idx)

! now do acceptance-rejection density and on the load window
        call rand_mt(temp)
        if(nmax*temp .lt. n .and. idx.gt.istart-1 .and. idx.lt.iend+1)then
            Nloaded = Nloaded+1

! create the particle


!       write(153,115) x,y,z,idx
115      format(' %%',4(1pe12.4))

! calculate vth at this location
            vth=sqrt(kb*T_end(idx)/m)

            p0%prev_coll = -1 ! no collisions have happened yet

! create the velocity of the particle
! fixed boundary speed at ufix
            call rand_mt(temp)
            vperp = vth*sqrt(-2.d0*log(temp))
            call rand_mt(temp)
            theta = twopi*temp
! chosen    p0%vx = vperp*cos(theta) + x/sqrt(x**2+z**2)*ufix
            p0%vx = vperp*cos(theta) + vr_end(idx)
            p0%vy = vperp*sin(theta)

            call rand_mt(temp)
            vperp = vth*sqrt(-2.d0*log(temp))
            call rand_mt(temp)
            theta = twopi*temp
! chosen    p0%vz =  vperp*sin(theta) + z/sqrt(x**2+z**2)*ufix
            p0%vz =  vperp*sin(theta) + vz_end(idx)

! create its position
            p0%x = x
            p0%y = 0.d0
            p0%z = z
            p0%r = x

! look at the particles being loaded
!           write(*,"( 6(1x,1pe12.4) )") x,z,p0%vx,p0%vy,p0%vz

            p0%flag=99 ! flag it as a particle born as a ghost particle

            p0%element=0  ! argon particle

! put the particle where it should go
            call cell_insert_particle(c, p0)

! we set the input trace element density to the input argon density
! to get good statistics. So here in the boundary we can just insert
! a trace particle every time we insert an argon particle and we will
! keep the densities proportional. If a different density is needed,
! use acceptance rejection here

!           p0%element=1 ! trace particle
!           call cell_insert_particle(c, p0)

!             write(311,101) mpi_rank,p0%x,p0%z,p0%vx,p0%vy,p0%vz
        end if
  101   format(i5,5(1x,1pe12.4))


    end do

! kluge
!     write(*,102) nmax, Ncandidates,Nloaded,mpi_rank,c%volume
102  format(' nmax, Ncandidates, Nloaded, mpi_rank: ',1(1x,1pe12.4),3(2x,i12),2x,1pe12.4)


     return

  end subroutine do_right_bc





!**********Electric Field Definitions:

    ! -- ElectricField_R --
    ! returns the magnitude of the electric field at a given (r,z) point
    real(8) function ElectricField_R( r, z )
        use constants_mod
        implicit none
        real(8) r
        real(8) z

        !ElectricField_R = -2.0*r*kb*3000.0/((0.001)**2)/(tr_q)  ! 1.6e-19 Coulombs, 3000 K, 0.001 meters.
        ! to make it faster, I just write it this way:
        !ElectricField_R = -517171.498076*r;

		ElectricField_R = 0.0;

    end function ElectricField_R


    ! -- ElectricField_Z --
    ! returns the magnitude of the electric field at a given (r,z) point
    real(8) function ElectricField_Z( r, z )
        implicit none
        real(8) r
        real(8) z

        ElectricField_Z = 0.0;

    end function ElectricField_Z




end module bounds_mod



