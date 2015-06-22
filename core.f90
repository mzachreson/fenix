



!----- constants_mod -----
! This fortran module defines all of the constants
module constants_mod
implicit none
  ! true constants
  real(8) pi, e0
  real(8) kb !boltzmann constant, in SI units

  !particle characteristics
  real(8) argon_diam        !"argon_diam" variable is NEVER USED!
  real(8) m !mass of argon atom, in kg

  real(8) tr_m
  real(8) tr_q    ! charge of trace particle -- in Coulombs

  real(8) aref, aref_tr, bref_tr, nu, nu_tr1,nu_tr2
  real(8) particle_code, trace_code !Variable for saving the name of the trace element.
                                    !Instructions for setting it are givine in initialize.f90
                                    !where if is filled.
 
  parameter (pi=3.14159265358979d0)
  parameter (argon_diam=142.d-12)
  parameter (m=6.634d-26)
  parameter (kb=1.381d-23)
  parameter (e0=8.854d-12)

  ! "aref" is used in "collide.f90"
  parameter(aref    = 7.925e-18)

  ! "nu" is used in "collide.f90"
  parameter(nu    = 0.22)


  !Now choose the trace element
! data trace_switch/2/ !1 - Barium; 2 - Calcium; 3 - Argon Ions
  !Trace collision data is loaded in initialize.f90, subroutine
  !initialize_trace_collide_data, go there to edit values
! data nanbu_switch/0/ !0 turns off nanbu collisions, 1 turns them on
                       !Nanbu only works if fluid_switch is also set to 1.
! data fluid_switch/1/ !switch that tells FENIX to track fluid properties
                       !(dens,temp,vr,vx,etc) these properties are stored in the
                       !cell type.


!Set trace
end module constants_mod





!----- simulation_mod -----
! This module stores general information about the simulation
module simulation_mod
  implicit none



  !the number of steps in the whole simulation
  integer nsteps

  !the current step in the simulation 
  integer cur_step

  ! tau is the timestep
  real(8) tau
  logical use_restart, single_proc_mode
  integer output_skip, restart_skip
  integer output_unit
  integer outputtrace_unit
  integer restart_unit
  integer tr_restart_unit
  integer adapt_unit
  integer partcount_unit

  !A few property switches
  integer trace_switch  ! Switch for changing the properties of the trace ion
  logical nanbu_switch !Switch to turn nanbu style trace collisions on and off
  logical fluid_switch !Switch for tracking fluid properties in each collision cell

  !This is the guesstimated temperature of the simulation
  !Used to figure out the size of ghost cells, etc. in initialization
  real(8) simT
  parameter(simT=5000.d0)
 
  ! N-effective
  ! or the number of real particles each simulated particle represents
  real(8) Nef
  ! These window parameters ignore the ghost cells
  real(8) rmin_global,rmax_global
  real(8) zmin_global,zmax_global
  ! These window parameters are per-processor, if applicable
  real(8) rmin_local,rmax_local
  real(8) zmin_local,zmax_local


  parameter (output_unit=250)
  parameter (outputtrace_unit=381)
  parameter (restart_unit=111)
  parameter (tr_restart_unit=628)
  parameter (adapt_unit=778)
  parameter (partcount_unit=253)

end module simulation_mod




!----- core_types_mod -----
! This fortran module defines all of the derived types needed by the core
module core_types_mod
implicit none
  
  type particle_type
    real(8) :: x,y,z
    real(8) :: r ! this variable is pre-computed for convenience; it is sqrt(x^2+y^2).  
                 ! Note that we need it square-rooted for the datacollect anyway, so we might 
                 ! as well just square root it right away
    real(8) :: vx,vy,vz
    integer :: prev_coll
	integer :: element ! 0=Argon, 1=TraceElement
    integer :: flag    ! This variable is so we can mark a particle and find it later
  end type particle_type
  
  type cell_type
    type(particle_type),pointer,dimension(:) :: ps
    integer :: num_parts,capacity,num_left
	integer :: partition ! "partition" points to the last normal argon particle in the array.
						 ! it partitions the normal particles from the trace particles.
						 ! From 1 to partition are normal particles, and from partition+1 to num_parts are trace particles.
						 ! if partition == num_parts, then there are no trace particles.
    
    real(8) min_r,max_r
    real(8) min_z,max_z
    real(8) min_th,max_th
    real(8) max_vrel2, tr_max_vrel2
    real(8) volume
    real(8) :: efield_r=0
    real(8) :: efield_z=0
    real(8) :: Te=0
    real(8) :: extra=0
    real(8) :: num_cand_remaining=0
    real(8) dens(50), vr_fluid(50), vz_fluid(50), Temp(50)
    real(8) densavg, vr_fluidavg, vz_fluidavg, tempavg
    real(8) :: aref_ambi
    integer mask
  end type cell_type

  type sampling_cell_type

    ! argon particles data
    real(8) sum_n
    real(8) sum_vx,sum_vy,sum_vz
    real(8) sum_vx2,sum_vy2,sum_vz2
    real(8) sum_vxvy,sum_vxvz,sum_vyvz

    ! trace particles data
    real(8) tr_sum_n
    real(8) tr_sum_vx, tr_sum_vy, tr_sum_vz
    real(8) tr_sum_vx2, tr_sum_vy2, tr_sum_vz2
    real(8) tr_sum_vxvy, tr_sum_vxvz, tr_sum_vyvz

    ! volume of the cell
	real(8) volume
  end type sampling_cell_type

  type mask_type
    integer code
    real(8) reflect !a sort of flag that says what kind of reflection to do
    real(8) s0z, s0r
    real(8) w1z,w1r,g1z,g1r ! w--direction of the metal, g--inward normal
    real(8) w2z,w2r,g2z,g2r
    real(8) T
  end type mask_type
  
  
  type point_type
    real(8) r, z
  end type point_type
  
  
  type polygon_type
    type(point_type), pointer, dimension(:) :: pts
    integer num_pts
    integer type_poly
    real(8) reflect_coeff, temperature
  end type polygon_type
  
  
end module core_types_mod






!----- parallel_info_mod -----
! This module stores information about the parallelization and has a few helper functions
! to know how the current process lines up with the sampling cells (i goes with r, j with z)
module parallel_info_mod
  use core_types_mod
  implicit none
  integer mpi_rank
  integer num_mpi_processes
  integer, allocatable, dimension(:) :: comp_reg_r
  integer, allocatable, dimension(:) :: comp_reg_z
  integer num_proc_r
  integer num_proc_z
  
  integer mpi_ierr
  type(cell_type) comm_cell_left
  type(cell_type) comm_cell_right
  type(cell_type) comm_cell_top
  type(cell_type) comm_cell_bottom
  character,allocatable,dimension(:) :: comm_buffer
  
contains

  ! Returns value specific to this mpi_rank
  integer function min_i()
    implicit none
    min_i = min_ip(mpi_rank)
  end function min_i
  
  integer function max_i()
    implicit none
!  write(*,*) 'function max_i input', mpi_rank
    max_i = max_ip(mpi_rank)
!  write(*,*) 'function max_i output', max_i
  end function max_i
  
  integer function min_j()
    implicit none
    min_j = min_jp(mpi_rank)
  end function min_j
  
  integer function max_j()
    implicit none
    max_j = max_jp(mpi_rank)
  end function max_j
  

  ! Returns value general to whatever mpi_rank
  integer function min_ip(proc)
    implicit none
    integer proc
  ! write(*,*) 'function min_ip input',proc,num_proc_z
    min_ip = comp_reg_r(proc/num_proc_z+1)+1
  ! write(*,*) 'function min_ip output',min_ip
  end function min_ip
  
  integer function max_ip(proc)
    implicit none
    integer proc
    max_ip = comp_reg_r(proc/num_proc_z+2)
  end function max_ip

  integer function min_jp(proc)
    implicit none
    integer proc
    min_jp = comp_reg_z(mod(proc,num_proc_z)+1)+1
  end function min_jp
  
  integer function max_jp(proc)
    implicit none
    integer proc
    max_jp = comp_reg_z(mod(proc,num_proc_z)+2)
  end function max_jp
  
end module parallel_info_mod






!----- core_data_mod -----
! This fortran module stores the global arrays belonging to the core module
! Note: Contains data formerly contained in core_data_mod, cells_mod, and polygons_mod
module core_data_mod
use core_types_mod
implicit none


  ! Masks
  type(mask_type),allocatable,dimension(:) :: masks
  ! max_mask and num_masks aren't necessarily the same. Example:  Masks: 1,2,3,5.  max_mask = 5.  num_masks = 4.
  integer max_mask, min_mask   
  integer num_masks
  
  ! Sampling Cells
  type(sampling_cell_type),allocatable,dimension(:,:) :: s_cells
  integer num_s_cells_r
  integer num_s_cells_z
  real(8) s_cells_dr
  real(8) s_cells_dz
  real(8) num_samples

  ! Global total number of collision cells
  integer num_cells_global

  ! This processor's collision cells
  type(cell_type), allocatable, dimension(:) :: cells
  integer num_cells
  integer starting_capacity
  parameter(starting_capacity=10)
  
  ! This processor's fine_cells array
  ! the first index is the r index, the second the z index
  integer, allocatable, dimension(:,:) :: fine_cells
  real(8) fine_cells_dr, fine_cells_dz
  integer num_fine_cells_r
  integer num_fine_cells_z

  ! Polygons
  type(polygon_type), allocatable, dimension(:) :: polys
  integer num_polys


end module core_data_mod





!----- helpers_mod -----
! This fortran module defines all of the low-level helper functions
module helpers_mod
implicit none
contains



  !--- is_inside ---
  ! This method determines whether a given position (r,z) is 'inside' a mask m
  logical function is_inside(r,z,m)
    use core_types_mod
    implicit none
    real(8) r,z
    type(mask_type)m
    real(8) pr, pz, inside1, inside2

    pr = r - m%s0r
    pz = z - m%s0z
    inside1 = pr*m%g1r+pz*m%g1z
    if(m%g2z .eq. 0 .and. m%g2r .eq. 0)then
      is_inside = inside1 .gt. 0.d0
      return
    end if
    inside2 = pr*m%g2r+pz*m%g2z
    if((m%w1r+m%w2r)*(m%g1r+m%g2r)+(m%w1z+m%w2z)*(m%g1z+m%g2z) .gt. 0)then
      is_inside = inside1 .gt. 0.d0 .and. inside2 .gt. 0.d0
    else
      is_inside = inside1 .gt. 0.d0 .or. inside2 .gt. 0.d0
    end if
    return
  end function is_inside

  
  !--- isnan ---
  ! returns true if the given number is NaN
  logical function isnan(num)
    implicit none
    real(8) num
    isnan = num .ne. num
    return
  end function isnan
  


  !--- maxwell_boltzmann ---
  ! This function creates a random velocity (just a single component) according to a particular temperature
  real(8) function maxwell_boltzmann(T)
    use constants_mod
    implicit none
    integer num_save
    parameter(num_save=256)
    real(8) T, temp, r, th
    integer i
    real(8),save :: norm_deviates(num_save)
    integer,save :: cur_dev

    if(cur_dev .eq. 0)then
      do i=1,num_save/2
        ! Box-Muller for random normal deviates
        call rand_mt(temp)
        r = sqrt(-2.d0*log(temp))
        call rand_mt(temp)
        th = 2.d0*pi*temp
        norm_deviates(2*i) = r*cos(th)
        norm_deviates(2*i-1) = r*sin(th)
      end do
      cur_dev = num_save
    end if
    ! Now the distribution of a single component of velocity is normally distributed,
    !    with variance k*T/m and mean 0 (this is the maxwell-boltzmann distribution)
    ! So to produce v, I just say v = sqrt(k*T/m) * N(0,1)
    !    where N(0,1) is a random, normally distributed 
    maxwell_boltzmann = sqrt(kb*T/m) * norm_deviates(cur_dev)
    cur_dev = cur_dev - 1
    return
  end function maxwell_boltzmann
  
  
  
  
  
  
  
  !--- create_particle ---
  ! This routine creates a particle at the given position
  ! it randomly assigns the particle's velocities according
  ! to the distribution defined by the mean vx, vy, and vz,
  ! and the temperature
  subroutine create_particle(p, T, vx, vy, vz, x, y, z)
    use core_types_mod
    implicit none
    type(particle_type) p
    real(8) T
    real(8) vx,vy,vz
    real(8) x,y,z

    p%prev_coll = -1 ! no collisions have happened yet

    ! create the velocity of the particle
    p%vx = vx + maxwell_boltzmann(T)
    p%vy = vy + maxwell_boltzmann(T)
    p%vz = vz + maxwell_boltzmann(T)

    ! put it where it belongs
    p%x = x
    p%y = y
    p%z = z
    p%r = sqrt(p%x**2 + p%y**2)

  end subroutine create_particle
  
  
  
  
  
  !--- cell_insert_particle ---
  ! Inserts a particle into a cell, growing the cell's capacity if needed
  subroutine cell_insert_particle(c, p)
    use core_types_mod
    use core_data_mod
    implicit none
    type(cell_type) c
    type(particle_type) p
    type(particle_type) temp
    type(particle_type),allocatable,dimension(:) :: new_parts
    integer new_capacity,i

    if(c%mask .ne. 0)then
      if(is_inside(p%r, p%z, masks(c%mask))) return
    end if

    if(c%capacity .le. c%num_parts+c%num_left)then
      !this is the case where there is not yet enough space
      new_capacity = int(1.5*c%capacity+1)
      allocate(new_parts(c%capacity))
      do i=1,c%capacity
        new_parts(i) = c%ps(i)
      end do
      deallocate(c%ps)
      c%ps => null()
      allocate(c%ps(new_capacity))
      !I've got to keep the right justification around
      c%ps(1:c%num_parts) = new_parts(1:c%num_parts)
      do i=c%capacity-c%num_left+1,c%capacity
        c%ps(i+new_capacity-c%capacity) = new_parts(i)
      end do
      c%capacity=new_capacity
      deallocate(new_parts)
    end if

    !now I know there is enough space
	if( p%element == 0 ) then
		! normal particle
		c%ps(c%num_parts+1) = c%ps(c%partition+1)
    	c%ps(c%partition+1) = p
    	c%num_parts = c%num_parts + 1
		c%partition = c%partition + 1
	else
		! trace particle
		c%ps(c%num_parts+1)=p
    	c%num_parts = c%num_parts + 1
	end if

  end subroutine cell_insert_particle
  


  !--- fill_flat_distrib ---
  ! fills a cell with a normal distribution, for use when beginning without a restart file.
  ! USED FOR DEBUGGING ONLY
  subroutine fill_flat_distrib( cl )
    use core_types_mod
    use simulation_mod
    implicit none

        type(cell_type) :: cl
    	integer :: num_cand_to_insert, i
    	type(particle_type) :: p0
    	real(8) :: temp  	
    	real(8) :: xx, yy, zz, rr
        real(8) :: delr2
    	
        ! this is just a guess at a number density.
        ! The programmer will have to tweak this if this isn't right when testing
        real(8) :: guess_number_density = 2.723d22

    	num_cand_to_insert = round_rand(guess_number_density * cl%volume / (Nef*100.0))
        delr2 = cl%max_r**2 - cl%min_r**2

		do i = 1,num_cand_to_insert
		
			! assigning the particle's position

            call rand_mt(temp)
            xx=sqrt(cl%min_r**2+temp*delr2)
            yy=0.d0

			call rand_mt(temp)
			zz = temp*(cl%max_z - cl%min_z) + cl%min_z
			rr = xx

            ! creating the particle
			call create_particle( p0, simT, 0.d0, 0.d0, 0.d0, xx, yy, zz )
            p0%element = 0
			
			! inserting the particle
			call cell_insert_particle(cl, p0)
				
		end do
   
  end subroutine fill_flat_distrib



  !--- rand_mt ---
  ! stores a random number in [0,1) in 'output'
  subroutine rand_mt(output)
    implicit none
    integer size_arr
    parameter(size_arr=(4*1024))
    real(8) output
    real(8),save,dimension(size_arr) :: saved_randoms
    !DEC$ ATTRIBUTES ALIGN:16 :: saved_randoms
    integer,save :: current_place

    if(current_place .eq. 0)then
      !call fill_array_close_open(saved_randoms(1), size_arr)
      call fill_array_close_open(saved_randoms, size_arr)
      current_place = size_arr
    end if
    output = 0.999999*saved_randoms(current_place)
    current_place = current_place - 1
  end subroutine rand_mt


  
  !--- fill_array_close_open---
  subroutine fill_array_close_open(arr, num_arr)
    implicit none
    real(8), dimension(:) :: arr
    integer num_arr, i

    do i=1,num_arr
      call random_number(arr(i))
    end do
  end subroutine fill_array_close_open



  !--- round_rand ---
  ! This algorithm just rounds either up or down randomly,
  ! in such a way as the mean of this function is mu.
  ! so if mu is 7.1, 10% of the time it rounds up and
  ! 90% of the time it rounds down.
  integer function round_rand(mu)
    implicit none
    real(8) mu
    integer int_val
    real(8) diff, temp

    int_val = int(mu)
    diff = mu - int_val
    call rand_mt(temp)
    if(temp .le. diff) then
      round_rand = int_val+1
    else
      round_rand = int_val
    end if
    return
  end function round_rand  





  ! "Sort" functions and Subroutines
  recursive subroutine sortZ(arr, n1, n2)
  	use core_types_mod
    implicit none
    type(particle_type),dimension(:) :: arr
    integer n1, n2
    integer k

    ! This is a mergesort, which terminates early with insertion sort
    if(n2-n1 .gt. 20)then
      k = (n1+n2)/2
      call sortZ(arr, n1, k)
      call sortZ(arr, k+1, n2)
      call mergeZ(arr, n1, n2, k)
    else if(n2-n1 .gt. 10)then
      k = (n1+n2)/2
      call isortZ(arr, n1, k)
      call isortZ(arr, k+1, n2)
      call mergeZ(arr, n1, n2, k)
    else
      call isortZ(arr, n1, n2)
    end if
  end subroutine sortZ

  
  
  

  

  ! insertion sort
  subroutine isortZ(arr, n1, n2)
  	use core_types_mod
    implicit none
    type(particle_type),dimension(:) :: arr
    integer n1, n2
    integer i,j
    logical go_back
    type(particle_type) temp

    

    do i=n1+1,n2
      temp = arr(i)
      j = i-1
      go_back = temp%z .lt. arr(j)%z
      do while(go_back)
        arr(j+1) = arr(j)
        j = j-1
        if(j .lt. n1)then
          go_back = .false.
        else
          go_back = temp%z .lt. arr(j)%z
        end if
      end do
      arr(j+1) = temp
    end do
  end subroutine isortZ

  

  ! merge
  subroutine mergeZ(arr, n1, n2, k)
  	use core_types_mod
    implicit none
    type(particle_type),dimension(:) :: arr
    integer n1, n2, k
    integer i,j,n
    type(particle_type) temparr(n1:n2)
    temparr(n1:n2) = arr(n1:n2)

    n = n1
    i = n1
    j = k+1
    !write(*,*)"n1,n2=",n1,n2
    !write(*,*)"k=",k
    do while(i .lt. k+1 .and. j .lt. n2+1)
      !write(*,*)"i,j=",i,j
      !write(*,*)"Comparing ", temparr(i)%z, temparr(j)%z
      if(temparr(i)%z .lt. temparr(j)%z)then
        arr(n) = temparr(i)
        !write(*,*)"Inserted ",temparr(i)%z
        i = i + 1
        n = n + 1
      else
        arr(n) = temparr(j)
        !write(*,*)"Inserted ",temparr(j)%z
        j = j + 1
        n = n + 1
      end if
    end do

    do while(i .lt. k+1)
      arr(n) = temparr(i)
      !write(*,*)"End inserted ",temparr(i)%z
      n = n + 1
      i = i + 1
    end do

    do while(j .lt. n2+1)
      arr(n) = temparr(j)
      !write(*,*)"End inserted ",temparr(j)%z
      n = n + 1
      j = j + 1
    end do

    

    if(n .ne. n2+1)then
      write(*,*)"Sorted wrong.", n, n2
      write(*,*)i, k
      write(*,*)j
      write(*,*)n1
      stop 101
    end if
  end subroutine mergeZ

  
  
  ! tests the 'sort' functions
  subroutine testsort
  	use core_types_mod

    implicit none
    type(particle_type),dimension(100) :: backup_arr
    type(particle_type),dimension(100) :: arr
    integer i,j
    real(8) temp

    do i=5,100
      ! create the arrays
      do j=1,i
        call rand_mt(temp)
        arr(j)%z = temp
        backup_arr(j)%z = temp
      end do
      ! sort them both, differently
      call isortZ(backup_arr, 1, i)
      call sortZ(arr, 1, i)
      ! check that backup_arr is nondecreasing (insertion sort test)
      do j=2,i
        if(backup_arr(j-1)%z .gt. backup_arr(j)%z)then
          write(*,*) "Bad insertion sorting!"
          write(*,*) "barr: ",backup_arr(:)%z
          stop 667
        end if
      end do

      ! check that they are the same
      do j=1,i
        if(arr(j)%z .ne. backup_arr(j)%z)then
          write(*,*) "Bad sorting!"
          write(*,*) "arr: ",arr(1:i)%z
          write(*,*) "barr: ",backup_arr(1:i)%z
          stop 668
        end if
      end do
    end do
    write(*,*) "Testsort complete."
  end subroutine testsort
  
end module helpers_mod









