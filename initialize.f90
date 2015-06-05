


module initialize_mod
  implicit none
	
  contains
  
  subroutine initialize_all()
    use bounds_mod
    use core_data_mod
    use simulation_mod
  	use helpers_mod
  	use parallel_info_mod
  	implicit none

    call basic_initialization

    call read_geometry
    call initialize_trace_collide_data  ! call before initialize_collision_cells so that tr_m is set

    call initialize_collision_cells
    call initialize_fine_cells
    call initialize_polys
    call initialize_masks
    
    call initialize_boundaries
    call initialize_output
    call initialize_communication
    call initialize_particles
 

    open(unit=partcount_unit,file='particle_count.txt',status='unknown')

    !write(*,*)mpi_rank,'Done initialization.'
    
  end subroutine initialize_all





  !----- basic_initialization
  ! this subroutine reads "input_data.txt" and sets up the MPI environment
  subroutine basic_initialization
    use parallel_info_mod
    use simulation_mod
    use core_data_mod
    implicit none
    integer in_u, i
	real(8) vth_no_need
	real(8) rmin_no_need, rmax_no_need, zmin_no_need, zmax_no_need
    parameter(in_u=200)
    ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'


    ! INITIALIZING MPI
    ! before anything else happens
    !*************************
    call MPI_INIT(mpi_ierr)
    if(mpi_ierr .ne. MPI_SUCCESS)  call MPI_ABORT(MPI_COMM_WORLD,0,mpi_ierr)

    ! find out how many MPI processes there are
    call MPI_COMM_SIZE(MPI_COMM_WORLD,num_mpi_processes,mpi_ierr)
    if(mpi_ierr .ne. MPI_SUCCESS)  call MPI_ABORT(MPI_COMM_WORLD,0,mpi_ierr)

    ! find out what rank this MPI process is
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_rank,mpi_ierr)
    if(mpi_ierr .ne. MPI_SUCCESS)  call MPI_ABORT(MPI_COMM_WORLD,0,mpi_ierr)
 
    ! Introductory WRITING
    if(mpi_rank .eq. 0)then
        write(*,*) '------- Fenix 10 (parallel version) -------'
        write(*,*) '-------- DSMC Simulation Software ---------'
        write(*,*) '-------- Brigham Young University ---------'
        write(*,*) ''
    end if
    ! READING "input_data.txt"
    !*************************
    if(mpi_rank == 0) write(*,*) ""
    if(mpi_rank == 0) write(*,*) "Reading input_data.txt:"
    if(mpi_rank == 0) write(*,*) '-------------------------------------------'
    if(mpi_rank == 0) write(*,*) ""

    open(unit=in_u,file='input_data.txt',status='old',action='read')
  
    read(in_u,*) Nef, vth_no_need
    if(mpi_rank == 0) then
       write(*,*) 'Nef, vth_no_need'
       write(*,*)  Nef, vth_no_need
    end if

    read(in_u,*) nsteps, tau
    if(mpi_rank == 0) then
       write(*,*) 'nsteps, tau'
       write(*,*) nsteps, tau
    end if

    read(in_u,*) use_restart, single_proc_mode
    if(mpi_rank == 0) then
       write(*,*) 'use_restart, single_proc_mode'
       write(*,*) use_restart, single_proc_mode
    end if

    read(in_u,*) output_skip, restart_skip
    if(mpi_rank == 0) then
       write(*,*) 'output_skip, restart_skip'
       write(*,*) output_skip, restart_skip
    end if

    read(in_u,*) fluid_switch, nanbu_switch
    if(mpi_rank == 0) then
       write(*,*) 'fluid_switch, nanbu_switch'
       write(*,*) fluid_switch, nanbu_switch
    end if

    read(in_u,*) trace_switch
    if(mpi_rank == 0) then
       write(*,*) 'trace_switch'
       write(*,*) trace_switch
    end if

    read(in_u,*) fine_cells_dr, fine_cells_dz
    if(mpi_rank == 0) then
       write(*,*) 'fine_cells_dr, fine_cells_dz'
       write(*,*) fine_cells_dr, fine_cells_dz
    end if

    read(in_u,*) s_cells_dr, s_cells_dz
    if(mpi_rank == 0) then
       write(*,*) 's_cells_dr, s_cells_dz'
       write(*,*) s_cells_dr, s_cells_dz
    end if


	! convert these to units of meters
	s_cells_dr = s_cells_dr * fine_cells_dr
	s_cells_dz = s_cells_dz * fine_cells_dz

	! we don't need this data here, but read it so we can get the data we need below it.
	read(in_u,*) rmin_no_need, rmax_no_need, zmin_no_need, zmax_no_need

	! MPI multi-processor information
    read(in_u,*) num_proc_r, num_proc_z
    if(mpi_rank == 0) then
       write(*,*) 'num_proc_r, num_proc_z'
       write(*,*) num_proc_r, num_proc_z
    end if


    ! single processor if block
    if( single_proc_mode ) then
        ! allocate the compute region arrays
        allocate(comp_reg_r(2))
        allocate(comp_reg_z(2))

        ! fill them in -- first one set to zero is correct. Each number in the array
        ! is the TOP cell of the section.  That makes zero the first one, being the "top" cell
        ! of the one just below the first cell.
        comp_reg_r(1) = 0
        comp_reg_z(1) = 0

        ! in single_proc_mode, we only care about the highest r boundary
        ! so we read though all the lines to the last one, and keep the last one.
        do i=2,num_proc_r+1
            ! only the last one is saved, in single_proc_mode
            read(in_u, *) comp_reg_r(2)
        end do

        ! in single_proc_mode, we only care about the highest z boundary
        ! so we read though all the lines to the last one, and keep the last one.
        do i=2,num_proc_z+1
            read(in_u, *) comp_reg_z(2)
        end do

        ! we are in single processor mode, so these equal 1
        num_proc_r = 1
        num_proc_z = 1

        ! in single processor mode, we don't want any more than about a couple million
        ! particles (or it would be too much for a processor to handle)
        ! to prevent that, all we can do is increase "Nef" by a factor of 10.
        Nef = Nef * 10

    ! end single processor if block


    ! multi processor if block
    else
        ! allocate the compute region arrays
        allocate(comp_reg_r(num_proc_r+1))
        allocate(comp_reg_z(num_proc_z+1))

        ! fill them in -- first one set to zero is correct. Each number in the array
        ! is the TOP cell of the section.  That makes zero the first one, being the "top" cell
        ! of the one just below the first cell.
        comp_reg_r(1) = 0
        comp_reg_z(1) = 0
    
        do i=2,num_proc_r+1

            read(in_u, *) comp_reg_r(i)
            if(mpi_rank == 0) then
               write(*, *) 'comp_reg_r(',i,')'
               write(*, *) comp_reg_r(i)
            end if

        end do

        do i=2,num_proc_z+1
            read(in_u, *) comp_reg_z(i)
            if(mpi_rank == 0) then
               write(*, *) 'comp_reg_z(',i,')'
               write(*, *) comp_reg_z(i)
            end if
        end do
    ! end multi processor if block
    end if

    ! closing the input file
    close(unit=in_u)



    ! PERFORMING SOME MPI CHECKS
    ! *************************
 
    ! Checking the MPI data against the input_data.txt info
    ! to make sure they are consistent.
    if(num_proc_r*num_proc_z .gt. num_mpi_processes) then
      write(*,*) 'run_input.txt provides for ', num_proc_r*num_proc_z, ' processes,'
      write(*,*) 'but the MPI environment only has ', num_mpi_processes, 'processes.'
      stop 669
    end if
    if(mpi_rank .ge. num_proc_r*num_proc_z)then
      write(*,*)'Too many processes.'
      stop 670
      !do while(.true.)
      !end do
    end if

  end subroutine basic_initialization


 



  !----- read_geometry -----
  ! This subroutine reads the geometry files, and creates the fine cells array
  subroutine read_geometry
    use constants_mod
    use core_data_mod
    use parallel_info_mod
    use simulation_mod
    implicit none
    integer param_u
    parameter(param_u=201)
    integer num_fine_cells_r_global
    integer num_fine_cells_z_global

  
    ! read the parameters file
    open(unit=param_u,file='geom/parameters.dat',status='old',action='read')
    read(param_u,*) num_cells_global
    read(param_u,*) num_fine_cells_r_global
    read(param_u,*) num_fine_cells_z_global
    read(param_u,*) num_polys
    read(param_u,*) num_masks
    read(param_u,*) fine_cells_dr
    read(param_u,*) fine_cells_dz
    rmin_global = 0
    read(param_u,*) rmax_global
    read(param_u,*) zmin_global
    read(param_u,*) zmax_global
    close(unit=param_u)

    ! Nef is now read in from "input_data.txt"
    ! no need to initialize it here.
  
    ! initialize the local window bounds
    rmin_local = rmin_global + (min_i()-1)*s_cells_dr
    rmax_local = rmin_global + max_i()*s_cells_dr
    zmin_local = zmin_global + (min_j()-1)*s_cells_dz
    zmax_local = zmin_global + max_j()*s_cells_dz
  
    ! initialize the number of sampling cells in each direction
    ! 'nint' rounds to the nearest integer
    num_s_cells_r = nint((rmax_global-rmin_global)/s_cells_dr)
    num_s_cells_z = nint((zmax_global-zmin_global)/s_cells_dz)
    ! check this
    if(comp_reg_r(num_proc_r+1) .ne. num_s_cells_r) then
      write(*,*) 'There are ', num_s_cells_r, ' sampling cells in r, but'
      write(*,*) comp_reg_r(num_proc_r+1), ' are assigned to a processor.'
      write(*,*) 'num_proc_r=',num_proc_r
      write(*,*) 'comp_reg_r=',comp_reg_r
      stop 671
    end if
    if(comp_reg_z(num_proc_z+1) .ne. num_s_cells_z) then
      write(*,*) 'There are ', num_s_cells_z, ' sampling cells in z, but'
      write(*,*) comp_reg_z(num_proc_z+1), ' are assigned to a processor.'
      stop 672
    end if
  
  end subroutine read_geometry






  !----- initialize_output -----
  ! This subroutine initializes the output.txt file, 
  ! reads in the "samplingcells.dat" file, and allocates the sampling cells
  subroutine initialize_output
    use parallel_info_mod
    use simulation_mod
    use core_data_mod
    use constants_mod
    implicit none
    integer i, j, indexi, indexj
	integer scl_unit
	parameter(scl_unit=392)

    ! only initialize the output file if rank=0
    if(mpi_rank .eq. 0)then

      if( single_proc_mode ) then
        open(unit=output_unit,file='output_testmode.txt',status='unknown')
        open(unit=outputtrace_unit,file='output_trace_testmode.txt',status='unknown')
      else
        open(unit=output_unit,file='output.txt',status='unknown')
        open(unit=outputtrace_unit,file='output_trace.txt',status='unknown')
      end if
      
      write(output_unit,"(2(1x,i5),2(1x,1pe12.4),1x,i6,3(1x,1pe12.4),1x,i4,4(1x,1pe12.4))") &
            num_s_cells_r, num_s_cells_z, &
            s_cells_dr,s_cells_dz, &
            output_skip, tau, m, particle_code, num_polys, &
            rmin_global, rmax_global, zmin_global, zmax_global

      write(outputtrace_unit,"(2(1x,i5),2(1x,1pe12.4),1x,i6,3(1x,1pe12.4),1x,i4,4(1x,1pe12.4))") &
            num_s_cells_r, num_s_cells_z, &
            s_cells_dr,s_cells_dz, &
            output_skip, tau, tr_m, trace_code, num_polys, &
            rmin_global, rmax_global, zmin_global, zmax_global
    end if

    ! normal output.txt file
    do i=1,num_polys
        write(output_unit,"(2(1x,i5),11(1pe12.4))") polys(i)%num_pts, polys(i)%type_poly, &
                polys(i)%reflect_coeff, polys(i)%temperature,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
        do j=1,polys(i)%num_pts
            write(output_unit,"(13(1pe12.4))") polys(i)%pts(j)%z, polys(i)%pts(j)%r, &
                0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
        end do
    end do

    ! trace particle output_trace.txt file
    do i=1,num_polys
        write(outputtrace_unit,"(2(1x,i5),11(1pe12.4))") polys(i)%num_pts, polys(i)%type_poly, &
                polys(i)%reflect_coeff, polys(i)%temperature,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
        do j=1,polys(i)%num_pts
            write(outputtrace_unit,"(13(1pe12.4))") polys(i)%pts(j)%z, polys(i)%pts(j)%r, &
                0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
        end do
    end do


    ! Here I actually allocate the sampling cells
    allocate(s_cells(num_s_cells_r,num_s_cells_z))


	! read in their volumes from "samplingcells.dat"
	open(unit=scl_unit,file='geom/samplingcells.dat',status='old',action='read')
	do i = 1,num_s_cells_r
		do j = 1,num_s_cells_z
			read(scl_unit,*) indexi, indexj, s_cells(i,j)%volume
           

			if( s_cells(i,j)%volume <= 0.d0 ) then
				! we'll just set it to something really small but non-zero
				s_cells(i,j)%volume = (fine_cells_dr**3)*1.d-15
			end if

			if( indexi .ne. i .or. indexj .ne. j ) then
				write(*,*) "Error.. sampling cells file written in wrong order, or is wrong dimensions"
				stop 921
			end if
		end do
	end do
	close(unit=scl_unit)

    ! Initialize the rest of their data (all but the volumes)
    call zero_s_cells
  end subroutine initialize_output





  !----- initialize_communication -----
  ! This subroutine initializes the communication cells and the communication buffer
  subroutine initialize_communication
    use parallel_info_mod
    use core_data_mod
    implicit none

    comm_cell_left%num_parts = 0
    comm_cell_left%partition = 0
    comm_cell_left%capacity = starting_capacity
    allocate(comm_cell_left%ps(starting_capacity))

    comm_cell_right%num_parts = 0
    comm_cell_right%partition = 0
    comm_cell_right%capacity = starting_capacity
    allocate(comm_cell_right%ps(starting_capacity))

    comm_cell_top%num_parts = 0
    comm_cell_top%partition = 0
    comm_cell_top%capacity = starting_capacity
    allocate(comm_cell_top%ps(starting_capacity))

    comm_cell_bottom%num_parts = 0
    comm_cell_bottom%partition = 0
    comm_cell_bottom%capacity = starting_capacity
    allocate(comm_cell_bottom%ps(starting_capacity))

    allocate(comm_buffer(8*8*1000000))
    call MPI_BUFFER_ATTACH(comm_buffer,8*8*1000000,mpi_ierr)
  end subroutine initialize_communication




  !----- initialize_particles -----
  ! This subroutine either initializes the collision cells by some algorithm or reads in the restart file
  subroutine initialize_particles
    use simulation_mod
    use helpers_mod
    use constants_mod
    use parallel_info_mod
    use core_data_mod
    use bounds_mod
    implicit none
    integer i
    real(8) r, z
  
    if(use_restart .and. .not. single_proc_mode )then
      call read_restart_file

    else
      do i=1,num_cells
        ! just use the boundary condition, for now
        !call do_left_bc(cells(i))
  
        ! fills each cell with a distribution
        ! of particles, just so we don't start from a vacuum.
        call fill_flat_distrib(cells(i))
      end do

    end if
  end subroutine initialize_particles






  !----- initialize_collision_cells -----
  ! This subroutine reads in the collision cells
  subroutine initialize_collision_cells
    use constants_mod
    use core_data_mod
    use helpers_mod
    use parallel_info_mod
    use simulation_mod
    implicit none
    integer collstuff_u
    parameter(collstuff_u=202)
    type(cell_type) c
    integer i,j


  
    ! first count how many cells there are
    ! also find out the maximum mask index is while reading in this file
    open(unit=collstuff_u,file='geom/collstuff.dat',status='old',action='read')
    j=0
    max_mask = 0
    min_mask = 0
    do i=1,num_cells_global
      read(collstuff_u,*) c%min_r, c%max_r, c%min_z, c%max_z, &
                          c%volume, c%mask
    
      if(c%mask .gt. max_mask)then
        max_mask = c%mask
      end if
    
      if(c%mask .lt. min_mask)then
        min_mask = c%mask
      end if

      if(nint(c%min_r/fine_cells_dr) .ge. nint(rmin_local/fine_cells_dr) .and. &
         nint(c%max_r/fine_cells_dr) .le. nint(rmax_local/fine_cells_dr) .and. &
         nint(c%min_z/fine_cells_dr) .ge. nint(zmin_local/fine_cells_dr) .and. &
         nint(c%max_z/fine_cells_dr) .le. nint(zmax_local/fine_cells_dr)) then
       
        j=j+1
      end if
    end do
    close(unit=collstuff_u)
  
    ! Allocate space
    num_cells = j
    allocate(cells(num_cells))
  
    ! Now read them all in
    open(unit=collstuff_u,file='geom/collstuff.dat',status='old',action='read')
    j=1
    do i=1,num_cells_global


      ! IF NO E-FIELD:
      !read(collstuff_u,*) c%min_r, c%max_r, c%min_z, c%max_z, &
      !                    c%volume, c%mask
      !c%efield_r = 0.0
      !c%efield_z = 0.0


      ! IF THERE IS AN E-FIELD DEFINED IN "collstuff.dat", THEN:
      read(collstuff_u,*) c%min_r, c%max_r, c%min_z, c%max_z, &
                          c%volume, c%mask, c%efield_r, c%efield_z, &
                          c%Te, c%extra
      if( isNaN(c%efield_r) ) c%efield_r = 0.0
      if( isNaN(c%efield_z) ) c%efield_z = 0.0


      if(nint(c%min_r/fine_cells_dr) .ge. nint(rmin_local/fine_cells_dr) .and. &
         nint(c%max_r/fine_cells_dr) .le. nint(rmax_local/fine_cells_dr) .and. &
         nint(c%min_z/fine_cells_dr) .ge. nint(zmin_local/fine_cells_dr) .and. &
         nint(c%max_z/fine_cells_dr) .le. nint(zmax_local/fine_cells_dr)) then
       
        c%min_th = -pi
        c%max_th = pi

        c%num_cand_remaining = 0.d0
        c%max_vrel2 = 3*kb*simT/m
        c%tr_max_vrel2 = 3*kb*simT/tr_m
        c%num_parts = 0
        c%partition = 0
        c%num_left = 0
        c%capacity = starting_capacity
        cells(j) = c
        allocate(cells(j)%ps(starting_capacity))
        j=j+1
      end if
    end do
    close(unit=collstuff_u)
  
  end subroutine initialize_collision_cells




  !----- initialize_fine_cells -----
  ! This subroutine creates the fine cells array
  subroutine initialize_fine_cells
    use core_data_mod
    use constants_mod
    use parallel_info_mod
    use simulation_mod
    implicit none
    integer k,i,j
    integer min_ic,max_ic
    integer min_jc,max_jc
    
    num_fine_cells_r = nint( (rmax_local-rmin_local)/fine_cells_dr )
    num_fine_cells_z = nint( (zmax_local-zmin_local)/fine_cells_dz )
    ! allocate the fine cells
    allocate(fine_cells(num_fine_cells_r, num_fine_cells_z))
    ! initialize them all to 0
    fine_cells = 0
    
    ! now fill in the fine cells
    do k=1,num_cells
      min_ic = nint( (cells(k)%min_r-rmin_local)/fine_cells_dr ) + 1
      max_ic = nint( (cells(k)%max_r-rmin_local)/fine_cells_dr )
      min_jc = nint( (cells(k)%min_z-zmin_local)/fine_cells_dz ) + 1
      max_jc = nint( (cells(k)%max_z-zmin_local)/fine_cells_dz )
      do i=min_ic,max_ic
        do j=min_jc,max_jc
          if(i .gt. num_fine_cells_r)then
            write(*,*)mpi_rank,k,min_ic,max_ic,num_fine_cells_r
          end if
          if(j .lt. 1)then
            write(*,*)mpi_rank,k,min_jc,max_jc,num_fine_cells_z
          end if
          if(fine_cells(i,j) .ne. 0)then
            write(*,*)mpi_rank,'Illegal overlap!',i,j,k
            stop 673
          end if
          fine_cells(i,j) = k
        end do
      end do
    end do
  
    ! now look for zeros in fine_cells
    do i=1,num_fine_cells_r
      do j=1,num_fine_cells_z
        if(fine_cells(i,j) .eq. 0)then
          write(*,*)mpi_rank,'Processor','could not create fine_cells array: some elements zero.'
          write(*,*)mpi_rank,'num_cells=',num_cells
          write(*,*)mpi_rank,i,j
          if(mpi_rank .lt. 15)then
            write(*,*)'Dumping cells for proc',mpi_rank
            write(*,*)mpi_rank,rmin_local,rmax_local
            write(*,*)mpi_rank,zmin_local,zmax_local
            do k=1,num_cells
              write(*,*)mpi_rank,k,cells(k)%min_r,cells(k)%max_r,cells(k)%min_z,cells(k)%max_z
            end do
          end if
          stop 674
        end if
      end do
    end do
  end subroutine initialize_fine_cells







  !----- initialize_polys -----
  ! This subroutine reads in the polygons
  subroutine initialize_polys
    use core_data_mod
    implicit none
    integer spoly_u
    parameter(spoly_u=203)
    integer i,j
  
    ! allocate space
    allocate(polys(num_polys))
  
    ! read in the polygons
    open(unit=spoly_u,file='geom/spoly.dat',status='old',action='read')
    do i=1,num_polys
      read(spoly_u,*) polys(i)%num_pts,polys(i)%type_poly,polys(i)%reflect_coeff,polys(i)%temperature
      ! allocate space for the points in the polygon
      allocate(polys(i)%pts(polys(i)%num_pts))
      do j=1,polys(i)%num_pts
        read(spoly_u,*) polys(i)%pts(j)%z, polys(i)%pts(j)%r
      end do
    end do
    close(unit=spoly_u)
  end subroutine initialize_polys




  !----- initialize_masks -----
  ! This subroutine reads in the mask data
  subroutine  initialize_masks
    use core_types_mod
    use core_data_mod
    implicit none
    integer maskstuff_u
    parameter(maskstuff_u=205)
    real(8) temp
    type(mask_type) msk
    integer i
  
    ! allocate the masks
    allocate(masks(min_mask:max_mask))
  
    ! read them in
    open(unit=maskstuff_u,file='geom/maskstuff.dat',status='old',action='read')
    do i=1,num_masks
       read(maskstuff_u,*) msk%code,msk%reflect, &
            msk%s0z,msk%s0r, &
            msk%w1z,msk%w1r,msk%g1z,msk%g1r, &
            msk%w2z,msk%w2r,msk%g2z,msk%g2r, &
            msk%T
       ! copy that mask into the right place
       temp = sqrt(msk%w1z**2 + msk%w1r**2)
       if(temp .gt. 0)then
          msk%w1z = msk%w1z/temp
          msk%w1r = msk%w1r/temp
       end if
       temp = sqrt(msk%g1z**2 + msk%g1r**2)
       if(temp .gt. 0)then
          msk%g1z = msk%g1z/temp
          msk%g1r = msk%g1r/temp
       end if
       temp = sqrt(msk%w2z**2 + msk%w2r**2)
       if(temp .gt. 0)then
          msk%w2z = msk%w2z/temp
          msk%w2r = msk%w2r/temp
       end if
       temp = sqrt(msk%g2z**2 + msk%g2r**2)
       if(temp .gt. 0)then
          msk%g2z = msk%g2z/temp
          msk%g2r = msk%g2r/temp
       end if
       masks(msk%code) = msk
    end do
    close(unit=maskstuff_u)
  end subroutine  initialize_masks




  !----- zero_s_cells -----
  ! Initializes or re-initializes the sampling cells to contain no data
  ! (except it leaves the volumes alone)
  ! The number of samples is set to 0
  subroutine zero_s_cells
    use core_data_mod
    implicit none
    integer i,j
  
    ! zero the number of samples
    num_samples = 0
    ! zero out the data
    do i=1,num_s_cells_r
      do j=1,num_s_cells_z
        ! normal particle data
        s_cells(i,j)%sum_n = 0
        s_cells(i,j)%sum_vx = 0
        s_cells(i,j)%sum_vy = 0
        s_cells(i,j)%sum_vz = 0
        s_cells(i,j)%sum_vx2 = 0
        s_cells(i,j)%sum_vy2 = 0
        s_cells(i,j)%sum_vz2 = 0
        s_cells(i,j)%sum_vxvy = 0
        s_cells(i,j)%sum_vxvz = 0
        s_cells(i,j)%sum_vyvz = 0

        ! trace particle data
        s_cells(i,j)%tr_sum_n = 0
        s_cells(i,j)%tr_sum_vx = 0
        s_cells(i,j)%tr_sum_vy = 0
        s_cells(i,j)%tr_sum_vz = 0
        s_cells(i,j)%tr_sum_vx2 = 0
        s_cells(i,j)%tr_sum_vy2 = 0
        s_cells(i,j)%tr_sum_vz2 = 0
        s_cells(i,j)%tr_sum_vxvy = 0
        s_cells(i,j)%tr_sum_vxvz = 0
        s_cells(i,j)%tr_sum_vyvz = 0

      end do
    end do
  end subroutine zero_s_cells





!----- read_restart_file -----
subroutine read_restart_file
  use parallel_info_mod
  use core_types_mod
  use core_data_mod
  use simulation_mod
  implicit none
  integer i,j,k,arg_then_tr
  integer onepercent,received_count,size_particle
  real(8) Nef_rstfle, Nef_rstfle_tr,first_r
  integer num_to_receive
  real(8) r, z
  integer fc_i, fc_j, total_num_parts
  type(particle_type) p
  character(10) ctime
  integer restart_flag

  ! arrays that we send with MPI:
  type(particle_type), dimension(1) :: p_buffer
  type(particle_type), allocatable, dimension(:) :: recv_buffer
  integer,             dimension(1) :: flag_buffer
  integer,             dimension(1) :: flag_recv

  ! array with which we make a "map" of the multiprocessor setup
  integer, allocatable, dimension(:) :: ireg
  integer, allocatable, dimension(:) :: jreg
  integer iproc

  ! More integers used with MPI:
  integer bytes_to_send, bytes_to_receive
  logical keep_going

  ! MPI header file
  include 'mpif.h'
  integer status_struct(MPI_STATUS_SIZE)
  integer iostatus

 ! Set what style of restart file is being read
 restart_flag=1 !0 is the binary file, 1 is the ascii style file


  allocate(ireg(num_s_cells_r))
  allocate(jreg(num_s_cells_z))
  size_particle = sizeof(p_buffer)
  bytes_to_send = size_particle


  ! RANK 0 reads restart file
  if( mpi_rank == 0 ) then

	! INITIALIZE MAP SO THAT WE KNOW WHICH RANK PARTICLES BELONG TO

    ! RADIAL
    ! We make an array defined like this: processor_r = ireg(s_cell_value_r)
	! We run one r-processor region at a time, and assign the
    ! processor value to the index in the array which corresponds to a s_cell value in R.
    k=0
    do i=1,num_proc_r
        do j=comp_reg_r(i)+1,comp_reg_r(i+1)
            k=k+1
            ireg(k)=(i-1)
        end do
    end do

    ! AXIAL
    ! We make an array defined like this: processor_z = jreg(s_cell_value_z)
	! We run one z-processor region at a time, and assign the
    ! processor value to the index in the array which corresponds to a s_cell value in Z.
    k=0
    do i=1,num_proc_z
        do j=comp_reg_z(i)+1,comp_reg_z(i+1)
           k=k+1
           jreg(k)=(i-1)
           !write(*,*) k,jreg(k)
        end do
    end do
	! we are now ready to read the restart file

  
    write(*,*) '-------------------------------------------'
    write(*,*) mpi_rank,'Reading binary restart file...'

  
    open(unit=restart_unit, file='restart_file.bin', status='old', action='read',form='unformatted')


    !Start by reading in the number of particles in the file
    read(restart_unit) total_num_parts

    write(*,*) ' total_num_parts in restart file: ',total_num_parts

    !Now read in the argon particle weight to compare with current Nef
    read(restart_unit) Nef_rstfle

    write(*,*) 'Nef from restart file: ',Nef_rstfle


    !Stop if current particle weight doesn't match the old particle weight
    if(Nef_rstfle.ne.Nef) then
       write(*,*) 'Nef and Nef_rstfle do not match'
       stop 773
    end if

    !Define a variable for tracking read progress
    onepercent = (total_num_parts / 100) + 1

    do k=1,total_num_parts

       if(mod(k,onepercent) == 0) then
          call date_and_time(time=ctime)
          write(*,*) "Reading restart_file.bin: ", int((real(k)/real(total_num_parts))*100.0), " % complete. Time: ", ctime
          write(*,*)  "This reads the binary restart file.  Ascii files will cause errors."
       end if       

       ! reading the particles (argon and trace)
       read(restart_unit) p
 

       ! any NaN values are skipped
       if(isnan(p%r) .or. isnan(p%z) .or. isnan(p%x) .or. isnan(p%y) &
           .or. isnan(p%vx) .or. isnan(p%vy) .or. isnan(p%vz)) cycle

       ! Any particles that are out of bounds are skipped
       if(p%r.gt.rmax_global .or. p%z.lt.zmin_global .or. p%z.gt.zmax_global ) cycle


       ! assign the particle to the first (and only) spot in an array 
       ! I do this because I think MPI_BSEND only works with arrays. (I think)
       p_buffer(1) = p

       ! here is where we either insert the particle if it's in our region, 
       ! or send it somewhere else if it isn't

       ! figure out which processor the particle belongs to
       i=p%r/s_cells_dr+1
       j=p%z/s_cells_dz+1
       iproc = jreg(j) + ireg(i)*num_proc_z

       if( iproc < 0 ) then
           write(*,*) "ERROR! iproc is less than zero!"
           stop 63912
       end if

       ! send it to the correct processor
       if( iproc .ne. 0 ) then
          ! in this section we are sending the particle to a processor other
          ! than processor 0

          ! tell the processor that there's still stuff coming
          flag_buffer(1) = 1 ! set this flag to 1
          ! the second argument of MPI_SEND is the number of elements of flag_buffer
          ! iproc tells which processor to send to
          ! the remaining arguments are standard
          ! this just sends the flag integer in flag_buffer
          call MPI_SEND(flag_buffer, 1, MPI_INTEGER, iproc, 0, MPI_COMM_WORLD, mpi_ierr)
          ! send the particle - this is where the actual particle data gets sent
          call MPI_BSEND(p_buffer, bytes_to_send, MPI_BYTE, iproc, 0, MPI_COMM_WORLD, mpi_ierr)

       else
          ! If it belongs to processor rank 0, then we call the function below.

          call insert_into_correct_cell(p)

          ! we clear these because the "insert_into_correct_cell" function will put
          ! the particle into the comm cells if it doesn't belong there, but we don't
          ! want to do that in this case. It shouldn't happen, though.  Maybe only very rarely.
          comm_cell_left%num_parts = 0
          comm_cell_right%num_parts = 0
          comm_cell_top%num_parts = 0
          comm_cell_bottom%num_parts = 0

          comm_cell_left%partition = 0
          comm_cell_right%partition = 0
          comm_cell_top%partition = 0
          comm_cell_bottom%partition = 0
       end if	
 
        !end of particle read loop
    end do

    ! sending flag that we're done sending particles
    flag_buffer(1) = 0
    do i = 1,num_mpi_processes-1
       call MPI_SEND(flag_buffer, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, mpi_ierr)
    end do

    !close file
    close(unit=restart_unit)
    write(*,*) mpi_rank,'Done reading binary restart file.'
 



  !end "SEND" portion of If-block




  ! RECEIVE
  else
 
    ! FIRST,
    !   RECEIVE ARGON PARTICLES
    ! Second,
    !   RECEIVE TRACE PARTICLES
        keep_going = .true.
        received_count = 0
        do while( keep_going .eqv. .true. )

            ! Probing a "flag" receive from rank 0 - this just waits for a send to come in
            call MPI_PROBE(0, 0, MPI_COMM_WORLD, status_struct, mpi_ierr)
        
            ! receiving "flag" integer
            call MPI_RECV(flag_recv, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status_struct, mpi_ierr)

            if( flag_recv(1) == 1 ) then

                ! Probing a Particle receive from rank 0
                call MPI_PROBE(0, 0, MPI_COMM_WORLD, status_struct, mpi_ierr)
                ! receive the # of bytes to receive
                call MPI_GET_COUNT(status_struct,MPI_BYTE,bytes_to_receive, mpi_ierr)

                !perform a check
                num_to_receive = bytes_to_receive/size_particle
                if(num_to_receive*size_particle .ne. bytes_to_receive)then
                    write(*,*) 'Discrepancy in read_restart_file subroutine.',num_to_receive,size_particle,bytes_to_receive
                    stop 668
                end if
                ! allocate space
                allocate(recv_buffer(num_to_receive))

                ! receive them
                call MPI_RECV(recv_buffer, bytes_to_receive, MPI_BYTE, 0, 0, MPI_COMM_WORLD, status_struct, mpi_ierr)

                do i = 1,num_to_receive
                    ! we now interpret these received bytes as a particle array and put it in p
                    p = recv_buffer(i)
                    call insert_into_correct_cell(p)
                end do

                deallocate(recv_buffer)

                ! we clear these because the "insert_into_correct_cell" function will put
                ! the particle into the comm cells if it doesn't belong there, but we don't
                ! want to do that in this case.  It shouldn't happen, though. Maybe only very rarely.
                comm_cell_left%num_parts = 0
                comm_cell_right%num_parts = 0
                comm_cell_top%num_parts = 0
                comm_cell_bottom%num_parts = 0

                comm_cell_left%partition = 0
                comm_cell_right%partition = 0
                comm_cell_top%partition = 0
                comm_cell_bottom%partition = 0

                ! iterating the number of particles received
                received_count = received_count + 1
                !if( mod(received_count,100000) == 0 ) then
                !    write(*,*) mpi_rank, "successfully received this many:", received_count
                !end if

            else
                keep_going = .false.
            end if

        ! end adding particles loop
        end do


  ! end of "RECEIVE" else-block
  end if

end subroutine read_restart_file



!----- write_restart_file -----
subroutine write_restart_file
  use parallel_info_mod
  use core_data_mod
  use simulation_mod
  implicit none
  integer i,k
  integer find_whole_num_parts !a function I call
  integer num_parts_tot(1), num_parts_arg(1)
  integer, allocatable, dimension(:) :: recv_array1, recv_array2
  integer ierr
  character, dimension(2) :: char_buf !just a place to receive the message 'OK'
  ! To use MPI (Message Passing Interface) I have to include this header file
  include 'mpif.h'
  integer status_struct(MPI_STATUS_SIZE)
  character num_buf*12 ! a place to put the processor number

! write the processor rank into the character variable num_buf
  write(num_buf,*)mpi_rank
! left adjust it
  num_buf = adjustl(num_buf)

  ! Part 1: Writing particles themselves to files

! open the files in the directory restart_files, each one tagged with num_buf, the
! mpi_rank. Do this for background gas and trace gas
  open(unit=restart_unit,file='restart_files/file'//trim(num_buf)//'.bin',status='unknown',form='unformatted')


! call the routine that writes the particle data into each of the files that
! were just created
  call write_restart_data(restart_unit)

! close all of these files
  close(unit=restart_unit)



  ! Part 2: Figure out total numbers of particles in each restart file

! allocate the receive arrays for doing the gather from all of the
! processors so that we know how many particles there are, both
! background and trace
  allocate(recv_array1(num_mpi_processes))


! each processor will now sum over its cells to get the total number of
! particles, background+trace
  num_parts_tot(1) = sum(cells(:)%num_parts) 
              
! each processor will now sum over its cells to get the total number of
! background particles only
  num_parts_arg(1) = sum(cells(:)%partition)

! the number of trace particles is the num_parts_tot(1)-num_parts_arg(1)

! gather num_parts_tot up and have them all collected into recv_array1
! arguments 1-3: variable to gathered, # of elements, variable type
! arguments 4-6: array to collect into, # of elements being passed to it, variable type
! argument 7: processor # of the receiving processor
! arguments 8-9: mpi commumication stuff
  call MPI_GATHER(num_parts_tot, 1, MPI_INTEGER, &
                  recv_array1, 1, MPI_INTEGER, &
                  0, MPI_COMM_WORLD, mpi_ierr)


  if(mpi_rank .eq. 0)then

! save the boundary condition files into files with new names
! so that the run can be restarted from this new restart file
! if the run crashes or runs out of time


   call system("cp top_bound.txt top_bound_save.txt")
   call system("cp end_bound.txt end_bound_save.txt")



! open the file into which the background gas particles will be gathered
! using a system cat command. Begin by writing the number of particles and Nef
    open(unit=231,file='restart_sum_argon.bin',status='unknown',form='unformatted')
    write(231) sum(recv_array1)
    write(231) Nef
    close(unit=231)

 
       ! Build the restart files in /restart_files with cat; them move them to the working directory.
       ! This is so that if the run dies (wall-time death, for instance) during the long cat process
       ! the restart file in the working directory will not be ruined.

! cat the 2-line header file with the processor particle files in restart_files and build
! the restart_file.txt restart file in restart_files
       call system("cat restart_sum_argon.bin restart_files/file*.bin >  restart_files/restart_file.bin")
! once it has been built down there, move it up to the execution directory
       call system("mv  restart_files/restart_file.bin restart_file.bin")
    end if

 


   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! deallocate the receive arrays
  deallocate(recv_array1)



! clean up the processor files in restart_files, left over from the previous run
  if(mpi_rank.eq.0) then
    call system("rm restart_files/file*.bin")
  end if


end subroutine write_restart_file




!----- write_restart_data -----
! this subroutine writes one process' particles out to the restart file
subroutine write_restart_data(restart_unit)
  use core_data_mod
  implicit none
  integer restart_unit, tr_restart_unit
  integer i,k
  type(particle_type) p

  ! all particles
  do i=1,num_cells
    do k=1,cells(i)%num_parts
      p = cells(i)%ps(k)
      write(restart_unit) p
    end do
  end do

 
end subroutine write_restart_data







!----- insert_into_correct_cell -----
! Used by "read_restart_file"
! This subroutine takes a particle and inserts it into the correct cell
! it calls cell_insert_particle, which is in helpers_mod in the core module
subroutine insert_into_correct_cell(p)
  use core_data_mod
  use parallel_info_mod
  use helpers_mod
  use constants_mod
  use simulation_mod
  implicit none
  type(particle_type) p
  integer fc_i,fc_j,i,sc_j,sc_i

  ! is the particle in a simulation region at all?
  if( p%r .lt. rmin_global .or. p%r .gt. rmax_global .or. &
      p%z .lt. zmin_global .or. p%z .gt. zmax_global )then
    return
  else
    ! is the particle is my simulation region?
    sc_i = 1 + (p%r - rmin_global)/s_cells_dr
    sc_j = 1 + (p%z - zmin_global)/s_cells_dz

    if(sc_j .lt. min_j())then
      if(sc_j .ge. 1) then
        call cell_insert_particle( comm_cell_left, p )
      end if
    else if(sc_j .gt. max_j())then
      if(sc_j .le. num_s_cells_z)then
        call cell_insert_particle( comm_cell_right, p )
      end if
    else if(sc_i .lt. min_i())then
      if(sc_i .ge. 1) then
        call cell_insert_particle( comm_cell_bottom, p )
      end if
    else if(sc_i .gt. max_i())then
      if(sc_i .le. num_s_cells_r)then
        call cell_insert_particle( comm_cell_top, p )
      end if
    else

      ! The particle is in my simulation region: insert it
      fc_i = 1 + (p%r - rmin_local)/fine_cells_dr
      fc_j = 1 + (p%z - zmin_local)/fine_cells_dz
      i = fine_cells(fc_i,fc_j)

      call cell_insert_particle( cells(i), p )
    end if
  end if

end subroutine insert_into_correct_cell

subroutine initialize_trace_collide_data

!This subroutine initializes the collision parameters for the trace particles
!Trace switch is set in core.f90
use constants_mod
use simulation_mod
!First load the argon id code:
particle_code=1.1180E+00
!The name of the element is encoded in a real number to make it backwards
!compatable with old post processing software.
!The first two digits save the first letter of the Element
!They are encoded as 11=A;12=B; up to 36=Z; 
!(The decimal point is ignored)
!Starting at 11 avoids rounding problems when the number is printed
!The 3rd and forth digits hold the second letter of the element.
! 00 means there is no second letter, then the rest of the alphabet
! is encoded as 01 = a, 02=b, up to 26=z;
!A few examples:
  ! for argon, the element is Ar, A=11, r=18,
  ! so the first four digits would be 1.118
  !Carbon - (C) C=13, and there is no second digit,
  !So C -> 1.300
!The fifth digit is set to 1 for ions and 0 for neutrals
!The second digit of the exponent saves whether the particle 
!was in the background flow (0) or is analyte (1)
!for example: Neutral argon background flow is 1.1180E+00
!whereas Argon ion analyte flow is 1.1181E+01

  if(trace_switch.eq.1)then !Barium

     tr_q = 1.60217646e-19 ! Charge in C
     tr_m = 2.280d-25 !Mass in kg                

     !Collision Parameters
     aref_tr = 4.61d-14
     nu_tr1 = 0.823
     bref_tr=8.52e-20
     nu_tr2=0d0
     trace_code=1.2011E+01

  elseif(trace_switch.eq.2)then !Calcium

     tr_q = 1.60217646e-19 ! Charge in C
     tr_m = 6.642e-26 !Mass in kg                

     !Collision Parameters
     aref_tr = 8.575d-14
     nu_tr1 = 0.842
     bref_tr=7.29e-20
     nu_tr2=0d0
     trace_code=1.3011E+01

  elseif(trace_switch.eq.3)then !Argon

     tr_q = 1.60217646e-19 ! Charge in C
     tr_m = 6.642e-26 !Mass in kg                

     !Collision Parameters
     aref_tr =  1.84d-17
     nu_tr1 =   0.162d0
     bref_tr =  1.10d-30
     nu_tr2=   -1.57d0
     trace_code=1.1181E+01
  elseif(trace_switch.eq.4)then !Lithium
     tr_q = 1.60217646e-19 ! Charge in C
     tr_m = 1.152e-26 !Mass in kg                

     !Collision Parameters
     aref_tr =8.163e-15 
     nu_tr1 = 0.604
     bref_tr = -2.333e-20
     nu_tr2 = 0
     trace_code=2.2091E+01

  elseif(trace_switch.eq.5)then !Mercury
     tr_q = 1.60217646e-19 ! Charge in C
     tr_m = 1.152e-26 !Mass in kg                

     !Collision Parameters
     aref_tr =2.451e-16
     nu_tr1 = 0.3389
     bref_tr = -2.39e-19
     nu_tr2 = 0
     trace_code=1.8071E+01

  end if



end subroutine initialize_trace_collide_data


end module initialize_mod






