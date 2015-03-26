






!----- fenix10 -----
! This is the main program
program fenix10

  ! from "core.f90":
  use constants_mod
  use simulation_mod
  use parallel_info_mod
  use core_data_mod
  ! from "initialize.f90":
  use initialize_mod
  ! from "bounds.f90":
  use bounds_mod
  ! from "collide.f90":
  use collide_mod
  !from "move.f90":
  use move_mod
  ! from "reloadcells.f90":
  use reloadcells_mod
  ! from "dataio.f90":
  use dataio_mod

  implicit none
  include 'mpif.h'
  
  character(10) ctime
  integer hour,minute,profileflag,err,i,j,writeflag
  integer maxnum,minnum
  integer debug_flag
  
  real(8) radial_avg_parts
  real(8) tic, toc, sec, all_tic, all_toc
  character num_buf*12 ! a place to put the processor number


! set profileflag to 1 if you want timing data on the subroutines called
! in the timestep loop below
! set profileflag to 2 if you want a log of the ghost cell
! particles created every time step
  data profileflag/0/
!set debug_flag to 1 to make fenix report when it enters and exits each
!subroutine along with the timestep number
  data debug_flag /0/


! open a file into which axis particles will be dumped for
! post-processing to analyze distribution functions



!   if(mpi_rank.eq.0) then
!     open(unit=311,file='dump.txt',status='unknown')
!   end if

!kluge: Open a file to dump grid densities for balancing number of particles on each processor
!   if(mpi_rank.eq.0) then
!     open(unit=747,file='cell_dens.txt',status='unknown')
!   end if

!Open A file to dump the ghost particles on.
   if(profileflag.eq.2) then
     open(unit=310,file='ghostparts.txt',status='unknown')
   end if

  call initialize_all()

  cur_step = 0
  call write_particles_per_processor()
  if( mpi_rank == 0 ) write(*,*) "Initial Standard Out:"
  call update_std_out()

!Time how long the program takes once it's loaded

     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      all_tic=3600*hour+60*minute+sec
     end if
  
!*********** top of timestep loop ********************
  do cur_step=1,nsteps


    if(mod(cur_step,output_skip).eq.0)then
        maxnum=0
        minnum=0
        radial_avg_parts=0
        do i=1,num_cells
           maxnum=max(maxnum,cells(i)%partition)

            if(cells(i)%min_r.lt.1.d-15) then
                ! compute the average number in each processor for the cells near r=0
                radial_avg_parts=radial_avg_parts+cells(i)%partition
                minnum=minnum+1
            end if
        end do
     radial_avg_parts=radial_avg_parts/real(minnum)

  !  write(*," ('Max in a cell, Avg. particles on axis, processor',i5,1x,1pe12.4,1x,i5) ") &
  !           maxnum,radial_avg_parts,mpi_rank
   end if


   !******************** do_bc *** boundary conditions
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if
   
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering do_bc, current step = ',cur_step
!   if(cur_step.gt.150) write(*,*) 'Entering do_bc, current step =',cur_step,'MPI_rank',mpi_rank

    call do_bc
    !When in debug mode, stop all processors after each subroutine
    if(debug_flag.eq.1) call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting do_bc, current step = ',cur_step


     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'do_bc ',mpi_rank,toc-tic,cur_step
     end if
   !******************** do_bc *** boundary conditions


   !******************** collide *** collide particles
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if

    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering collide, current step = ',cur_step
 !  if(cur_step.gt.150) write(*,*) 'Entering collide, current step =',cur_step,'MPI_rank',mpi_rank
    call collide

    if(debug_flag.eq.1) call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting collide, current step = ',cur_step


     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'collide ',mpi_rank,toc-tic,cur_step
     end if
   !******************** collide *** collide particles


   !******************** move *** move particles
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if

    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering move, current step = ',cur_step
 !  if(cur_step.gt.150) write(*,*) 'Entering move, current step =',cur_step,'MPI_rank',mpi_rank
    call move

  


   !if(mpi_rank.eq.0) then
   !  write(*,*) gc(1)%partition,gc(1)%num_parts-gc(1)%partition,' $$$'
  !end if

     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'move ',mpi_rank,toc-tic,cur_step
     end if
    if(debug_flag.eq.1) call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting move, current step = ',cur_step
   !******************** move *** move particles


   !******************** reload_cells *** reload cells with particles
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if

    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering reload_cells, current step = ',cur_step
   !if(cur_step.gt.150) write(*,*) 'Entering reload_cells, current step =',cur_step,'MPI_rank',mpi_rank
    call reload_cells

   !if(cur_step.gt.150) write(*,*) 'Exiting reload_cells, current step =',cur_step,'MPI_rank',mpi_rank
    

     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'reload_cells ',mpi_rank,toc-tic,cur_step
     end if
    if(debug_flag.eq.1) call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting reload_cells, current step = ',cur_step
   !******************** reload_cells *** reload cells with particles


   !******************** send_recv_parts *** move particles between processors
     if(profileflag.eq.1) then
	!call MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr) !Added to testing actual
	  ! time spent  in send_recv_parts subroutine
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if

    ! call it twice to handle the corner neighbors
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering send_recv_parts, current step = ',cur_step
    call send_recv_parts
    call send_recv_parts

    if(debug_flag.eq.1) call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting send_recv_parts, current step = ',cur_step

    !Now keep the processors on the same timestep
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'send_recv_parts ',mpi_rank,toc-tic,cur_step
     end if
   !******************** send_recv_parts *** move particles between processors

   !******************** collect_data *** collect data from all processors into output arrays
     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      tic=3600*hour+60*minute+sec
     end if

    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Entering collect_data, current step = ',cur_step
 !  if(cur_step.gt.150) write(*,*) 'Entering collect_data, current step =',cur_step,'MPI_rank',mpi_rank
    call collect_data

    if(debug_flag.eq.1.and.mpi_rank.eq.0) write(*,*) 'Exiting collect_data, current step = ',cur_step
  

     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3,i5)") 'collect_data ',mpi_rank,toc-tic,cur_step
     end if
   !******************** collect_data *** collect data from all processors into output arrays


    ! collected data written intermittently
    if(mod(cur_step,output_skip).eq.0)then
!     write(*,*) ' Going into write_collected_data ',mpi_rank
      call write_collected_data
!     write(*,*) ' Coming out of write_collected_data ',mpi_rank
!     call update_bc_data  ! turn off boundary cell updating
      ! now erase the sampling cells, the data is no longer needed
      call zero_s_cells
    end if

    ! Standard Out updated intermittently
    if(mod(cur_step,10).eq.0)then
      call update_std_out()
    end if

    ! writes the number of particles per processor to the standard out, intermittently
    if(mod(cur_step,500).eq.0)then
      call write_particles_per_processor()
    end if

    ! restart file re-written intermittently
    if(mod(cur_step,restart_skip).eq.0)then
     write(*,*) mpi_rank, "Writing Restart File..."
     call write_restart_file()
     writeflag=1
    else
     writeflag=0
    end if

  end do
!*********** end of timestep loop **********************

  ! the restart_file is re-written one last time at the end of the run
  ! unless you just did it

  writeflag=1 ! don't do it

  if(writeflag.eq.0) then
     call write_restart_file
  end if
!End program timer

     if(profileflag.eq.1) then
      call date_and_time(time=ctime)
      read(ctime,"(i2,i2,f6.3)") hour,minute,sec
      all_toc=3600*hour+60*minute+sec
      write(*,"(a15,i5,f8.3)") 'Total time ',mpi_rank,all_toc-all_tic
     end if

  if( mpi_rank == 0 ) write(*,*) "End of program reached"

  ! I still have to finalize the MPI environment
  call MPI_FINALIZE(mpi_ierr)
  if(mpi_ierr .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,0,mpi_ierr)
  
end program fenix10







