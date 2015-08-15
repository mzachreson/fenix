






module dataio_mod
  implicit none
  
  contains

  !----- collect_data -----
  ! This subroutine collects the data over this processors' sampling cells
  subroutine collect_data
    use constants_mod
    use core_data_mod
    use simulation_mod
    implicit none
    integer i, k, s_i, s_j
    type(particle_type) p
  
    do i=1,num_cells

      ! data for normal argon particles
      do k=1,cells(i)%partition
        p = cells(i)%ps(k)
        s_i = 1 + (p%r - rmin_global)/s_cells_dr
        s_j = 1 + (p%z - zmin_global)/s_cells_dz
      
        s_cells(s_i,s_j)%sum_vx = s_cells(s_i,s_j)%sum_vx + p%vx
        s_cells(s_i,s_j)%sum_vy = s_cells(s_i,s_j)%sum_vy + p%vy
        s_cells(s_i,s_j)%sum_vz = s_cells(s_i,s_j)%sum_vz + p%vz
        s_cells(s_i,s_j)%sum_vx2 = s_cells(s_i,s_j)%sum_vx2 + p%vx**2
        s_cells(s_i,s_j)%sum_vy2 = s_cells(s_i,s_j)%sum_vy2 + p%vy**2
        s_cells(s_i,s_j)%sum_vz2 = s_cells(s_i,s_j)%sum_vz2 + p%vz**2
        s_cells(s_i,s_j)%sum_vxvy = s_cells(s_i,s_j)%sum_vxvy + p%vx*p%vy
        s_cells(s_i,s_j)%sum_vxvz = s_cells(s_i,s_j)%sum_vxvz + p%vx*p%vz
        s_cells(s_i,s_j)%sum_vyvz = s_cells(s_i,s_j)%sum_vyvz + p%vy*p%vz
        s_cells(s_i,s_j)%sum_n = s_cells(s_i,s_j)%sum_n + 1.d0
      end do

      ! data for trace particles
      if( cells(i)%partition < cells(i)%num_parts ) then
        do k=cells(i)%partition+1, cells(i)%num_parts
            p = cells(i)%ps(k)
            s_i = 1 + (p%r - rmin_global)/s_cells_dr
            s_j = 1 + (p%z - zmin_global)/s_cells_dz
      
            s_cells(s_i,s_j)%tr_sum_vx = s_cells(s_i,s_j)%tr_sum_vx + p%vx
            s_cells(s_i,s_j)%tr_sum_vy = s_cells(s_i,s_j)%tr_sum_vy + p%vy
            s_cells(s_i,s_j)%tr_sum_vz = s_cells(s_i,s_j)%tr_sum_vz + p%vz
            s_cells(s_i,s_j)%tr_sum_vx2 = s_cells(s_i,s_j)%tr_sum_vx2 + p%vx**2
            s_cells(s_i,s_j)%tr_sum_vy2 = s_cells(s_i,s_j)%tr_sum_vy2 + p%vy**2
            s_cells(s_i,s_j)%tr_sum_vz2 = s_cells(s_i,s_j)%tr_sum_vz2 + p%vz**2
            s_cells(s_i,s_j)%tr_sum_vxvy = s_cells(s_i,s_j)%tr_sum_vxvy + p%vx*p%vy
            s_cells(s_i,s_j)%tr_sum_vxvz = s_cells(s_i,s_j)%tr_sum_vxvz + p%vx*p%vz
            s_cells(s_i,s_j)%tr_sum_vyvz = s_cells(s_i,s_j)%tr_sum_vyvz + p%vy*p%vz
            s_cells(s_i,s_j)%tr_sum_n = s_cells(s_i,s_j)%tr_sum_n + 1.d0
        end do
      end if

    end do

    num_samples = num_samples + 1
  end subroutine collect_data





  !----- write_collected_data
  ! This is an attempt to clean up data output by changing the form of data
  !on each processor to a summable form and then using MPI_
  subroutine write_collected_data(logdens,ambipolar_switch)
    use constants_mod
    use simulation_mod
    use parallel_info_mod
    use core_data_mod
    implicit none

    real(8) r, z, V


    integer i,j,k,l,n
    real(8) logdens(num_s_cells_r,num_s_cells_z)
    real(8),allocatable,dimension(:,:,:) :: particle_info
    real(8),allocatable,dimension(:,:,:) :: particle_info_all

! ambipolar field
    integer ambipolar_switch
    real(8),allocatable,dimension(:,:) :: newlogdens
    real(8),allocatable,dimension(:,:) :: Er
    real(8),allocatable,dimension(:,:) :: Ez
    real(8),allocatable,dimension(:,:) :: Te  ! electron temperature in eV
    real(8),allocatable,dimension(:,:) :: work
    real(8) dr,dz,under,chk
    real(8) x(100),a11,a12,a22,b1,b2,denom,a,b
    integer stride
    ! counters for controlling ambipolar field calculation
    integer ncall,nsteady,FlagAmbipolar,FirstTime,di
    real(8) Npart,Nprev

    data ncall/0/, nsteady/0/,Nprev/0/,Npart/0/,FlagAmbipolar/0/,FirstTime/0/

  ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'

!allocating arrays to dump averages into
! everybody needs particle_info because they send their data in it
    allocate(particle_info(20,num_s_cells_r,num_s_cells_z))

! allocate ambipolar field log(density) array

! these are only needed by processor 0, but why not
    allocate(newlogdens(num_s_cells_r,num_s_cells_z))
    allocate(Te(num_s_cells_r,num_s_cells_z))
    allocate(work(num_s_cells_r,num_s_cells_z))



!everyone needs these
    allocate(Er(num_s_cells_r,num_s_cells_z))
    allocate(Ez(num_s_cells_r,num_s_cells_z))

! only mpi_rank=0 needs the array where the sum goes
!   if(mpi_rank.eq.0) then
       allocate(particle_info_all(20,num_s_cells_r,num_s_cells_z))
 !     write(*,*) 'particle_info_all allocated'
!   end if

    !Make all the processors catch up to avoid problems with changing
    !particle_info_all
 !  write(*,*) 'Before particle_info barrier', mpi_rank
    call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
 !  write(*,*) 'After particle_info barrier', mpi_rank

!kluge: Open a file to dump grid densities for balancing number of particles on each processor
!   open(unit='747',file='cell_dens.txt',status='unknown')


  do j=1,num_s_cells_z
     do i=1,num_s_cells_r
        
       ! normal particles averaging
       if (s_cells(i,j)%sum_n /= 0) then
		!Needs comments on what the averages mean
          particle_info(1,i,j) = s_cells(i,j)%sum_n
          particle_info(8,i,j) = s_cells(i,j)%sum_vx
          particle_info(9,i,j) = s_cells(i,j)%sum_vy
          particle_info(10,i,j) = s_cells(i,j)%sum_vz
          particle_info(2,i,j) = s_cells(i,j)%sum_vx2
          particle_info(3,i,j) = s_cells(i,j)%sum_vy2
          particle_info(4,i,j) = s_cells(i,j)%sum_vz2
          particle_info(5,i,j) = s_cells(i,j)%sum_vxvy
          particle_info(6,i,j) = s_cells(i,j)%sum_vxvz
          particle_info(7,i,j)= s_cells(i,j)%sum_vyvz
       else
          do k=1,10
	     particle_info(k,i,j)=0
          end do
       end if


       ! trace particles averaging
       if (s_cells(i,j)%tr_sum_n /= 0) then
          particle_info(11,i,j) = s_cells(i,j)%tr_sum_n
          particle_info(18,i,j) = s_cells(i,j)%tr_sum_vx
          particle_info(19,i,j) = s_cells(i,j)%tr_sum_vy
          particle_info(20,i,j) = s_cells(i,j)%tr_sum_vz
          particle_info(12,i,j) = s_cells(i,j)%tr_sum_vx2
          particle_info(13,i,j) = s_cells(i,j)%tr_sum_vy2
          particle_info(14,i,j) = s_cells(i,j)%tr_sum_vz2
          particle_info(15,i,j) = s_cells(i,j)%tr_sum_vxvy
          particle_info(16,i,j) = s_cells(i,j)%tr_sum_vxvz
          particle_info(17,i,j) = s_cells(i,j)%tr_sum_vyvz
       else
          do k=11,20
	     particle_info(k,i,j)=0
          end do
       end if
     end do
  end do  

 ! write(*,*) 'Calling MPI_REDUCE on particle_info_all',mpi_rank
!Bringing all of the info onto the first processor using MPI_SUM as the
!operator so that everybody is summed into particle_info_all on Processor 0
   call MPI_REDUCE(particle_info,particle_info_all,20*num_s_cells_r*num_s_cells_z,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
!Writing the data to output.txt and output_trace.txt
 ! write(*,*) 'Writing to file',mpi_rank

  ! Now, only on Processor 0, write the summed and averaged data to the argon and trace output files
  ! We also use this information in update_bc_data, so first thing we will load the summed data
  ! from all of the processors into the cell data for processor 0. This will mess these data up, but only
  ! for long enough for us to (1) call update_bc_data, then (2) zero out all of the sampling cell data
  ! on all of the processors.  In main.f90 the order must be:  (1) write the collected data (this routine)
  ! (2) update the boundary conditions (3) zero the sampling cell data

!************* top of mpi_rank=0 if block
  if (mpi_rank == 0) then

     ! increment the ncall counter, number of times that write_collected_data has been called
     ncall = ncall+1
     write(*,*) 'ncall in write_collected_data: ',ncall

     ! zero the Npart counter for counting how many trace simulation particles have
     ! been accumulated for this call to write_collected_data, and save the previous
     ! value in Nprev
     Nprev = Npart
     Npart = 0

     ! find the current Npart
     do j=1,num_s_cells_z
     do i=1,num_s_cells_r    
         ! add up the trace particles since the last call to write_collected_data
          Npart = Npart + particle_info_all(11,i,j)
     end do
     end do

     ! Compare Nprev and Npart. If they are close, start averaging logdens
     ! and computing the ambipolar field

     chk = (Npart-Nprev)*2.d0/(Npart+Nprev)
  ! kluge
     write(*,153) Npart,Nprev,chk
153  format(' Npart, Nprev, chk: ',3(1pe12.4))

     if(abs(chk).lt.0.003) then
       FlagAmbipolar=1 ! tell the code that it's time to start averaging logdens
                       ! and computing Eambipolar. The ambipolar field will be
                       ! crude at first, then get better. Or, you could wait
                       ! a while for logdens to accumulate before starting to
                       ! compute Eambipolar. nsteady will count how many times
                       ! logdens has been accumulated
     end if

     ! with the ambipolar flag turned on, it's time to start incrementing the
     ! steady state counter
     if(FlagAmbipolar.eq.1) then
       nsteady = nsteady+1
     end if

   ! kluge
   ! write(*,*) 'nsteady, FlagAmbipolar: ',nsteady,FlagAmbipolar

     do j=1,num_s_cells_z
     do i=1,num_s_cells_r    

      ! load the summed data from all processors into the cell data for processor 0
      ! for use by update_bc_data, argon information
        s_cells(i,j)%sum_n = particle_info_all(1,i,j)
        s_cells(i,j)%sum_vx = particle_info_all(8,i,j)
        s_cells(i,j)%sum_vy = particle_info_all(9,i,j)
        s_cells(i,j)%sum_vz = particle_info_all(10,i,j)
        s_cells(i,j)%sum_vx2 = particle_info_all(2,i,j)
        s_cells(i,j)%sum_vy2 = particle_info_all(3,i,j)
        s_cells(i,j)%sum_vz2 = particle_info_all(4,i,j)
        s_cells(i,j)%sum_vxvy = particle_info_all(5,i,j)
        s_cells(i,j)%sum_vxvz = particle_info_all(6,i,j)
        s_cells(i,j)%sum_vyvz = particle_info_all(7,i,j)

          !Loading Geometry info
          r = (i*s_cells_dr + (i-1)*s_cells_dr)/2.d0
          z = (j*s_cells_dz + (j-1)*s_cells_dz)/2.d0

          ! "s_cells(i,j)%volume" is a corrected volume, taking into account
          ! intersection with the metal.  It was read in from "samplingcells.dat"
          V = s_cells(i,j)%volume

          !Averaging cell data and writing it for Argon

          ! Turn particle_info_all(1,i,j) into argon density and put it in particle_info(1,i,j)
          if (particle_info_all(1,i,j)/=0) then
             particle_info(1,i,j) = Nef * particle_info_all(1,i,j) / num_samples / V

             ! average the other quantities from particle_info_all and put them in particle_info
             do k=2,10
                particle_info(k,i,j)=particle_info_all(k,i,j) / particle_info_all(1,i,j)
             end do
          else
             do k=1,10
                particle_info(k,i,j)=0
             end do
          end if

          write(output_unit,"(13(1pe12.4))") (particle_info(k,i,j),k=1,10),r,z,V
          !Kluge: writing out number of particles in each cell
          !write(747," ( 'num parts, r grid, z grid',1pe12.4,2(i5) ) ") particle_info_all(1,i,j),i,j




          !Averaging cell data and writing it for trace particles

          ! turn particle_info_all(11,i,j) into trace density and put it in particle_info(11,i,j)
          if (particle_info_all(11,i,j)/=0) then
             particle_info(11,i,j) = Nef * particle_info_all(11,i,j) / num_samples / V

 ! load the amibpolar field array logdens, non-zero
            if(nsteady.ge.1) then
             newlogdens(i,j)=log(particle_info(11,i,j))  ! use the trace density in position 11
            end if

             ! average the other quantities from particle_info_all and put them in particle_info
             do k=12,20
                particle_info(k,i,j)=particle_info_all(k,i,j) / particle_info_all(11,i,j)
             end do
          else
             do k=11,10
                particle_info(k,i,j)=0
             end do


 ! load the ambipolar field array logdens, zero since there is no density if control is here
            if(nsteady.ge.1) then
             newlogdens(i,j)=0.d0
            end if

          end if

          ! write the data to the trace output file from particle_info
          write(outputtrace_unit,"(13(1pe12.4))") (particle_info(k,i,j),k=11,20),r,z,V
     end do
     end do


 ! log(n) has bad statistics at r=0 sometime, fix them by fitting to r.^2, linear
 !*********************************************************************



          do j=1,num_s_cells_z

              di=5  ! correct the first di-1 points, if needed

              do  k=1,di-1

                  i=di-k

                  do n=1,di
                   x(n)=i+n-.5  ! grid points beyond i, di of them
                  end do

                 ! least squares fit to these points, a + b*(i-.5)^2
                  a11=0.d0
                  a12=0.d0
                  a22=0.d0
                  b1=0.d0
                  b2=0.d0


                  do  n=1,di
                     a11=a11+1.d0
                     a12=a12+x(n)**2
                     a22=a22+x(n)**4
                     b1=b1+newlogdens(i+n,j)
                     b2=b2+newlogdens(i+n,j)*x(n)**2
                  end do

                  ! solve for a and b
                  denom=a11*a22-a12**2
                  a=(a22*b1-a12*b2)/denom
                  b=(a11*b2-a12*b1)/denom

                  chk=abs(1-(a+b*(i-.5)**2)/newlogdens(i,j))  ! check the error on log(n) at i


                  ! if log(n(i)) isn't close to the fit to its neighbors, replace it with the fit
                  if(chk.gt.0.01) then
                      newlogdens(i,j)=a+ b*(i-.5)**2
                  end if

              end do

          end do



 !******* end of log(n) correction *****************
          ! Perform a running average of logdens and compute the ambipolar electric field
          ! as soon as steady state is reached (nsteady starts incrementing)

 !***top of ambipolar if block ******************************
     if(nsteady.ge.1.and.ambipolar_switch.eq.1) then

      ! if this is the first time, just set logdens to the newlogdens that was just computed
      if(nsteady.eq.1) then
       logdens = newlogdens   ! takes care of the un-initialized logdens the first time, just being careful
      else ! average the new one with the previous average if nsteady>1

           ! logdens should get smoother as the run proceeds, if it gets saved from call to call
        logdens = ( real(nsteady-1)*logdens + newlogdens )/real(nsteady)

      end if

!     write(*,*) '%%% logdens check ',logdens(30,75),newlogdens(30,75)

! build a reasonable electron temperature function Te(r,z)

      do i=1,num_s_cells_r
      do j=1,num_s_cells_z

         r = rmin_global + (i-.5d0)*s_cells_dr
         z = zmin_global + (j-.5d0)*s_cells_dz

         Te(i,j) = 2000.d0/( 1+(r**2+z**2)/1d-3**2 ) + 4000.d0 ! Te in K, spherical
         Te(i,j) = Te(i,j)/11609.d0  ! convert to eV

      end do
      end do

! Use the routine Ecalc to calculate the ambipolar electric field.

! set the stride for fitting to log(n) in Ecalc
      stride=6

      call Ecalc(num_s_cells_r,num_s_cells_z,s_cells_dr,s_cells_dz,logdens,Te,Er,Ez,work,stride)

! note that Er and Ez are only loaded on processor 0 at this point


!kluge
write(*,*) 'Proc 0 value of Er(30,75) ',Er(30,75)

    else

      ! not steady yet, or ambipolar_switch=0, zero the ambipolar electric field
      do i=1,num_s_cells_r
      do j=1,num_s_cells_z
        Er(i,j) = 0.d0
        Ez(i,j) = 0.d0
      end do
      end do

    end if
!*****bottom of ambipolar if block *************************


  end if  ! bottom of mpi_rank=0 if block
!************* bottom of mpi_rank=0 if block

! let everyone catch up, just in cast it's needed
    call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

! now send Er and Ez to every processor on their own version of the sampling grid
    call MPI_BCAST(Er,num_s_cells_r*num_s_cells_z,MPI_DOUBLE_PRECISION, &
                  0,MPI_COMM_WORLD,mpi_ierr)


    call MPI_BCAST(Ez,num_s_cells_r*num_s_cells_z,MPI_DOUBLE_PRECISION, &
                  0,MPI_COMM_WORLD,mpi_ierr)


! now run over the collision cells and load the electric fields from the
! corresponding sampling cell

do n=1,num_cells

   ! find the center of the cell
   r=.5*(cells(i)%min_r+cells(i)%max_r )
   z=.5*(cells(i)%min_z+cells(i)%max_z )

   ! now figure out the (i,j) indices in the sampling cells that this
   ! collision cell is in

   i=(r-rmin_global)/s_cells_dr + 1
   j=(z-zmin_global)/s_cells_dz + 1

! assign each collision cell the E-values of the sampling cell they are in. A bit crude.
! also underrelax so that we get some extra averaging
   under = 1.d0  ! if 1, don't underrelax
   cells(i)%efield_r = Er(i,j)*under + (1.d0-under)*cells(i)%efield_r
   cells(i)%efield_z = Ez(i,j)*under + (1.d0-under)*cells(i)%efield_z

   ! don't put electric field in cells that are completely inside metal


    if(cells(i)%mask.eq.-1) then
      cells(i)%efield_r = 0.d0
      cells(i)%efield_z = 0.d0
    end if


end do




     deallocate(particle_info_all)
     deallocate(particle_info)
     deallocate(newlogdens)
     deallocate(Te)
     deallocate(work)

end subroutine write_collected_data


  !----- update_std_out -----
  ! Tell the user what's going on
  subroutine update_std_out()
    use parallel_info_mod
    use simulation_mod
    use core_data_mod
    implicit none
    integer num_parts(1)
    integer num_trace_parts(1)
    integer, allocatable, dimension(:) :: recv_array
    integer, allocatable, dimension(:) :: recv_array2
    character(10) ctime
    ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'
    integer status_struct(MPI_STATUS_SIZE)

    allocate(recv_array(num_mpi_processes))
    allocate(recv_array2(num_mpi_processes))

    ! tell process 0 how many NORMAL particles there are ("partition" is the number of normal particles)
    num_parts(1) = sum(cells(:)%partition)
    num_trace_parts(1) = sum(cells(:)%num_parts - cells(:)%partition)

!   Check to see who has reached this point
!   write(*,*) 'Argon parts - ', num_parts, 'Trace parts - ', num_trace_parts, &
!                             'mpi rank', mpi_rank
    call MPI_GATHER(num_parts, 1, MPI_INTEGER, &
                    recv_array, 1, MPI_INTEGER, &
                    0, MPI_COMM_WORLD, mpi_ierr)

    call MPI_GATHER(num_trace_parts, 1, MPI_INTEGER, &
                    recv_array2, 1, MPI_INTEGER, &
                    0, MPI_COMM_WORLD, mpi_ierr)

    ! only actually write to stdout if rank 0
    if(mpi_rank .eq. 0)then
      write(*,*)'Finished ',cur_step,' of ',nsteps,' steps.'
      call date_and_time(time=ctime)
      write(*,*)'Current Time: ', ctime
      write(*,*)'Num Normal Particles:',sum(recv_array)
      write(*,*)'Num Trace Particles:',sum(recv_array2)

      ! writing to the "particle_count.txt" file -- allows for an
      ! easy read into MATLAB to track the number of particles
      ! in the simulation
      write(partcount_unit,"(I14,3x,I14,3x,I14)") cur_step, sum(recv_array), sum(recv_array2)
    end if

    ! deallocating the receive arrays
    deallocate(recv_array)
    deallocate(recv_array2)

  end subroutine update_std_out



  subroutine write_particles_per_processor()
    use parallel_info_mod
    use simulation_mod
    use core_data_mod
    implicit none

    integer maxparts,i,imax
    real(8) totparts


    ! this is the array (on each processor) for helping the user better
    ! balance the processors in the axial index j
    real(8), allocatable, dimension(:) :: Narray,Narray_all,dens
    integer, allocatable, dimension(:) :: top,bottom,newtop
    integer j,n,m,k,nr_block,nz_block,Nproc,Nz
    real(8) Ni,Ntotal,Ntarget,f

    ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'

    allocate(Narray(num_mpi_processes))

    if(mpi_rank.eq.0) then

       allocate(Narray_all(num_mpi_processes))
       allocate(top(num_mpi_processes))
       allocate(bottom(num_mpi_processes))
       allocate(dens(num_mpi_processes+1))
       allocate(newtop(num_mpi_processes))

       do i=1,num_mpi_processes
         Narray_all(i)=0.d0
         top(i)=0
         bottom(i)=0
         dens(i)=0
         newtop(i)=0

       end do


    end if

    do i=1,num_mpi_processes
      Narray(i)=0.d0
    end do


    totparts=0.d0
    maxparts=0

    do i=1,num_cells
      totparts=totparts + cells(i)%num_parts
      if(cells(i)%num_parts.gt.maxparts) then
         maxparts=cells(i)%num_parts
         imax=i
      end if
    end do


     write(*,"('Processor, Npart, Maxpart, big cell index: ',i4,2x,1pe12.4,2(2x,i10))") mpi_rank, totparts, maxparts,imax



     ! load Narray for this processor with its number of particles
     Narray(mpi_rank+1)=totparts

! let everybody catch up
    call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)


   ! add all of the Narray's from all the processors into Narray_all on processor 0
   call MPI_REDUCE(Narray,Narray_all,num_mpi_processes,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)




   ! Processor 0 will now use the information in the array Narray_all to give the user a suggestion about how
   ! to better arrange the processor boundaries in the axial direction (not tested for more than 1 radial segment)
   ! Fix this!

   !**************  processor 0 **********************
   if(mpi_rank.eq.0) then


   ! read input_data.txt and pull out the top indices in j

   open(unit=83,file='input_data.txt',status='old')


      ! skip 9 lines
      do i=1,9
         read(83,*)
      end do

      ! read the processor block information
      read(83,*) nr_block,nz_block

      ! read the radial block information, and skip it
      do i=1,nr_block
       read(83,*)
      end do

      do i=1,nz_block
        read(83,*) top(i)
      end do

      close(unit=83)


      ! this code is borrowed from ProcessorBalance.m



      ! find the total number of particles
      Ntotal=sum(Narray_all)
      Nproc = num_mpi_processes

      bottom(2:Nproc) = top(1:Nproc-1)+1;  ! bottom index of each processor
      bottom(1)=1;
      Nz=top(Nproc) !  number of z (j) sampling cells

      ! for each processor, calculate the particle density in units
      ! of particles/sampling cell
      do i=1,Nproc
         dens(i)=Narray_all(i)/(top(i)-bottom(i)) !not really the density (denominator too small by 1)
                                                  !but the algorithm seems to work better this way
      end do
      ! extend it by one using the last density
      dens(Nproc+1) = dens(Nproc)


      ! find the target particles/processor
      Ntarget = Ntotal/Nproc;

      ! add particles until you exceed Ntarget, building a new top
      ! array as you go


      ! run over the processors and figure out, for each one, what its
      ! top j-index should be
      j=0  ! start the j-index
      Ni=0  ! start the approximate particles in a processor array

      do n=1,Nproc

          do while(Ni.le.Ntarget)
              j=j+1
              Ni=Ni+dens(n)
          end do

          ! we have gone over the limit; assign newtop and reset Ni

          ! be careful not to always err on the long side - you counted all
          ! of the particles in the last j as being in the previous processor

          f=(Ni-Ntarget)/Ntarget
          newtop(n)=j
          Ni=-f*dens(n+1)


      end do


      ! if newtop(Nproc) isn't Nz, stretch the indices

      f=real(Nz)/real(newtop(Nproc))

      write(*,*) '****************** suggested axial processor divisions ***********'
      do i=1,Nproc
         newtop(i) = f*newtop(i)+0.5d0 ! round
         write(*,*) newtop(i)
      end do
      write(*,*) '****************** suggested axial processor divisions ***********'



   end if
   !**************  processor 0 **********************

    call  MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)


    deallocate(Narray)

    if(mpi_rank.eq.0) then

       deallocate(Narray_all)
       deallocate(top)
       deallocate(bottom)
       deallocate(dens)
       deallocate(newtop)

    end if

    
  end subroutine write_particles_per_processor


subroutine Ecalc(nr,nz,dr,dz,logn,Te,Er,Ez,work,stride)

! this routine interpolates in (r,z) in integer form (i,j)
! to find the radial and axial derivatives of log(n) to
! determine the ambipolar electric field:

! E = -Te(eV) grad(log(n))

use ieee_arithmetic


implicit none
integer nr,nz,i,j,m,n,k
integer stride,imin,imax,jmin,jmax
integer msk(100)
real(8) dr,dz,y(100),fmax,x(100)
real(8) logn(nr,nz),Er(nr,nz),Ez(nr,nz),Te(nr,nz),work(nr,nz)
real(8) S0,Sx,Sxx,Sxxx,Sxxxx,R0,Rx,Rxx,b,c,slope,denom,denom2
real(8) x0,z0
integer nzero


fmax=0.d0
do i=1,nr
do j=1,nz
   fmax=max(fmax,logn(i,j))
end do
end do



! run over the grid
do j=1,nz
   do i=1,nr

   ! don't find the field at points where log(n)=0
   if(abs(logn(i,j)).lt.1.d-5) then
      Er(i,j)=0.d0
      Ez(i,j)=0.d0
      cycle
   end if

   ! Start radial gradient section


    ! choose imin and imax on the radial grid for doing the fit
      imin=i-stride
      imax=i+stride


      ! set the value of x at which to evaluate the derivative
      x0 = (i-0.5d0)**2  ! normal values, center of the data set
                         ! if the data set is shifted, this is still
                         ! the x-value where the derivative should be
                         ! evaluated

    ! shift imin and imax upward if imin is below 1
      if(imin.lt.1) then
         imin=1
         imax=1+2*stride
      end if

    ! shift imin and imax downward if imax>nr
      if(imax.gt.nr) then
         imin=nr-2*stride
         imax=nr
      end if


      ! now run over the points in the fitting set and load x=r^2/dr^2 and y=log(n)
      do m=imin,imax

         x(m-imin+1)=(m-0.5)**2  ! the x-values to use, squared almost integers
         y(m-imin+1)=logn(m,j)  ! set a string of log(n) values to fit to

         msk(m-imin+1)=0 ! set a mask value for the data being fitted: 0, density is too low
         if(y(m-imin+1)/fmax.gt.1.e-3)  msk(m-imin+1)=1 ! include this point in the fit

      end do


      ! if points in msk are zero at either end then we have hit an inner conductor and we should shift the data values
      nzero = 1+2*stride-sum(msk) ! number of zero msk points



      ! bottom shift
      if(msk(1).eq.0) then
        do m=imin+nzero,imax+nzero
         x(m-imin-nzero+1)=(m-0.5)**2
         y(m-imin-nzero+1)=logn(m,j) ! shift the data
        end do
      end if

      ! top shift
      if(msk(1+2*stride).eq.0) then
        do m=imin-nzero,imax-nzero
         x(m-imin+nzero+1)=(m-.5)**2
         y(m-imin+nzero+1)=logn(m,j) ! shift the data
        end do
      end if

      ! do the least-square sums
      S0=0.d0
      Sx=0.d0
      Sxx=0.d0
      Sxxx=0.d0
      Sxxxx=0.d0
      R0=0.d0
      Rx=0.d0
      Rxx=0.d0

      do m=1,1+2*stride
         S0=S0+1.d0
         Sx=Sx+x(m)
         Sxx=Sxx+x(m)**2
         Sxxx=Sxxx+x(m)**3
         Sxxxx=Sxxxx+x(m)**4
         R0=R0+y(m)
         Rx=Rx+y(m)*x(m)
         Rxx=Rxx+y(m)*x(m)**2
      end do


      ! find the least-squares slope at the radial point i by solving the
      ! least-squares set of equations for b and c
      !

      ! a + b*Sx + c*Sxx = R0
      ! a*Sx + b*Sxx + c*Sxxx = R0x
      ! a*Sxx + b*Sxxx + c*Sxxxx = R0xx

         denom = Sxx*(S0*Sxxxx+2.d0*Sxxx*Sx-Sxx**2)-S0*Sxxx**2-Sx**2*Sxxxx

         b = -(R0*(Sx*Sxxxx-Sxxx*Sxx)+Rx*(Sxx**2-S0*Sxxxx)+Rxx*(S0*Sxxx-Sx*Sxx))/denom



         c = -(R0*(Sxx**2-Sxxx*Sx)+Rx*(S0*Sxxx-Sxx*Sx)+Rxx*(Sx**2-S0*Sxx))/denom


! kluge check out a case

! if(i.eq.1  .and.j.eq.78) then
!    write(*,*) 'x ',(x(k),k=1,1+2*stride)
!    write(*,*) 'y ',(y(k),k=1,1+2*stride)
!    write(*,*) 'msk ',(msk(k),k=1,1+2*stride)
!    write(*,*) 'b,c ',b,c
!    stop 666
! end if


      ! convert to radial derivative
       slope = 2.d0*b*sqrt(x0)/dr  + 4.d0*c*x0**1.5/dr

      Er(i,j) = -Te(i,j)*slope
! trouble at bad spots occasionally - fix them
      if(ieee_is_nan(Er(i,j)).or.abs(Er(i,j)).gt.1.d15) then
        Er(i,j)=0.d0
      end if




!   write(*,103) i,slope,logn(i,j)
 103 format(' row: ',i4,'  slope, log(n) ' ,2(1x,1pe12.4))



   ! End radial gradient section




   ! Start axial gradient section


    ! choose jmin and jmax on the radial grid for doing the fit
      jmin=j-stride
      jmax=j+stride


      ! set the value of z at which to evaluate the derivative
      z0 = (j-0.5d0)     ! normal values, center of the data set;
                         ! if the data set is shifted, this is still
                         ! the x-value where the derivative should be
                         ! evaluated

    ! shift jmin and jmax upward if jmin is below 1
      if(jmin.lt.1) then
         jmin=1
         jmax=1+2*stride
      end if

    ! shift jmin and jmax downward if jmax>nz
      if(jmax.gt.nz) then
         jmin=nz-2*stride
         jmax=nz
      end if


      ! now run over the points in the fitting set and load x=z/dz and y=log(n)
      do m=jmin,jmax

         x(m-jmin+1)=(m-0.5)    ! the x-values to use, half-integers
         y(m-jmin+1)=logn(i,m)  ! set a string of log(n) values to fit to

         msk(m-jmin+1)=0 ! set a mask value for the data being fitted: 0, density is too low
         if(y(m-jmin+1)/fmax.gt.1.e-3)  msk(m-jmin+1)=1 ! include this point in the fit

      end do


      ! if points in msk are zero at either end then we have hit an inner conductor and we should shift the data values
      nzero = 1+2*stride-sum(msk) ! number of zero msk points


      ! bottom shift
      if(msk(1).eq.0) then
        do m=jmin+nzero,jmax+nzero
         x(m-jmin-nzero+1)=(m-0.5)
         y(m-jmin-nzero+1)=logn(i,m) ! shift the data
        end do
      end if

      ! top shift
      if(msk(1+2*stride).eq.0) then
        do m=jmin-nzero,jmax-nzero
         x(m-jmin+nzero+1)=(m-.5)
         y(m-jmin+nzero+1)=logn(i,m) ! shift the data
        end do
      end if


      ! do the least-square sums
      S0=0.d0
      Sx=0.d0
      Sxx=0.d0
      Sxxx=0.d0
      Sxxxx=0.d0
      R0=0.d0
      Rx=0.d0
      Rxx=0.d0

      do m=1,1+2*stride
         S0=S0+1.d0
         Sx=Sx+x(m)
         Sxx=Sxx+x(m)**2
         Sxxx=Sxxx+x(m)**3
         Sxxxx=Sxxxx+x(m)**4
         R0=R0+y(m)
         Rx=Rx+y(m)*x(m)
         Rxx=Rxx+y(m)*x(m)**2
      end do

      ! find the least-squares slope at the axial point j by solving the
      ! least-squares set of equations for b and c
      !

      ! a + b*Sx + c*Sxx = R0
      ! a*Sx + b*Sxx + c*Sxxx = R0x
      ! a*Sxx + b*Sxxx + c*Sxxxx = R0xx

         denom = Sxx*(S0*Sxxxx+2.d0*Sxxx*Sx-Sxx**2)-S0*Sxxx**2-Sx**2*Sxxxx

         b = -(R0*(Sx*Sxxxx-Sxxx*Sxx)+Rx*(Sxx**2-S0*Sxxxx)+Rxx*(S0*Sxxx-Sx*Sxx))/denom



         c = -(R0*(Sxx**2-Sxxx*Sx)+Rx*(S0*Sxxx-Sxx*Sx)+Rxx*(Sx**2-S0*Sxx))/denom


      ! convert to axial derivative
       slope = b/dz  + 2.d0*c*z0/dz

      Ez(i,j) = -Te(i,j)*slope
! trouble at bad spots occasionally - fix them
      if(ieee_is_nan(Ez(i,j)).or.abs(Ez(i,j)).gt.1.d15) then
        Ez(i,j)=0.d0
      end if

! testing
! if(i.eq.166) then
!      write(*,*) 'x ',x
!      write(*,*) 'y ',y
!      write(*,*) 'msk ',msk
!      write(*,*) 'jmin,jmax',jmin,jmax
!      write(*,*) 'Ez ',Ez(i,j)
!      stop 666
! end if




!   write(*,104) j,slope,logn(i,j),c,b
 104 format(' j: ',i4,'  slope, log(n) c b',    4(1x,1pe12.4))



   ! End axial gradient section

! kluge
if(abs(Ez(i,j)).gt.1e15) then
   write(*,*) '$$',i,j,logn(i,j),Ez(i,j)
end if




   end do
!  write(*,*) 'End of i= ',i
end do
! end of loop over grid

! Because we are differentiating noisy data this is a bit rough
! run over the interior points and do a smoothing pass

do i=2,nr-1
do j=2,nz-1

   work(i,j) = 0.25d0*Er(i,j) + 0.125d0*(Er(i,j+1)+Er(i,j-1)+Er(i+1,j)+Er(i-1,j) ) &
                              + 0.0625*(Er(i+1,j+1)+Er(i+1,j-1)+Er(i-1,j-1)+Er(i-1,j+1) )
end do
end do


do i=2,nr-1
do j=2,nz-1

   Er(i,j) = work(i,j)

end do
end do


do i=2,nr-1
do j=2,nz-1

   work(i,j) = 0.25d0*Ez(i,j) + 0.125d0*(Ez(i,j+1)+Ez(i,j-1)+Ez(i+1,j)+Ez(i-1,j) ) &
                              + 0.0625*(Ez(i+1,j+1)+Ez(i+1,j-1)+Ez(i-1,j-1)+Ez(i-1,j+1) )
end do
end do


do i=2,nr-1
do j=2,nz-1

   Ez(i,j) = work(i,j)

end do
end do

! dump the ambipolar field to a file for use later
! in loading the ambipolar field components into collstuff.dat

open(unit=511,file='Eambipolar.txt',status='unknown')
do j=1,nz
do i=1,nr

   write(511,105) i,j,Er(i,j),Ez(i,j),logn(i,j),Te(i,j)
105 format(2(i5),4(1pe16.8))

end do
end do

close(unit=511,status='keep')




return

end subroutine Ecalc
  



  subroutine fixlogdens(newlogdens,num_s_cells_r,num_s_cells_z)
  implicit none
  integer num_s_cells_r,num_s_cells_z
  real(8) newlogdens(num_s_cells_r,num_s_cells_z)
  integer di,i,j,k,m,n
  real(8) x(100),a11,a12,a22,b1,b2,denom,a,b,chk

 ! log(n) has bad statistics at r=0 sometime, fix them by fitting to r.^2, linear

          do j=1,num_s_cells_z

              di=5  ! correct the first di-1 points, if needed

              do  k=1,di-1

                  i=di-k

                  do n=1,di
                   x(n)=i+n-.5  ! grid points beyond i, di of them
                  end do

                 ! least squares fit to these points, a + b*(i-.5)^2
                  a11=0.d0
                  a12=0.d0
                  a22=0.d0
                  b1=0.d0
                  b2=0.d0


                  do  n=1,di
                     a11=a11+1.d0
                     a12=a12+x(n)**2
                     a22=a22+x(n)**4
                     b1=b1+newlogdens(i+n,j)
                     b2=b2+newlogdens(i+n,j)*x(n)**2
                  end do

                  ! solve for a and b
                  denom=a11*a22-a12**2
                  a=(a22*b1-a12*b2)/denom
                  b=(a11*b2-a12*b1)/denom

                  chk=abs(1-(a+b*(i-.5)**2)/newlogdens(i,j))  ! check the error on log(n) at i


                  ! if log(n(i)) isn't close to the fit to its neighbors, replace it with the fit
                  if(chk.gt.0.01) then
                      newlogdens(i,j)=a+ b*(i-.5)**2
                  end if

              end do

          end do


  end subroutine fixlogdens

  
  
end module dataio_mod





