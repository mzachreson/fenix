






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
  ! This subroutine collects the data from all of the processors, and writes it to a file
  subroutine write_collected_data_old
    use constants_mod
    use simulation_mod
    use parallel_info_mod
    use core_data_mod
    implicit none

    real(8) avx,avy,avz
    real(8) avx2,avy2,avz2
    real(8) avxvy,avxvz,avyvz
    real(8) an

    real(8) tr_avx,tr_avy,tr_avz
    real(8) tr_avx2,tr_avy2,tr_avz2
    real(8) tr_avxvy,tr_avxvz,tr_avyvz
    real(8) tr_an

    real(8) r, z, V

    integer i,j,j_c,k,l, count_s_cell, s_i, s_j
    integer diffi,diffj
    real(8),allocatable,dimension(:) :: send_buf
    real(8),allocatable,dimension(:) :: recv_buf
    integer size_buf
    integer recv_counts(num_mpi_processes)
    integer displs(num_mpi_processes)
  ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'
    size_buf = 20*(max_i()-min_i()+1)*(max_j()-min_j()+1)
    allocate(send_buf(size_buf))

    k = 1
    do j=min_j(),max_j()
      do i=min_i(),max_i()
        ! normal argon particles data
        send_buf(20*k-19) = s_cells(i,j)%sum_vx
        send_buf(20*k-18) = s_cells(i,j)%sum_vy
        send_buf(20*k-17) = s_cells(i,j)%sum_vz
        send_buf(20*k-16) = s_cells(i,j)%sum_vx2
        send_buf(20*k-15) = s_cells(i,j)%sum_vy2
        send_buf(20*k-14) = s_cells(i,j)%sum_vz2
        send_buf(20*k-13) = s_cells(i,j)%sum_vxvy
        send_buf(20*k-12) = s_cells(i,j)%sum_vxvz
        send_buf(20*k-11) = s_cells(i,j)%sum_vyvz
        send_buf(20*k-10) = s_cells(i,j)%sum_n

        ! trace particle data
        send_buf(20*k-9) = s_cells(i,j)%tr_sum_vx
        send_buf(20*k-8) = s_cells(i,j)%tr_sum_vy
        send_buf(20*k-7) = s_cells(i,j)%tr_sum_vz
        send_buf(20*k-6) = s_cells(i,j)%tr_sum_vx2
        send_buf(20*k-5) = s_cells(i,j)%tr_sum_vy2
        send_buf(20*k-4) = s_cells(i,j)%tr_sum_vz2
        send_buf(20*k-3) = s_cells(i,j)%tr_sum_vxvy
        send_buf(20*k-2) = s_cells(i,j)%tr_sum_vxvz
        send_buf(20*k-1) = s_cells(i,j)%tr_sum_vyvz
        send_buf(20*k-0) = s_cells(i,j)%tr_sum_n
        k=k+1
      end do
    end do

    ! allocate a receive buffer
    allocate(recv_buf(20*num_s_cells_r*num_s_cells_z))

    ! set up the amount of gathering from each process
    do i=0,num_mpi_processes-1
      recv_counts(i+1) = 20*(max_ip(i)-min_ip(i)+1)*(max_jp(i)-min_jp(i)+1)
      if(i.eq.0)then
        displs(i+1) = 0
      else
        displs(i+1) = displs(i)+recv_counts(i)
      end if
    end do

    ! gather all of the collected data
    call MPI_ALLGATHERV(send_buf,size_buf,MPI_DOUBLE_PRECISION, &
                  recv_buf,recv_counts,displs,MPI_DOUBLE_PRECISION, &
                  MPI_COMM_WORLD,mpi_ierr)

    ! free up the send buffer
    deallocate(send_buf)

    ! collect the data into the sampling cells (all processes)
    do j_c=0,num_proc_z-1
      diffj = max_jp(j_c)-min_jp(j_c)
      do j=0,diffj
        do k=0,num_proc_r-1
          diffi = max_ip(num_proc_z*k) - min_ip(num_proc_z*k)
          do l=0,diffi
            i = 1 + (displs(num_proc_z*k+j_c+1)/20) + (1+diffi)*j + l
            s_i = min_ip(num_proc_z*k) + l
            s_j = min_jp(j_c) + j
            !if(mpi_rank .eq. 0)then
            !  write(*,*)i,s_i,s_j
            !end if

            ! normal argon particles data
            s_cells(s_i,s_j)%sum_vx = recv_buf(20*i-19)
            s_cells(s_i,s_j)%sum_vy = recv_buf(20*i-18)
            s_cells(s_i,s_j)%sum_vz = recv_buf(20*i-17)
            s_cells(s_i,s_j)%sum_vx2 = recv_buf(20*i-16)
            s_cells(s_i,s_j)%sum_vy2 = recv_buf(20*i-15)
            s_cells(s_i,s_j)%sum_vz2 = recv_buf(20*i-14)
            s_cells(s_i,s_j)%sum_vxvy = recv_buf(20*i-13)
            s_cells(s_i,s_j)%sum_vxvz = recv_buf(20*i-12)
            s_cells(s_i,s_j)%sum_vyvz = recv_buf(20*i-11)
            s_cells(s_i,s_j)%sum_n = recv_buf(20*i-10)

            ! trace particle data
            s_cells(s_i,s_j)%tr_sum_vx = recv_buf(20*i-9)
            s_cells(s_i,s_j)%tr_sum_vy = recv_buf(20*i-8)
            s_cells(s_i,s_j)%tr_sum_vz = recv_buf(20*i-7)
            s_cells(s_i,s_j)%tr_sum_vx2 = recv_buf(20*i-6)
            s_cells(s_i,s_j)%tr_sum_vy2 = recv_buf(20*i-5)
            s_cells(s_i,s_j)%tr_sum_vz2 = recv_buf(20*i-4)
            s_cells(s_i,s_j)%tr_sum_vxvy = recv_buf(20*i-3)
            s_cells(s_i,s_j)%tr_sum_vxvz = recv_buf(20*i-2)
            s_cells(s_i,s_j)%tr_sum_vyvz = recv_buf(20*i-1)
            s_cells(s_i,s_j)%tr_sum_n = recv_buf(20*i-0)

          end do
        end do
      end do
    end do

    ! deallocate the receive buffer
    deallocate(recv_buf)

    !only write the data if you're the root
    if(mpi_rank .eq. 0)then
      do j=1,num_s_cells_z
        do i=1,num_s_cells_r

          ! normal particles averaging
          avx = s_cells(i,j)%sum_vx / s_cells(i,j)%sum_n
          avy = s_cells(i,j)%sum_vy / s_cells(i,j)%sum_n
          avz = s_cells(i,j)%sum_vz / s_cells(i,j)%sum_n
          avx2 = s_cells(i,j)%sum_vx2 / s_cells(i,j)%sum_n
          avy2 = s_cells(i,j)%sum_vy2 / s_cells(i,j)%sum_n
          avz2 = s_cells(i,j)%sum_vz2 / s_cells(i,j)%sum_n
          avxvy = s_cells(i,j)%sum_vxvy/s_cells(i,j)%sum_n
          avxvz = s_cells(i,j)%sum_vxvz/s_cells(i,j)%sum_n
          avyvz = s_cells(i,j)%sum_vyvz/s_cells(i,j)%sum_n

          ! trace particles averaging
          tr_avx = s_cells(i,j)%tr_sum_vx / s_cells(i,j)%tr_sum_n
          tr_avy = s_cells(i,j)%tr_sum_vy / s_cells(i,j)%tr_sum_n
          tr_avz = s_cells(i,j)%tr_sum_vz / s_cells(i,j)%tr_sum_n
          tr_avx2 = s_cells(i,j)%tr_sum_vx2 / s_cells(i,j)%tr_sum_n
          tr_avy2 = s_cells(i,j)%tr_sum_vy2 / s_cells(i,j)%tr_sum_n
          tr_avz2 = s_cells(i,j)%tr_sum_vz2 / s_cells(i,j)%tr_sum_n
          tr_avxvy = s_cells(i,j)%tr_sum_vxvy/s_cells(i,j)%tr_sum_n
          tr_avxvz = s_cells(i,j)%tr_sum_vxvz/s_cells(i,j)%tr_sum_n
          tr_avyvz = s_cells(i,j)%tr_sum_vyvz/s_cells(i,j)%tr_sum_n

          r = (i*s_cells_dr + (i-1)*s_cells_dr)/2.d0
          z = (j*s_cells_dz + (j-1)*s_cells_dz)/2.d0

          !V = (s_cells_dz*pi*((i*s_cells_dr)**2 - ((i-1)*s_cells_dr)**2))

          ! "s_cells(i,j)%volume" is a corrected volume, taking into account
          ! intersection with the metal.  It was read in from "samplingcells.dat"
          V = s_cells(i,j)%volume

          an = Nef * s_cells(i,j)%sum_n / num_samples / V
          tr_an = Nef * s_cells(i,j)%tr_sum_n / num_samples / V

          if(s_cells(i,j)%sum_n .eq. 0)then
            avx = 0
            avy = 0
            avz = 0
            avx2 = 0
            avy2 = 0
            avz2 = 0
            avxvy = 0
            avxvz = 0
            avyvz = 0
          end if

          if(s_cells(i,j)%tr_sum_n .eq. 0)then
            tr_avx = 0
            tr_avy = 0
            tr_avz = 0
            tr_avx2 = 0
            tr_avy2 = 0
            tr_avz2 = 0
            tr_avxvy = 0
            tr_avxvz = 0
            tr_avyvz = 0
          end if


          write(output_unit,"(13(1pe12.4))") &
                          an, &
                          avx2, avy2, avz2, &
                          avxvy, avxvz, avyvz, &
                          avx, avy, avz, &
                          r, z, V

          write(outputtrace_unit,"(13(1pe12.4))") &
                          tr_an, &
                          tr_avx2, tr_avy2, tr_avz2, &
                          tr_avxvy, tr_avxvz, tr_avyvz, &
                          tr_avx, tr_avy, tr_avz, &
                          r, z, V


        end do
      end do
    end if
  end subroutine write_collected_data_old


  !----- write_collected_data2
  ! This is an attempt to clean up data output by changing the form of data 
  !on each processor to a summable form and then using MPI_
  subroutine write_collected_data
    use constants_mod
    use simulation_mod
    use parallel_info_mod
    use core_data_mod
    implicit none

    real(8) r, z, V

    integer i,j,k,l
    real(8),allocatable,dimension(:,:,:) :: particle_info
    real(8),allocatable,dimension(:,:,:) :: particle_info_all
  ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'

!allocating arrays to dump averages into
! everybody needs particle_info because they send their data in it
    allocate(particle_info(20,num_s_cells_r,num_s_cells_z))    

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
!Bringing all of the info onto the first processor
   call MPI_REDUCE(particle_info,particle_info_all,20*num_s_cells_r*num_s_cells_z,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
!Writing the data to output.txt and output_trace.txt
 ! write(*,*) 'Writing to file',mpi_rank
  if (mpi_rank == 0) then
     do j=1,num_s_cells_z
     do i=1,num_s_cells_r    
          !Loading Geometry info
          r = (i*s_cells_dr + (i-1)*s_cells_dr)/2.d0
          z = (j*s_cells_dz + (j-1)*s_cells_dz)/2.d0

          ! "s_cells(i,j)%volume" is a corrected volume, taking into account
          ! intersection with the metal.  It was read in from "samplingcells.dat"
          V = s_cells(i,j)%volume

          !Averaging cell data and writing it for Argon
          if (particle_info_all(1,i,j)/=0) then
             particle_info(1,i,j) = Nef * particle_info_all(1,i,j) / num_samples / V
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
          if (particle_info_all(11,i,j)/=0) then
             particle_info(11,i,j) = Nef * particle_info_all(11,i,j) / num_samples / V
             do k=12,20
                particle_info(k,i,j)=particle_info_all(k,i,j) / particle_info_all(11,i,j)
             end do
          else
             do k=11,10
                particle_info(k,i,j)=0
             end do
          end if
          write(outputtrace_unit,"(13(1pe12.4))") (particle_info(k,i,j),k=11,20),r,z,V
     end do
     end do


  end if

     deallocate(particle_info_all)
deallocate(particle_info)

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

    integer totparts,maxparts,i,imax

    totparts=0
    maxparts=0

    do i=1,num_cells
      totparts=totparts + cells(i)%num_parts
      if(cells(i)%num_parts.gt.maxparts) then
         maxparts=cells(i)%num_parts
         imax=i
      end if
    end do


!    write(*,"('Processor, Npart, Maxpart, big cell index: ',4(2x,I12))") mpi_rank, totparts, maxparts,imax

    
  end subroutine write_particles_per_processor
  
  
  
  
end module dataio_mod





