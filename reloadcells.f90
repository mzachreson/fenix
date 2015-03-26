





module reloadcells_mod
  implicit none
	
  contains




  !----- reload_cells -----
  ! This subroutine reloads all of the cells
  subroutine reload_cells
    use core_data_mod
    use bounds_mod
    use parallel_info_mod
    use simulation_mod
    implicit none
    integer i
  
    ! zero out the communication cells
    comm_cell_left%num_parts = 0
    comm_cell_right%num_parts = 0
    comm_cell_top%num_parts = 0
    comm_cell_bottom%num_parts = 0

    comm_cell_left%partition = 0
    comm_cell_right%partition = 0
    comm_cell_top%partition = 0
    comm_cell_bottom%partition = 0
  
    !reload all the normal cells
    do i=1,num_cells

      call reload_cell(cells(i))
    end do

    !reload the ghost cells
    do i=1,num_gc
      call reload_ghost_cell(gc(i))
    end do
  end subroutine reload_cells




  !----- reload_cell -----
  ! This subroutine reloads a single collision cell
  subroutine reload_cell(cl)
  	use initialize_mod !for reference to "insert_into_correct_cell"
    use core_types_mod
    implicit none
    type(cell_type) cl
    type(particle_type) p

    do while (cl%num_left.gt.0)

      ! Load the particle into variable "p" and adjust "num_left" appropriately
      p = cl%ps(cl%capacity-cl%num_left+1)
      cl%num_left = cl%num_left - 1

      ! insert_into_correct_cell also puts stuff into the communication cells, as needed
      call insert_into_correct_cell(p)

    end do
  end subroutine reload_cell




  !----- reload_ghost_cell -----
  ! This subroutine loads all particles which have moved out of this ghost cell
  ! into neighboring cells
  subroutine reload_ghost_cell(cl)
    use core_types_mod
    use core_data_mod
    use parallel_info_mod
    use helpers_mod
    use constants_mod
    use simulation_mod
    implicit none
    type(cell_type) cl
    type(particle_type) p
    integer sc_i,sc_j,fc_i,fc_j,i

    !write(*,*)mpi_rank,'Welcome to ghost cell reload.',cl%num_left,cl%num_parts

    ! Cycle through all particles on the far right of the array (right of the index "capacity - num_left")
    ! These are the particles which need to be loaded into a different cell than they are in now.
    do while (cl%num_left.gt.0)

      ! Load the particle into variable "p" and adjust "num_left" appropriately
      p = cl%ps(cl%capacity-cl%num_left+1)
      cl%num_left = cl%num_left - 1
    
      ! ignore the ones that are not in the simulation region
      sc_i = 1 + (p%r - rmin_global)/s_cells_dr
      sc_j = 1 + (p%z - zmin_global)/s_cells_dz
      if(sc_j .lt. min_j() .or. &
         sc_j .gt. max_j()) cycle
      if(sc_i .lt. min_i() .or. &
         sc_i .gt. max_i()) cycle

      ! continue
      fc_i = 1 + (p%r - rmin_local)/fine_cells_dr
      fc_j = 1 + (p%z - zmin_local)/fine_cells_dz
      if( fc_i .ge. 1 .and. fc_i .le. num_fine_cells_r .and. &
          fc_j .ge. 1 .and. fc_j .le. num_fine_cells_z )then
        i = fine_cells(fc_i,fc_j)
        if(cells(i)%mask .eq. 0)then
          call cell_insert_particle( cells(i), p )
          !write(*,*)mpi_rank,"insert successful."
        else
          if(.not.is_inside(p%r, p%z, masks(cells(i)%mask)))then
            call cell_insert_particle( cells(i), p )
            !write(*,*)mpi_rank,"insert successful."
          end if
        end if
      else
        !write(*,*)mpi_rank,"Messed-up insert."
      end if
    end do
  end subroutine reload_ghost_cell



  !----- send_recv_parts -----
  ! This subroutine sends and receives particles from other processors
  subroutine send_recv_parts
    use parallel_info_mod
    use core_data_mod
    use initialize_mod !for reference to "insert_into_correct_cell"
    implicit none
    integer num_to_send
    integer bytes_to_send
    integer dest, i, j
    type(particle_type), allocatable, dimension(:) :: recv_buffer
    integer bytes_to_receive,num_to_receive
    type(particle_type) p
  ! To use MPI (Message Passing Interface) I have to include this header file
    include 'mpif.h'
    integer status_struct(MPI_STATUS_SIZE)
    integer, save :: size_particle !the size of a particle when in an array
    integer neighbors(4)

    !write(*,*)mpi_rank,"In send_recv_parts."
    if(size_particle .eq. 0)then
	   ! The number of bytes a single particle takes in memory.
       size_particle = sizeof(comm_cell_left%ps(:))/comm_cell_left%capacity
    end if
  
    ! send the data to the left
    if(mod(mpi_rank,num_proc_z) .ne. 0)then
      num_to_send = comm_cell_left%num_parts
      bytes_to_send = sizeof(comm_cell_left%ps(1:num_to_send))

      ! send the data
      call MPI_BSEND(comm_cell_left%ps(1:num_to_send), bytes_to_send, MPI_BYTE, &
                     mpi_rank-1, 0, MPI_COMM_WORLD, mpi_ierr)
      comm_cell_left%num_parts = 0
      comm_cell_left%partition = 0
      neighbors(1) = mpi_rank-1
    else
      neighbors(1) = -1
    end if

    ! send the data to the right
    if(mod(mpi_rank+1,num_proc_z) .ne. 0)then
      num_to_send = comm_cell_right%num_parts
      bytes_to_send = sizeof(comm_cell_right%ps(1:num_to_send))

      ! send the data
      call MPI_BSEND(comm_cell_right%ps(1:num_to_send), bytes_to_send, MPI_BYTE, &
                     mpi_rank+1, 0, MPI_COMM_WORLD, mpi_ierr)
      comm_cell_right%num_parts = 0
      comm_cell_right%partition = 0
      neighbors(2) = mpi_rank+1
    else
      neighbors(2) = -1
    end if
  
    ! send the data to the top
    if(mpi_rank .lt. (num_proc_r-1)*num_proc_z)then
      num_to_send = comm_cell_top%num_parts
      bytes_to_send = sizeof(comm_cell_top%ps(1:num_to_send))

      ! send the data
      call MPI_BSEND(comm_cell_top%ps(1:num_to_send), bytes_to_send, MPI_BYTE, &
                     mpi_rank+num_proc_z, 0, MPI_COMM_WORLD, mpi_ierr)
      comm_cell_top%num_parts = 0
      comm_cell_top%partition = 0
      neighbors(3) = mpi_rank+num_proc_z
    else
      neighbors(3) = -1
    end if
  
    ! send the data to the bottom
    if(mpi_rank .ge. num_proc_z)then
      num_to_send = comm_cell_bottom%num_parts
      bytes_to_send = sizeof(comm_cell_bottom%ps(1:num_to_send))

      ! send the data
      call MPI_BSEND(comm_cell_bottom%ps(1:num_to_send), bytes_to_send, &
                     MPI_BYTE, mpi_rank-num_proc_z, 0,MPI_COMM_WORLD, mpi_ierr)
      comm_cell_bottom%num_parts = 0
      comm_cell_bottom%partition = 0
      neighbors(4) = mpi_rank-num_proc_z
    else
      neighbors(4) = -1
    end if

    !write(*,*)mpi_rank,"About to get particles."


    ! now I can receive from the neighbors
    do i=1,4
      if(neighbors(i) .eq. -1) cycle

      ! blocking probe call
     ! write(*,*)mpi_rank,"About to probe particles."
     ! write(*,*) 'processor, i, neighbors(i)',mpi_rank,i,neighbors(i)
      call MPI_PROBE(neighbors(i), 0, MPI_COMM_WORLD,status_struct, mpi_ierr)
     !write(*,*) 'status_struct,mpi_ierr,processor',status_struct,mpi_ierr,mpi_rank
      !write(*,*)mpi_rank,"Probed."
      ! find out how many to receive
      call MPI_GET_COUNT(status_struct,MPI_BYTE,bytes_to_receive, mpi_ierr)
      ! allocate space
      num_to_receive = bytes_to_receive/size_particle
      if(num_to_receive*size_particle .ne. bytes_to_receive)then
        write(*,*) 'Discrepancy.',num_to_receive,size_particle,bytes_to_receive
        stop
      end if
      allocate(recv_buffer(num_to_receive))
      ! receive them
      call MPI_RECV(recv_buffer,bytes_to_receive,MPI_BYTE, neighbors(i), 0,MPI_COMM_WORLD, status_struct, mpi_ierr)
  
      ! now put each one in the cell it belongs to
      do j=1,num_to_receive
        call insert_into_correct_cell(recv_buffer(j))
      end do

      !throw away the buffer
      deallocate(recv_buffer)
    end do
  end subroutine send_recv_parts
  



end module reloadcells_mod




