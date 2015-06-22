
module inflow_mod
    use constants_mod
    use simulation_mod
    implicit none

    integer :: nr, nz
    real(8) :: drblk, dzblk, zleft, Nphys, Nsim, rbottom, zshift
    real(8) :: facr,facz
    real(8) :: w, a, b, n_not, floor !! used for calculating "ndens_tr_grid".
    real(8), pointer, dimension(:,:) ::     ndens_grid
    real(8), pointer, dimension(:,:) ::  ndens_tr_grid
    real(8), pointer, dimension(:,:) ::         T_grid
    real(8), pointer, dimension(:,:) :: Tparallel_grid
    real(8), pointer, dimension(:,:) ::     Tperp_grid
    real(8), pointer, dimension(:,:) ::        vx_grid
    real(8), pointer, dimension(:,:) ::        vy_grid
    real(8), pointer, dimension(:,:) ::        vz_grid
    real(8), pointer, dimension(:,:) ::         r_grid
    real(8), pointer, dimension(:,:) ::         z_grid

    integer barf_unit
    parameter(barf_unit=131)

    contains

    subroutine init_inflow(mpi_rank)
        implicit none

        integer i,j,mpi_rank

        real(8) navg,iavg

        ! importing data from which "barf file" data is generated

        ! open "grid.txt"
        open(unit=barf_unit,file="inflowdata/grid.txt",action="read")

        ! read in grid information, cell-edge grid for bilinear interpolation
        read(barf_unit,*) rbottom
        read(barf_unit,*) zleft
        read(barf_unit,*) nr
        read(barf_unit,*) nz
        read(barf_unit,*) drblk
        read(barf_unit,*) dzblk

        ! zshift is the width of the load region in z
        zshift = (nz-1)*dzblk

        close(unit=barf_unit)

        ! allocate arrays
        allocate(     ndens_grid(nr, nz) )
        allocate(  ndens_tr_grid(nr, nz) )
        allocate(         T_grid(nr, nz) )
        allocate( Tparallel_grid(nr, nz) )
        allocate(     Tperp_grid(nr, nz) )
        allocate(        vx_grid(nr, nz) )
        allocate(        vy_grid(nr, nz) )
        allocate(        vz_grid(nr, nz) )
        allocate(         r_grid(nr, nz) )
        allocate(         z_grid(nr, nz) )

 
        ! ngrid.txt
        open(unit=barf_unit,file="inflowdata/ngrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (ndens_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! Tgrid.txt
        open(unit=barf_unit,file="inflowdata/Tgrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (T_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! Tparallelgrid.txt
        open(unit=barf_unit,file="inflowdata/Tparallelgrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (Tparallel_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! Tperpgrid.txt
        open(unit=barf_unit,file="inflowdata/Tperpgrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (Tperp_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! vxgrid.txt
        open(unit=barf_unit,file="inflowdata/vxgrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (vx_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! vygrid.txt
        open(unit=barf_unit,file="inflowdata/vygrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (vy_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! vzgrid.txt
        open(unit=barf_unit,file="inflowdata/vzgrid.txt",action="read")
        do i = 1,nr
            read(barf_unit,*) (vz_grid(i,j),j=1,nz)
        end do
        close(unit=barf_unit)

        ! r_grid, z_grid
        do j = 1,nz
            do i = 1,nr
                ! load the cell edge grid points
                !
                !
                r_grid(i,j) = i*drblk - drblk
                z_grid(i,j) = (j*dzblk - dzblk) - zshift + zmin_global
            end do
        end do

        Nphys = 0
        navg=0.d0
        iavg=0.d0
        ! perform a trapezoid integration over the load region
        do i = 1,nr

            if(i.eq.1.or.i.eq.nr) then
               facr=.5d0
            else
               facr=1.d0
            end if

            do j = 1,nz


               if(j.eq.1.or.j.eq.nz) then
                  facz=.5d0
               else
                  facz=1.d0
               end if

                Nphys = Nphys + facr*facz*ndens_grid(i,j)*drblk*dzblk*r_grid(i,j)*2.d0*pi  !I'm taking wedge_width off this, and replacing it with 2*pi

                iavg = iavg + i  ! average the density over the load grid, cylindrical weighting
                navg = navg + ndens_grid(i,j)*i
            end do
        end do

        navg = navg/iavg

        Nsim = Nphys/Nef

        ! ndens_tr_grid
        w = 0.15d-3
        a = 0.45d0
        b = 0.8d0
        floor = .1d0

        ! set the central trace density to the average of the argon density over
        ! ndens_grid

       	n_not = 5*navg    ! set the trace density artificially high because we
                        ! treat them as a trace species, i.e., no self collisions
                        ! only collisions with argon. This means that the trace
                        ! density is arbitrary, and we choose it high to get good
                        ! statistics
        do i = 1,nr
            do j = 1,nz
!               ndens_tr_grid(i,j) = n_not*( exp( -a*(r_grid(i,j)/w)**2 - b*(r_grid(i,j)/w)**4 ) + floor )/(1.d0+floor)
! decent fit to Haibin's data at 10 mm downstream
                ndens_tr_grid(i,j) = n_not*( exp(-3.0d6*r_grid(i,j)**1.7) + .20 )/1.20
            end do
        end do


    end subroutine init_inflow



    ! -- fill_inflow_cell --
    ! builds a particle in a random location in the given cell dimension, and 
    ! gives the particle characteristics according to the "inflowdata" files
    !  ** This replaces the "barf file" functionality.
    subroutine fill_inflow_cell(c)
        use core_types_mod
        use simulation_mod
        use constants_mod
        use helpers_mod
        implicit none

        type(cell_type) c
        type(particle_type) p0
        integer i,j,k
        integer Nparts,Nloaded
        real(8) Nparts_real, frac, rnd, phi, vth, dblint
        real(8) topdens, intrp, xtry, ztry, r1, r2, r12, facr2, fr, fz
        real(8) intrpvx, intrpvy, intrpvz

        frac = 0.0

        !open(unit=1234,file='inflowstat.txt',status='unknown')

        ! the load grid is cell-edge, so there is one less cell in each dimension
        ! than there are points

        do i = 1,nr-1
            do j = 1,nz-1

                r1=r_grid(i,j)
                r12=r1**2
                r2=r1+drblk
                facr2=r2**2-r1**2

                ! load the maximum corner density
                topdens = max( ndens_grid(i,j),ndens_grid(i+1,j),ndens_grid(i,j+1),ndens_grid(i+1,j+1) )


                ! figure out how many simulation particles to put into this section of the cell

                dblint = 0.25*( ndens_grid(i,j)*r_grid(i,j)+ndens_grid(i+1,j)*r_grid(i+1,j) + &
                               ndens_grid(i,j+1)*r_grid(i,j+1)+ndens_grid(i+1,j+1)*r_grid(i+1,j+1) )

                !Nparts = (wedge_width*dblint*drblk*dzblk) / Nphys * Nsim + frac
                Nparts = (2.d0*pi*dblint*drblk*dzblk) / Nphys * Nsim + frac

                frac = (Nparts - int(Nparts))
                Nparts = int(Nparts)

                ! Putting the particles into this section of the cell according to the 
                ! distribution specified in the inflow data files.

                Nloaded=0
                ! write(*,*) 'i,j,Nparts',i,j,Nparts
                do while ( Nloaded<Nparts )

                    ! putting in the particles

                    ! propose a particle position

                    call rand_mt(fr)
                    xtry = sqrt(r12 + fr*facr2) ! r^2 weighting
                    call rand_mt(fz)
                    ztry = z_grid(i,j) + fz*dzblk	

                    ! interpolate the density at this position
                    intrp = fr*fz*ndens_grid(i+1,j+1)+fr*(1.d0-fz)*ndens_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*ndens_grid(i,j)+(1.d0-fr)*fz*ndens_grid(i,j+1)

                    ! use acceptance-rejection to match the particle load to the desired density
                    call rand_mt(rnd)
                    if(intrp < rnd*topdens) cycle

                    ! if you get to here you are loading a particle
                    Nloaded=Nloaded+1
                    
                    p0%x = xtry
                    p0%y = 0.
                    p0%z = ztry
                    p0%r = sqrt(p0%x**2+p0%z**2)

                    ! load the velocities

                    ! interpolate the temperature at this position
                    intrp = fr*fz*T_grid(i+1,j+1)+fr*(1.d0-fz)*T_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*T_grid(i,j)+(1.d0-fr)*fz*T_grid(i,j+1)

                    vth = sqrt(kb*intrp/m)


                    ! interpolate velocities at this position
                    intrpvx = fr*fz*vx_grid(i+1,j+1)+fr*(1.d0-fz)*vx_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vx_grid(i,j)+(1.d0-fr)*fz*vx_grid(i,j+1)
                    intrpvy = fr*fz*vy_grid(i+1,j+1)+fr*(1.d0-fz)*vy_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vy_grid(i,j)+(1.d0-fr)*fz*vy_grid(i,j+1)
                    intrpvz = fr*fz*vz_grid(i+1,j+1)+fr*(1.d0-fz)*vz_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vz_grid(i,j)+(1.d0-fr)*fz*vz_grid(i,j+1)

                    call rand_mt(rnd)
                    call rand_mt(phi)
                    phi = 2.d0*pi*phi
                    p0%vx = vth*sqrt(-2.d0*log(rnd))*cos(phi) + intrpvx
                    p0%vy = vth*sqrt(-2.d0*log(rnd))*sin(phi) + intrpvy

                    call rand_mt(rnd)
                    call rand_mt(phi)
                    phi = 2.d0*pi*phi
                    p0%vz = vth*sqrt(-2.d0*log(rnd))*cos(phi) + intrpvz

                    p0%element = 0

                    call cell_insert_particle(c,p0)

                    !write(*,'(a3,6(1pe12.4))') '## ',p0%x, p0%y, p0%z, p0%vx, p0%vy, p0%vz

            
                end do
            end do
        end do

        !close(unit=1234)

    end subroutine fill_inflow_cell



    ! -- fill_inflow_cell_trace
    ! builds a particle in a random location in the given cell dimension, and 
    ! gives the particle characteristics according to the formula used in
    ! init_inflow
    subroutine fill_inflow_cell_trace(c)
        use core_types_mod
        use simulation_mod
        use constants_mod
        use helpers_mod
        implicit none

        type(cell_type) c
        type(particle_type) p0
        integer i,j,k
        integer Nparts,Nloaded
        real(8) Nparts_real, frac, rnd, phi, vth, dblint
        real(8) topdens, intrp, xtry, ztry, r1, r2, r12, facr2, fr, fz
        real(8) intrpvx, intrpvy, intrpvz

        frac = 0.0

        !open(unit=1235,file='inflowstat_trace.txt',status='unknown')

        ! the load grid is cell-edge, so there is one less cell in each dimension
        ! than there are points

        do i = 1,nr-1
            do j = 1,nz-1

                r1=r_grid(i,j)
                r12=r1**2
                r2=r1+drblk
                facr2=r2**2-r1**2

                ! load the maximum corner density
                topdens = max( ndens_tr_grid(i,j),ndens_tr_grid(i+1,j),ndens_tr_grid(i,j+1),ndens_tr_grid(i+1,j+1) )


                ! figure out how many simulation particles to put into this section of the cell

                dblint = 0.25*( ndens_tr_grid(i,j)*r_grid(i,j)+ndens_tr_grid(i+1,j)*r_grid(i+1,j) + &
                               ndens_tr_grid(i,j+1)*r_grid(i,j+1)+ndens_tr_grid(i+1,j+1)*r_grid(i+1,j+1) )

                !Nparts = (wedge_width*dblint*drblk*dzblk) / Nphys * Nsim + frac
                Nparts = (2.d0*pi*dblint*drblk*dzblk) / Nphys * Nsim + frac

                frac = (Nparts - int(Nparts))
                Nparts = int(Nparts)

                ! Putting the particles into this section of the cell according to the 
                ! distribution specified in the inflow data files.

                Nloaded=0
                ! write(*,*) 'i,j,Nparts',i,j,Nparts
                do while ( Nloaded<Nparts )

                    ! putting in the particles

                    ! propose a particle position

                    call rand_mt(fr)
                    xtry = sqrt(r12 + fr*facr2) ! r^2 weighting
                    call rand_mt(fz)
                    ztry = z_grid(i,j) + fz*dzblk	

                    ! interpolate the density at this position
                    intrp = fr*fz*ndens_tr_grid(i+1,j+1)+fr*(1.d0-fz)*ndens_tr_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*ndens_tr_grid(i,j)+(1.d0-fr)*fz*ndens_tr_grid(i,j+1)

                    ! use acceptance-rejection to match the particle load to the desired density
                    call rand_mt(rnd)
                    if(intrp < rnd*topdens) cycle

                    ! if you get to here you are loading a particle
                    Nloaded=Nloaded+1
                    
                    p0%x = xtry
                    p0%y = 0.
                    p0%z = ztry
                    p0%r = sqrt(p0%x**2+p0%z**2)

                    ! load the velocities using the argon values - assuming thermal
                    ! equilibrium and trace flow velocity locked to argon flow velocity

                    ! interpolate the temperature at this position
                    intrp = fr*fz*T_grid(i+1,j+1)+fr*(1.d0-fz)*T_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*T_grid(i,j)+(1.d0-fr)*fz*T_grid(i,j+1)

                    vth = sqrt(kb*intrp/tr_m)


                    ! interpolate velocities at this position
                    intrpvx = fr*fz*vx_grid(i+1,j+1)+fr*(1.d0-fz)*vx_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vx_grid(i,j)+(1.d0-fr)*fz*vx_grid(i,j+1)
                    intrpvy = fr*fz*vy_grid(i+1,j+1)+fr*(1.d0-fz)*vy_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vy_grid(i,j)+(1.d0-fr)*fz*vy_grid(i,j+1)
                    intrpvz = fr*fz*vz_grid(i+1,j+1)+fr*(1.d0-fz)*vz_grid(i+1,j)+ &
                                (1.d0-fr)*(1.d0-fz)*vz_grid(i,j)+(1.d0-fr)*fz*vz_grid(i,j+1)

                    call rand_mt(rnd)
                    call rand_mt(phi)
                    phi = 2.d0*pi*phi
                    p0%vx = vth*sqrt(-2.d0*log(rnd))*cos(phi) + intrpvx
                    p0%vy = vth*sqrt(-2.d0*log(rnd))*sin(phi) + intrpvy

                    call rand_mt(rnd)
                    call rand_mt(phi)
                    phi = 2.d0*pi*phi
                    p0%vz = vth*sqrt(-2.d0*log(rnd))*cos(phi) + intrpvz

                    p0%element = 1

                    call cell_insert_particle(c,p0)

                    !write(*,'(a3,6(1pe12.4))') '## ',p0%x, p0%y, p0%z, p0%vx, p0%vy, p0%vz

            
                end do
            end do
        end do

        !close(unit=1235)

    end subroutine fill_inflow_cell_trace





end module inflow_mod


