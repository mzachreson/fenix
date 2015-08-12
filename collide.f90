


!----- collide_mod -----
! This fortran module defines the routine for collisions
module collide_mod
implicit none
contains

  !----- collide -----
  ! this subroutine makes all of the collision cells collide their particles
  subroutine collide
    use core_data_mod
    use simulation_mod
    use constants_mod
    use parallel_info_mod
!kluge   use parallel_info_mod
    use parallel_info_mod
    implicit none
    integer i
    real(8) tot_coll,tot_cand
! kluge
!   integer iflag
!   iflag=0


! watch candidates vs. collisions
    tot_coll=0.d0
    tot_cand=0.d0

    do i=1,num_cells

! kluge look for cells that have z=3.3.89e-4 in them
!     if(cells(i)%min_z.lt.3.89e-4.and.cells(i)%max_z.gt.3.89e-4) then
!        write(*,*) mpi_rank,i,cells(i)%min_r,cells(i)%max_r
!        iflag = 1
!     end if

      ! do normal particle collisions
      call cell_collide_subdivide(cells(i),tot_coll,tot_cand)
!     write(*,*) '#$# ',mpi_rank,tot_cand,tot_coll


      ! do trace particle collisions

      call trace_collide(cells(i),i)
      ! do trace particle self collisions
      if(nanbu_switch) call trace_self_collide(cells(i))
    end do


!kluge
!      write(*,101) mpi_rank,tot_cand/tot_coll,tot_cand,tot_coll
101   format(' Collide:  tot_cand/tot_coll: ',i5,3(1pe12.4))


  end subroutine collide






  subroutine cell_collide_subdivide(cl,tot_coll,tot_cand)
    use helpers_mod
    use core_types_mod
    use constants_mod
    use simulation_mod

    implicit none
    type(cell_type) cl
    integer num_fracs, j, kmin, kmax
    real(8) zmin, zmax, frac, dn,sigmavr,vexp,sigmavrmax_local

!kluge
    integer num_parts,jj,k,i,num_coll,Npairs,num_cand
    real(8) tot_coll,tot_cand
    real(8) vrel2,avgsigmavr,coll_target


  ! if this is the first time there are 2 or more particles in this cell,
  ! run over the collision pairs and find the maximum value of sigma*vr
    vexp = 0.5d0-nu


    if(cl%firstflag.eq.1.and.cl%partition.ge.2) then
	cl%firstflag=0
        cl%sigmavrmax=0.d0
        do i=1,cl%partition
         do j=i+1,cl%partition

           vrel2 =  (cl%ps(j)%vx-cl%ps(i)%vx)**2 + &
                    (cl%ps(j)%vy-cl%ps(i)%vy)**2 + &
                    (cl%ps(j)%vz-cl%ps(i)%vz)**2
           sigmavr = aref*vrel2**vexp
           cl%sigmavrmax = max(cl%sigmavrmax,sigmavr)

         end do
       end do
   end if

!test  !kluge: average sigma*vr by running over all of the pairs
!test  ! in just one collision cell chosen by finding the one with a particular point in it
!test
!test      vexp=0.5d0-nu
!test
!test        num_parts=cl%partition
!test
!test        num_coll=0
!test        if(cl%min_r.lt.1.e-3.and.cl%max_r.gt.1.d-3.and.cl%min_z.lt.1.d-3.and.cl%max_z.gt.1.d-3 ) then
!test          Npairs=cl%partition*(cl%partition-1)/2
!test
!test          if(Npairs.gt.0) then
!test              avgsigmavr=0.d0
!test              do jj=1,num_parts
!test              do i=jj+1,num_parts
!test                       k = i
!test                       j = jj
!test
!test                       vrel2 =  (cl%ps(j)%vx-cl%ps(k)%vx)**2 + &
!test                       (cl%ps(j)%vy-cl%ps(k)%vy)**2 + &
!test                       (cl%ps(j)%vz-cl%ps(k)%vz)**2
!test                     ! write(*,*) 'j,k,vrel2',j,k,vrel2
!test                       avgsigmavr = avgsigmavr + aref*vrel2**vexp
!test
!test              end do
!test              end do
!test
!test            ! write(*,*) 'aref_tr,nu_tr1,bref_tr,nu_tr2', aref_tr,nu_tr1,bref_tr,nu_tr2
!test
!test              avgsigmavr = avgsigmavr/Npairs
!test              write(*,403) cl%partition,avgsigmavr,cl%volume,Nef
!test 403         format(' Nparts,avgsigmavr,cl%volume,Nef',i5,3(1x,1pe12.4))
!test              coll_target = 0.5d0*Nef*cl%partition*(cl%partition-1)/cl%volume*avgsigmavr*tau
!test              write(*,*) 'Test cell, collision target: ',coll_target
!test
!test
!test          end if
!test
!test        end if
!test  !kluge

    ! zero out all of the previous collision partners
    cl%ps(:)%prev_coll = 0

    ! control growth of sigmavrmax
    cl%sigmavrmax = cl%sigmavrmax *0.99

    ! set sigmavrmax_local
    sigmavrmax_local = cl%sigmavrmax

    if(cl%partition .le. 1) return

    ! I want to make sure that each sub-collision cell has at least 7 particles
    ! Bird says 7 is great; 20 is better (correct number of collisions per cell)
    num_fracs = cl%partition / 20.d0

    if(num_fracs .le. 1)then
      call cell_collide(cl, 1.d0, 1, cl%partition,num_coll,num_cand,sigmavrmax_local)
!testif(cl%min_r.lt.1.e-3.and.cl%max_r.gt.1.d-3.and.cl%min_z.lt.1.d-3.and.cl%max_z.gt.1.d-3 ) then
!test      tot_coll=tot_coll+num_coll
!test      tot_cand=tot_cand+num_cand
!testend if
    else
      call sortZ(cl%ps, 1, cl%partition)
      dn = dble(cl%partition)/dble(num_fracs)
      do j=1,num_fracs
        kmin = nint((j-1)*dn)+1
        kmax = nint(j*dn)
        ! find zmin
        if(j .eq. 1)then
          zmin = cl%min_z
        else
          zmin = (cl%ps(kmin-1)%z + cl%ps(kmin)%z)/2.d0
        end if
        ! find zmax
        if(j .eq. num_fracs)then
          zmax = cl%max_z
        else
          zmax = (cl%ps(kmax)%z + cl%ps(kmax+1)%z)/2.d0
        end if
        ! collide!
        frac = (zmax-zmin)/(cl%max_z - cl%min_z)
        call cell_collide(cl, frac, kmin, kmax,num_coll,num_cand,sigmavrmax_local)

      if(cl%min_r.lt.1.e-3.and.cl%max_r.gt.1.d-3.and.cl%min_z.lt.1.d-3.and.cl%max_z.gt.1.d-3 ) then
              tot_coll=tot_coll+num_coll
              tot_cand=tot_cand+num_cand
      end if
      end do
    end if


    ! update the cell value of sigmavrmax here after all of the subdivisions of the cell
    ! have been collided
    cl%sigmavrmax = sigmavrmax_local



!test            if(cl%min_r.lt.1.e-3.and.cl%max_r.gt.1.d-3.and.cl%min_z.lt.1.d-3.and.cl%max_z.gt.1.d-3 ) then
!test             write(*,*) 'No. of argon collisions and target: ',tot_coll,coll_target
!test            end if



  end subroutine cell_collide_subdivide

  !--- cell_collide (VSS) ---
  ! This subroutine collides the molecules in a given cell
  ! according to the variable soft sphere collider
  subroutine cell_collide(cl, frac, n1, n2,num_coll,num_cand,sigmavrmax_local)  ! collides argon with argon
    use simulation_mod
    use core_types_mod
    use constants_mod
    use helpers_mod
    implicit none
    type(cell_type) cl
    real(8) frac
    integer n1, n2
    integer num_parts

    real(8) rnum_candidates
    integer num_candidates,num_cand
    integer i,j,k,l,jj
    real(8) dist2, min_dist2
    real(8) temp,vrel2,vrel
    real(8) phi,cos_theta,sin_theta,cos_phi,sin_phi,dvz_sin_phi,vrel_cos_phi
    real(8) dvx,dvy,dvz
    real(8) vrx,vry,vrz
    real(8) vperp,vperpi
    real(8) vrx_cm,vry_cm,vrz_cm
    real(8) cmx,cmy,cmz
    integer nearswitch,num_coll
    real(8) alpha,vexp,rexp
    real(8) sigmavr,Npairs,avgsigmavr,coll_target,sigmavrmax_local


    data alpha/1.66/
    vexp=0.5d0-nu
    rexp=1.d0/alpha

    num_parts = n2-n1+1



    num_coll=0 ! zero this before processing collisions

    ! compute number of candidates
    num_candidates = round_rand(Nef*tau*0.5d0*cl%sigmavrmax*num_parts*(num_parts-1.d0)/(frac*cl%volume)) !&
                    !+ cl%num_cand_remaining
    !num_candidates = int(rnum_candidates)

    !cl%num_cand_remaining = rnum_candidates - num_candidates

    ! main collide loop
    num_coll=0
    num_cand=num_candidates
    do i=1,num_candidates
      ! select an atom
      call rand_mt(temp)
      j = num_parts*temp*.9999999+n1

      ! select another atom, different from the first
      k=0
      min_dist2=1.d99
      do l=n1,n2
        if(l .ne. j .and. l .ne. cl%ps(j)%prev_coll)then
          dist2 = (cl%ps(l)%x - cl%ps(j)%x)**2
          if(dist2 .lt. min_dist2)then
            dist2 = dist2 + (cl%ps(l)%y - cl%ps(j)%y)**2 &
                          + (cl%ps(l)%z - cl%ps(j)%z)**2
            if(dist2 .lt. min_dist2)then
              min_dist2 = dist2
              k = l
            end if
          end if
        end if
      end do
      if(k .eq. 0)then
        cycle ! if there are no possible partners, skip this candidate.
      end if

      !if(j.gt.n2 .or. k.gt.n2)then
      !  write(*,*) 'bad random numbers in cell collider'
      !  stop
      !end if

      ! acceptance-rejection
      vrel2 =  (cl%ps(j)%vx-cl%ps(k)%vx)**2 + &
               (cl%ps(j)%vy-cl%ps(k)%vy)**2 + &
               (cl%ps(j)%vz-cl%ps(k)%vz)**2

! get sigma*vr and push the maximum up if necessary
      sigmavr = aref*vrel2**vexp
      sigmavrmax_local = max(sigmavr,sigmavrmax_local)

      call rand_mt(temp)
      ! if rejecting this candidate, call 'cycle'
      if(sigmavr  .le. temp*cl%sigmavrmax)then
        cycle
      end if
      ! this candidate is accepted

      num_coll=num_coll+1
      !update the collision partners
      cl%ps(j)%prev_coll = k
      cl%ps(k)%prev_coll = j


      ! compute the scattering angles in the CM frame
      call rand_mt(temp)
      phi = 2*pi*temp
      cos_phi=cos(phi)
      sin_phi=sin(phi)

      call rand_mt(temp)
      cos_theta = 2.d0*temp**rexp - 1.d0
      sin_theta = sqrt(1-cos_theta**2)

      if(isnan(sin_theta))then
        write(*,*) 'sin_theta is NaN in collider'
        stop 172
      end if

      !compute center of mass velocity
      cmx = (cl%ps(j)%vx + cl%ps(k)%vx)*0.5d0
      cmy = (cl%ps(j)%vy + cl%ps(k)%vy)*0.5d0
      cmz = (cl%ps(j)%vz + cl%ps(k)%vz)*0.5d0

      ! compute the relative velocity
      dvx = cl%ps(j)%vx - cl%ps(k)%vx
      dvy = cl%ps(j)%vy - cl%ps(k)%vy
      dvz = cl%ps(j)%vz - cl%ps(k)%vz
      vrel = sqrt(vrel2)

      dvz_sin_phi = dvz*sin_phi
      vrel_cos_phi = vrel*cos_phi

      !compute the new relative velocity components
      !with z in the direction of the original relative velocity vector
      vrx_cm = vrel*sin_theta*cos_phi
      vry_cm = vrel*sin_theta*sin_phi
      vrz_cm = vrel*cos_theta

      ! now transform this scattered velocity back to the frame of the simulation

      vperp = sqrt(dvx**2 + dvy**2)
      vperpi = 1.d0/vperp
      if(vperp .gt. 0.d0)then
        vrx = cos_theta*dvx + sin_theta*vperpi*(vrel_cos_phi*dvy-dvz_sin_phi*dvx)
        vry = cos_theta*dvy - sin_theta*vperpi*(vrel_cos_phi*dvx+dvz_sin_phi*dvy)
      else
        vrx = sin_theta*dvz*cos_phi
        vry = sin_theta*dvz*sin_phi
      end if
      vrz = cos_theta*dvz + sin_theta*sin_phi*vperp

!      if(isnan(vrx) .or. isnan(vry) .or. isnan(vrz))then
!        write(*,*)'collision caused NaN.'
!        write(*,*)'vperp=',vperp
!        stop 117
!      end if

      !now set the velocities
      cl%ps(j)%vx = cmx - vrx*0.5d0
      cl%ps(j)%vy = cmy - vry*0.5d0
      cl%ps(j)%vz = cmz - vrz*0.5d0
      cl%ps(k)%vx = cmx + vrx*0.5d0
      cl%ps(k)%vy = cmy + vry*0.5d0
      cl%ps(k)%vz = cmz + vrz*0.5d0

      if(isnan(cl%ps(j)%vx) .or. isnan(cl%ps(j)%vy) .or. isnan(cl%ps(j)%vz) )then
        write(*,*)'Found Nan.'
        write(*,*)cmx,cmy
        write(*,*)vrx,vry
        write(*,*)dvx,dvy
        stop 123
      end if

    ! end of main collide loop
    end do



 ! kluge watch a cell at z=1 mm, r=1 mm - see how many collisions there are
!test if(cl%min_r.lt.1d-3.and.cl%max_r.gt.1.e-3.and.cl%min_z.lt.1d-3.and.cl%max_z.gt.1.d-3) then
!test    write(*,303) num_parts,num_candidates,num_coll
!test    write(*,304) cl%sigmavrmax,frac,cl%volume
!test 303 format(' $new 303 ',3(1x,i5))
!test 304 format(' $new 304 ',3(1x,1pe12.4))
!test end if

    return
  end subroutine cell_collide




    ! collides all the trace-element particles
    ! which are in the ps array from index "partition+1" to "num_parts"
    subroutine trace_collide(cl,icell)  ! trace particles collide with argon, argon is unchanged
        use core_types_mod
        use simulation_mod
        use constants_mod
        use helpers_mod
    use parallel_info_mod

        implicit none
        integer icell,num_coll
        real(8) avgsigmavr,Npairs,coll_target,sigmavrmax_local

        type(cell_type) :: cl
        integer :: num_tr_parts, num_norm_parts

        real(8) rnum_candidates
        integer num_candidates
        integer i,j,k,l
        real(8) dist2, min_dist2
        real(8) temp,vrel2,vrel
        real(8) phi,cos_theta,sin_theta
        real(8) dvx,dvy,dvz
        real(8) vrx,vry,vrz
        real(8) vperp
        real(8) vrx_cm,vry_cm,vrz_cm
        real(8) cmx,cmy,cmz
        integer nearswitch
        real(8) alpha,vexp,rexp
        real(8) sigmavrmax,sigmavr,sigmamax,sigma
        real(8) vxold,vyold,vzold

! set the local value of sigmavrmax for the cell
        sigmavrmax_local=cl%tr_sigmavrmax

        num_coll=0

        num_tr_parts = cl%num_parts - cl%partition
        num_norm_parts = cl%partition


  ! if this is the first time there are 2 or more particles in this cell,
  ! run over the collision pairs and find the maximum value of sigma*vr
    if(cl%firstflag_tr.eq.1.and.num_tr_parts.gt.1.and.num_norm_parts.gt.1) then
	cl%firstflag_tr=0
        cl%tr_sigmavrmax=0.d0
             do j=1,num_tr_parts
             do i=1,num_norm_parts

                k = num_norm_parts + j
                vrel2 =  (cl%ps(i)%vx-cl%ps(k)%vx)**2 + &
                (cl%ps(i)%vy-cl%ps(k)%vy)**2 + &
                (cl%ps(i)%vz-cl%ps(k)%vz)**2
                cl%tr_sigmavrmax = max(cl%tr_sigmavrmax, aref_tr*vrel2**(.5d0-nu_tr1)+bref_tr*vrel2**(.5d0-nu_tr2) )

             end do
             end do


    end if




        data alpha/1.66/
        rexp=1.d0/alpha





!test   !kluge: average sigma*vr by running over all of the pairs
!test   ! in just one collision cell chosen with the if statement below
!test
!test         num_coll=0
!test         if(cl%min_r.lt.0.39d-3.and.cl%max_r.gt.0.39d-3.and.cl%min_z.lt.0.2d-3.and.cl%max_z.gt.0.2d-3) then
!test           Npairs=num_tr_parts*num_norm_parts
!test   !       write(*,*) 'num_tr_parts,num_norm_parts,Npairs:',num_tr_parts,num_norm_parts,Npairs
!test   !       write(*,*) 'rmin,rmax,zmin,zmax: ',cl%min_r,cl%max_r,cl%min_z,cl%max_z
!test   !       write(*,*) 'aref_ambi: ',cl%aref_ambi
!test
!test           if(Npairs.gt.0) then
!test               avgsigmavr=0.d0
!test               do j=1,num_tr_parts
!test               do i=1,num_norm_parts
!test                        k = num_norm_parts + j
!test                        vrel2 =  (cl%ps(j)%vx-cl%ps(k)%vx)**2 + &
!test                        (cl%ps(j)%vy-cl%ps(k)%vy)**2 + &
!test                        (cl%ps(j)%vz-cl%ps(k)%vz)**2
!test                      ! write(*,*) 'j,k,vrel2',j,k,vrel2
!test                        avgsigmavr = avgsigmavr +  aref_tr*vrel2**(.5d0-nu_tr1)+bref_tr*vrel2**(.5d0-nu_tr2)
!test
!test               end do
!test               end do
!test
!test             ! write(*,*) 'aref_tr,nu_tr1,bref_tr,nu_tr2', aref_tr,nu_tr1,bref_tr,nu_tr2
!test
!test               avgsigmavr = avgsigmavr/Npairs
!test             ! write(*,*) 'avgsigmavr,cl%volume,Nef',avgsigmavr,cl%volume,Nef
!test               coll_target = Nef*num_norm_parts/cl%volume*avgsigmavr*tau*num_tr_parts
!test             ! write(*,*) 'Test cell, collision target: ',coll_target
!test
!test
!test           end if
!test
!test         end if
!test   !kluge: average sigma*vr by running over all of the pairs

        ! too few particles to collide, so return
        if( num_norm_parts <= 1 .or. num_tr_parts <= 1 ) return

        ! control growth of sigmavrmax
        cl%tr_sigmavrmax = cl%tr_sigmavrmax *0.99

        ! compute number of candidates


! round_rand method, rounding routine allows small numbers to round up
        num_candidates = round_rand(Nef*tau*cl%tr_sigmavrmax*num_norm_parts*num_tr_parts/(cl%volume))

       ! kluge
          if(num_candidates.gt.10000) then
            write(*,*) 'Big candidates: ',num_candidates
            write(*,*) 'r and z',cl%max_r,cl%max_z
            write(*,*) 'processor ',mpi_rank
            write(*,*) 'cl%tr_sigmavrmax',cl%tr_sigmavrmax
          end if

! fraction method
!        rnum_candidates = (Nef*tau*cl%tr_sigmavrmax*num_norm_parts*num_tr_parts/(cl%volume)) + cl%num_cand_remaining
!        num_candidates = int(rnum_candidates)
!        cl%num_cand_remaining = rnum_candidates - num_candidates


         ! Questions about the above equation:
        !  1. Do we use "num_norm_parts" or "num_tr_parts", or both? Answer: normal for density, trace for the the number to collide
        !  2. Will we use a different version of "Nef", for the trace particles? Will we also use the normal particles' Nef?


        ! start of main collide loop


        do i = 1,num_candidates

            ! select an argon (normal) particle
            call rand_mt(temp)
            j = num_norm_parts*temp*0.999999 + 1

            ! select a trace atom
            call rand_mt(temp)
            k = num_norm_parts + num_tr_parts*temp*0.999999 + 1

            if(    j < 1             .or. j > cl%partition   &
              .or. k <= cl%partition .or. k > cl%num_parts   )then
              write(*,*) 'bad random numbers in cell TRACE collider'
              stop 298
            end if

            ! acceptance-rejection
            vrel2 =  (cl%ps(j)%vx-cl%ps(k)%vx)**2 + &
                     (cl%ps(j)%vy-cl%ps(k)%vy)**2 + &
                     (cl%ps(j)%vz-cl%ps(k)%vz)**2

! compute sigmavr and adjust the ceiling on tr_sigmavrmax
            !Two nu cross section, but avoid enormous cross sections if vr turns out to be too small
            !and sigma is singular at small vr

            sigma = (aref_tr/vrel2**nu_tr1+bref_tr/vrel2**nu_tr2)
            sigmamax = 1d-15
            sigma = min(sigma,sigmamax)

            sigmavr = sigma*sqrt(vrel2);
            sigmavrmax_local = max(sigmavr, sigmavrmax_local)

            call rand_mt(temp)
            ! if rejecting this candidate, call 'cycle'
            if(sigmavr .le. temp*cl%tr_sigmavrmax) then
                cycle
            end if

            !kluge
            num_coll=num_coll+1

            !give new speeds
            call rand_mt(temp)
            phi = 2*pi*temp
            call rand_mt(temp)
            cos_theta = 2.d0*temp**rexp - 1.d0
            sin_theta = sqrt(1-cos_theta**2)

            if(isnan(sin_theta))then
                write(*,*) 'sin_theta is NaN in collider'
                stop 137
            end if

            !compute center of mass speed
            cmx = (m*cl%ps(j)%vx + tr_m*cl%ps(k)%vx)/(m+tr_m)
            cmy = (m*cl%ps(j)%vy + tr_m*cl%ps(k)%vy)/(m+tr_m)
            cmz = (m*cl%ps(j)%vz + tr_m*cl%ps(k)%vz)/(m+tr_m)

            dvx = cl%ps(j)%vx - cl%ps(k)%vx
            dvy = cl%ps(j)%vy - cl%ps(k)%vy
            dvz = cl%ps(j)%vz - cl%ps(k)%vz
            vrel = sqrt(vrel2)
            !compute relative velocity in center-of-mass frame
            ! I let the z-axis be the relative velocity vector
            vrx_cm = vrel*sin_theta*cos(phi)
            vry_cm = vrel*sin_theta*sin(phi)
            vrz_cm = vrel*cos_theta
            ! now change coordinates back into world frame
            ! now back to real-world frame
            vperp = sqrt(dvx**2 + dvy**2)
            if(vperp .gt. 0.d0)then
                vrx = cos_theta*dvx + sin_theta/vperp*(vrel*dvy*cos(phi)-dvz*dvx*sin(phi))
                vry = cos_theta*dvy - sin_theta/vperp*(vrel*dvx*cos(phi)+dvz*dvy*sin(phi))
            else
                vrx = sin_theta*dvz*cos(phi)
                vry = sin_theta*dvz*sin(phi)
            end if
            vrz = cos_theta*dvz + sin_theta*sin(phi)*vperp
            if(isnan(vrx) .or. isnan(vry) .or. isnan(vrz))then
                write(*,*)'collision caused NaN.'
                write(*,*)'vperp=',vperp
                stop 143
            end if

            !now set the velocities
            cl%ps(k)%vx = cmx - m*vrx/(m+tr_m)
            cl%ps(k)%vy = cmy - m*vry/(m+tr_m)
            cl%ps(k)%vz = cmz - m*vrz/(m+tr_m)


            if(isnan(cl%ps(j)%vx) .or. isnan(cl%ps(j)%vy) .or. isnan(cl%ps(j)%vz) )then
                write(*,*)'Found Nan.'
                write(*,*)cmx,cmy
                write(*,*)vrx,vry
                write(*,*)dvx,dvy
                stop 161
            end if

        ! end of main collide loop
        end do

!kluge watch for num_candidates/num_coll too big

if(real(num_candidates)/real(num_coll).gt.10..and.num_coll.ne.0) then
   write(*,*) 'r,z,ratio',cl%min_r,cl%min_z,real(num_candidates)/real(num_coll)
end if


 ! kluge inspect the collision testing cell
!test       if(cl%min_r.lt.0.39d-3.and.cl%max_r.gt.0.39d-3.and.cl%min_z.lt.0.2d-3.and.cl%max_z.gt.0.2d-3) then
!test        write(*,*) 'No. of trace collisions and target: ',num_coll,coll_target
!test        write(*,*) 'candidates: ',num_candidates
!test       end if
 ! kluge inspect the collision testing cell


    ! update the cell value of sigmavrmax
    cl%tr_sigmavrmax = sigmavrmax_local

    end subroutine trace_collide


  subroutine trace_self_collide(cl)
    use helpers_mod
    use core_types_mod
    use constants_mod
    use simulation_mod
    use parallel_info_mod
    implicit none
    type(cell_type) cl
    integer :: num_parts, num_ions, num_electrons, num_part_types, Ncoll
    real(8) :: vthe, me
    real(8) :: dt, A, vrel(3), h(3), vperp, vrelmag, Aab, lD, mu
    real(8) :: taun !Nanbu's tau, not to be confused with the simulation timestep, tau
    real(8) :: vxe, vye, vze
    real(8) :: loglambda, psi, coschi, sinchi, exp2A, ndens
    integer tmpint1, tmpint2, i,j,k,l
    integer, allocatable, dimension(:) :: coll_pairs ! variable for generating random collision pairs
    real(8) tmp, tmp1, tmp2



    !Fix unloaded variables if this is the first time step
    if(cur_step.eq.1)then
       cl%densavg=0.0d0
       cl%tempavg=0.0d0
    end if

    !make sure there is an ion to collide
    num_ions=cl%num_parts-cl%partition
    if(num_ions.lt.1) return

    !Skip if temp and dens are not loaded
    if(cl%densavg.lt.1d-6.or.cl%tempavg.lt.1d-6) return


    !Find electron thermal velocity for random thermal electrons
    me=9.11d-31
    vthe=sqrt(kb*cl%Te/me)
    ndens=2d0*cl%densavg*1e3/Nef !gives the number density of all charged species, assuming ni=ne

   !Now do collisions
   !Make random collision pairs list

   allocate(coll_pairs(num_ions))

   do i=1,num_ions
      if (i.gt.num_ions) write(*,*) 'Coll_pairs on line 457'
      coll_pairs(i)=i+cl%partition
   end do

   !Randomize lists using the Knuth shuffle
      do i=1,num_ions-1
         l=num_ions+1-i
         call rand_mt(tmp)
         k=floor(l*tmp+1)
         if (l.gt.num_ions) write(*,*) 'Coll_pairs on line 467'

         if (k.gt.num_ions) write(*,*) 'Coll_pairs on line 468', k,l,i,tmp
         tmpint1=coll_pairs(l)
         coll_pairs(l)=coll_pairs(k)
         coll_pairs(k)=tmpint1

      end do


   ! Now do electron-ion collisions
   Ncoll=ceiling(num_ions/2d0)  !Number of interspecies collisions, rounded up
                                ! So that single ions collide with an electron




   lD=sqrt(e0*kb*cl%tempavg/cl%densavg/tr_q**2) !Debye length
   mu=tr_m*me/(tr_m+me)


       loglambda=log(4*pi*e0*mu*((3*kb*(cl%tempavg/tr_m+cl%Te/me)))*lD/tr_q**2)


       Aab=ndens*.25/pi*(tr_q**2/e0/mu)**2*loglambda


      if (Ncoll.gt.num_ions) write(*,*) 'Ncoll too big'
   do i=1,Ncoll
      k=coll_pairs(i) !ion number
      ! Now create a thermal electron to collide with
      call rand_mt(tmp1)
      call rand_mt(tmp2)
      tmp2=2.d0*pi*tmp2
      tmp=vthe*sqrt(-2.d0*log(tmp1))
      vxe=tmp*cos(tmp2)
      vye=tmp*sin(tmp2)
      call rand_mt(tmp2)
      tmp2=2*pi*tmp2
      vze=tmp*cos(tmp2)
 !  write(15,'(3(1pe12.4),I5)') ps(i)%vx, ps(i)%vy, ps(i)%vz, 1

      vrel(1)=cl%ps(k)%vx-vxe
      vrel(2)=cl%ps(k)%vy-vye
      vrel(3)=cl%ps(k)%vz-vze
      vperp=sqrt(vrel(2)**2+vrel(3)**2)
      vrelmag=sqrt(vperp**2+vrel(1)**2)

!write(*,*) 'Relative velocity', vrel(1), vrel(2), vrel(3), vperp, vrelmag
!write(*,*) 'Loop info', i,k,l
      taun=Aab*tau/vrelmag**3
      tmp=exp(-taun)
      A=1.d0/(1.d0-tmp)+(-2.1939+2.6742*tmp-1.2038*tmp**2+1.4485*tmp**3)*exp(-.7856d0/(1.d0-tmp))


! write(*,*) taun,loglambda,A,' %%'
      call rand_mt(tmp1)
      psi=tmp1*2.d0*pi


      h(1)=vperp*cos(psi)
      h(2)=-(vrel(1)*vrel(2)*cos(psi)+vrelmag*vrel(3)*sin(psi))/vperp
      h(3)=-(vrel(1)*vrel(3)*cos(psi)-vrelmag*vrel(2)*sin(psi))/vperp

! Do the Nanbu inversion to find cos(chi)
      call rand_mt(tmp1)
      exp2A = exp(-2.d0*A)
      coschi = 1.d0 + 1.d0/A*log(tmp1*(1-exp2A)+exp2A)
      sinchi = sqrt(1.d0-coschi**2)




      !set post collision velocities
      cl%ps(k)%vx=cl%ps(k)%vx+mu/tr_m*(vrel(1)*(1.d0-coschi)+h(1)*sinchi)
      cl%ps(k)%vy=cl%ps(k)%vy+mu/tr_m*(vrel(2)*(1.d0-coschi)+h(2)*sinchi)
      cl%ps(k)%vz=cl%ps(k)%vz+mu/tr_m*(vrel(3)*(1.d0-coschi)+h(3)*sinchi)

   if(cl%ps(k)%vy.ne.cl%ps(k)%vy)then
      write(*,*) 'Ions Broke - colliding with electrons'
      write(*,*) l,k, cl%ps(k)%vx, cl%ps(k)%vy, cl%ps(k)%vz
      if(taun.ne.taun)then
      write(*,*) 'passed'
      end if
      stop 666
   end if
    end do


!Now do ion-ion collisions
    !Skip if only one ion remains
    if(num_ions-Ncoll.le.1) return


    !Redefine reduced mass
       mu=.5*tr_m
       loglambda=log(4*pi*e0*mu*((6*kb*(cl%tempavg/tr_m)))*lD/tr_m**2)
   Aab=ndens*.25/pi*(tr_m**2/e0/mu)**2*loglambda


   do i=1,INT((num_ions-Ncoll)/2)
      if (2*i+Ncoll.gt.num_ions) write(*,*) 'Coll_pairs on line 569'

      k=coll_pairs(2*i+Ncoll)
      l=coll_pairs(2*i+Ncoll-1)
      vrel(1)=cl%ps(k)%vx-cl%ps(l)%vx
      vrel(2)=cl%ps(k)%vy-cl%ps(l)%vy
      vrel(3)=cl%ps(k)%vz-cl%ps(l)%vz
      vperp=sqrt(vrel(2)**2+vrel(3)**2)
      vrelmag=sqrt(vperp**2+vrel(1)**2)

!write(*,*) 'Relative velocity', vrel(1), vrel(2), vrel(3), vperp, vrelmag

      taun=Aab*tau/vrelmag**3
      tmp=exp(-taun)
      A=1.d0/(1.d0-tmp)+(-2.1939+2.6742*tmp-1.2038*tmp**2+1.4485*tmp**3)*exp(-.7856d0/(1.d0-tmp))


! write(*,*) taun,loglambda,A,' %%'
      call rand_mt(tmp1)
      psi=tmp1*2.d0*pi


      h(1)=vperp*cos(psi)
      h(2)=-(vrel(1)*vrel(2)*cos(psi)+vrelmag*vrel(3)*sin(psi))/vperp
      h(3)=-(vrel(1)*vrel(3)*cos(psi)-vrelmag*vrel(2)*sin(psi))/vperp

! Do the Nanbu inversion to find cos(chi)
      call rand_mt(tmp1)
      exp2A = exp(-2.d0*A)
      coschi = 1.d0 + 1.d0/A*log(tmp1*(1-exp2A)+exp2A)
      sinchi = sqrt(1.d0-coschi**2)




      !set post collision velocities
      cl%ps(l)%vx=cl%ps(l)%vx-mu/tr_m*(vrel(1)*(1.d0-coschi)+h(1)*sinchi)
      cl%ps(l)%vy=cl%ps(l)%vy-mu/tr_m*(vrel(2)*(1.d0-coschi)+h(2)*sinchi)
      cl%ps(l)%vz=cl%ps(l)%vz-mu/tr_m*(vrel(3)*(1.d0-coschi)+h(3)*sinchi)
      cl%ps(k)%vx=cl%ps(k)%vx+mu/tr_m*(vrel(1)*(1.d0-coschi)+h(1)*sinchi)
      cl%ps(k)%vy=cl%ps(k)%vy+mu/tr_m*(vrel(2)*(1.d0-coschi)+h(2)*sinchi)
      cl%ps(k)%vz=cl%ps(k)%vz+mu/tr_m*(vrel(3)*(1.d0-coschi)+h(3)*sinchi)

   if(cl%ps(l)%vy.ne.cl%ps(l)%vy)then
      write(*,*) 'Ions Broke- colliding with ions'
      write(*,*) 'First bits', A, coschi, sinchi, mpi_rank
      write(*,*) 'vrel', mpi_rank, vrel(1), vrel(2), vrel(3)
      write(*,*) 'Prelim', Aab, tau, cl%densavg, taun, loglambda, cl%tempavg
      stop 666
   end if

   end do

   deallocate(coll_pairs)



  end subroutine trace_self_collide






end module collide_mod













