


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
!kluge   use parallel_info_mod
    implicit none
    integer i
! kluge
!   integer iflag
!   iflag=0


    do i=1,num_cells

! kluge look for cells that have z=3.3.89e-4 in them
!     if(cells(i)%min_z.lt.3.89e-4.and.cells(i)%max_z.gt.3.89e-4) then
!        write(*,*) mpi_rank,i,cells(i)%min_r,cells(i)%max_r
!        iflag = 1
!     end if

      ! do normal particle collisions
      call cell_collide_subdivide(cells(i))
      ! do trace particle collisions
      call trace_collide(cells(i),i)
      ! do trace particle self collisions
      if(nanbu_switch) call trace_self_collide(cells(i))
    end do

! kluge
!     if(iflag.eq.1) stop 666


  end subroutine collide






  subroutine cell_collide_subdivide(cl)
    use helpers_mod
    use core_types_mod

    implicit none
    type(cell_type) cl
    integer num_fracs, j, kmin, kmax
    real(8) zmin, zmax, frac, dn

    ! zero out all of the previous collision partners
    cl%ps(:)%prev_coll = 0

    ! control growth of max_vrel2
    cl%max_vrel2 = cl%max_vrel2 * .999

    if(cl%partition .le. 1) return

    ! I want to make sure that each sub-collision cell has at least 7 particles
    num_fracs = cl%partition / 7

    if(num_fracs .le. 1)then
      call cell_collide(cl, 1.d0, 1, cl%partition)
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
        call cell_collide(cl, frac, kmin, kmax)
      end do
    end if
  end subroutine cell_collide_subdivide




  !--- cell_collide (VSS) ---
  ! This subroutine collides the molecules in a given cell
  ! according to the variable soft sphere collider
  subroutine cell_collide(cl, frac, n1, n2)  ! collides argon with argon
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
    integer num_candidates
    integer i,j,k,l
    real(8) dist2, min_dist2
    real(8) temp,vrel2,vrel
    real(8) phi,cos_theta,sin_theta,cos_phi,sin_phi,dvz_sin_phi,vrel_cos_phi
    real(8) dvx,dvy,dvz
    real(8) vrx,vry,vrz
    real(8) vperp,vperpi
    real(8) vrx_cm,vry_cm,vrz_cm
    real(8) cmx,cmy,cmz
    integer nearswitch
    real(8) alpha,vexp,rexp
    real(8) sigmavrmax,sigmachk
    data alpha/1.66/
    vexp=0.5d0-nu
    rexp=1.d0/alpha

    num_parts = n2-n1+1

    ! compute number of candidates
    sigmavrmax = aref*cl%max_vrel2**vexp
    num_candidates = round_rand(Nef*tau*0.5d0*sigmavrmax*num_parts*(num_parts-1.d0)/(frac*cl%volume)) !&
                    !+ cl%num_cand_remaining
    !num_candidates = int(rnum_candidates)

    !cl%num_cand_remaining = rnum_candidates - num_candidates

    ! main collide loop
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

      cl%max_vrel2 = max(vrel2, cl%max_vrel2)
      sigmachk = aref*vrel2**vexp
      call rand_mt(temp)
      ! if rejecting this candidate, call 'cycle'
      if(sigmachk .le. temp*sigmavrmax)then
        cycle
      end if
      ! this candidate is accepted
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
      !num_collided = num_collided + 1

      if(isnan(cl%ps(j)%vx) .or. isnan(cl%ps(j)%vy) .or. isnan(cl%ps(j)%vz) )then
        write(*,*)'Found Nan.'
        write(*,*)cmx,cmy
        write(*,*)vrx,vry
        write(*,*)dvx,dvy
        stop 123
      end if

    ! end of main collide loop
    end do
    return
  end subroutine cell_collide




    ! collides all the trace-element particles
    ! which are in the ps array from index "partition+1" to "num_parts"
    subroutine trace_collide(cl,icell)  ! trace particles collide with argon, argon is unchanged
        use core_types_mod
        use simulation_mod
        use constants_mod
        use helpers_mod
! kluge
    use parallel_info_mod

        implicit none

! kluge
        integer icell

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
        real(8) sigmavrmax,sigmachk
        data alpha/1.66/
  !     vexp=0.5d0-nu_tr
        rexp=1.d0/alpha

        if(cur_step.eq.1) cl%aref_ambi=2.d0  !Set Te=Ti on the first step.
                                           !aref_ambi is later calculated in
                                           !move.f90

        num_tr_parts = cl%num_parts - cl%partition
        num_norm_parts = cl%partition

        ! too few particles to collide, so return
        if( num_norm_parts <= 1 .or. num_tr_parts <= 1 ) return

        ! control growth of tr_max_vrel2
        cl%tr_max_vrel2 = cl%tr_max_vrel2 * .999

        ! compute number of candidates
       !sigmavrmax = aref_tr/cl%aref_ambi*cl%tr_max_vrel2**vexp !old cross
       !section
       !New two nu cross section:
       sigmavrmax = (aref_tr*cl%tr_max_vrel2**(0.5d0-nu_tr1)+bref_tr*cl%tr_max_vrel2**(0.5d0-nu_tr2))/cl%aref_ambi

        num_candidates = round_rand(Nef*tau*sigmavrmax*num_norm_parts*(num_tr_parts-0.d0)/(cl%volume))

        !Check for errors in aref_ambi implementation
        if(sigmavrmax.ne.sigmavrmax)then
           write(*,*) 'broken cross-section', cl%aref_ambi, cur_step
           stop 1234
        end if

! kluge
!        if(mpi_rank.eq.0.and.icell.eq.501950) then
!            write(*,*) 'num_norm_parts  num_tr_parts  cl%tr_max_vrel2'
!            write(*,*)  num_norm_parts, num_tr_parts, cl%tr_max_vrel2
!            write(*,*) 'sigmavrmax  num_candidates  Nef'
!            write(*,*)  sigmavrmax , num_candidates,  Nef
!            write(*,*) 'tau  cl%volume'
!            write(*,*)  tau, cl%volume
!!           stop 666
!        end if



         ! Questions about the above equation:
        !  1. Do we use "num_norm_parts" or "num_tr_parts", or both?
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

            cl%tr_max_vrel2 = max(vrel2, cl%tr_max_vrel2)

           ! sigmachk = aref_tr/cl%aref_ambi*vrel2**vexp !old code
            !Two nu cross section
            sigmachk = (aref_tr*vrel2**(0.5d0-nu_tr1)+bref_tr*vrel2**(0.5d0-nu_tr2))/cl%aref_ambi

            call rand_mt(temp)
            ! if rejecting this candidate, call 'cycle'
            if(sigmachk .le. temp*sigmavrmax)then
                cycle
            end if

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













