






!----- move_mod -----
! This fortran module defines the routine for moving, including interactions with walls
module move_mod
    use helpers_mod
    implicit none
    contains


    !----- move -----
    ! this subroutine makes all of the collision cells move their particles
    subroutine move
        use core_data_mod
        use bounds_mod
        use simulation_mod
        use parallel_info_mod
        implicit none
        real(8) dN,dE,Nnet,Enet
        integer i,j

        Nnet=0.d0
        Enet=0.d0
        dN=0.d0
        dE=0.d0
  
        !write(*,*) "starting cell move:"
        do i=1,num_cells
            call cell_move(cells(i), tau,dN,dE)
!           Nnet=Nnet+dN
!           Enet=Enet+dE
        end do




        !write(*,*) "starting ghost cell move: (num_gc:", num_gc, "):"
        do i=1,num_gc
            !write(*,*) "ghost cell number: ", i

            ! watch i=1, mpi_rank=0
            if(mpi_rank.eq.0.and.i.eq.1) then
             dN=-10.d0
            else
             dN=0.d0
            end if

            call cell_move(gc(i), tau,dN,dE)


!           Nnet=Nnet+dN
!           Enet=Enet+dE
        end do

!      if(mpi_rank.eq.0) then
!       stop 666
!      end if


!       if(Nnet.gt.0.d0) then
!          write(*,*) 'Nnet, Enet: ',mpi_rank,Nnet,Enet
!       end if


    end subroutine move
  






    !--- cell_move ---
    ! Moves all the particles in the cell, performing interactions with the metals 
    subroutine cell_move(cl, tau,dN,dE)
        use core_types_mod
        use core_data_mod
        use helpers_mod
        use constants_mod
        use bounds_mod
        use parallel_info_mod
        use simulation_mod
	    implicit none

        !inputted values	
        type(cell_type) cl

        type(particle_type) p0
        real(8) tau

        ! we will know the details of a mask like this: masks(cl%mask)
        type(mask_type) msk
        type(particle_type) temp
        real(8) th1, th2, u1, u2
        real(8) answer,dr,dz,taubar,dperp1,dperp2,dperp,vtot2
        real(8) a,b,c,e,f,d,t1,t2,myhit,chk1,chk2,uh1,uh2
        real(8) thit,uhit,zhit,rhit,vx,vy,vz,di,x,y,vxp,vyp
        real(8) grhit,gzhit,wrhit,wzhit,vn,vtmp,offset,nx,ny,nz,px,py,pz
        real(8) E1,E2
        integer i,j,k,count,klook,flag,count1,count2
        real(8) tmprand, tmprand2
        real(8) vthwall, vg, vw, v_2D, theta
        real(8) Nin,Nout,Ein,Eout,zold,znew,zchkleft,dE,dN
        real(8) rnd,vth,phi,r0,z0,n0r,n0z
        real(8) flux, densavg, vbar, vavg, vrfluid, vzfluid, vadd, vbarr, vbarz
        real(8) tempavg, vr_fluidavg,vz_fluidavg, num_avg

        integer numtrace,numargon,ichk,num_dens
        integer :: part_count=0, cur_proc=0

        !This section tracks the fluid properties in the cell and stores them
        !over several timesteps.  The number of timesteps that it saves is set
        !in core.f90, in the cell_type.  Turn this section on when needed by
        !setting fluid_switch to 1 in core.f90, in the constants_mod.  This
        !section must be on in order for Nanbu collsions to work
        if(fluid_switch) then
           !Read in the number of timesteps that will be stored.    
           num_dens=size(cl%dens)

           !initialize the variables
           if(cur_step.eq.1)then
              do k=1,num_dens
                 cl%dens(k)=0d0
                 cl%vr_fluid(k)=0d0
                 cl%vz_fluid(k)=0d0
                 cl%temp(k)=0d0
              end do
           end if
        
   
           !Now save old densities
           do k=1,num_dens-1
               cl%dens(k)=cl%dens(k+1)
               cl%vr_fluid(k)=cl%vr_fluid(k+1)
               cl%vz_fluid(k)=cl%vz_fluid(k+1)
               cl%temp(k)=cl%temp(k+1)
           end do
           !only find values when there is a particle in the cell.
           if(cl%num_parts-cl%partition.gt.0)then
              ! Now load the last element of the array with the current density
              cl%dens(num_dens)=dble((cl%num_parts-cl%partition)*Nef*1e-3)/cl%volume       
              cl%vr_fluid(num_dens)=sum(cl%ps(cl%partition+1:cl%num_parts)%vx)/(cl%num_parts-cl%partition)
        
              cl%vz_fluid(num_dens)=sum(cl%ps(cl%partition+1:cl%num_parts)%vz)/(cl%num_parts-cl%partition)
           end if    
           !If the 
           num_avg=dble(min(num_dens,cur_step))
           cl%vr_fluidavg=sum(cl%vr_fluid)/num_avg
           cl%vz_fluidavg=sum(cl%vz_fluid)/num_avg
           cl%densavg=sum(cl%dens)/num_avg

           !Now find the mean of v**2, then scale to temperature.
           cl%temp(num_dens)=sum((cl%ps(cl%partition+1:cl%num_parts)%vx-cl%vr_fluidavg)**2+(cl%ps(cl%partition+1:cl%num_parts)%vz-cl%vz_fluidavg)**2)/num_avg
           !Now scale into temperature
           cl%temp(num_dens)=tr_m*cl%temp(num_dens)/3d0/kb

           !Now take the average and compute the new cross section scale
           cl%tempavg=sum(cl%temp)/num_avg
           !Now avoid INF and NAN by making a zero temp just very small
           if(cl%tempavg.lt.1) cl%tempavg=1d0

           if(cl%tempavg.ne.cl%tempavg)then
               write(*,*) '##', num_avg,  cl%num_parts-cl%partition, cl%vz_fluid(num_dens), cl%vr_fluid(num_dens)
               write(*,*) 'TEMP', cl%temp
               write(*,*) 'VR', cl%vr_fluid
               write(*,*) 'VZ', cl%vz_fluid
               stop 666
           end if
        end if    
! variables for tracking particles and energy entering and leaving

!        zchkleft = 1.038d-4*.5
!        Nin=0.d0
!        Nout=0.d0
!        Ein=0.d0
!        Eout=0.d0


!        write(*,*) cl%num_parts-cl%partition,cl%min_r,cl%min_z,mpi_rank


        ! set the size of the offset which will be used to move a particle that hits
        ! the wall a little bit into the cell before computing the next bounce
        offset = (cl%max_z-cl%min_z)*1.d-4


        ! start of move if block
        ! if the mask code is 0, just move them




!*********start mask=0, ballistic move and return ***************
        if(cl%mask.eq.0) then


            !****Start Ballistic move Loop over particles
            k = 1

!           if(cl%ps(cl%num_parts)%element.eq.1.and.mpi_rank.eq.0) then
!               write(*,*) 'Trace particle found',cl%partition,cl%num_parts-cl%partition
!               count1=0
!               count2=0
!               do i=1,cl%num_parts
!                 if(cl%ps(i)%element.eq.0) count1=count1+1
!                 if(cl%ps(i)%element.eq.1) count2=count2+1
!               end do
!               write(*,*) 'Normal and trace:',count1,count2
!               write(*,*) 'Trace move info',tr_q,tau,tr_m,cl%efield_r,cl%efield_z

!               stop 666
!           end if



            do while(k .le. cl%num_parts)


! special: dump particles
!        if(cur_step.eq.8.and.mpi_rank.eq.0) then

! special: dump particles


		if( cl%ps(k)%element == 0 ) then
			! normal particle, ballistic move
                	cl%ps(k)%x = cl%ps(k)%x + cl%ps(k)%vx*tau
                	cl%ps(k)%y = cl%ps(k)%y + cl%ps(k)%vy*tau

                        zold = cl%ps(k)%z
                	cl%ps(k)%z = cl%ps(k)%z + cl%ps(k)%vz*tau
                        znew=cl%ps(k)%z


		else
					! trace particle, special move (1 of 2)

                    ! step 1: modify the velocity as a result of the electric field
                    ! vnew_vector = vold_vector + (tr_q*tau/tr_m)*ElectricField_R(position_old)

                    cl%ps(k)%vx = cl%ps(k)%vx + (tr_q*tau/tr_m)*cl%efield_r
                    cl%ps(k)%vz = cl%ps(k)%vz + (tr_q*tau/tr_m)*cl%efield_z

                    ! step 2: Move the particle with the new velocity
                	cl%ps(k)%x = cl%ps(k)%x + cl%ps(k)%vx*tau
                	cl%ps(k)%y = cl%ps(k)%y + cl%ps(k)%vy*tau
                	cl%ps(k)%z = cl%ps(k)%z + cl%ps(k)%vz*tau
		end if

! rotate the particle to theta=0 (x-z plane) where x=d and y=0
! also set ps%r=d

                d=sqrt(cl%ps(k)%x**2+cl%ps(k)%y**2)
                di=1.d0/d

                vx=cl%ps(k)%vx
                vy=cl%ps(k)%vy
                cl%ps(k)%vx=di*(cl%ps(k)%x*vx+cl%ps(k)%y*vy)
                cl%ps(k)%vy=di*(cl%ps(k)%x*vy-cl%ps(k)%y*vx)
                cl%ps(k)%x=d
                cl%ps(k)%y=0.d0
                cl%ps(k)%r=d


!                        zchkleft=0.d0
!                        ! watch particles coming in from the left
!                        if( sign(1.d0,zchkleft-zold).ne.sign(1.d0,zchkleft-znew) ) then

!                        if(cl%ps(k)%element.eq.1) then
!                            write(*,*) 'Trace coming in',cl%ps(k)%z,cl%ps(k)%r,cl%ps(k)%vz,cl%ps(k)%vx
!                        end if

                            ! entering particle
!                           if(cl%ps(k)%vz.gt.0.d0) then
!                              Nin = Nin + 1.d0
!                              Ein = Ein + .5d0*m*(cl%ps(k)%vx**2 + cl%ps(k)%vy**2 + cl%ps(k)%vz**2 )
!                           else
!                              Nout = Nout + 1.d0
!                              Eout = Eout + .5d0*m*(cl%ps(k)%vx**2 + cl%ps(k)%vy**2 + cl%ps(k)%vz**2 )
!                           end if
!                        end if
 !                       ! watch particles coming in from the left

                ! Checking if particles have moved out of the region, 
                ! and if so, swapping them out
                ! This is because particles remaining in the cell are lined in the array from left 
                !     to right, the max being "num_parts", while particles which
                !     need to be transferred to other cells (in module "reload_cells")
                !     are lined in the array from right to left, the max being "num_left"
                if(  cl%ps(k)%r .lt. cl%min_r .or. &
                    cl%ps(k)%r .gt. cl%max_r .or. &
                    cl%ps(k)%z .lt. cl%min_z .or. &
                    cl%ps(k)%z .gt. cl%max_z ) &
                    then

                    !swap the current particle out
					if( cl%ps(k)%element == 0 ) then
						! normal particle
                    	temp = cl%ps(k)
                    	cl%ps(k) = cl%ps(cl%partition)
						cl%ps(cl%partition) = cl%ps(cl%num_parts)
                    	cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                        cl%partition = cl%partition - 1
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					else
!                                     write(*,*) 'Moved trace to another cell',cl%min_r
						! trace particle
						temp = cl%ps(k)
						cl%ps(k) = cl%ps(cl%num_parts)
						cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					end if
                else
                    k = k + 1 !Increment counter (instance 1 of 1 in this do loop)
                end if
            end do
            !***end ballistic loop over particles
            return
        end if
!***********end mask=0, ballistic *******************************

!************** start mask<0, inside metal**********************
!Empty all cells inside of the metal
        if(cl%mask.lt.0)then
           cl%num_parts = 0 
           cl%partition = 0 
           cl%num_left = 0   
           return
        end if
               
!**************** end mask<0, inside metal*********************

! if we get to here the mask code isn't zero
! unload the mask code
        msk=masks(cl%mask)

! this is a value which is needed to calculate thermal reflection.
! It is the thermal velocity of a particle at the temperature of a wall
        vthwall = sqrt( kb * msk%T / m )
 

!********** start of loop over particles **********************
        k = 1
        do while(k .le. cl%num_parts)

            zold = cl%ps(k)%z

            !****FINDING OUT IF IT'S CLOSE ENOUGH TO HIT A WALL

            ! find the perpendicular distance to the boundary(ies)
            dperp1=(cl%ps(k)%r-msk%s0r)*msk%g1r + (cl%ps(k)%z-msk%s0z)*msk%g1z
            dperp1=abs(dperp1)

            if( abs(msk%w2r) .gt. 1.0d-9 .or. abs(msk%w2z) .gt. 1.0d-9 ) then
                dperp2=(cl%ps(k)%r-msk%s0r)*msk%g2r + (cl%ps(k)%z-msk%s0z)*msk%g2z
                dperp2=abs(dperp2)
            else
                dperp2=dperp1
            end if
            dperp=min(dperp1,dperp2)

            ! find the square of the particle speed
            vtot2=cl%ps(k)%vx**2+cl%ps(k)%vy**2+cl%ps(k)%vz**2



            ! if the particle can't go far enough in time tau to reach the wall
            ! along the shortest path it can't reach it at all - just move the
            ! particle
            if( vtot2*tau**2 < dperp**2 ) then

				! move the particle ballistically
					
				if( cl%ps(k)%element == 0 ) then
					! normal particle, ballistic move
                    cl%ps(k)%x = cl%ps(k)%x + cl%ps(k)%vx*tau
                    cl%ps(k)%y = cl%ps(k)%y + cl%ps(k)%vy*tau

                    cl%ps(k)%z = cl%ps(k)%z + cl%ps(k)%vz*tau
                    znew = cl%ps(k)%z


				else
					! trace particle, special move (2 of 2)

                    ! step 1: modify the velocity as a result of the electric field
                    ! vnew_vector = vold_vector + (tr_q*tau/tr_m)*ElectricField_R(position_old)

                    cl%ps(k)%vx = cl%ps(k)%vx + (tr_q*tau/tr_m)*cl%efield_r
                    cl%ps(k)%vz = cl%ps(k)%vz + (tr_q*tau/tr_m)*cl%efield_z

                    ! step 2: Move the particle with the new velocity
                	cl%ps(k)%x = cl%ps(k)%x + cl%ps(k)%vx*tau
                	cl%ps(k)%y = cl%ps(k)%y + cl%ps(k)%vy*tau

                	cl%ps(k)%z = cl%ps(k)%z + cl%ps(k)%vz*tau
                	znew = cl%ps(k)%z
				end if

                ! rotate the particle to theta=0 (x-z plane) where x=d and y=0
                ! also set ps%r=d

                d=sqrt(cl%ps(k)%x**2+cl%ps(k)%y**2)
                di=1.d0/d

                vx=cl%ps(k)%vx
                vy=cl%ps(k)%vy
                cl%ps(k)%vx=di*(cl%ps(k)%x*vx+cl%ps(k)%y*vy)
                cl%ps(k)%vy=di*(cl%ps(k)%x*vy-cl%ps(k)%y*vx)
                cl%ps(k)%x=d
                cl%ps(k)%y=0.d0
                cl%ps(k)%r=d

                ! this particle can't hit the wall, so we just moved it ballistically
            

                ! Checking if particle has moved out of the region, 
                ! and if so, swapping it out
                ! This is because particles remaining in the cell are lined in the array from left 
                !     to right, the max being "num_parts", while particles which
                !     need to be transferred to other cells (in module "reload_cells")
                !     are lined in the array from right to left, the max being "num_left"
                if(  cl%ps(k)%r .lt. cl%min_r .or. &
                    cl%ps(k)%r .gt. cl%max_r .or. &
                    cl%ps(k)%z .lt. cl%min_z .or. &
                    cl%ps(k)%z .gt. cl%max_z ) &
                    then

                    !swap the current particle out
					if( cl%ps(k)%element == 0 ) then
						! normal particle
                    	temp = cl%ps(k)
                    	cl%ps(k) = cl%ps(cl%partition)
						cl%ps(cl%partition) = cl%ps(cl%num_parts)
                    	cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                        cl%partition = cl%partition - 1
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					else
						! trace particle
						temp = cl%ps(k)
						cl%ps(k) = cl%ps(cl%num_parts)
						cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					end if
                else
                    k = k + 1 !Increment counter (instance 1 of 2 in this do loop)
                end if

                cycle
		    end if   
				
				
!******** top of multiple bounce while loop*************

			! set the total time of this step (handles multiple bounces)
            taubar=tau
            count=0
            do while (taubar.gt.0.d0)
                count=count+1
                if(count.ge.200) stop 678
                

! do the hard work of seeing if we hit the wall
! ********** time to hit segment 1 **************************************

                
                if(abs(msk%w1z) .le. 1.0d-9 .and. abs(msk%w1r) .le. 1.0d-9 ) then
                    write(*,*) "mask doesn't exist! w1r and w1z both zero"
                    stop 277

                ! vertical line is special
                else if(abs(msk%w1z) .le. 1.0d-9) then

                    t1 = (msk%s0z-cl%ps(k)%z)/cl%ps(k)%vz  ! much less round-off error
                    u1 = (sqrt((cl%ps(k)%x+t1*cl%ps(k)%vx)**2+ &
                            (cl%ps(k)%y+t1*cl%ps(k)%vy)**2) - msk%s0r)/msk%w1r
                    t2 = 1.d99 !just mark tau2 and s2 as infeasible
                    u2 = -1

                    ! check to make sure that u is positive and that the crossing time is positive
                    if(t1.lt.0.d0.or.u1.lt.0.d0) t1=1.d99

                    th1=t1
                    uh1=u1

                ! general case: quadratic
                else

                        a = msk%w1z**2*(cl%ps(k)%vx**2+cl%ps(k)%vy**2) - &
                        (msk%w1r*cl%ps(k)%vz)**2
                        e = msk%w1z*msk%s0r-msk%w1r*(msk%s0z-cl%ps(k)%z)
                        f = msk%w1z*cl%ps(k)%x
                        b=2*(msk%w1z*cl%ps(k)%vx*f-msk%w1r*cl%ps(k)%vz*e)
                        c=(f+e)*(f-e)

                        d=b**2-4*a*c

                        if(d.lt.0.d0) then
                            th1=1.d99
                        else
                            d=sqrt(d)
                            a=0.5d0/a
                            t1=(-b+d)*a
                            u1 = (cl%ps(k)%z+cl%ps(k)%vz*t1-msk%s0z)/msk%w1z
                            chk1 = msk%s0r+msk%w1r*u1
                            t2=(-b-d)*a
                            u2 = (cl%ps(k)%z+cl%ps(k)%vz*t2-msk%s0z)/msk%w1z
                            chk2 = msk%s0r+msk%w1r*u2


                            if(t1.le.0.d0.or.chk1.lt.0.d0.or.u1.lt.0.d0) t1=1.d99
                            if(t2.le.0.d0.or.chk2.lt.0.d0.or.u2.lt.0.d0) t2=1.d99

                            if(t1.lt.t2) then
                                th1=t1
                                uh1=u1
                            else
                                th1=t2
                                uh1=u2
                            end if

                        end if

                end if

                ! set thit and the point on the segment where the crossing occurs
                ! if there is a valid crossing
                ! also load grhit and gzhit with the components of g, the inward pointing
                ! unit vector, at the contact point
                if(th1.lt.1.d99.and.th1.lt.taubar) then
                        thit=th1
                        uhit=uh1
                        zhit=msk%s0z+uhit*msk%w1z
                        rhit=abs(msk%s0r+uhit*msk%w1r)
                        grhit=msk%g1r
                        gzhit=msk%g1z
                        wrhit=msk%w1r
                        wzhit=msk%w1z
                else
                        thit=1.d99
                        zhit=0.d0
                        rhit=0.d0
                        grhit=0.d0
                        gzhit=0.d0
                        wrhit=0.d0
                        wzhit=0.d0
                end if

! ********** end   time to hit segment 1 **************************************


! ********** begin time to hit segment 2 **************************************

                ! segment 2 exists if the w2 unit vector is nonzero
                if(abs(msk%w2r)+abs(msk%w2z).ge.1.0d-9) then

                    ! vertical line is special
                    if(abs(msk%w2z) .le. 1.0d-9)then

                        t1 = (msk%s0z-cl%ps(k)%z)/cl%ps(k)%vz  ! much less round-off error
                        u1 = (sqrt((cl%ps(k)%x+t1*cl%ps(k)%vx)**2+ &
                            (cl%ps(k)%y+t1*cl%ps(k)%vy)**2) - msk%s0r)/msk%w2r
                        t2 = 1.d99 !just mark tau2 and s2 as infeasible
                        u2 = -1

                        ! check to make sure that u is positive
                        if(t1.lt.0.d0.or.u1.lt.0.d0) t1=1.d99

                        th2=t1
                        uh2=u1

                    ! general case: quadratic
                    else

                        a = msk%w2z**2*(cl%ps(k)%vx**2+cl%ps(k)%vy**2) - &
                            (msk%w2r*cl%ps(k)%vz)**2
                        e = msk%w2z*msk%s0r-msk%w2r*(msk%s0z-cl%ps(k)%z)
                        f = msk%w2z*cl%ps(k)%x
                        b=2*(msk%w2z*cl%ps(k)%vx*f-msk%w2r*cl%ps(k)%vz*e)
                        c=(f+e)*(f-e)

                        d=b**2-4*a*c

                        if(d.lt.0.d0) then
                            th2=1.d99
                        else
                            d=sqrt(d)
                            a=0.5d0/a
                            t1=(-b+d)*a
                            u1 = (cl%ps(k)%z+cl%ps(k)%vz*t1-msk%s0z)/msk%w2z
                            chk1 = msk%s0r+msk%w2r*u1
                            t2=(-b-d)*a
                            u2 = (cl%ps(k)%z+cl%ps(k)%vz*t2-msk%s0z)/msk%w2z
                            chk2 = msk%s0r+msk%w2r*u2


                            if(t1.le.0.d0.or.chk1.lt.0.d0.or.u1.lt.0.d0) t1=1.d99
                            if(t2.le.0.d0.or.chk2.lt.0.d0.or.u2.lt.0.d0) t2=1.d99

                            if(t1.lt.t2) then
                                th2=t1
                                uh2=u1
                            else
                                th2=t2
                                uh2=u2
                            end if

                        end if

                    end if

                    ! if segment 2 exists we need to choose the soonest crossing between segment 1 and 2


                    ! set thit and the point on the segment where the crossing occurs
                    ! if there is a valid crossing and if it happens sooner than that
                    ! of segment 1
                    ! also find the inward pointing g vector at the contact point
                    if(th2.lt.th1.and.th2.lt.taubar) then

                        thit=th2
                        uhit=uh2
                        zhit=msk%s0z+uhit*msk%w2z
                        rhit=abs(msk%s0r+uhit*msk%w2r)
                        grhit=msk%g2r
                        gzhit=msk%g2z
                        wrhit=msk%w2r
                        wzhit=msk%w2z
                    end if

                end if
! ********** end time to hit segment 2 **************************************


!*********** move the particles
                ! if thit=1.d99 use taubar up with a ballistic move
                if(thit.gt.1.d98) then

                    cl%ps(k)%x = cl%ps(k)%x + cl%ps(k)%vx*taubar
                    cl%ps(k)%y = cl%ps(k)%y + cl%ps(k)%vy*taubar

                    cl%ps(k)%z = cl%ps(k)%z + cl%ps(k)%vz*taubar
                    znew = cl%ps(k)%z



                    ! rotate the particle to theta=0 (x-z plane)

                    d=sqrt(cl%ps(k)%x**2+cl%ps(k)%y**2)
                    di=1.d0/d

                    vx=cl%ps(k)%vx
                    vy=cl%ps(k)%vy
                    cl%ps(k)%vx=di*(cl%ps(k)%x*vx+cl%ps(k)%y*vy)
                    cl%ps(k)%vy=di*(cl%ps(k)%x*vy-cl%ps(k)%y*vx)
                    cl%ps(k)%x=d
                    cl%ps(k)%y=0.d0
                    cl%ps(k)%r=d


                    ! set taubar to zero because it's used up
                    taubar=0.d0
                    
                else

                    count = count + 1  ! counting how many bounces the particle had




!********begin SPECULAR/THERMALIZED REFLECTION***********
    ! We already know the location of the collision with the wall (rhit, zhit, from above)
    ! Now we just need to calculate the reflection off the wall

                    call rand_mt(tmprand)
                    !tmprand = 0.d0 ! BUG FINDING -- always do specular, just to see if this changes anything
                    !tmprand = 1.0 ! BUG FINDING -- always do thermalized, just to see if this changes anything
                    if( tmprand <= msk%reflect ) then

                        ! put the particle at the segment contact point, reset its velocity,
                        ! and decrement taubar

                        ! move the particle in the x-y plane for rotation purposes
                        x = cl%ps(k)%x + cl%ps(k)%vx*thit
                        y = cl%ps(k)%y + cl%ps(k)%vy*thit

                        ! now rotate the velocities
                        vxp = cl%ps(k)%vx
                        vyp = cl%ps(k)%vy
                        vz = cl%ps(k)%vz

                        d = sqrt(x**2+y**2)
                        di = 1.d0/d

                        vx = di*(x*vxp+y*vyp)
                        vy = di*(x*vyp-y*vxp)

                        ! specular reflection
                        ! do the specular reflection using the vector formula, but make it work
                        ! for particles coming from inside the metal as well - for these particles
                        ! kill the normal component and then sent the particle off with |vn|
                        ! away from the wall
                        ! v(new) = v(old) - (vnormal+abs(vnormal))*g(unit normal into wall)

                        vn = grhit*vx+gzhit*vz

                        vtmp=vn+abs(vn)
                        cl%ps(k)%vx = vx-vtmp*grhit
                        cl%ps(k)%vz = vz-vtmp*gzhit
                        cl%ps(k)%vy = vy

                    else
                        !thermalized reflection

                        ! find a random angle, theta
                        call rand_mt(tmprand2)
                        theta = 2*pi*tmprand2

                        ! sample a 2D maxwellian speed distribution
                        call rand_mt(tmprand2)
                        v_2D = vthwall*sqrt(-2.0d0 * log(tmprand2))
                        
                        !since vy is always parallel to the wall, we can set vy right now
                        cl%ps(k)%vy = v_2D*cos(theta)

                        ! vw is referring to the other parallel direction, parallel to the surface, and perpendicular to vy
                        ! vw is a scalar
                        vw = v_2D * sin(theta)

                        ! vg is referring the normal distribution, which is a 2D maxwellian speed distribution (not obvious)
                        ! vg is a scalar
                        call rand_mt(tmprand2)
                        vg = vthwall*sqrt(-2.0d0 * log(tmprand2)) ! again but with a new random number
                        
                        !  vector velocity v = (vw * w_unitvector) - (vg * g_unitvector)
                        !  vx = dot(v,r_unitvector)
                        !  vz = dot(v,z_unitvector)
                        cl%ps(k)%vz = -vg*gzhit + vw*wzhit
                        cl%ps(k)%vx = -vg*grhit + vw*wrhit
                        
                    end if
! ********** end Specular/Thermalized reflection

                    ! put the particle on the wall, and then just a little bit away
                    ! in the direction of -g, so that thit won't come back zero next time

                    cl%ps(k)%x = abs(rhit-grhit*offset)
                    cl%ps(k)%y = 0.d0
                    cl%ps(k)%z = zhit-gzhit*offset
                    znew = cl%ps(k)%z
                    cl%ps(k)%r = cl%ps(k)%x

                    ! decrement taubar
                    taubar=taubar-thit

                end if

            end do
!****** bottom of multiple bounce while loop *********

!****************REMOVE TRACE PARTICLES THAT HIT THE WALL TO FORM SHEATH******
                 if(cl%ps(k)%element.eq.1.and.count.gt.0) then
                   !Kick the particle outside of the simulation region so that
                   !it will be removed
                      cl%ps(k)%r = 2e3
                      cl%ps(k)%z = 2e3
                    
                 end if
!********END OF TRACE PARTICLE REMOVAL************************

!                        ! watch particles coming in from the left
!                        if( sign(1.d0,zchkleft-zold).ne.sign(1.d0,zchkleft-znew) ) then
!
!                           ! entering particle
!                           if(cl%ps(k)%vz.gt.0.d0) then
!                              Nin = Nin + 1.d0
!                              Ein = Ein + .5d0*m*( cl%ps(k)%vx**2 + cl%ps(k)%vy**2 + cl%ps(k)%vz**2 )
!                           else
!                              Nout = Nout + 1.d0
!                              Eout = Eout + .5d0*m*( cl%ps(k)%vx**2 + cl%ps(k)%vy**2 + cl%ps(k)%vz**2 )
!                           end if
!                        end if
!                        ! watch particles coming in from the left




             !******** Checking for particles exiting the region,
             !         and if so, swapping them out
             ! This is because particles remaining in the cell are lined in the array from left
             !     to right, the max being "num_parts", while particles which
            !     need to be transferred to other cells (in module "reload_cells")
            !     are lined in the array from right to left, the max being "num_left"
            if(  cl%ps(k)%r .lt. cl%min_r .or. &
                cl%ps(k)%r .gt. cl%max_r .or. &
                cl%ps(k)%z .lt. cl%min_z .or. &
                cl%ps(k)%z .gt. cl%max_z ) &
                then

                    !swap the current particle out
					if( cl%ps(k)%element == 0 ) then
						! normal particle
                    	temp = cl%ps(k)
                    	cl%ps(k) = cl%ps(cl%partition)
						cl%ps(cl%partition) = cl%ps(cl%num_parts)
                    	cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                        cl%partition = cl%partition - 1
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					else
						! trace particle
						temp = cl%ps(k)
						cl%ps(k) = cl%ps(cl%num_parts)
						cl%ps(cl%capacity-cl%num_left) = temp
                    	!update the counts
                    	cl%num_parts = cl%num_parts - 1
                    	cl%num_left = cl%num_left + 1
					end if
            else
                k = k + 1 !Increment counter (instance 2 of 2 in this do loop)
            end if
        end do
!****** end of loop over particles **************

!       dN=(Nin-Nout)*Nef
!       dE=Ein-Eout


    end subroutine cell_move




 
end module move_mod








