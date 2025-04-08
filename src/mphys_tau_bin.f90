!
! This Module contains the code required to call 
! Tel-aviv university (TAU) bin microphysics subroutines
! within the 1-D framework
!
!
module mphys_tau_bin

  Use parameters, only : num_h_moments, num_h_bins, &
       nspecies, mom_names, h_names, mom_units, max_char_len, &
       num_aero_moments,num_aero_bins, aero_mom_init
  Use column_variables
  Use physconst, only : p0, this_r_on_cp=>r_on_cp, pi
  Use mphys_tau_bin_declare
  
  Use switches_bin
  Use module_mp_tau_bin 
  Use module_bin_init

 implicit none

  real :: q_lem(JMINP:JMAXP, KKP, NQP)
  real :: th_lem(JMINP:JMAXP, KKP)
  real :: sq_lem(JMINP:JMAXP, KKP, NQP)
  real :: sth_lem(JMINP:JMAXP, KKP)
  real :: w_lem(JMINP:JMAXP, KKP)

  integer i

  !Logical switches 
  logical :: l_ice=.False.
  logical :: micro_unset=.True.

 type qindex
     integer :: ispecies ! Species index
     integer :: imoment  ! moment index
  end type qindex
     
  type(qindex), allocatable :: qindices(:)
  integer :: nqs ! number of qindices (possibly different from nqp)


contains

  subroutine set_micro

    rdt=1./dt

    ! vapour
    
    iq=1
    iqv=iq

    if (num_h_bins(1) >= 11) then
      IMICROBIN=1
      IRAINBIN=1
       
      call qcount(iqss, iq)   ! advected supersat.

!      call qcount(iql, iq)   ! total water for diag
      
      if (num_aero_moments(1) >= 1) then
         do i = 1,ln2
            call qcount(IAERO_BIN(i),iq) ! aerosol bins
         enddo
      end if
      do i = 1,lk
        call qcount(ICDKG_BIN(i),iq) ! cloud mass bins
      enddo
      do i = 1,lk
        call qcount(ICDNC_BIN(i),iq) ! cloud number bins
      enddo
      
      nqs=aero_bin+ln2+lk+lk+1

      allocate(qindices(nqs))      

      ! Set qindices for later use
      if (num_aero_moments(1) >= 1) then
         do i = 1,ln2
            qindices(IAERO_BIN(i))%ispecies=1
            qindices(IAERO_BIN(i))%imoment=1
         enddo
      end if
      if (num_h_bins(1) >= 4)then ! cloud mass
         do i = 1,lk
            qindices(ICDKG_BIN(i))%ispecies=1
            qindices(ICDKG_BIN(i))%imoment=1
            qindices(ICDNC_BIN(i))%ispecies=1
            qindices(ICDNC_BIN(i))%imoment=2
         enddo
      end if

    end if ! bin model selected
  
  end subroutine set_micro


 
  subroutine qcount(var, count)
    integer, intent(out) :: var
    integer, intent(inout) ::  count

    count=count+1
    var=count
  end subroutine qcount

  subroutine mphys_tau_bin_interface

    ! Warm bin scheme
    integer ::j, k, iq


    ! Set up input arrays...
    rprefrcp(2:kkp)=exner(1:kkp-1,nx) ! I think this is upside-down in LEM
    ! AH - 04/03/10, line below leads to divide by 0 
    !      corrected by setting array to 2:kkp. Problem highlighted
    !      by Theotonio Pauliquevis
    ! prefrcp(:)=1./rprefrcp(:)
    prefrcp(2:kkp)=1./rprefrcp(2:kkp)

    prefn(2:kkp)=p0*exner(1:kkp-1,nx)**(1./this_r_on_cp)

    dzn(2:kkp)=dz_half(1:kkp-1)

    rhon(2:kkp)=rho(1:kkp-1)
    rdz_on_rhon(2:kkp)=1./(dz(1:kkp-1)*rhon(2:kkp))
    ! Reference temperature (this is fixed in lem, but 
    ! shouldn't make a difference for microphysics if we 
    ! just set it to be the current profile (i.e. th'=0)
    tref(2:kkp)=theta(1:kkp-1,nx)*exner(1:kkp-1,nx)

    ! Set up microphysics species
    if (micro_unset)then ! See later call to microset
       call set_micro
    end if


    do j=jminp,jmaxp
       q_lem (j, 2:kkp, iqv) = qv(1:kkp-1, j)
       q_lem (j, 2:kkp, iqss) = ss(1:kkp-1, j)

       do iq=1,ln2
         ih=qindices(IAERO_BIN(iq))%ispecies
         imom=qindices(IAERO_BIN(iq))%imoment
         do k=1,nz-1
             q_lem (j, k+1, IAERO_BIN(iq)) = aerosol(k,j,ih)%moments(iq,imom)
         end do
       enddo
       do iq=1,lk
         ! mass bins
         ih=qindices(ICDKG_BIN(iq))%ispecies
         imom=qindices(ICDKG_BIN(iq))%imoment
         do k=1,nz-1
            q_lem (j, k+1, ICDKG_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
         end do
         ! number bins
         ih=qindices(ICDNC_BIN(iq))%ispecies
         imom=qindices(ICDNC_BIN(iq))%imoment
         do k=1,nz-1
            q_lem (j, k+1, ICDNC_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
         end do
       enddo
       th_lem (j, :) = 0.0
       w_lem(j,2:kkp)=w_half(1:kkp-1,j)
    end do

    if (micro_unset)then
        call bin_init !initialises the cloud bin categories
        call data     !reads in and sets the coll-coal kernal

       DO IQ = 1,LN2
         ih=qindices(IAERO_BIN(iq))%ispecies
         imom=qindices(IAERO_BIN(iq))%imoment
         DO K=1,nz-1
           DO J = JMINP , JMAXP
             CCNORIG(J,K+1,IQ) = aerosol(k,j,ih)%moments(iq,imom)
           ENDDO
         ENDDO
       ENDDO

       DO K = 1, KKP
         DO J = JMINP, JMAXP
           TOTCCNORIG(J,K) = 0.0
           DO IQ = 1, LN2
             TOTCCNORIG(J,K) = TOTCCNORIG(J,K) + CCNORIG(J,K,IQ)
           ENDDO
         ENDDO
       ENDDO
       DO IQ = 1, Ln2
          CCNORIGTOT(IQ) = 0.0
          CCNORIGAVG(IQ) = 0.0
           DO K = 1, KKP
             DO J = JMINP, JMAXP
               CCNORIGTOT(IQ) = CCNORIGTOT(IQ) + CCNORIG(J,K,IQ)
             ENDDO
           ENDDO
           CCNORIGAVG(IQ) = CCNORIGTOT(IQ)/(JJP*KKP)
       ENDDO    
        micro_unset=.False.
    end if

    ! This bit doesn't yet use the switches on advection...
    do j=jminp,jmaxp
       sth_lem(j,2:kkp)=dtheta_adv(1:kkp-1,j)+dtheta_div(1:kkp-1,j)
       sq_lem(j,2:kkp,iqv)=dqv_adv(1:kkp-1,j)+dqv_div(1:kkp-1,j)

       sq_lem(j,2:kkp,iqss)=dss_adv(1:kkp-1,j)+dss_div(1:kkp-1,j)
       
       do iq=1,ln2
          ih=qindices(iaero_bin(iq))%ispecies
          imom=qindices(iaero_bin(iq))%imoment
          do k=1,nz-1
             sq_lem(j,k+1,iaero_bin(iq))=(daerosol_adv(k,j,ih)%moments(iq,imom) &
                  + daerosol_div(k,j,ih)%moments(iq,imom))
          end do
       enddo
       do iq=1,lk
          ih=qindices(icdkg_bin(iq))%ispecies
          imom=qindices(icdkg_bin(iq))%imoment
          do k=1,nz-1
             sq_lem(j,k+1,icdkg_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                  + dhydrometeors_div(k,j,ih)%moments(iq,imom)
          end do
          ih=qindices(icdnc_bin(iq))%ispecies
          imom=qindices(icdnc_bin(iq))%imoment
          do k=1,nz-1
             sq_lem(j,k+1,icdnc_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                  + dhydrometeors_div(k,j,ih)%moments(iq,imom)
          end do
       end do
     end do
! test if the transport has moved mass and number around
     DO K = 2, nz
        DO J = jminp,jmaxp
           DO IQ = 1, LK
              CALL ADVECTcheck(J,K,iq,DT,Q_lem(J,K,ICDKG_BIN(iq)),              &
    &                 Q_lem(J,K,ICDNC_BIN(iq)),SQ_lem(J,K,ICDKG_BIN(iq)),  &
    &                 SQ_lem(J,K,ICDNC_BIN(iq)))
           ENDDO
        ENDDO
     ENDDO

     call tau_bin(1, th_lem, q_lem, sth_lem, sq_lem, dt, rdt )


     do j=jminp,jmaxp
        sth_lem(j,2:kkp)=sth_lem(j,2:kkp)-(dtheta_adv(1:kkp-1,j)+dtheta_div(1:kkp-1,j))
        sq_lem(j,2:kkp,iqv)=sq_lem(j,2:kkp,iqv)-(dqv_adv(1:kkp-1,j)+dqv_div(1:kkp-1,j))

        do iq=1,ln2
           ih=qindices(iaero_bin(iq))%ispecies
           imom=qindices(iaero_bin(iq))%imoment
           do k=1,nz-1
              sq_lem(j,k+1,iaero_bin(iq))=sq_lem(j,k+1,iaero_bin(iq))       &
                   - (daerosol_adv(k,j,ih)%moments(iq,imom)                &
                   + daerosol_div(k,j,ih)%moments(iq,imom))
           end do
        enddo
        
        do iq=1,lk
           ih=qindices(icdkg_bin(iq))%ispecies
           imom=qindices(icdkg_bin(iq))%imoment
           do k=1,nz-1
              sq_lem(j,k+1,icdkg_bin(iq))= sq_lem(j,k+1,icdkg_bin(iq))      &
                   - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                   + dhydrometeors_div(k,j,ih)%moments(iq,imom))
           end do
           ih=qindices(icdnc_bin(iq))%ispecies
           imom=qindices(icdnc_bin(iq))%imoment
           do k=1,nz-1
              sq_lem(j,k+1,icdnc_bin(iq))= sq_lem(j,k+1,icdnc_bin(iq))      &
                   - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                   + dhydrometeors_div(k,j,ih)%moments(iq,imom))
           end do
        end do
        
        ! For now set no microphysics on the bottom level - this would be 
       ! better done by having a subterranian level 0 in column variables
        sq_lem(j,1,1)=0
        
        dtheta_mphys(1:kkp-1,j)=sth_lem(j,2:kkp)

        dqv_mphys(1:kkp-1,j)=sq_lem(j,2:kkp,iqv)
        
       !
       ! update supersaturation field here (not in step fields)
       !
        ss(1:kkp-1,j) = q_lem(j,2:kkp,iqss)
    
        do iq=1,ln2
           ih=qindices(iaero_bin(iq))%ispecies
           imom=qindices(iaero_bin(iq))%imoment
           do k=1,nz-1
              daerosol_mphys(k,j,ih)%moments(iq,imom) =              &
                   sq_lem(j,k+1,iaero_bin(iq))
           end do
        enddo

        do iq=1,lk
           ih=qindices(icdkg_bin(iq))%ispecies
           imom=qindices(icdkg_bin(iq))%imoment
           do k=1,nz-1
              dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =         &
                   sq_lem(j,k+1,icdkg_bin(iq))
           end do
           ih=qindices(icdnc_bin(iq))%ispecies
           imom=qindices(icdnc_bin(iq))%imoment
           do k=1,nz-1
              dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =        &
                   sq_lem(j,k+1,icdnc_bin(iq))
           end do
        end do
     end do
     
   end subroutine mphys_tau_bin_interface

!DECK ADVcheck
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE  ADVECTcheck(J,K,IQ,DT,ZQmass,ZQnum,Sourcemass,&
           Sourcenum)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
      IMPLICIT NONE
      
!CALL PRAMETR
!CALL RI
!CALL XD      
!CALL GRID1     
!local variable
      REAL :: QmassFLD, ZQmass, Qmass
      REAL :: Qnumfield, ZQnum, Qnum
      REAL :: Sourcemass
      REAL :: Sourcenum
      REAL :: AVG, AVGinit
      REAL :: DT,RDT
!loop counters
      INTEGER J, K, IQ      
          
      RDT = 1.0/DT

!First calculate the new field      
      QmassFLD=(ZQmass+(DT*Sourcemass))
      Qnumfield = (ZQnum+(DT*Sourcenum))
!Change units to microphys units (just for consistency
 
      QmassFLD = (QmassFLD*rhon(k))/1.e3
      Qnumfield = (Qnumfield*rhon(k))/1.e6
      Sourcemass = (Sourcemass*rhon(k))/1.e3
      Sourcenum = (Sourcenum*rhon(k))/1.e6

!calculate the average particle size for the bin

      IF(Qnumfield  >  0.0)THEN
        AVG =QmassFLD/Qnumfield
        IF(AVG >  (2.*X_BIN(IQ))) THEN
          sourcenum=((QmassFLD/(2.*X_BIN(IQ)-1.e-20))-                        &
     &                             ((ZQnum*rhon(k))/1.e6))*RDT
        ENDIF
        IF(AVG <  X_BIN(IQ).AND.AVG >  0.0) THEN
          sourcenum=((QmassFLD/(X_BIN(IQ)+1.e-20))-                           &
     &                             ((ZQnum*rhon(k))/1.e6))*RDT
        ENDIF

      ENDIF
!
      Sourcemass = (Sourcemass*1.e3)/rhon(k)
      Sourcenum  = (Sourcenum*1.e6)/rhon(k)
      
!do positivity check after normalisation as changing SQ by normalisation 
!can lead to negatives      
      IF((ZQmass+(DT*Sourcemass) <  0.0).or.                            &
     &                    (ZQnum+(DT*Sourcenum) <  0.0)) THEN
        Sourcemass = -0.99999 * ZQmass/DT
        Sourcenum =  -0.99999 * ZQnum/DT
      ENDIF

      IF (ABS(Sourcemass) <  1.e-20.OR. ABS(Sourcenum) <  1.e-20) THEN
          Sourcemass = 0.0
          Sourcenum = 0.0
      ENDIF

      END subroutine ADVECTcheck



end module
