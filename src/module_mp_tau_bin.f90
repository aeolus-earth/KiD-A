!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module contains the Tel-Aviv university (TAU) warm,
! size-resolved cloud microphysics scheme described in Tzivion et al,
! JAS (1987,JAS), Feingold et al (1988, JAS), and Tzivion et al (1989,
! JAS).
!
! In this version of the TAU microphysics the cloud drop size
! distribution is divided into 34 bins with a radii range of 1.56 to
! 3200 microns and mass doubling from one bin to the next. The
! method of moments (Tzivion et al. 1987, JAS) is used to solve for
! mass and number concentration in each size bin that result from
! diffusional growth (Tzivion et al 1989, JAS), collision-coalescence
! and collisional breakup (Tzivion et al, 1987 and Feingold et al,
! 1989, JAS). Sedimentation is performed with a first-order upwind scheme.
! Aerosol are represented by a single prognostic variable that is
! assumed to be ammonium sulfate with a log-normal distribution (Stevens
! et al 1996, JAS).
!
! The numerical methods and code in this module have been used in
! a variety of 2-D and 3-D dynamical frameworks to investigate a number
! of cloud microphysical problems. For example, drizzle production in marine
! Sc (Feingold et al, 1996), the dynamic and microphysical
! details of non-precipitating and precipitating marine Sc
! (Stevens et al,JAS, 1996 & 1998), the effect of drizzle on cloud optical
! depth and susceptibility (Feingold et al, JGR, 1997),
! the role of giant CCN in marine Sc,
! (Feingold et al, JAS, 1999), the role of giant CCN in cumulus (Yin et al,
! Atmospheric Research, 2000), turbulence, condensation and liquid water
! transport in non-precipitating marine Sc (Wang et al, JAS, 2003) and
! aerosol-cloud interactions (Feingold et al, GRL, 2005; Jiang et al, JGR,
! 2006; Xue and Feingold,JAS,2006; Xue et al, JAS, 2008; Hill et al, JAS,
! 2009)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module module_mp_tau_bin

  Use parameters, only : num_h_moments, num_h_bins, &
       nspecies, mom_names, h_names, mom_units, max_char_len, &
       num_aero_moments,num_aero_bins, aero_mom_init
  Use column_variables
  Use physconst, only : p0, this_r_on_cp=>r_on_cp, pi
  Use mphys_tau_bin_declare
  Use diagnostics, only: save_dg, i_dgtime, save_binData
  Use common_physics, only: qsaturation

  Use switches, only: l_sediment, mphys_var, l_fix_aerosols
  Use switches_bin

  IMPLICIT NONE
  
  !generic 1d and 2d KiD arrays for diagnostics
  real, dimension(KKP) :: field ! field for 1-D KiD diagnostics
  real, dimension(KKP,JMINP:JMAXP):: field_2d ! field for 1-D KiD diagnostics

  contains

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE tau_bin(I,TH,Q,STH,SQ,DT,RDT)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!      
      IMPLICIT NONE
!     
!variables for fixed radiative flux
!CALL RAD_MOD      
       REAL ::                                                          &
     &  RADFAC                                                          &
                  ! factor to change from flux to heating rate
     & ,QL_1D                                                           &
                  ! q_l
     & ,SWF_NET                                                         &
                  ! Net SW flux (W/m2)
     & ,SWF_DN                                                          &
                  ! Downward SW flux (W/m2)
     & ,SWF_UP    ! Upward SW flux (W/m2)
      REAL, DIMENSION(jjp) ::                                           &
     &  Z_BL_TOP  ! BL top used in radiation (m)
      REAL, DIMENSION(JJP,KKP) :: QL, QV
      REAL, DIMENSION(jjp,kkp) :: fnt_lw, sth_lw
!     variables for fixed cooling rate
      REAL :: cool_kday,sth_lwmax,qttol
      INTEGER,DIMENSION(JMINP:JMAXP) ::  K_BL_TOP 
! Subprogram arguments
! IN  
      INTEGER, INTENT(IN) :: I
      REAL, INTENT(IN) ::                                               &
     &    DT                                                            &
                   ! timestep 
     &   ,RDT      ! 1/timestep
      REAL, INTENT(IN), DIMENSION(JMINP:JMAXP,KKP) ::                       &
     &    TH       ! potential temperature perturbation
      REAL, INTENT(IN), DIMENSION(JMINP:JMAXP,KKP,NQP) ::                   &
     &    Q        ! moisture fields
      
! INOUT    
      REAL, INTENT(INOUT), DIMENSION(JMINP:JMAXP,KKP) ::                    &
     &    STH      ! potential temperature tendency 
      REAL, INTENT(INOUT), DIMENSION(JMINP:JMAXP,KKP,NQP) ::                &
     &    SQ        ! moisture fields' tendency      
!
      integer ijj,ikk
!
!END of RADIATION DECLARATIONS 
!
! Local variables
!
      
      REAL,DIMENSION(JMINP:JMAXP,KKP) :: DQLDT,DNQLDT,QLOLD,QLNEW           &
     &                             , NQLOLD,NQLNEW, RAINOLD,RAINNOOLD   &
     &                             ,RAINNEW,RAINNONEW,DRAINDT,DRAINNODT &
     &                             , QSATMIX                            &
                                             !saturation mixing ratio
     &                             , RH, RH2                            &
                                             !relative humidity 
     &                             , DTHDT                              &
                                           ! microphysical tendency 
                                           ! in potential temperature
     &                             , DQVDT                              &
                                           ! microphysical tendency 
                                           ! in water vapour
     &                             , TBASE                              &
                                           !value of T around which 
                                           !saturation functions base
     &                             ,QVOLD,QTOT 
!      
      REAL,DIMENSION(JMINP:JMAXP,KKP,LN2) :: CCNfrac_loc 
      REAL :: THOLD                                                     &
                    !potential temp pertubation
     & ,      THNEW                                                     &
                    !potential temp pertubation after bin microphysics
     & ,      PMB                                                       &
                    !reference pressure in millibars
     & ,      vapour                                                    &
                     !total water vapour for the whole grid
     & ,      water                                                     &
                    ! total water for the whole grid
     & , u0, T_OLD, totloss, dccndt,TH2T,qln_tot
      REAL :: TOTCCNNUC, TOTEVAP,totevap2,TOTCCNREG,fein_ccnfrac
      REAL, DIMENSION(JMINP:JMAXP, KKP) :: cell_ccnreg
      REAL, DIMENSION(JMINP:JMAXP, KKP) :: t_0, qst_0
      REAL :: cloudrop
      REAL, DIMENSION(JMINP:JMAXP, KKP) :: ds_0_temp ! supersat @ beginning of the timestep
      REAL :: ds_0      ! dynamic term for supersaturation total trans
      REAL :: QmassFLD,QnumFLD
      real :: t
!
!local integers
!     
      INTEGER, DIMENSION(JMINP:JMAXP) :: KQLINV
      INTEGER :: KQLMAX
      INTEGER J,K,IQ, N
      character(2) :: str2


!Functions
!      
!      REAL,EXTERNAL :: QSATURATION                                      &
                                   !specific humidity for this case
!     &  , RELHUM  !function for relative humidty 
!      
! 
      IF(IMICROBIN == 1.AND.IRAINP == 0) THEN

!Calculate the domain total bin resolved CCN number for regen
!when Ln2 >  3
 
        TOTCCNNUC = 0.0
        DO IQ = 1, Ln2
           CCNNEWTOT(IQ) = 0.0
           CCNNEWAVG(IQ) = 0.0
           DO K = 2, KKP
              DO J = JMINP, JMAXP
                 CCNNEWTOT(IQ) =  CCNNEWTOT(IQ) +(Q(J,K,IAERO_BIN(IQ))+    &
                      (SQ(J,K,IAERO_BIN(IQ))*DT))
              ENDDO
           ENDDO
           TOTCCNNUC = TOTCCNNUC + (CCNORIGTOT(IQ) - CCNNEWTOT(IQ))
        ENDDO
        
     ENDIF

     eps = epsilon(1.0d0)
     
     do K = 2, KKP
        do J = JMINP,JMAXP
           DO IQ=1,LK
!AMKORIG and ANKORIG for use in positivity check of mass and number
!change due to bin microphysics
              AMKORIG(J,K,IQ)=Q(J,K,ICDKG_BIN(IQ)) + &
                   (SQ(J,K,ICDKG_BIN(IQ))*DT)
              ANKORIG(J,K,IQ)=Q(J,K,ICDNC_BIN(IQ)) + &
                   (SQ(J,K,ICDNC_BIN(IQ))*DT)
              if (AMKORIG(J,K,IQ) < eps .or.         &
                   ANKORIG(J,K,IQ) < eps .or.        &
                   AMKORIG(J,K,IQ) > ANKORIG(J,K,IQ)) then
                 AMKORIG(J,K,IQ) = 0.0
                 ANKORIG(J,K,IQ) = 0.0
              endif
           ENDDO
        ENDDO
     ENDDO
       
     if (l_sediment) then
           
        CALL BIN_SEDIMENT(I,DT,AMKORIG,ANKORIG,Q,SQ,RDT)
        
     endif                    ! sedimentation calculation 

       
     totevap = 0.0
     totevap2 = 0.0
     totccnreg = 0.0
     DS_0 = 0.0
     QLOLD(:,:) = 0.0
     NQLOLD(:,:) = 0.0
     QLNEW(:,:) = 0.0
     NQLNEW(:,:) = 0.0
 
     DO K=2,KKP
        PMB = 0.01*PREFN(K)
        DO J=JMINP,JMAXP
! 1. set base values
!
!!a) calc ss at beginning of DT for Get_forcing
           T_0(j,k) = TREF(K) + TH(j,k)*RPREFRCP(K)
           QST_0(j,k) = QSATURATION(t_0(j,k),PMB)
           DS_0 = q(j,k,IQSS) + (SQ(J,K,IQSS)*DT)
            
!b) Set base values after dynamics for use in microphys            
           THOLD=TH(J,K) + (STH(J,K)*DT)
           QVOLD(J,K)=Q(J,K,IQV) + (SQ(J,K,IQV)*DT)
           DO IQ = 1, LK
              QLOLD(J,K)=QLOLD(J,K)+(Q(J,K,ICDKG_BIN(IQ))               &
     &        +(SQ(J,K,ICDKG_BIN(IQ))*DT))
              NQLOLD(J,K)=NQLOLD(J,K)+(Q(J,K,ICDNC_BIN(IQ))               &
     &        +(SQ(J,K,ICDNC_BIN(IQ))*DT))
           ENDDO
            
           TBASE(J,K)=TREF(K) + THOLD*RPREFRCP(K)
!
! 2. calculate saturation functions
!
           QSATPW(J,K)  = QSATURATION(TBASE(J,K),PMB)
!
! 3. call cloudbin to calculate the bin microphysics (i.e.
!   nucleation, condensation, evaporation, collection, breakup)
!          
            CALL CLOUDBIN(I,J,K,Q,SQ,AMKORIG,ANKORIG,QSATPW,RH,TBASE,TREF,    &
     &                  DQVDT(J,K),DT,RDT,PMB,QVOLD(J,K),QLOLD(J,K),    &
     &                  totevap,totccnreg,DS_0)

!
! 4. calculate the change in theta due to bin microphysics            
!
            THNEW=(TBASE(J,K) - TREF(K))*PREFRCP(K)
            DTHDT(J,K)=(THNEW-THOLD)*RDT
            totevap = totevap + CDNCEVAP(J,K)
        ENDDO
        
      ENDDO      
!
!regeneration for single prognostic CCN variable 
      if (.not. l_fix_aerosols) then        
        DO IQ = 1, LN2
          CCNNEWTOT(IQ) = 0.0
          CCNNEWAVG(IQ) = 0.0
          DO J = JMINP, JMAXP
            DO K = 2, KKP
              CCNNEWTOT(IQ) =  CCNNEWTOT(IQ) + (Q(J,K,IAERO_BIN(IQ))+   &
     &          (SQ(J,K,IAERO_BIN(IQ))*DT))
            ENDDO
          ENDDO
          TOTCCNNUC = TOTCCNNUC + (CCNORIGTOT(IQ) - CCNNEWTOT(IQ))
          CCNNEWAVG(IQ) = CCNNEWTOT(IQ)/(JJP*(KKP))
        ENDDO
        DO IQ = 1, LN2
          DO J = JMINP, JMAXP
            DO K = 2, KKP
              IF(CCNNEWAVG(IQ) <  CCNORIGAVG(IQ)                        &
     &           .AND.CDNCEVAP(J,K) >  0.0) THEN
                DCCNDT = CDNCEVAP(J,K)*RDT
                SQ(J,K,IAERO_BIN(IQ)) = SQ(J,K,IAERO_BIN(IQ)) + DCCNDT
                dqn_reg(J,K) = DCCNDT
                CDNCEVAP(J,K) = 0.0
              ELSE
                dqn_reg(j,k) = 0.0
                CDNCEVAP(J,K) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO       
        if (jjp == 1) then 
           field(:) = dqn_reg(jjp,:)
           call save_dg(field,'ccn_reg', i_dgtime, & 
                units='#/kg/s',dim='z')
        else 
           do  K = 2,KKP 
              do J = JMINP, JMAXP
                 field_2d(k, j) = dqn_reg(j,k)
              enddo
           enddo
           call save_dg(field_2d(1:kkp,1:jjp),'ccn_reg', i_dgtime, & 
                units='#/kg/s',dim='z,x')
        endif
      endif

      if (jjp == 1) then
         field(:) = dqn_act(jjp,:)
         call save_dg(field,'ccn_act', i_dgtime, &
              units='#/kg/s',dim='z')
      else
         do  K = 2,KKP 
            do J = JMINP, JMAXP
               field_2d(k, j) = dqn_act(j,k)
            enddo
         enddo
         call save_dg(field_2d(1:kkp,1:jjp), 'ccn_act', i_dgtime, units, dim='z,x')
      endif
         
      DO K = 2,KKP 
         DO J = JMINP,JMAXP
            STH(J,K) = STH(J,K) + DTHDT(J,K)
            SQ(J,K,IQV) = SQ(J,K,IQV)+DQVDT(J,K)
            DO IQ=1,LK
!Call microcheck to update the source fields for bin resolved mass
!and number, and check that small numbers are not causing erroneous 
!values of mass and number that lead to numerical instabilities
      
               IF((Q(J,K,ICDKG_BIN(IQ))+(SQ(J,K,ICDKG_BIN(IQ))*DT))      &
                < 0.0 &
                .or.(Q(J,K,ICDNC_BIN(IQ))+(SQ(J,K,ICDNC_BIN(IQ))*DT))      &
                < 0.0) THEN
                  QLNEW(J,K)=QLNEW(J,K)
                  NQLNEW(J,K)=NQLNEW(J,K)
               ELSE
                  QLNEW(J,K)=QLNEW(J,K)+(Q(J,K,ICDKG_BIN(IQ))               &
                                  +(SQ(J,K,ICDKG_BIN(IQ))*DT))
                  NQLNEW(J,K)=NQLNEW(J,K)+(Q(J,K,ICDNC_BIN(IQ))             &
                                  +(SQ(J,K,ICDNC_BIN(IQ))*DT))
               ENDIF
            ENDDO
            DQLDT(J,K)=(QLNEW(J,K)-QLOLD(J,K))*RDT
            DNQLDT(J,K)=(NQLNEW(J,K)-NQLOLD(J,K))*RDT
         ENDDO
      ENDDO
! 6. calculate effective radius for radiation 
        CALL REFFCALC(Q,SQ,DT,RDT)
 !      
! 7. update diagnostic fields if necessary
!              
       DO K = 2, KKP
         DO J = JMINP, JMAXP
            IF(l_dodgs)THEN
              dth_dt(j,k) = DTHDT(J,K)
              dq_dt(j,k,iqv) = DQVDT(J,K)
            ENDIF
          ENDDO
        ENDDO
!
      END subroutine TAU_BIN

!***********************************************************************
      SUBROUTINE CLOUDBIN(I,J,K,Q,SQ,AMKORIG,ANKORIG,QST,RH,TBASE,TREF, &
           DQVDT,DT,RDT,PMB,QVOLD,QLOLD,totevap,totccnreg,DS_0)
!***********************************************************************
      IMPLICIT NONE
!This routine is the main loop for bin microphysics, it calculates
!the changes in mass and nmber concentrations of CCN and cloud 
!droplets and writes the CCN number into array of Q(J,K,IAERO_BIN),
!Cloud drop mass into Q(J,K,ICDKG_BIN) and cloud drop number into 
!Q(J,K,ICDNC_BIN).
      
      REAL,DIMENSION(JMINP:JMAXP,KKP,NQP)::                                 &
     & SQ                                                               &
          !Source term of the water fields
     &, Q
      REAL,INTENT(INOUT),DIMENSION(JMINP:JMAXP,KKP,LK) ::                      &
      AMKORIG,ANKORIG
      REAL,DIMENSION(JMINP:JMAXP,KKP) ::                                    &
     & QS                                                               &
          !saturation mixing ratio of water calced in Qsaturation
     & , QST,QST2                                                       &
                  !specific humidity calculated in Spechum
     & , RH, RH2 !relative humidity calculated in relhum
      INTEGER J,K,L,I !loop counter 
      
      REAL,PARAMETER :: RRV=461 !gas constant of moist air J/kg/K
      
      REAL,PARAMETER :: AL=597  !latent heat of evaporation and cond
      REAL,PARAMETER :: T00=273.15 !freezing point of water in K,from
                                   !the LES
      
      REAL,PARAMETER :: CPBIN=0.24 !specific heat of water vapour
      REAL,PARAMETER :: AAR=0.24E-3

!      REAL, PARAMETER :: DG1 = 0.10E-4
!      REAL, PARAMETER :: SG1 = 1.5

!      REAL :: DG1 ! mean aerosol diameter in cm (not m!) = 0.10E-4 !
!      REAL :: SG1 ! = 1.5
      
      REAL,DIMENSION(JJP,KKP) :: AMY
      REAL,DIMENSION (JMINP:JMAXP,KKP):: TBASE!temp calced in Tau_bin
                                !using tref+th(1/exner function). passed into 
                                !all other routines as TEMP
      REAL, DIMENSION(JMINP:JMAXP,KKP) :: CDNCEVAP1, CCNTOT, CCNTOTNEW
      REAL :: CN1
      REAL,DIMENSION(KKP)::TREF !reference temperature from LES
      
      REAL DT  !timestep from tau_bin
      REAL RDT !1/timestep from tau_bin
      REAL DM !the condensated mass calced in EVAP and COND routines
      real dm_cloud_evap ! the evaporated mass
      real dm_rain_evap ! the evaporated mass
      REAL DDDD!supersat percentage

      REAL DCCNDT!tendency in aerosol due to bin microphysics
      REAL DAMKDT!tendency in bin resolved mass cloud due to microphysics
      REAL DANKDT!tendency in bin resolved number cloud due to microhysics
      REAL QLNEW!post bin microphysics total liquid water kg/kg
      REAL NQLNEW!post bin microphysics total liquid water number conc /kg
      REAL QVNEW!vapour mixing ratio after the nucleation, before code
                !proceeds to evap or cond, after evap or cond, this var
                !QV + change due to nucleation and cond or evap. It is
                !used to calc DQVDT
      REAL QVOLD!vapor mixing ratio before bin micro
      REAL DQLDT!change in liquid water as a result of bin microphys
      REAL DNQLDT!change in number concentration of liquid water
      REAL DQVDT!change in water vapour due to microphysics
      REAL PMB   ! pressure in mb
      REAL DMCCDT!change in mass over due to nucleation
      REAL SQmassbef,sqnumbef
      REAL QLOLD!liquid water mass before microphysics
      REAL NQLOLD!liquid water number concentration before microphysics
      REAL,DIMENSION(JMINP:JMAXP,KKP,LK)::BIN_AVG, BIN_AVGold
      REAL,DIMENSION(JMINP:JMAXP,KKP,LK)::AMK_init,ANK_init,                &
     &                                AMK_res, ANK_res,                 &
     &                                AMK_sub, ANK_sub
      REAL AVG, masstemp, notemp, ANKwrong, ANKoldwrong, evapbin
!     code added 1/5/08 - AH - changing declarations of get_forcing
!     variables to arrays, for debugging+clarity purposes      
      REAL,DIMENSION(JMINP:JMAXP,KKP) :: TAU                                &
                                         !microphys forcing
     &                ,VSW                                              &
     &                ,EN                                               &
                            !estimated ETA for end of timestep
     &                ,EA   !average ETA for the timestep
      REAL :: QST_nuc,ds_force,temp_nuc,mpos
      !added 28/04/08 for new inhom mixing method
      REAL :: ds_0      ! dynamic term for supersaturation total trans
      REAL,DIMENSION(JMINP:JMAXP,KKP) :: tau_dum
      REAL ::  delm
!      
      REAL,DIMENSION(JMINP:JMAXP,KKP) :: NCC,MCC,AN1_bef_adv,AN1_aft_adv
      REAL totnuc, totcdnc, totevap, totccnreg,DtraceDT

      real, dimension(kkp) :: auto_mass, auto_num, auto_con_mass, d_rmass
      real :: rmass_tot_orig, rmass_tot_new
      real :: t_test, dtcalc
      
      LOGICAL L_ERROR
      INTEGER ISTOPU
      Integer LC !the smallest bin that will nucleate, used in subnuc
                                !and regeneration
      INTEGER :: LT, IT, inhom_evap
      INTEGER :: ccn_pos, cloud_pos, count, loop_count
      character(2) :: str2
       
!      REAL,EXTERNAL::QSATURATION,XACT
! Set values of bins before microphysics, these values do not include the
! effect of dynamics from this timestep. I do not think this is needed 
! as all I want is a tendency due ti microphysics. The tendency for aerosol
! is calculated in SUBNUC, as this is where nucleation is calculated

      if (jjp == 1) then
         call save_dg(k,rhon(k),'density', & 
              i_dgtime, units='kg/m3',dim='z')
      endif
         

      QVNEW = QVOLD
      DO L=1,LN2
        CCNOLD(J,K,L)=Q(J,K,IAERO_BIN(L)) + (SQ(J,K,IAERO_BIN(L))*DT)
        CCN(J,K,L) = CCNOLD(J,K,L) !this line is required so that 
        !CCN does not equal zero if DS < 0, i.e. for calc of aerosol
                                !tendency 
      ENDDO

      DO L=1,LK
        AMKOLD(J,K,L)=AMKORIG(J,K,L)
        ANKOLD(J,K,L)=ANKORIG(J,K,L)        
      ENDDO

!****************************************************************
!          DS() IS THE SUPERSATURATED QUANTITY
!****************************************************************
      DS(J,K)=QVNEW-QST(J,K)
      DUS(J,K)=DS(J,K)/QST(J,K)
      QST_nuc = QST(J,K) !used in activation calculation
      TEMP_NUC = TBASE(J,K) 
      DS_FORCE = 0.0
      
      IF (DT.GT.4.0) THEN
         LT = CEILING(DT/1.0)
         DTcalc = DT/REAL(LT)
      ELSE
         LT=1
         DTcalc=DT
      ENDIF
     
      CDNCEVAP(j,k) = 0.0
      am1_diag(j,k) = 0.0
      an1_diag(j,k) = 0.0
     
!AH 0410 - moving all unit conversion outside
!        the calls for get_force, cond and evap_new
!        to minimise precision errors, also don't change
!        amkold or ankold
      DO L=1,LK
         IF(AMKOLD(J,K,L) > eps .AND.       &
              ANKOLD(J,K,L) > eps )THEN
            am0(l)=(amkold(j,k,l)*Rhon(k))/1.E3 
            an0(l)=(ankold(j,k,l)*Rhon(k))/1.E6
            AMN(L)=am0(L)/an0(L)
         ELSE
            AMN(L)=0.
            am0(l)=0.
            an0(l)=0.
         ENDIF
         am1_diag(j,k) = am1_diag(j,k)+am0(l)
         an1_diag(j,k) = an1_diag(j,k)+an0(l)
      ENDDO 
      
      DO it = 1,lt

         IF(IINHOM_mix == 0) then
            VSW(J,K) = MIN(1.,-1.E+10*MAX(DS(j,k),DS_0))
            CALL GET_FORCING(J,K,QST,TBASE,DS_0,DS,                      &
                 AM0,AN0,DTcalc,TAU,EN,EA,VSW(J,K))

            inhom_evap = 0
                
         ENDIF
         
         DM=0.0
         dm_cloud_evap = 0.0
         dm_rain_evap = 0.0
         NCC(J,K)=0.0
         MCC(J,K)=0.0
         DUS1(J,K)=0.0
         tevap_bb(j,k) = 0.0
         tevap_squires(j,k) = 0.0
         t_mix(j,k) = 0.0
         t_mix_approx(J,K) = 0.0
         da_no(J,K) = 0.0
         da_no_rat(J,K) = 0.0
         r_int(j,k) = 0.0
         r_bar(j,k) = 0.0
         ssat(j,k) = 0.0
         evap_a(j,K) = 0.0

         ssat(J,K) = EA(J,K)/QST(J,K)

         DO L=1,LK
            DIEMC(L)=0.0
            DIENC(L)=0.0
            DIEMD(L)=0.0
            DIEND(L)=0.0
            ANK(J,k,L) = 0.0
            AMK(j,k,l) = 0.0
         ENDDO         
         IF (IINHOM_MIX == 0) then          
            DS_force = tau(j,k)
            tau_dum(j,k) = tau(j,k)
         endif

         IF (DS_force >  eps) THEN        
!*****************************************************************        
!        CONDENSATION
!     COND RECEIVES MKOLD,NKOLD RETURNS MK,NK
!****************************************************************
            AN1OLD(J,K) = 0.0
            AM1OLD(j,k) = 0.0
            DO L=1,LK
               AN1OLD(J,K) =AN1OLD(J,K)+AN0(l) 
               AM1OLD(j,K) = AM1OLD(j,K)+am0(l)
            ENDDO

!***************************************************************      
            CALL COND_new(J,K,DM,TBASE,QST,RH,Q,SQ,DT,RDT,PMB,QVNEW,      &
                 TAU_dum,it,LT)
!***************************************************************      

            CALL REBIN(J,K)
            
            AN1(J,K) = 0.0
            AM1(j,k) = 0.0
            DO L=1,LK 
               AN1(J,K) = AN1(J,K) + (ANK(J,K,L))
               AM1(J,K) = AM1(j,k) + (AMK(j,k,l))
            ENDDO
               
            IF(ABS((AN1(J,K))-(AN1OLD(J,K))) >  1.) THEN
               print *, 'cond conservation prob after rebin'
               print *, 'AN1OLD', AN1OLD(j,k),k,j
               print *, 'AN1', AN1(j,k),k,j
            ENDIF

       ! new code to calc the change in total mass and num from 
       ! rain dropsize bins due to condensation. This is 
       ! sort-of equivalent to the autoconversion in the Bulk scheme   
            auto_mass(k) = 0.0
            auto_num(k) = 0.0
            DO L = 16, LK
                  
               auto_mass(k) = auto_mass(k) + (AMK(J,K,L) - AM0(L))
               auto_num(k) = auto_num(k) + (ANK(J,K,L) - AN0(L))
               
            enddo
               
            if (auto_mass(k) < 0.0 .or. auto_num(k) < 0.0 ) then
               auto_mass(k) = 0.0
               auto_num(k) = 0.0
            endif
            
            if (jjp == 1) then
               call save_dg(k,(auto_mass(k)*1.e3/rhon(k))/dt,'cond_c_r_mass', & 
                    i_dgtime, units='kg/kg/s',dim='z')    
               call save_dg(k,(auto_num(k)*1.e6/rhon(k))/dt,'cond_c_r_num', &
                    i_dgtime, units='#/kg/s',dim='z')
            else
               call save_dg(k,j,(auto_mass(k)*1.e3/rhon(k))/dt,'cond_c_r_mass', & 
                    i_dgtime, units='kg/kg/s',dim='z,x')    
               call save_dg(k,j,(auto_num(k)*1.e6/rhon(k))/dt,'cond_c_r_num', &
                    i_dgtime, units='#/kg/s',dim='z,x')
            endif

         ELSEIF(DS_force < -eps) then !DS_force < 0.0
!****************************************************************
!                     EVAPORATION
!     EVAP RECEIVES MKD,NKD RETURNS MK,NK
!***********************************************************
            AN1OLD(J,K) = 0.0
               
            DO L=1,LK
               AN1OLD(J,K) =AN1OLD(J,K)+AN0(l) 
            ENDDO
               
            IF (AN1OLD(J,K) > eps) then
               
               CALL EVAP_new(I,J,K,DM,TBASE,QST,RH,Q,SQ,DT,RDT,PMB,      &
                    QVNEW,TAU_dum,EA,it,LT,inhom_evap,delm)
               ! 
               CALL REBIN(J,K)
               !  
               AN1(J,K) = 0.0
               DO L=1,LK 
                  AN1(J,K) = AN1(J,K) + (ANK(J,K,L))        
               ENDDO
               
            ELSE
               DO L=1,LK
                  AMK(J,K,L)=AM0(L)
                  ANK(J,K,L)=AN0(L)
               ENDDO
            ENDIF
         ELSE
            
            DO L=1,LK
               AMK(J,K,L)=am0(L)
               ANK(J,K,L)=an0(L)
            ENDDO
            
         ENDIF
                  
            !update fields 
         an1old(j,k) = 0.0
         an1(j,k) = 0.0
         DO L=1,LK
            ! n0 and m0 in cjs units for next lt iteration if
            ! needed
            AN0(L)=ank(j,k,L)
            AM0(L)=amk(j,k,L)
            IF(AM0(L) >  eps .AND.                       &
                 AN0(L) > eps)THEN
               AMN(L)=AM0(L)/AN0(L)
            ELSE
               AMN(L)=0.
            ENDIF
!
! AH 0410 - convert mass and number back to kg/kg and #/kg 
!           respectively for source term calc and thermodynamic
!           updates
!             
            amk(j,k,l)=amk(j,k,l)*1.e3/rhon(k)      
            ank(j,k,l)=ank(j,k,l)*1.e6/rhon(k)     
            if(amk(j,k,l) < eps.or.                 &
                 ank(j,k,l) < eps) then
               amk(j,k,l) = 0.0
               ank(j,k,l) = 0.0
            endif
            an1old(j,k) = an1old(j,k) + ankorig(j,k,l)
            an1(j,k) = an1(j,k) + ank(j,k,l)

         ENDDO

! AH 0410 - work out total mass and number change due to cond
!           and evap
         dm = 0.0
         do l = 1,lk
            dm = dm + (amk(j,k,l) - amkorig(j,k,l))
         enddo
! AH 0410 - Update thermodynamic field once evap and cond
!           finished
         if ( it .ne. lt )  DS_0 = EN(j,k)
         TBASE(J,K)=TBASE(J,K)+AL*DM/CPBIN
         QVNEW = QVNEW - DM      
         QST(J,K)=QSATURATION(TBASE(J,K),PMB) 
         DS(J,K)=QVNEW-QST(J,K)

      ENDDO          !end of iteration over LT  
      
!subtract total number after evaporation, from total number before evap
!this value tells the number of droplets that have evaporated completely  
!
      IF (AN1(J,K) >  AN1OLD(J,K)) THEN
         CDNCEVAP(J,K) = 0.0
      ELSE
         CDNCEVAP(J,K) = (AN1OLD(J,K) - (AN1(J,K)))
      ENDIF
      
      if (jjp > 1) then
         !    diagnostics for total cond/evap rate liquid water (change in mass)
         call save_dg(k,j,dm/dt,'dm_ce', i_dgtime, &
              units='kg/kg/s',dim='z,x')
         !    cond/evap rate of cloud 
         do l = 1,lk_cloud
            dm_cloud_evap = dm_cloud_evap + (amk(j,k,l) - amkorig(j,k,l))
         enddo
         call save_dg(k,j,dm_cloud_evap/dt,'dm_cloud_ce', i_dgtime, &
              units='kg/kg/s',dim='z,x')
         !     cond/evap rate of rain 
         do l = lk_cloud+1, lk
            dm_rain_evap = dm_rain_evap + (amk(j,k,l) - amkorig(j,k,l))
         enddo
         call save_dg(k,j,dm_rain_evap/dt,'dm_rain_ce', i_dgtime, &
              units='kg/kg/s',dim='z,x')
      else
         !    diagnostics for total cond/evap rate liquid water (change in mass)
         call save_dg(k,dm/dt,'dm_ce', i_dgtime, &
              units='kg/kg/s',dim='z')
         !    cond/evap rate of cloud 
         do l = 1,lk_cloud
            dm_cloud_evap = dm_cloud_evap + (amk(j,k,l) - amkorig(j,k,l))
         enddo
         call save_dg(k,dm_cloud_evap/dt,'dm_cloud', i_dgtime, &
              units='kg/kg/s',dim='z')
         !     cond/evap rate of rain 
         do l = lk_cloud+1, lk
            dm_rain_evap = dm_rain_evap + (amk(j,k,l) - amkorig(j,k,l))
         enddo
         call save_dg(k,dm_rain_evap/dt,'dm_rain', i_dgtime, &
              units='kg/kg/s',dim='z')        
      endif
!               
            

! do activation after cond/evap, using updated supersat for timestep
       EA(J,K) = DS(J,K)

      IF(EA(J,K) >  0.0) THEN
 
        AN1OLD(J,K) = 0.0
        AM1OLD(J,K) = 0.0
        CCNTOT(J,K) = 0.0
        DO L = 1, LK
          AN1OLD(J,K)  = AN1OLD(J,K) +ANK(J,K,L)
          AM1OLD(J,K)  = AM1OLD(J,K) +AMK(J,K,L)
        ENDDO

        DO L = 1, LN2
          CCNTOT(J,K) =  CCNTOT(J,K) + CCNOLD(J,K,L)
        ENDDO
   
        AN1(J,K) = 0.0
        AM1(J,K) = 0.0
        amkcc(1) = 0.0
        ankcc(1) = 0.0
        
        DO L = 1, LK
          AN1(J,K) = AN1(J,K) + ANK(J,K,L)
        ENDDO

        AM1(J,K) = AMK(J,K,1)
        DDDD=(EA(J,K)/QST(J,K))*100
        
        ccn_pos = 1
        cloud_pos = 1      
         
        CN1 = CCN(j,k,ccn_pos)/rhon(k) ! convert to /kg
    
        if (.not. l_fix_aerosols) then
        
           AN2(J,K) = MAX(0.,                                            &
     &       XACT(tbase(j,k),DDDD,DG1,SG1,dcrit(j,k))*                  & 
              CN1)

           !adjust CCN to show the removal of CCN by activation
            CCN(j,k,ccn_pos) = (CN1 - AN2(J,K))*rhon(k) ! convert to /m^3

        else
           AN2(J,K) = MAX(0.,                                            &
     &       XACT(tbase(j,k),DDDD,DG1,SG1,dcrit(j,k))*                  &
              CN1-AN1(J,K))
!     &       CCN(J,K,ccn_pos)-AN1(J,K))
           
        endif

        dqn_act(J,K) = AN2(J,K)/dt
        
        ANK(J,K,1) = ANK(J,K,1) + AN2(J,K)
        ! the 0.25 factor attempts to accomodate for the fact
        ! that not all CCN will grow to the size of the 1st bin
        AMK(J,K,1) = AMK(J,K,1) + AN2(J,K)*XK(1)*0.25
        AMK(J,K,1) = MAX(ANK(J,K,1)*XK(1)*1.01, AMK(J,K,1))
       
        MCC(J,K) = AMK(J,K,1) - AM1(J,K)
        
        AM1(j,k) = 0.0
        DO L = 1, LK
          AM1(J,K) = AM1(J,K) + AMK(J,K,L)
        ENDDO
            
!calculate the new termodynamic fields following activation          
        TBASE(J,K)=TBASE(J,K)+AL/CPBIN*MCC(J,K)
        QVNEW = QVNEW - MCC(J,K)
        QST(J,K)=QSATURATION(TBASE(J,K),PMB)                        
     
      ENDIF                          ! end of do activation

      AN1(j,k) = 0.0
      DO L = 1, LK
        AN1(J,K) = AN1(J,K) + ANK(J,K,L)
      ENDDO

! AH 0410 - to update supersat for advection. The update of ss is 
!           performed here not in stepfields (code commented out in
!           stepfields)
!
      DS(J,K) =  QVNEW - QST(J,K)
      q(j,k,iqss) = DS(j,k)

!*************************************************************
!        UPDATING CHANGE IN MASS CAUSED BY MICROPHYSICS
!             AFTER  CONDENSATION-EVAPORATION
!*************************************************************
      DO L=1,LK
         DIEMC(L)=AMK(J,K,L)-AMKORIG(J,K,L)!AMKORIG is the mass prior to 
                !bin micro, therefore DIEMC (old variable name) is 
                !mass resulting from bin micro
         DIENC(L)=ANK(J,K,L)-ANKORIG(J,K,L)!ANKORIG is the number prior to 
                !bin micro, therefore DIENC (old variable name) is 
                !number resulting from bin micro
      ENDDO

      if (IRAINBIN == 1.AND.IMICROBIN == 1) then
!
!************************************************************
!            COLLECTION + BREAKUP
!************************************************************
!
! AH 0410 - The collision-coalescence code uses the original format 
!           of the TAU model (so is an older version than that used 
!           in RAMS). AMKORIG and ANKORIG, which are the same as 
!           AMKOLD and ANKOLD) are used to initialise calculations in
!           SXY, SCONC and BREAK. Unit conversions are performed 
!           in the respective subroutines (this is not neccesary and 
!           should be moved out so that AN0 and AM0 are used
! 
          rmass_cw(k) = 0.0

          loop_count = 0
          count = 0
          DO L=1,LK
            AMK(J,K,L)=AMKORIG(J,K,L)
            ANK(J,K,L)=ANKORIG(J,K,L)
          ENDDO

          AM1(J,K)=0.0
          AN1(J,K)=0.0
          DO L=1,LK
            AM1(J,K)=AM1(J,K)+AMKORIG(J,K,L)
            AN1(J,K)=AN1(J,K)+ANKORIG(J,K,L)
          ENDDO


          IF(AM1(J,K) >  lk*eps .AND. AN1(J,K) > lk*eps ) THEN
 
            if (l_coll_coal) then
               

               CALL SXY(J,K,DT)
               CALL SCONC(J,K,DT)
             
               if (l_break) then
                  CALL BREAK(J,K,DT)
               endif
 
            endif
      
!*************************************************************
!        UPDATING CHANGE IN MASS CAUSED BY MICROPHYSICS
!             AFTER  COLLECTION & BREAKUP
!*************************************************************
            DO L=1,LK
               DIEMD(L)=AMK(J,K,L)-AMKORIG(J,K,L)
               DIEND(L)=ANK(J,K,L)-ANKORIG(J,K,L)
            ENDDO

            d_rmass(k) = 0.0
            
            do l = 16, lk                 
               d_rmass(k) = d_rmass(k) + DIEMD(L)
            enddo
            
            if (jjp == 1) then
               call save_dg(k,d_rmass(k)/dt,'drain_tot', i_dgtime, &
                    units='kg/kg/s',dim='z')
            else
               call save_dg(k,j,d_rmass(k)/dt,'drain_tot', i_dgtime, &
                    units='kg/kg/s',dim='z,x')
            endif
          ELSE
            DO L=1,LK
             DIEMD(L) = 0.0
             DIEND(L)= 0.0
           ENDDO 
         ENDIF        
      ENDIF    !if IRAINBIN == 1.AND.IMICROBIN == 1

!***********************************************************
!     UPDATING MASS & CONCENTRATION AFTER MICROPHYSICS
!***********************************************************
!      AM1(J,K)=0.0
!      AN1(J,K)=0.0
      DO L=1,LK
        IF(IRAINBIN == 1.AND.IMICROBIN == 1) THEN
        !add microphysical change resulting from  nuc,cond/evap,
        !collision-coalscence, breakup and sedimentation
          AMK(J,K,L)=DIEMC(L)+AMKORIG(J,K,L)+DIEMD(L)+QL_SED(J,K,L)
          ANK(J,K,L)=DIENC(L)+ANKORIG(J,K,L)+DIEND(L)+QLN_SED(J,K,L)
        ELSE!only add microphysical change resulting from  nuc,cond/evap
          AMK(J,K,L)=DIEMC(L)+AMKORIG(J,K,L)
          ANK(J,K,L)=DIENC(L)+ANKORIG(J,K,L)
        ENDIF
      ENDDO
!   
!     Calculate the change in "rain bin" mass and number
!     (i.e. bin greater than  100 microns diameter)
      if (j == 0) then
        rmass_tot_orig = 0.0
        rmass_tot_new = 0.0
      !  d_rmass(k) = 0.0
        auto_con_mass(k) = 0.0

        do l=16,lk
          rmass_tot_orig = rmass_tot_orig + AMKORIG(J,K,L)
          rmass_tot_new = rmass_tot_new + AMK(J,K,L)
      enddo

     ! d_rmass(k) = rmass_tot_new - rmass_tot_orig    

      auto_con_mass(k) =  d_rmass(k) - rmass_cw(k)

      if (jjp == 1) then
         call save_dg(k,auto_con_mass(k)/dt,'drain_auto', i_dgtime, &
              units='kg/kg/s',dim='z') 
         
         call save_dg(k,rmass_cw(k)/dt,'drain_rain', i_dgtime, &
              units='kg/kg/s',dim='z') 
      else
         call save_dg(k,j,auto_con_mass(k)/dt,'drain_auto', i_dgtime, &
              units='kg/kg/s',dim='z,x') 
         
         call save_dg(k,j,rmass_cw(k)/dt,'drain_rain', i_dgtime, &
              units='kg/kg/s',dim='z,x') 

      endif
       
      endif
!**********************************************************
!     Update mass and conc, tendencies due to micro and end
!**********************************************************
!First update change in aerosol (if aerosol is not fixed)
       if (.not. l_fix_aerosols) then
         DO L = 1,LN2
            DCCNDT=(CCN(J,K,L)-CCNOLD(J,K,L))*RDT
            IF ((Q(J,K,IAERO_BIN(L))+(SQ(J,K,IAERO_BIN(L))+DCCNDT)*DT)     &
     &         <  0.0) THEN
               DCCNDT=-0.9999*CCNOLD(J,K,L)*RDT
            ENDIF
            SQ(J,K,IAERO_BIN(L))=SQ(J,K,IAERO_BIN(L))+DCCNDT
         ENDDO
       endif

       DO L=1,LK         
!AH 0410 - Calculate the change in mass and number due to microphysics i.e. 
!          nucleation, cond, evap, collection and sedimentation
!
          IF(IRAINBIN == 1.AND.IMICROBIN == 1) THEN
           ! sum the effects of cond/evap+coll/coal+sedi
           ! on mass and number
             DAMKDT=(DIEMC(L)+DIEMD(L)+QL_SED(J,K,L)) * RDT
             DANKDT=(DIENC(L)+DIEND(L)+QLN_SED(J,K,L)) * RDT
          ELSE
             ! only use the term from cond/evap (and nuc)
             DAMKDT = DIEMC(L) * RDT
             DANKDT = DIENC(L) * RDT
          ENDIF

!AH 0410 - Call microcheck to update the source fields for bin resolved mass
!          and number with DAMKDT and DANKDT
!
          CALL MICROcheck(J,K,L,DT,AMKORIG(J,K,L), &
               ANKORIG(j,k,l),SQ(J,K,ICDKG_BIN(L)), &
               SQ(J,K,ICDNC_BIN(L)),DAMKDT,DANKDT,RDT,Q(J,K,ICDKG_BIN(L)) &
               ,Q(J,K,ICDNC_BIN(L)))

       ENDDO
      
      dD(1)=xkk1(1)*2.0*1.e6
      do l = 2, lk 
         dD(l)=(xkk1(l)-xkk1(l-1))*2.0*1.e6
      end do

      call save_binData((xkk1(:)*2.0*1.e6),'bins_D_upper', &
           units='microns', longname='upper bound of droplet diameter')
      call save_binData(xk(:),'bins_mass_upper', &
           units='kg', longname='upper bound of droplet mass')
      call save_binData(xkmean(:),'bins', &
           units='kg', longname='mean droplet mass')
      call save_binData(((xkmean(:)*6./3141.59)**(1./3.)*1.e6), &
           'bins_D', units='microns' &
           , longname='mean droplet diameter')
      call save_binData(dD(:), 'dD', units='microns' &
           , longname='width of bin')

! Finally calc the microphys change in QV
      DQVDT=(QVNEW-QVOLD)*RDT
!        
      END subroutine CLOUDBIN

!*******************************************************************
      SUBROUTINE SXY(KR,II,DT)
!*******************************************************************
! * *
! * SUBROUTINE SXY CALCULATES THE EVOLUTION OF THE SPECTRUM FOR
! * GOLOVINS KERNEL (X+Y) USING THE METHOD OF MOMENTS
! * THE CONNECTION BETWEEN MOMENTS CAN EITHER BE TAKEN AS CONSTANT
! * OR CALCULATED IN EACH CATEGORY IN EACH TIME STEP (AP1(K)).
! * NORMALIZATION OF THE APPROXIMATED DISTRIBUTION FUNCTION
! * CAN EITHER BE DONE WITH MOMENTS MK AND NK (PSI2(K) AND F2(K))
! * OR WITH MOMENTS MK AND ZK (PSI (K) AND F(K))
! * CATEGORIES 1-33 WITH 1-15
! * * * * * *
      IMPLICIT NONE

      INTEGER KR,II!loop counters for J and K respect!, different because
      !an arse to change the whole code (lazy but in a rush!)
      INTEGER K,L !loop counter for cloud bins, confused yet - I am! I
      !think more than one person wrote this code
      INTEGER I,LRK
      REAL :: DT   ! timestep 
      REAL :: AVG, sum, sun
!      
      DO L=1,LK
         CM(L)=AMKORIG(KR,II,L)*RHON(II)/1.E+03
         CN(L)=ANKORIG(KR,II,L)*RHON(II)/1.E+06

         IF(CM(L) >  eps .AND.CN(L) > eps)THEN
            SCM(L)=CM(L)/CN(L)
        !check to see whether SCM (mean single particle mass in a 
        !bin) is larger than twice maximum bin if so not used in
        !collision - adrian 01/06/05
            IF(SCM(L) >  X_BIN(L+1)) THEN
               SCM(L) = X_BIN(L+1) !upper bin boundary
            ENDIF
            IF (SCM(L) <  X_BIN(L)) THEN
               SCM(L) = X_BIN(L) !lower bin boundary
            ENDIF
         ELSE
            SCM(L)=0.
         ENDIF
      ENDDO ! end loop over LK bins
      
      DO K=1,LK
        IF (SCM(K) <  X_BIN(K).OR.SCM(K) >  (2*X_BIN(K)).OR.SCM(K) == 0.0) THEN
          AP=1.0
        ELSE
          AP=0.5*(1.+3.*(X_BIN(K)/SCM(K))-2*((X_BIN(K)/SCM(K))**2.))
        ENDIF
!        

        AZ0(K)=AP*SCM(K)*CM(K)
        AM3(K)=(AP**3)*(SCM(K)**2)*CM(K)
        AM4(K)=(AP**6)*(SCM(K)**3)*CM(K)
        PSI2(K)=2./X_BIN(K)*(CM(K)/X_BIN(K)-CN(K))
        F2(K)=2./X_BIN(K)*(2.*CN(K)-CM(K)/X_BIN(K))
!
        IF (SCM(K) <  X_BIN(K)) THEN
           F2(K)=2.*CN(K)/X_BIN(K)
           PSI2(K)=0.0
        ENDIF
        !        
        IF (SCM(K) >  X_BIN(K+1)) THEN
           PSI2(K)=2.*CN(K)/X_BIN(K)
           F2(K)=0.0
        ENDIF
!        
        SM1=0.
        SM2=0.
        SM3=0.
        SM4=0.
        SM5=0.
        SN1=0.
        SN2=0.
        SN3=0.
        SN4=0.
!        
       DO I=K,LK ! do bin K to bin LK
          !
          IF(K >= 16) EXIT ! exit I = K, LK
          !
          SM5=SM5+KBAR(I,K)*(AZ0(K)*CN(I)+CM(K)*CM(I))
          SN4=SN4+KBAR(I,K)*(CN(K)*CM(I)+CM(K)*CN(I))
       ENDDO ! enddo bin K to bin LK


        IF(K == 1)THEN
          CM(K)=CM(K)-SM5*DT
          CN(K)=CN(K)-SN4*DT
          !
          CYCLE ! execute next iteration of K = 1, LK
          !
        ENDIF
        IF(K-1 < 16) then
          SM3=KBAR(K-1,K-1)*(AZ0(K-1)*CN(K-1)+CM(K-1)**2) 
        ENDIF
        DO I=1,K-1
           !
           IF(I >= 16) EXIT ! exit K = 1, LK 
           !
            SM2=SM2+KBAR(K,I)*(4.*X2(K)*PSI2(K)*CM(I)+                  &
     &      X_BIN(K)/2.*(4.*PSI2(K)+F2(K))*AZ0(I)-(2.*PSI2(K)-F2(K))*AM3(I) &
     &      -1./(2.*X_BIN(K))*(PSI2(K)-F2(K))*AM4(I)+PSI2(K)*AM3(I))

            SN3=SN3+KBAR(K,I)*(2.*X_BIN(K)*PSI2(K)*CM(I)+0.5*F2(K)*AZ0(I)   &
     &      -1./(2.*X_BIN(K))*(PSI2(K)-F2(K))*AM3(I))
 
        ENDDO
        DO I=1,K-1
           !
           IF(I >= 16) EXIT ! exit I = 1, K-1
           !
           SM4=SM4+KBAR(K,I)*(CN(K)*AZ0(I)+CM(K)*CM(I))
        ENDDO
        !
        IF(K-1 < 16)THEN
          SN2=KBAR(K-1,K-1)*CN(K-1)*CM(K-1)
        ENDIF
        !
        IF(K == 2)THEN
          CM(K)=CM(K)+(SM4+SM3-SM5-SM2)*DT
          CN(K)=CN(K)+(SN2-SN3-SN4)*DT
          !
          CYCLE ! execute next iteration of K = 1, LK
          !
        ENDIF
        DO I=1,K-2
           !
           IF(I >= 16) EXIT ! exit loop I = 1, K-2
           !
            SM1=SM1+KBAR(K-1,I)*(4.*X2(K-1)*PSI2(K-1)*CM(I)+            &
     &      X_BIN(K-1)/2.*(4.*PSI2(K-1)+F2(K-1))*AZ0(I)-                    &
     &      (2.*PSI2(K-1)-F2(K-1))*AM3(I)+PSI2(K-1)*AM3(I)              &
     &      -1./(2.*X_BIN(K-1))*(PSI2(K-1)-F2(K-1))*AM4(I))
            SN1=SN1+KBAR(K-1,I)*(2.*X_BIN(K-1)*PSI2(K-1)*CM(I)+0.5*F2(K-1)  &
     &      *AZ0(I)-1./(2.*X_BIN(K-1))*(PSI2(K-1)-F2(K-1))*AM3(I))
!
         ENDDO
         
         ! update mass and number fields with terms from collision
         !
         IF (K == LK) THEN
            CM(LK)=CM(LK)+(SM1+SM3+SM4)*DT
            CN(LK)=CN(LK)+(SN1+SN2)*DT
         ELSE
            CM(K)=CM(K)+(SM1-SM2+SM3+SM4-SM5)*DT
            CN(K)=CN(K)+(SN1+SN2-SN3-SN4)*DT
         ENDIF
      ENDDO ! enddo loop K = 1, LK


      ! convert mass and number to kg/kg and #/kg
      !
      DO L=1,LK
         CM(L)=CM(L)/RHON(II)*1.E+03
         CN(L)=CN(L)/RHON(II)*1.E+06
      ENDDO
      LRK=15
!      LRK=LK
      DO L=1,LRK
         ANK(KR,II,L)=CN(L)
         AMK(KR,II,L)=CM(L)
      ENDDO

      ! save diagnostics of for accretion rates (pracw, cloud mass and 
      ! rracw, cloud number)
      !
      SUM=0.
      SUN=0.
      DO L=LRK+1,LK
         SUM=SUM+CM(L) - AMKORIG(KR,II,L)
         SUN=SUN+CN(L) - ANKORIG(KR,II,L)
      end DO
         
      if (jjp == 1) then 
         call save_dg(ii, sum/dt,'pracw', i_dgtime, &
              units='kg/kg/s',dim='z')
         call save_dg(ii,sun/dt,'rracw', i_dgtime, &
              units='#/kg/s',dim='z')
      else
         call save_dg(ii, kr, sum/dt,'pracw', i_dgtime, &
              units='kg/kg/s',dim='z,x')
         call save_dg(ii, kr,sun/dt,'rracw', i_dgtime, &
              units='#/kg/s',dim='z,x')
      endif
      rmass_cw(ii) = SUM
      rnum_cw(ii) = SUN
      
      END subroutine SXY

!****************************************************************
      SUBROUTINE SCONC(KR,II,DT)
! ***********************************************************************
! * EVOLUTION DUE TO COLLECTION USING A TWO MOMENT METHOD
! * COLLECTION KERNEL IS CONSTANT IN A CATEGORY
! * LOW AND LIST COLLECTION EFFICIENCIES ARE USED FOR THE LARGER DROPS
! ***********************************************************************
! * *
      IMPLICIT NONE

      INTEGER KR, II !J and K respect!
      INTEGER I !loop counter
      INTEGER L,K!loop counter for the cloud bins
      REAL :: DT !timestep from the 1d framework
      REAL :: AVG,masstemp,notemp, sum, sun

!      
      DO  L=16,LK
        SM0(L-15)=AMKORIG(KR,II,L)*RHON(II)/1.E+03
        SN0(L-15)=ANKORIG(KR,II,L)*RHON(II)/1.E+06
      ENDDO
      
      DO L=1,LK_big
         SM1=0.
         SM2=0.
         SM3=0.
         SM4=0.
         SM5=0.
         SN1=0.
         SN2=0.
         SN3=0.
         SN4=0.
         CA(L)=(2./(AX(L)**2))*(SM0(L)-AX(L)*SN0(L))
         CB(L)=(2.*SM0(L)-3.*AX(L)*SN0(L))/(AX(L)**3)
         CC(L)=2.*AX(L)*CA(L)
         CD(L)=(3.*SM0(L)-4.*AX(L)*SN0(L))/(AX(L)**2)
         
         IF(SM0(L) > eps .AND. SN0(L) > eps)THEN
            SMC(L)=SM0(L)/SN0(L)
            IF(SMC(L) >  AX(L+1))  THEN
               SMC(L) = (AX(L+1)) !upper bin boundary
            ENDIF
            IF(SMC(L) <  AX(L))  THEN
               SMC(L) = AX(L) !lower bin boundary
            ENDIF
         ELSE
            SMC(L)=0.
         ENDIF
!c
         IF (SMC(L) <  AX(L).OR.SMC(L) >  (2.*AX(L)).OR.SMC(L) == 0.0)   &
     &     THEN
            AP=1.0
         ELSE
            AP=0.5*(1.+3.*(AX(L)/SMC(L))-2*((AX(L)/SMC(L))**2.))
         ENDIF
         !
         !AH 0410 - original threshold was 1.e-20 for SN0 and no 
         !          threshold for mass (SM0). This can lead to overflows
         !          with ifort 11 and an AIX compiler. Thus, hardcoded
         !          threshold to be 1.e-15 and added in mass threshold
         !          This is not ideal - I will work on a more dynamic, 
         !          generic fix
         !
         IF(SN0(L) > eps .and. SM0(L) > eps)THEN
            AZ0(L)=AP*SM0(L)**2/SN0(L)
            AM3(L)=AP**3*SM0(L)**3/(SN0(L)**2)
         ELSE
            AZ0(L)=0.
            AM3(L)=0.
         ENDIF
         DO I=L,LK_big-1
            SN3=SN3+SN0(I)*ACON(I,L)
         ENDDO
         SM4=SN3*SM0(L)
         IF(L == 1)THEN
            SN(L)=SN0(L)-SN0(L)*SN3*DT
            SM(L)=SM0(L)-SM4*DT
            ! 
            CYCLE ! execute the next iteration of L=1,LK_big
            !
         ENDIF
         SN4=ACON(L-1,L-1)*0.5*SN0(L-1)**2
         SM5=ACON(L-1,L-1)*SN0(L-1)*SM0(L-1)
         DO I=1,L-1
            IF(SMC(L) >  AX(L))THEN
               SN2=SN2+(CA(L)*SM0(I)-CB(L)*AZ0(I))*ACON(L,I)
               SM2=SM2+((CA(L)-CD(L))*AZ0(I)+CC(L)*SM0(I)-CB(L)*AM3(I))      &
     &             *ACON(L,I)
            ENDIF
         ENDDO
         DO I=1,L-1
            SM3=SM3+SM0(I)*ACON(L,I)
         ENDDO
         SM3=SM3*SN0(L)
         SN3=SN3*SN0(L)
         IF(L == 2)THEN
            SN(L)=SN0(L)+(SN4-SN3-SN2)*DT
            SM(L)=SM0(L)+(SM3-SM4-SM2+SM5)*DT
            ! 
            CYCLE ! execute the next iteration of L=1,LK_big
            !
         ENDIF
         DO I=1,L-2
            IF(SMC(L-1) >  AX(L-1))THEN
               SN1=SN1+(CA(L-1)*SM0(I)-CB(L-1)*AZ0(I))*ACON(L-1,I)
               SM1=SM1+((CA(L-1)-CD(L-1))*AZ0(I)+CC(L-1)*SM0(I)-CB(L-1)    &
     &             *AM3(I))*ACON(L-1,I)
            ENDIF
         ENDDO
         !
         ! update the mass and number fields
         !
         IF (L == (LK-15)) THEN
            SM(LK-15)=SM0(LK-15)+(SM1+SM3+SM5)*DT
            SN(LK-15)=SN0(LK-15)+(SN1+SN4-SN3)*DT
         ELSE
            SM(L)=SM0(L)+(SM1-SM2+SM3-SM4+SM5)*DT
            SN(L)=SN0(L)+(SN1-SN2-SN3+SN4)*DT
         ENDIF
      ENDDO ! end do loop L=1,LK_big
      !
      ! store rain collecting rain mass rate (pracr) and 
      ! rain collection rain number rate (rracr)
      !
      !if (kr==1)then
         SUM=0.
         SUN=0.
         DO L=1,LK_big
            SUM=SUM+SM(L)*1.E+03/RHON(II) - AMKORIG(KR,II,L+15)
            SUN=SUN+SN(L)*1.E+06/RHON(II) - ANKORIG(KR,II,L+15)
         end DO
         
         if (jjp == 1) then
            call save_dg(ii, sum/dt,'pracr', i_dgtime, &
                 units='kg/kg/s',dim='z')
            call save_dg(ii, sun/dt,'rracr', i_dgtime, &
                 units='#/kg/s',dim='z')
         else
            call save_dg(ii,kr, sum/dt,'pracr', i_dgtime, &
                 units='kg/kg/s',dim='z,x')
            call save_dg(ii,kr, sun/dt,'rracr', i_dgtime, &
                 units='#/kg/s',dim='z,x')
         endif
      !end if
!**********************************************************************
!      IF BREAKUP IS NOT INVOKED THEN USE THE NEXT FOUR LINES
!      OTHERWISE PUT COMMENTS ON THEM
!*********************************************************************
      if(.not. l_break) then
         DO L=16,LK
            AMK(KR,II,L)=SM(L-15)*1.E+03/RHON(II)+(CM(L)                    &
     &                                  -AMKORIG(KR,II,L))
            ANK(KR,II,L)=SN(L-15)*1.E+06/RHON(II)+(CN(L)                    &
     &                                  -ANKORIG(KR,II,L))
         ENDDO
      endif

      END subroutine SCONC

!**********************************************************************
      SUBROUTINE BREAK(KR,II,DT)
!**********************************************************************
      IMPLICIT NONE

      INTEGER KR, II !J and K loop counters from 1d-mod respect!
      INTEGER K, J, I 
      INTEGER L
      REAL :: DT !timestep 
      REAL :: AVG, sum, sun
      SUM=0.
      SUN=0.
      DO K=1,LK_big
         SM1=0.
         SM2=0.
         SM3=0.
         SN1=0.
         SN2=0.
         SN3=0.

         IF (SMC(K) > eps ) THEN
            IF((SMC(K) >  2.*AX(K)).AND.(SM0(K) <                         &
     &         eps .or. SN0(K) < eps)) THEN
               SMC(K) = (2*AX(K)) !upper bin boundary
            ENDIF
            IF((SMC(K) <  AX(K)).AND.(SM0(K) <                            &
     &         eps .or.SN0(K) < eps)) THEN 
               SMC(K) = AX(K) !lower bin boundary
            ENDIF
         ELSE
            SMC(K)=0.
         ENDIF
!
         DO I=9,LK_big-1
            IF(K >  (I+1)) CYCLE ! execute next iteration of I=9,LK_big-1 
            IF(SMC(I) <  AX(I)) CYCLE ! execute next iteration of I=9,LK_big-1 
            DO J=6,I-1
               IF(SMC(J) <  AX(J)) CYCLE ! execute next iteration of J=6,I-1
               IF(SMC(I) >  AX(I+1))THEN
                  SM1=SM1+(SN0(I)*SM0(J)+SN0(I)*AX(I+1)*SN0(J))                 &
     &                 *PLL(I,J,K)*KIJ(I,J)
               ELSEIF(SMC(J) >  AX(J+1))THEN
                  SM1=SM1+(SN0(I)*SN0(J)*AX(J+1)+SM0(I)*SN0(J))                 &
     &                *PLL(I,J,K)*KIJ(I,J)
               ELSE
                  SM1=SM1+(SN0(I)*SM0(J)+SM0(I)*SN0(J))                         &
     &                *PLL(I,J,K)*KIJ(I,J)
               ENDIF
            ENDDO ! end of loop J=6,I-1
         ENDDO ! end of loop I=9,LK_big-1

         SN1=SM1*AX(K)
         SM1=SM1*AX(K)**2*1.5
         !
         ! AH 0410 - this loop looks a bit suspicious to me, check 
         !           against original code
         DO I=9,9
            IF(SMC(I) <  AX(I)) CYCLE 
            IF(K >  (I+1)) CYCLE
            IF(SMC(I) >  AX(I+1))THEN
               SM2=SM2+SN0(I)*AX(I+1)*SN0(I)*KIJ(I,I)*PLL(I,I,K)
            ELSE
               SM2=SM2+SM0(I)*SN0(I)*KIJ(I,I)*PLL(I,I,K)
            ENDIF
         ENDDO
         !
         SN2=SM2*AX(K)
         SM2=SM2*1.5*AX(K)**2
         IF ((K >= 6) .and. (K /= LK_big) .and. (SMC(K) >=  AX(K))) THEN
            !
            DO J=6,lk_big-1
               IF(SMC(J) <  AX(J)) CYCLE ! execute next iteration of J=6,lk_big-1
               IF(J >  K)THEN
                  IF(J <  9) CYCLE ! execute next iteration of J=6,lk_big-1 
                  SN3=SN3+SN0(J)*PLM(J,K)*KIJ(J,K)
               ELSEIF(K >= J)THEN
                  IF(K <  9) CYCLE ! execute next iteration of J=6,lk_big-1
                  SN3=SN3+SN0(J)*PLM(K,J)*KIJ(K,J)
               ENDIF
            ENDDO ! end of loop J=6,lk_big-1
            !
            IF(SMC(K) >  AX(K+1))THEN
               SM3=1.5*SN3*SN0(K)*AX(K+1)
            ELSE
               SM3=1.5*SN3*SM0(K)
            ENDIF
            SN3=1.5*SN3*SN0(K)
         ENDIF 
         
         if (kr==0)then
            DO L=1,LK_big
               SUM=SUM+((SM1+SM2-SM3)*DT)*1.E+03/RHON(II) - AMKORIG(KR,II,L+15)
            end DO
            
            if (jjp == 1) then
               call save_dg(ii,sum/DT,'prbrk', i_dgtime, &
                    units='kg/kg/s',dim='z')
            else
               call save_dg(ii,kr,sum/DT,'prbrk', i_dgtime, &
                    units='kg/kg/s',dim='z,x')               
            endif
            rmass_cw(ii) = rmass_cw(ii)+SUM
         end if
         
         SM(K)=SM(K)+(SM1+SM2-SM3)*DT
         SN(K)=SN(K)+(SN1+SN2-SN3)*DT
         !
      ENDDO ! end of loop  K=1,LK_big (main loop for entire routine)
      
      DO L=16,LK
         AMK(KR,II,L)=SM(L-15)*1.E+03/RHON(II) +(CM(L)                 &
     &                -AMKORIG(KR,II,L))
          ANK(KR,II,L)=SN(L-15)*1.E+06/RHON(II) +(CN(L)                 &
     &                -ANKORIG(KR,II,L))
      ENDDO

      END subroutine BREAK

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE BIN_SEDIMENT(I,DT,AMKORIG,ANKORIG,Q,SQ,RDT)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      
      IMPLICIT NONE

! This routine calculates the tendency for each bin due to sedimentation
! Unlike earlier versions, these tendencies (QL_ and QLN_SED) are saved 
! and passed to cloudbin using common block SEDIM so they can be added 
! to the other microphys tends at the end of cloudbin (in MICROcheck)
      
      INTEGER :: I !x grid 
      REAL :: DT !timestep 
      REAL :: RDT !reciprocal timestep
      REAL, DIMENSION(JMINP:JMAXP,KKP,NQP):: Q !moisture fields 
      REAL,INTENT(IN),DIMENSION(JMINP:JMAXP,KKP,LK)::AMKORIG,ANKORIG
      REAL, DIMENSION(JMINP:JMAXP,KKP,NQP):: SQ!moisture tendency
      REAL, DIMENSION(JMINP:JMAXP,KKP,LK) :: QLORIG
      REAL, DIMENSION(JMINP:JMAXP,KKP,LK) :: QLNORIG
      REAL :: DM
      REAL :: DN
      REAL, DIMENSION(KKP) :: DMTOT
      REAL, DIMENSION(KKP) :: DNTOT
      INTEGER,PARAMETER :: IORD=1
      INTEGER J,K,L,ITER
      REAL,DIMENSION(JMINP:JMAXP,KKP,LK)::binmass,binnum
      REAL, DIMENSION(JMINP:JMAXP,KKP) :: SED_FLUX, sed_field,fztemp
      REAL sqmassbef, sqnumbef,conserve_mass, conserve_no
      REAL, DIMENSION(KKP) :: mass_before,num_before,mass_after,&
                               num_after
     ! REAL, DIMENSION(KKP,LK) :: vt
      REAL :: DIAM
      
       character(2) :: str2   
 
     dmtot(:) = 0.0
     dntot(:) = 0.0
     vt(:,:) = 0.0 
     mass_before(:) = 0.0
     num_before(:) = 0.0
     mass_after(:) = 0.0
     num_after(:) = 0.0
     bin_pcpt(:) = 0.0   
      DO L=1,LK
        write(str2,'(i2.2)')l
        DO J=JMINP,JMAXP
          DO K=2,KKP
            QLORIG(J,K,L) = AMKORIG(J,K,L) 
            QLNORIG(J,K,L) = ANKORIG(J,K,L)
            IF(QLORIG(J,K,L) > eps.AND.             &
     &         QLNORIG(J,K,L) > eps)THEN
                  AMS(J,K)=QLORIG(J,K,L)/QLNORIG(J,K,L)
                  VT(J,K)=-ALP(L)*(AMS(J,K)*1000.)**BET(L)
            ELSE
              AMS(J,K)=0.
              VT(J,K)=0.
            ENDIF
!*************************************************************
!  DIVIDE BY 100 TO GET FALL VELOCITIES OF M/SEC
!*************************************************************
            VBAR(J,K)=VT(J,K)/100.
          ENDDO
        ENDDO
        
        DO K = 2, KKP-1
           DO J = JMINP, JMAXP
!1) WORK OUT THE MASS SEDIMENTATION FOR EACH BIN            
!calculate the value of the q field (q_sed) following sedimentation            
              QL_SED(J,K,L) = QLORIG(J,K,L) - DT*(RHON(K+1)* QLORIG(J,K+1,L)*    &
     &               VBAR(J,K+1) - RHON(K)* QLORIG(J,K,L)*                &
     &               VBAR(J,K))/(RHON(K)*DZN(K))
!work out the change in q due to sedimentation 
              DM = QL_SED(J,K,L) - QLORIG(J,K,L)
              QL_SED(J,K,L) = DM
              IF (K == 2) THEN
                PUDDLE(J,1,I) = PUDDLE(j,1,i) +                           &
     &           (ABS((RHON(K)* QLORIG(J,K,L)*VBAR(J,K))))*DT
                BIN_PCPT(J) = BIN_PCPT(J) +                           &
     &           (ABS((RHON(K)*QLORIG(J,K,L)*VBAR(J,K))))
              ENDIF

              IF (j .EQ. jjp) THEN
                dmtot(k) = dmtot(k)+dm
              endif
!
!2) WORK OUT THE NUMBER SEDIMENTATION FOR EACH BIN
!calculate the value of the q field (q_sed) following sedimentation            
              QLN_SED(J,K,L) = QLNORIG(J,K,L) - DT*(RHON(K+1)* QLNORIG(J,K+1,L)* &
     &               VBAR(J,K+1) - RHON(K)* QLNORIG(J,K,L)*               &
     &               VBAR(J,K))/(RHON(K)*DZN(K))
!work out the change in q due to sedimentation 
              DN = QLN_SED(J,K,L) - QLNORIG(J,K,L)
              QLN_SED(J,K,L) = DN
!add the DM divided by timestep (DT) to the source term for Q
!              SQ(J,K,ICDNC_BIN(L)) = SQ(J,K,ICDNC_BIN(L)) + DN*RDT 
              IF (j .EQ. jjp) THEN
                 dntot(k) = dntot(k)+dn
              endif           
            ENDDO
          ENDDO
          DO J = JMINP, JMAXP
             QL_SED(J,KKP,L) = QLORIG(J,KKP,L) - DT*(-RHON(KKP)* QLORIG(J,KKP,L)*                &
     &               VBAR(J,KKP))/(RHON(KKP)*DZN(KKP))
             DM = QL_SED(J,KKP,L) - QLORIG(J,KKP,L)
             QL_SED(J,KKP,L) = DM
             
             QLN_SED(J,KKP,L) = QLNORIG(J,KKP,L) - DT*(-RHON(KKP)* QLNORIG(J,KKP,L)*                &
     &               VBAR(J,KKP))/(RHON(KKP)*DZN(KKP))
             DN = QLN_SED(J,KKP,L) - QLNORIG(J,KKP,L)
             QLN_SED(J,KKP,L) = DN
          ENDDO
             
!deal with K=1
          DO J = JMINP, JMAXP
            SQ(J,1,ICDKG_BIN(L)) = 0.0
            SQ(J,1,ICDNC_BIN(L)) = 0.0
         ENDDO
      ENDDO
      
      name='surface_ppt_for_warm_bin'
      units='kg/kg ms-1'
      value=sum(bin_pcpt(1:JJP))/JJP

      call save_dg(value, name, i_dgtime, units,dim='time')
      
      if (jjp > 1) then
         call save_dg(bin_pcpt, name, i_dgtime, units,dim='time')
      endif

      END subroutine BIN_SEDIMENT
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
      SUBROUTINE REFFCALC (Q,SQ,DT,RDT) 
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      IMPLICIT NONE
      
!This routine uses the multi-moment method to calc rain rate 
!and cloud
!      
!arguments
      INTEGER :: I
      REAL :: DT !timestep 
      REAL :: RDT !reciprocal timestep
      REAL, DIMENSION(JMINP:JMAXP,KKP,NQP)::Q,                              &
     &                                 SQ!moisture tendency      
!      
!Local variables
      REAL :: PK
      REAL :: XI_AVE
      REAL :: AELO
      REAL :: RRM
      REAL :: SUR0M
      REAL :: SRR0H
      REAL :: SWR
      REAL :: DMJ
      REAL :: DMASJ  
!Local loop counters
      INTEGER J, K, L
      
      RRM = DYY*JJP
      DO L=1,LK
        DO K=2,KKP
          DO J=1,JMAXP
            AMK(J,K,L) = Q(J,K,ICDKG_BIN(L)) + (SQ(J,K,ICDKG_BIN(L))*DT)
            ANK(J,K,L) = Q(J,K,ICDNC_BIN(L)) + (SQ(J,K,ICDNC_BIN(L))*DT)
            IF (AMK(J,K,L) <  0.0.or.ANK(j,k,L) <  0.0) THEN
              AMK(J,K,L) = 0.0
              ANK(J,K,L) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!***********************************************************
!     UPDATING OVERALL MASS AND CONC
!***********************************************************
      DO J=1,JJP
         DO K=2,KKP
            AM1(J,K)=0.0
            AN1(J,K)=0.0
            ztotal(J,K)=0.0
            reff(J,K)=0.0
            DO L=1,LK
               AM1(J,K)=AM1(J,K)+AMK(J,K,L)
               AN1(J,K)=AN1(J,K)+ANK(J,K,L)
               reff(J,K)=reff(J,K)+amk(J,K,L)**(2./3.)*ank(J,K,L)        &
     &                   **(1./3.)
            enddo
            if(reff(J,K) > eps .and.AM1(J,K) >= 1.e-6) then
               reff(J,K)=(3./(4.*pi*1000.))**(1./3.)*AM1(J,K)/XI23         &
     &                   /reff(J,K)
            else
               reff(J,K)=0.0
            endif
         enddo
      enddo
      if (jjp == 1) then
         field(:) = reff(jjp,:)*1.e6
         call save_dg(field,'r_eff', i_dgtime, units='microns',dim='z')
      else
         DO J=1,JJP
            DO K=2,KKP
               field_2d(k, j) = reff(j,k)*1.e6
            enddo
         enddo
         call save_dg(field_2d,'r_eff', i_dgtime, units='microns',dim='z')
      endif

      END subroutine REFFCALC
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE REBIN(J,K)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!Purpose: test routine - this checks whether the particles lie 
!         outside their corresponding bin boundaries. If they 
!         do the mass and number are transfered to the correct
!         bin, and the number (ANK) and mass (AMK) are updated 
!         for the next microphysical steps
      IMPLICIT NONE
      
!CALL PRAMETR
!      
!CALL SD
!      
!CALL RI      
      REAL, DIMENSION(JMINP:JMAXP, KKP, LK) ::                              &
     & AVG,                                                             &
             !mean droplet size
     & DMK,                                                             &
             !change in mass due to transfer 
     & DNK   !change in number due to transfer
      REAL DUM
      INTEGER :: J,K,L !loop counters
      INTEGER :: NEWBIN!the bin that mass and no.are tranferred
                       !to
      
!1) zero DMK and DNK
      DO L = 1, LK
        DMK(J,K,L) = 0.0
        DNK(J,K,L) = 0.0
      ENDDO
!
!2) calculate  the mean droplet size
      DO L = 1, LK
        IF (ANK(J,K,L) > 0.0) THEN
          AVG(J,K,L) = (AMK(J,K,L)/ANK(J,K,L))
        ELSE
          AVG(J,K,L) = 0.0
        ENDIF
      ENDDO
!      
!3a) check to see whether the mean droplet size lies
!   within its corresponding bin boundaries
      DO L = 1, LK
        IF (AVG(J,K,L) >  0.0) THEN
          IF (AVG(J,K,L) <  XK_gr(L).OR.AVG(J,K,L) >  (2*XK_gr(L))) THEN
            IF (AVG(J,K,L) <= 0.0) THEN
              NEWBIN = 1
            ELSE
              DUM = 1.0 + (LOG(AVG(J,K,L)/XK_gr(1))/LOG(2.0))
              NEWBIN = INT(DUM)
              IF (NEWBIN >  LK)  NEWBIN = LK
              IF (NEWBIN <= 0.0) NEWBIN = 1
            ENDIF
!3b) calculate the changes in bin mass and number            
            DNK(J,K,NEWBIN)=DNK(J,K,NEWBIN) + ANK(J,K,L)
            DMK(J,K,NEWBIN)=DMK(J,K,NEWBIN) + AMK(J,K,L)
            DNK(J,K,L)=DNK(J,K,L) - ANK(J,K,L)
            DMK(J,K,L)=DMK(J,K,L) - AMK(J,K,L)
          ENDIF
        ENDIF
      ENDDO
!      
!4) calculate the new numbers and mass. This the change (DN and DM)
!   added to ANK and AMK. The sourceterms are calced at the end
!   of microphysics, they are not updated here.
      DO L=1, LK
        ANK(J,K,L)=ANK(J,K,L)+DNK(J,K,L)
        AMK(J,K,L)=AMK(J,K,L)+DMK(J,K,L)
      ENDDO
      
      END subroutine REBIN
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE  MICROcheck(J,K,IQ,DT,Qmass,Qnum,Sourcemass, &
           Sourcenum,DAMKDT,DANKDT,RDT,Qmass_orig,Qnum_orig)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
      IMPLICIT NONE
      
!local variable
      REAL :: AVG, avginit, avgorig
      REAL :: QmassFLD, Qmass,Qmass_orig
      REAL :: Qnumfield, Qnum,Qnum_orig
      REAL :: Sourcemass!SQ for mass passed in and out of routine
      REAL :: Sourcenum!SQ for number passed in and out of routine
      REAL :: dmdt!change in mass with time
      REAL :: dndt!change in number with time
      REAL :: DAMKDT
      REAL :: DANKDT
      REAL :: DT
      REAL :: RDT !1/timestep from tau
!loop counters
      INTEGER J, K, IQ, ISTOPU
!     
      QmassFLD = Qmass_orig+((Sourcemass+DAMKDT)*DT)
      Qnumfield = Qnum_orig+((Sourcenum+DANKDT)*DT)
!        
!positivity check for negative values of mass. 
!If mass or number goes less than zero after micro
!assume all mass is lost for that bin, set SQ to just above q value

      IF (QmassFLD.LT.0.0.or.Qnumfield.lt.0.0) THEN
        DAMKDT = -0.99999*Qmass*RDT
        DANKDT = -0.99999*Qnum*RDT
      ENDIF
!
      Sourcemass = Sourcemass + DAMKDT
      Sourcenum = Sourcenum + DANKDT

      IF (ABS(Sourcemass).LT.eps .OR.                 &
           ABS(Sourcenum).LT.eps) THEN
         Sourcemass = 0.0
         Sourcenum = 0.0
       ENDIF
       
      END subroutine MICROCHECK

!----------------------------------------------------------------------
! XACT: Computes the fraction of aerosol in a log normal dist.
! activated at a given temperature (K) and supersaturation (%)
! 
      REAL FUNCTION XACT(TEMP,S,DG,SG,DCRIT)

      IMPLICIT NONE

      REAL :: TEMP,S,DG,SG,A,DCRIT,GAMMA
      REAL, PARAMETER :: BETA = 0.651232
!      REAL, EXTERNAL :: TOPLOG
!
! Calculations of critical dry ccn radius associated with equilibrium 
! nucleation.
!
      GAMMA = LOG(S/100.+1.)
      A = (3.298E-05 - 6.7172E-08*(TEMP-273.15))/TEMP
      DCRIT =2.*((4.*A*A*A)/(27.*GAMMA*GAMMA*BETA))**(0.3333333333)
!
! Calculation of number larger than critical radius
!
      XACT=TOPLOG(DCRIT,DG,SG,0.)

      END function XACT
!
!----------------------------------------------------------------------
! TOPLOG: Finds the moments associated with the upper portion of the
! lognormal function specified by the parameters dg and sg
!
      FUNCTION TOPLOG(DMX,DG,SG,ORD)

      IMPLICIT NONE

      REAL :: TOPLOG,DMX,DG,SG,ORD
      REAL :: LNS,YMX,EMX
      REAL, PARAMETER :: C0=0.7071067814

      LNS = LOG((SG))
      YMX = (LOG((DMX)/(DG))/LNS - (ORD)*LNS)*C0
      EMX = ERF(YMX)

      TOPLOG = ((DG)**(ORD)                                         &
     &             *EXP(.50*(ORD*ORD)*LNS*LNS)                          &
     &             *(1.0-EMX)*.50)

      END function TOPLOG

!----------------------------------------------------------------------
! FUNCTION ERF:  This subroutine uses the polynomial fit to the error
! function given in Abramowitz and Stegun 7.1.28
!
      FUNCTION ERF(X)

      IMPLICIT NONE
      REAL :: X,ERF
      REAL :: Z,Y
      REAL, PARAMETER :: A1=0.0705230784
      REAL, PARAMETER :: A2=0.0422820123
      REAL, PARAMETER :: A3=0.0092705272
      REAL, PARAMETER :: A4=0.0001520143
      REAL, PARAMETER :: A5=0.0002765672
      REAL, PARAMETER :: A6=0.0000430638

      Y=ABS(X)
      Z=(1.D0+Y*(A1+Y*(A2+Y*(A3+Y*(A4+Y*(A5+Y*A6))))))
      Z=Z*Z
      Z=Z*Z
      Z=Z*Z
      Z=Z*Z

      ERF=(Y - Y/Z)/(X+eps)
      ERF=MIN( 1.0,ERF)
      ERF=MAX(-1.0,ERF)

      END function ERF
   
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE COND_NEW(J,K,DM,TBASE,QST,RH,Q,SQ,DT,RDT               &
     &                    ,PMB,QVNEW,TAU,it,lt)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
!This is a newer version of the COND routine from the orig.
!cloud model supplied by Yin Yan. Although the orig routine 
!works OK, it failed to run when code was optimised to -O2
!and compiled on pgf90 6.1. Therefore I have updated the code
!using newer cond and evap routines supplied by Graham 
!Feingold - (Adrian Hill, October 2006)
!
      IMPLICIT NONE

      INTEGER J,K,L,Jl !Jl is an inner loop of l
      INTEGER JK,LT,IT !components of the big loop, not sure if loop is req
                    !but is doing no harm except taling time
      INTEGER I,KK,JJ !loop counters
      INTEGER IK                                                        &
                !lower category boundary, eq to KK within the loop
     &,IJ                                                               &
          !upper category boundary, eq to JJ within the loop 
     &,NJ                                                               &
         !number of categories between the upper and lower bounds IJ-IK 
     &,AH 
      REAL,DIMENSION(JMINP:JMAXP,KKP,NQP):: Q,SQ 
  
      REAL,DIMENSION(JMINP:JMAXP,KKP) ::                                    &
     & QS                                                               &
          !saturation mixing ratio of water calced in Qsaturation
     & , QST                                                            &
             !specific humidity calculated in Spechum
     & , RH !relative humidity
      
      INTEGER,PARAMETER :: RRV=461 !gas constant of moist air J/kg/K
      
      INTEGER,PARAMETER :: AL=597  !latent heat of evaporation and cond
      
      REAL,PARAMETER :: T00=273.15 !freezing point of water in K,from
                                   !the LES
      
      REAL,PARAMETER :: CPBIN=0.24 !specific heat of water vapour
      
      REAL,PARAMETER :: AAR=0.24E-3
      
      REAL,DIMENSION (JMINP:JMAXP,KKP):: TBASE !temp calced in Lathe
                   !using tref+th(1/exner function).
      REAL DT !time calced in the LES
      REAL RDT !1/timestep
!      REAL DTmic !microphysical timestep if DT > 4
      REAL DM  ! mass of condensated water
      REAL DMM !mass of condensed water at the end of the subroutine
      
      REAL pres,R0 ! pressure and density
      REAL F,FF,FFF,HH,HHH,GG,GGG,SM1,SN1,SUMDS !not sure, but these are calced
      !and used in the big calculation of AMK and ANK
!      REAL, EXTERNAL :: AA,CPT,FG
      REAL :: QVNEW                                                     &
                    !vapour mixing ratio after subnu and cond
     &       ,QVOLD                                                     &
                    !vapour mixing ratio before bin microphysics
     &       ,TBOLDcon
      REAL :: PMB !pressure in mb
      REAL :: DMDT,SQbef
      REAL :: DTcalc
      REAL,DIMENSION(LK) :: DPHASEM
!New variables for from Grahams code
      REAL :: sdm !total change in condensed mass
      REAL :: sddm !itegrated mass over all bin (I think!)
                  !sddm=sddm+an0(l)*amn(l)**(1./3.)
      REAL :: delm !theoretical mass change in mass
      REAL, DIMENSION(LK) :: s1, s2, s3, s4
                  !variables for calculation of change in mass and
                  !number due to condensation 
      REAL :: fk, psik, ex2, fx1, fy, fz, fxn, fac
      REAL,DIMENSION(JMINP:JMAXP,KKP)::TAU

      Real :: t_test
      
      INTEGER :: i0                                                     &
                    !the bin in which yk(l) is located
     &          ,in                                                     &
                    !the bin where zk(l) is located NB.zk > yk, in>i0
     &          ,idiff                                                  &
                       ! i0 - in
     &          ,i1                                                     &
                    ! i0 + 1
     &  ,j1
      
! 
!      REAL,EXTERNAL::QSATURATION, FGV, F1, F2
       
       
       r0= RHON(K)*1.e-3
       pres = PREFN(K)*10
       ex2 = 0.6666666667    
! 
! for high supersaturations or large DT reduce the microphysical time step. 
! the microphysical DTcalc, can be
! changed allowing multiple condensation iterations, this is done for
! supersaturations greater that 4% as specified below, or timestep >  2s.
      do l=1,LK
        dphasem(l)=0.0         
!         sm0(l)=(smk0(l))
!         xk(l)=(x_bin(l))
      enddo
!         xk(lk+1)=dble(x_bin(lk+1))
! 
! begin the timestep iterations (usually only one) for condensation 
! calculations.  first we calculate the function f which is related to
! the integral of s(t)dt in the reference (appendix b).  necessary for
! the calculation of f is the functions a(p,t) and s2(p,t) where the
! former is a thermodynamic function and the latter is the growth rate
! of a single droplet.  both are given in appendix b of the reference.
! if f is below a certain threshold, there is no condensation we exit 
! from the routine.
!
      dds=ds(j,k)

      sddm=0.0

      do l=1,lk
         sddm=sddm+an0(l)*amn(l)**(1./3.)
      enddo
!
!tau is calced in get_forcing
      f = tau(j,k)
!
      if(f <  1.0e-20)then
         do l=1,lk
            amk(j,k,l)=am0(l)
            ank(j,k,l)=an0(l)
            !               smk(j,k,l)=sm0(l)
         enddo
         return
      endif
!
! delm is the theoretical value of mass change, as per equation (33)
!
      delm=f*1.50*xi*sddm
!
! preperation of some standard functions of the bin number for 
! use in the condensation calculations.  functions on determined
! depending on whether or not the average bin mass is within its
! bin range

 
      do i=1,lk
         if(amn(i) <  x_bin(i))then
            ff=1.0
         elseif(amn(i) >  x_bin(i+1))then
            ff=bin_res
         else
            ff=amn(i)/x_bin(i)
         endif
         fk=2.0*an0(i)/(x_bin(i)*(bin_res-1.0))*                               &
              (1.0 - 1./(bin_res-1.0) * (ff-1.0) )
         psik=2.*an0(i)/(x_bin(i)* (bin_res-1.0)*(bin_res-1.0))*               &
              (ff-1.0) 

! These for solute
         fac=1.0/(x_bin(i+1)-x_bin(i))
         s1(i)=fac *                                                       &
              (x_bin(i+1)* fk - x_bin(i) * psik )

! For solute

         s2(i)=0.50*fac * (psik - fk)
         s3(i)=fac* (bin_res * x_bin(i)*x_bin(i) * (fk-psik) )

         s4(i)=fac*                                                        &
              ( x_bin(i) * (bin_res*psik -fk ) )

      enddo

!
! begin iterating through the bins for the actual condensation
! calculations
!

      do l=1,lk

         zk(l)=(max((x_bin(l+1)**(ex2) - f),0.0))
         zk(l)=zk(l)*sqrt(zk(l))
         if(zk(l) >  x_bin(1))then
            yk(l)=max(x_bin(1),(max((x_bin(l)**(ex2) - f),0.0)**(1.50)))
            if(abs(yk(l)-x_bin(l)) <  1.0e-20)yk(l)=x_bin(l)
            if(abs(zk(l)-x_bin(l+1)) <  1.0e-20)zk(l)=x_bin(l+1)
            !
! locate the bin in which yk(l) is located, designate as i0
! 
            do i=1,lk
               if(yk(l) <  x_bin(i+1))then
                  i0=i 
                  exit
               endif
            enddo
! 
! locate the bin in which zk(l) is located; for efficiency, start
! at i0 (since zk>yk), futhermore designate as in
! 
!               do i=i0,lk+1
            do i=i0,lk
               if(zk(l) <  x_bin(i+1))then
                  in=i
                  exit
               endif
            enddo

            idiff=in-i0
            i1=i0+1
! 
! below the parameters fx1,fy, fz,fxn,gg, and ggg are parameters
! appearing in the analytic solution to integrals in the reference
! equations 12 and 13.  furthermore functions f1 and func2 are defined
! because the form of this computation is repeatedly called and they
! simplify the code
! 
            fx1 =sqrt((x_bin(i1)**ex2 + f))
            fy=sqrt(( yk(l)**ex2 + f))
            fz =sqrt(( zk(l)**ex2 + f))
! 
! solutions to eqs. 8 and 9 (#3) for various values of nj.
!    idiff=0 implies yk and zk are located in the same bin.
!    idiff=1 implies zk is in a bin one higher than yk.
!    idiff>1 implies zk is in a bin more than one higher than yk. 
! the latter case should never appear given mass doubling between bins
!
            if(idiff == 0)then
             ank(j,k,l)=(zk(l)-yk(l))*(s1(i0)+s2(i0)*(zk(l)+yk(l))) 
             amk(j,k,l)= s4(in)*(f1(zk(l),fz,f)-f1(yk(l),fy,f))         &
                  + s3(in)*(func2(f,fz)-func2(f,fy))

!             smk(l)=(zk(l)-yk(l))*(s1s(i0)+s2s(i0)*(zk(l)+yk(l))) 

! Note: The following doesn't work if the bin originally has no solute in it (sm0(l)=0)!
!            dn=(ank(j,k,l)-an0(l))
!            smk(l)=sm0(l)*(1.+dn/(an0(l)+1.d-30) )
               
          elseif(idiff == 1)then
             fxn=sqrt(x_bin(in)**(ex2)+f)
             ank(j,k,l)=(x_bin(i1)-yk(l))*(s1(i0)+s2(i0)*(x_bin(i1)+yk(l)))      &
                  +(zk(l) -x_bin(in))*(s1(in)+s2(in)*(zk(l)+x_bin(in)))
!            smk(l)=(x_bin(i1)-yk(l))*(s1s(i0)+s2s(i0)*(x_bin(i1)+yk(l)))
!     %           +(zk(l) -x_bin(in))*(s1s(in)+s2s(in)*(zk(l)+x_bin(in)))


!            dn=(ank(j,k,l)-an0(l))
!            smk(l)=sm0(l)*(1.+dn/(an0(l)+1.d-30) )

             amk(j,k,l)=s4(in)*f1(zk(l),fz,f)-s4(i0)*f1(yk(l),fy,f)   &
                  +s3(in)*func2(f,fz)      -s3(i0)*func2(f,fy)          &
                  + (s4(i0)-s4(in))*f1(x_bin(in),fx1,f)               &
                  + (s3(i0)-s3(in))*func2(f,fx1) 

          elseif(idiff >= 2)then
             fxn=sqrt(x_bin(in)**(ex2)+f)
             sn1=0.
             sm1=0.
!                  soln=0.
             do jl=i1,in-1
                j1=jl+1
                sn1=sn1+an0(jl)

!                     soln=soln+sm0(jl)

                gg=sqrt(x_bin(j1)**(ex2)+f)
                ggg=sqrt(x_bin(jl)**(ex2)+f)
                sm1=s4(jl)*(f1(x_bin(j1),gg,f) - f1(x_bin(jl),ggg,f))      &
                     +s3(jl)*(func2(f,gg)-func2(f,ggg)) + sm1
             enddo
             ank(j,k,l)=(x_bin(i1)-yk(l))*(s1(i0)+s2(i0)*(x_bin(i1)+yk(l)))       &
                  +(zk(l) -x_bin(in))*(s1(in)+s2(in)*(zk(l)+x_bin(in)))           &
                  +sn1
!           smk(l)=(x_bin(i1)-yk(l))*(s1s(i0)+s2s(i0)*(x_bin(i1)+yk(l)))
!     %          +(zk(l) -x_bin(in))*(s1s(in)+s2s(in)*(zk(l)+x_bin(in))) 
!     %            +soln


!            dn=(ank(j,k,l)-an0(l))
!            smk(l)=sm0(l)*(1.+dn/(an0(l)+1.d-30) )
             amk(j,k,l)= s4(in)*(f1(zk(l),fz,f) - f1(x_bin(in),fxn,f)) &
                  + s4(i0)*(f1(x_bin(i1),fx1,f)- f1(yk(l),fy,f))      &
                  + s3(in)*(func2(f,fz)-func2(f,fxn))                   &
                  + s3(i0)*(func2(f,fx1)-func2(f,fy))                   &
                  + sm1

          endif
       else
          ank(j,k,l)=0.0
          amk(j,k,l)=0.0
       endif
    enddo
  
    end subroutine COND_NEW

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      FUNCTION AA(P,J,K,temp,rst)  
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      IMPLICIT NONE
      
!CALL PRAMETR
      REAL :: AA
!Intent (IN)       
      INTEGER,INTENT(IN):: J,K
      REAL,INTENT(IN):: TEMP 
      REAL,INTENT(IN) :: RST
      REAL,INTENT(IN) :: P

      REAL,PARAMETER :: AL = 597.0

      REAL,PARAMETER :: CPBIN = 0.24
! * *
! * FUNCTION REQD. FOR CALCULATING THE CHANGE IN SATURATION
! * MIXING RATIO, ACCORDING TO SOONG (1974) (A6)
!     * *
      AA=1.+4098.17*AL*RST/(CPBIN*(TEMP-35.86)**2.0)

      END function AA

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      FUNCTION CPT(P,J,K,temp,rst)
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      IMPLICIT NONE

!CALL PRAMETR
!CALL PARAMS
      INTEGER,INTENT(IN):: J,K
      REAL,INTENT(IN) :: TEMP
      REAL,INTENT(IN) :: RST
      
      REAL,INTENT(IN) :: P
      
      REAL :: CPT      

      REAL,PARAMETER :: EPS_cpt=0.622
      REAL :: E_cpt
      REAL :: DEN

      E_cpt=P/EPS_cpt*RST
      DEN=1./(5.58E+10/(TEMP**2)+2.04E+04*TEMP*1000./E_cpt)
      CPT=4.*PI*(3./(4.*PI))**(1./3.)*DEN/RST

      END function CPT
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      function f1(a,b,c)
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      IMPLICIT NONE    

      real :: a,b,c,d,e,ex1,f1
      
      ex1 = 0.3333333333

      d=a**ex1
      e=b*d

      f1=0.5*b*b*b*b*e - 0.125*c*(a*b + c*(2.5*e+1.5*c*log(b+d)))

      end function f1
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      function func2(a,b)
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! just a functional form commonly encountered in the condensation 
! routine, resulting from the analytical solution of the integrals
! for condensation.  this function assumes a is positive as it should
! be for condensation.
      IMPLICIT NONE

      REAL :: a,b,c,ext1,func2
      
      ext1 = 0.3333333333

      c=sqrt(a)
      func2=1.5*(2*b*(a+(b*b)*ext1)+a*c*log((b-c)/(b+c)))

      end function func2

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE EVAP_NEW(igrid,J,K,DM,TBASE,QST,RH,Q,SQ,DT,RDT         &
     &                    ,PMB,QVNEW,TAU,EA,it,LT,inhom_evap,delm)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
!This is a newer version of the EVAP routine from the orig.
!cloud model supplied by Yin Yan. Although the orig routine 
!works OK, it failed to run when code was optimised to -O2
!and compiled on pgf90 6.1. Therefore I have updated the code
!using newer cond and evap routines supplied by Graham 
!Feingold - (Adrian Hill, October 2006)

      IMPLICIT NONE

      INTEGER J,K,L,AH, IQ, JL,igrid
      INTEGER JK,LT !components of the big loop, not sure if loop is req
                    !but is doing no harm except taling time
      INTEGER I,KK,JJ !loop counters
      INTEGER IK                                                        &
                !lower category boundary, eq to KK within the loop
     &,IJ                                                               &
          !upper category boundary, eq to JJ within the loop 
     &,NJ                                                               &
         !number of categories between the upper and lower bounds IJ-IK
     &  ,KP
      INTEGER :: inhom_evap
      REAL,DIMENSION(JMINP:JMAXP,KKP,NQP):: Q,SQ 
      REAL,DIMENSION(JMINP:JMAXP,KKP) ::                                    &
     & QS                                                               &
          !saturation mixing ratio of water calced in Qsaturation
     & , QST                                                            &
              !specific humidity calculated in Spechum
     &  , RH                    ! relative humidity
      REAL, DIMENSION(JMINP:JMAXP,KKP):: TAU,EA,tau_sub,tau_res
      INTEGER,PARAMETER :: RRV=461 !gas constant of moist air J/kg/K
      INTEGER,PARAMETER :: AL=597  !latent heat of evaporation and cond
      REAL,PARAMETER :: T00=273.15 !freezing point of water in K,from
                                   !the LES
      REAL,PARAMETER :: CPBIN=0.24
      REAL,PARAMETER :: AAR=0.24E-3
      REAL, DIMENSION(JMINP:JMAXP,KKP):: TBASE !temp calced in Tau
                 !using tref+th(1/exner function).
      REAL DT !timestep calced in the LES
      REAL RDT !1/timestep
!      REAL DTmic !microphysical timestep if DT > than 4
      REAL, DIMENSION(LK):: ZVEN

      REAL :: DM,dm_sub,dm_res    ! mass of condensated water
      REAL DMM !mass of condensed water at the end of the subroutine
      
      REAL pres,R0 ! pressure and density
      REAL XMAX
      REAL SUMX,SUMY,DDSD,AMMD
      REAL F,FF,FFF,HH,HHH,GG,GGG,SM1,SN1,SUMDS !not sure, but these are calced
      !and used in the big calculation of AMK and ANK
!      REAL, EXTERNAL :: AA,CPT,FGV
      REAL :: QVNEW                                                     &
                    !vapour mixing ratio after subnuc and evap
     &       ,QVOLD                                                     &
                    !vapour mixing ratio before bin microphysics
     &       ,SQbef
      REAL :: DMDT!tendency in cond mass due to evap
      REAL :: PMB,TBOLD
      REAL,DIMENSION(JMINP:JMAXP,KKP,LK)::AMKBEF !mass conc before evap
      REAL :: DTcalc, ex1, ex2, sldeni
!
!New declarations for this version of evap
!
      REAL, DIMENSION(LK) :: dphasem,dphasen,dphasem_sub                &
     &  ,dphasem_res
      
      REAL :: rofac
      REAL :: sdn !total change in CDNC
      REAL :: sdd, sddm, chi
       ! components for weighted average for scaling c(p,t)
      REAL :: zvenfac1, zvenfac2
       ! ventilation coefficients for evap
      REAL :: delm ! theoretical calculation of mass change
      REAL :: fk, psik, dn
      REAL :: fx0, fx1,fy, fz, fy1, fxn, fac
       !parameters appearing in the analytic solution to integrals
       !see Tzivion et al (1989)
      REAL :: xnevap !the total number of evaporated drops
      REAL, DIMENSION(LK) :: s1, s2, s3, s4
      REAL :: S
      
      REAL,PARAMETER :: UM = 1.72E-04 !dynamic viscosity of air
      real :: amk_tot,ank_tot
      INTEGER :: it                                                     &
                    !number of iterations for evap
     &          ,il1                                                    &
                    ! the first bin
     &          ,i0                                                     &
                    ! bin where yk is located
     &          ,in                                                     &
                    ! bin where zk is located
     &          ,idiff                                                  &
                       ! in - i0
     &          ,i1                                                     &
                    ! i0 + 1
     &  ,j1
      ! local declarations for inhomogeneous mixing
      REAL,DIMENSION(JMINP:JMAXP,KKP)::am1new,am1old,an1old
      ! 5/1/08 - local declarations for test mass and number
      !          calculated using tau_res and tau_sub
      REAL, DIMENSION(LK) :: AM0_ITER,AN0_ITER
      ! declarations for JR06 mixing in EVAP_MIX
      REAL,DIMENSION(JMINP:JMAXP,KKP)::                                     &
     &  L_evap                ! evaporation limiter as calced using JR06
      REAL, DIMENSION(JMINP:JMAXP,KKP,LK)::R_k,m_k !mean rad for each bin
      
      real :: t_test
      
!      REAL,EXTERNAL :: QSATURATION,f1, f3

!
! initialization and scaling of variables. dt is the evaporation time
! step, if t=dt then the evaporation time step is the same as the 
! dynamic timestep.  also we calculate the dynamic viscosity of air 
! (dv), and the scaling factor for density (r0fac) where sldeni is the
! inverse of the sea level density in cgs units.  the schmidt number
! for use in the ventilation coefficient calculations is also found.
!
      ex1 = 0.3333333333
      ex2 = 0.6666666667 
      sldeni = 784.3137255       

      r0= RHON(K)*1.e-3         !units for EVAP and COND are g/cm^3  
      pres = PREFN(K)*10 

      amk_tot = 0.0
      ank_tot = 0.0
      am1old(j,k) = 0.0
      an1old(j,k) = 0.0
      am1new(j,k) = 0.0

!     new routine call to calculate the evap and mix timescales.
!     can be used to calc evap limiters
!!$      IF (IEVAP_MIX_JR06 == 1) THEN
!!$        S = (QVNEW/QST(J,K)) - 1
!!$        R_bar(J,K) = 0.0
!!$        R_Int(J,K) = 0.0
!!$        DO L = 1, LK
!!$          IF(AMKOLD(J,K,L) > epsilon(AMKOLD(J,K,L)).AND.    &
!!$     &            ANKOLD(J,K,L) >  epsilon(ANKOLD(J,K,L)))THEN
!!$            m_k(J,K,L)=AMKOLD(J,K,L)/ANKOLD(J,K,L)
!!$          ELSE
!!$            m_k(J,K,L)=0.
!!$          ENDIF
!!$          an1old(j,k) = an1old(j,k)+ankold(j,k,l)
!!$          am1old(j,k) = am1old(j,k)+amkold(j,k,l)
!!$          R_k(J,K,L) = ((3.0*m_k(J,K,L))/(4.0*1000.0*PI))**(1.0/3.0)
!!$          R_Int(J,K) = R_Int(J,K)+ (R_k(J,K,L)*ankold(j,k,l))
!!$        ENDDO
!!$        IF (an1old(j,k) > epsilon(an1old(j,k))) THEN
!!$          R_bar(J,K) = R_Int(J,K)/an1old(j,k)
!!$!          R_bar(J,K) = ((3.0*am1old(j,k))/(4.0*PI*1000.0*an1old(j,k)))
!!$!     $               **(1./3.)
!!$        ENDIF
!!$       
!!$!        CALL EVAP_MIX(igrid,J,K,PREFN,TBASE,QVNEW,S,L_evap,QST,EA)
!!$      ENDIF
      
!     work out the mass before micro, which is used to work fractional
!     change in mass for inhom_mix parameterisation      
      IF(inhom_evap  ==  1) then
        am1old(J,K) = 0.0
        DO L=1,LK
          am0_iter(l) = am0(l)
          an0_iter(l) = an0(l)
          am1old(j,k) = am1old(j,k) + AM0(L)
        ENDDO
      endif
!
! zeroing total mass and number change scalars and initializing scratch
! arrays according to spectra given by dynamics.
!

      do l=1,lk
         dphasem(l)=0.0
      enddo
!
! beginning the timestep iterations, at each timestep we calculate the
! necessary thermodyanimc functions (see reference and comments in 
! condensation.)
!
     sdn=0.0
     dds=ds(j,k)
     sdd=0.0
     sddm=0.0
!
! calculations of ventilation coefficients using reynolds numbers
! scaled to appropriate density.
!

     do l=1,lk
! calculation of components of weighted average for scaling c(p,t)
!
        chi=an0(l)*amn(l)**ex1
        sddm=sddm+chi
     enddo

! tau is calced in get_forcing
     f = tau(j,k)
!
! delm is the theoretical value of mass change, as per equation (33)
!
     delm=f*1.50*xi*sddm

!     Conditional for mixing assumption, if inhom=0, run
!     normal evap scheme

     if(inhom_evap == 0) then  
!
! preperation of some standard functions of the bin number for 
! use in the condensation calculations.  functions on determined
! depending on whether or not the average bin mass is within its
! bin range
!
!
        do i=1,lk
           if(amn(i) <  x_bin(i))then
              ff=1.0
           elseif(amn(i) >  x_bin(i+1))then
              ff=bin_res
           else
              ff=amn(i)/x_bin(i)
           endif
           fk=2.0*an0(i)/(x_bin(i)*(bin_res-1.0))*                               &
                (1.0 - 1./(bin_res-1.0) * (ff-1.0) )
           psik=2.*an0(i)/(x_bin(i)* (bin_res-1.0)*(bin_res-1.0))*               &
                (ff-1.0) 
           fac=1.0/(x_bin(i+1)-x_bin(i))
           s1(i)=fac *                                                       &
                (x_bin(i+1)* fk - x_bin(i) * psik )
           s2(i)=0.50*fac * (psik - fk)

! For solute
           s3(i)=fac* (bin_res * x_bin(i)*x_bin(i) * (fk-psik) )
           s4(i)=fac*                                                        &
                ( x_bin(i) * (bin_res*psik -fk ) )



!       ff=amn(i)/x_bin(i)
!     fk=2.0*an0(i)/(x_bin(i)*(bin_res-1.0))*
!    +       (1.0 - 1./(bin_res-1.0) * (ff-1.0) )
!     psik=2.*an0(i)/(x_bin(i)* (bin_res-1.0)*(bin_res-1.0))*
!    +       (ff-1.0) 

        enddo
! begin iterating through the bins for the actual evaporation 
! calculations, so long as f is greater than a certain threshold.
!
           if(abs(f) >  1.0e-20)then 
              do l=1,lk
                 yk(l)=(x_bin(l)**(ex2)-f)**(1.50)
                 zk(l)=min((x_bin(l+1)**(ex2)-f)**1.50,x_bin(lk+1))
                 if(abs(yk(l)-x_bin(l)) <  1.0e-20)yk(l)=x_bin(l)
                 if(abs(zk(l)-x_bin(l+1)) <  1.0e-20)zk(l)=x_bin(l+1)
                 if(x_bin(lk+1) >  yk(l))then
!
! locate the bin in which yk(l) is located, designate as i0
! 
                    il1=0
                    i0=0
                    in=0
                    do i=1,lk
                       if(yk(l) <  x_bin(i+1))then
                          i0=i 
                          exit
                       endif
                    enddo

                    if(l == 1)il1=i0
! 
! locate the bin in which zk(l) is located; for efficiency, start
! at i0 (since zk>yk), futhermore designate as in
! 
!this loop goes out of bounds do i=i0,lk+1, hence change
                    do i=i0,lk 
                       if(zk(l) <  x_bin(i+1))then
                          in=i
                          exit
                       endif
                    enddo

                    idiff=in-i0
                    i1=i0+1
! 
! below the parameters fx1,fy, fz,fxn,gg, and ggg are parameters
! appearing in the analytic solution to integrals in the reference
! equations 12 and 13.  furthermore functions f1 and f3 are defined
! because the form of this computation is repeatedly called and they
! simplify the code
! 
!               fx0 = sqrt(( x_bin(1)**ex2 + f))
                    fx1 = sqrt(( x_bin(i1)**ex2 + f))
                    fy  = sqrt(( yk(l)**ex2 + f))
                    fz  = sqrt(( zk(l)**ex2 + f))
! 
                    fy1  = sqrt(( yk(1)**ex2 + f))
! solutions to eqs. 8 and 9 (#3) for various values of nj.
!    idiff=0 implies yk and zk are located in the same bin.
!    idiff=1 implies zk is in a bin one higher than yk.
!    idiff>1 implies zk is in a bin more than one higher than yk. 
! the latter case should never appear given mass doubling between bins
!
! * *
! * xnevap is the no. of evaporated drops
! * It equals the integral from x_bin(1) to yk(1) of n(m)dm
! * which is precisely that # that leaves the drop distribution
! * *
                    xnevap=0.
                    xnevap=                                                  &
                         (yk(1)-x_bin(1))*(s1(1)+s2(1)*(yk(1)+x_bin(1)))
!              xmevap=0.
!              xmevap=s4(1)*(f1(yk(1),fy1,f)-f1(x_bin(1),fx0,f))
!    %                  + s3(1)*(f3(f,fy1)-f3(f,fx0))

                    if(idiff == 0)then
                       ank(j,k,l)=(zk(l)-yk(l))*(s1(i0)+s2(i0)*(zk(l)+yk(l))&
                            )
!                  smk(l)=(zk(l)-yk(l))*(s1s(i0)+s2s(i0)*(zk(l)+yk(l)))
                       amk(j,k,l)= s4(in)*(f1(zk(l),fz,f)-f1(yk(l),fy,f))    &
                            + s3(in)*(f3(f,fz)-f3(f,fy))

                    elseif(idiff == 1)then
                       fxn=sqrt(x_bin(in)**(ex2)+f)
                       ank(j,k,l)=(x_bin(i1)-yk(l))*(s1(i0)+s2(i0)*(x_bin(i1)+yk(l)))&
                            +(zk(l) -x_bin(in))*(s1(in)+s2(in)*(zk(l) +x_bin(in)))
                       amk(j,k,l)= s4(in)*f1(zk(l),fz,f) + s3(in)*f3(f,fz)   &
                            - s4(i0)*f1(yk(l),fy,f) - s3(i0)*f3(f,fy)       &
                            + (s4(i0)-s4(in))*f1(x_bin(in),fx1,f)               &
                            + (s3(i0)-s3(in))*f3(f,fx1)


                    elseif(idiff >= 2)then
                       fxn=sqrt(x_bin(in)**(ex2)+f)
                       sn1=0.0
                       sm1=0.0

                       do jl=i1,in-1
                          j1=jl+1
                          sn1=sn1+an0(jl)
                          gg=sqrt(x_bin(j1)**(ex2)+f)
                          ggg=sqrt(x_bin(jl)**(ex2)+f)
                          sm1=s4(jl)*(f1(x_bin(j1),gg,f) - f1(x_bin(jl),ggg,f))      &
                               +s3(jl)*(f3(f,gg)-f3(f,ggg)) + sm1
                       enddo
                       ank(j,k,l)=(x_bin(i1)-yk(l))*(s1(i0)+s2(i0)*(x_bin(i1)+yk(l)))&
                            +(zk(l) -x_bin(in))*(s1(in)+s2(in)*(zk(l)+x_bin(in)))   &
                            +sn1
                       amk(j,k,l)= s4(in)*(f1(zk(l),fz,f) - f1(x_bin(in),fxn,f)) &
                            + s4(i0)*(f1(x_bin(i1),fx1,f)- f1(yk(l),fy,f))      &
                            + s3(in)*(f3(f,fz)-f3(f,fxn))                   &
                            + s3(i0)*(f3(f,fx1)-f3(f,fy))                   &
                            + sm1


                    endif
                 else
                    ank(j,k,l)=0.
                    amk(j,k,l)=0.
                 endif
              enddo
           else
              do l=1,lk
                 amk(j,k,l)=am0(l)
                 ank(j,k,l)=an0(l)
              enddo
           endif

        else ! inhom_evap = 1, do ext inhom mixing
  
!     calculate CDNC and mass assuming extreme inhomogeneous mixing.
!     This is achieved by assuming the bulk mass change (sum(amk)-sum(amkd))
!     calculated by the bin micro is correct. Use the old (pre-evap)
!     and the new (post-evap) bulk mass to calc % change in mass. Then
!     use this % to reduce all bin mass and number
!          

           delm=f*1.50*xi*sddm          
           am1new(j,k) = 0.0 
           am1new(j,k) = am1old(j,k) + delm

           IF(am1old(j,k) > eps.and.                &
                am1new(j,k) > eps) then
              
              DO l = 1, lk
                 amk(j,k,l) = (am1new(j,k)/am1old(j,k))*am0_iter(l)
                 ank(j,k,l) = (am1new(j,k)/am1old(j,k))*an0_iter(l)
              ENDDO
           ELSE
              DO l = 1, lk
                 amk(j,k,l) = 0.0
                 ank(j,k,l) = 0.0
              ENDDO
           ENDIF
        ENDIF ! endif to inhom_evap = 0
        
      end subroutine EVAP_new

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      function fgv(P,K,SM1,RHON,DTcalc)
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! fgv is a function required for evaluation of the evaporation and
! condensation integrals. see #3, appendix b, eq. 32.
! fgv = 2/3 * c(p,t) * integral ds dt
      IMPLICIT NONE
 
!CALL PRAMETR
!CALL SD
!CALL CON
!CALL CONEV
      INTEGER K
      REAL FGV
      REAL, DIMENSION(KKP)::RHON
      REAL SM1
      REAL P
      REAL R0
      REAL DTcalc
!      REAL,PARAMETER:: DTcalc=2 !current model timestep from the LES

      if(sm1 < eps)then
         fgv=0.
      else
         fgv=2.*dds*(1.0-exp(-xi*yyy*zzz*sm1*dtcalc))/                  &
     &          (3.0*xi*yyy*sm1)
      endif
      end function fgv

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      function f3(a,b)
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! just a functional form commonly encountered in the evaporation
! routine, resulting from the analytical solution of the integrals
! for evaporation.  this function assumes that a is negative, as it
! should be for evaporation.
!
      IMPLICIT NONE

      REAL :: a,b,c,f3

      if (a >  0.0) print *, 'tau prob', a
      c=sqrt(-1.*a)
      f3=3.*a*(b - c*atan(b/c)) + b*b*b 

      end function f3
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
! GET_FORCING:   This routine gets the integrated forcing for cnd/evp 
! in accordance with the analytic/numerical method.  It also returns
! various measures of the excess vapor amount for use in nucleation 
! and in diagnostic tests.
! 
      SUBROUTINE GET_FORCING(J,K,RST,TEMP,ETA_0,ETA_D,                  &
     &             AM_gf,AN_gf,DT,TAU,ETA_N,ETA_A,VSW)

      IMPLICIT NONE

!CALL PRAMETR
!CALL RW
!CALL GRID1

! Intent (IN) arguments      
      INTEGER, INTENT(IN) ::J,K
      INTEGER L
      REAL, INTENT(IN),DIMENSION(JMINP:JMAXP,KKP) :: TEMP !temp calced in Tau
                                !using tref+th(1/exner function).
      REAL,INTENT(IN),DIMENSION(JMINP:JMAXP,KKP) :: RST
      REAL,INTENT(IN),DIMENSION(JMINP:JMAXP,KKP) :: ETA_D
      REAL,INTENT(IN) :: ETA_0
!      REAL,INTENT(IN),DIMENSION(JMINP:JMAXP,KKP,LK)::AMKOLD
!      REAL,INTENT(IN),DIMENSION(JMINP:JMAXP,KKP,LK)::ANKOLD
      REAL,INTENT(IN) :: DT
      REAL,INTENT(IN) :: VSW
! Intent (OUT) arguments 
      REAL,INTENT(OUT),DIMENSION(JMINP:JMAXP,KKP):: TAU,ETA_A,ETA_N
!     Local variables
      REAL,DIMENSION(JMINP:JMAXP,KKP)::ETA_DT
      REAL,DIMENSION(LK) :: AM_GF,AN_GF,AMN_GF
      REAL :: CHI,SDDM,ZZZ,XI,EX1,EX2,SCR,L0,test,pres
      REAL :: FSW
!     functions
!      REAL,EXTERNAL :: CPT,AA,ETA_INT,ETA_NEW

! L0 is set to 0 as there is no gas kinetics in cond_new and evap_new
! If gas kinetics included, set l0 to 6.5E-4

      PARAMETER (L0=0.0)
      PARAMETER (XI=0.9840654,EX1=1./3.,EX2=2./3.)
      
      pres = PREFN(K)*10

      SDDM=0.
      DO L=1,LK
         SCR=(AM_GF(L)/(AN_GF(L)+eps))**EX1
         SDDM=SDDM+AN_GF(L)*(SCR*SCR/(SCR+L0+EPS))*MAX(1.,VNTF(L)*VSW)
      ENDDO

      ZZZ=cpt(pres,j,k,temp(j,k),rst(j,k))
!
      IF(SDDM <  eps*lk)THEN
         TAU(J,K)=0.
         ETA_DT=0.
         ETA_N(J,K)=ETA_D(J,K)
         RETURN
      ENDIF
      
!
! Eta at the next timestep and the integral of eta are estimated
!
      CHI   = -XI*ZZZ*SDDM*AA(pres,j,k,temp(j,k),RST(j,k))
     
      if (abs(CHI) < eps) Then
         TAU(J,K)=0.
         ETA_DT=0.
         ETA_N(J,K)=ETA_D(J,K)
         RETURN
      endif

      ETA_N(J,K) = ETA_NEW(ETA_0,ETA_D(J,K),CHI,DT)
      ETA_DT(J,K) = ETA_INT(ETA_0,ETA_D(J,K),CHI,DT)
       
!
! Here the integral for eta over the timestep is bounded by its values
! at the endpoints, erroneous values of eta_dt sometimes are computed 
! for vanishingly small amounts of mass in the presence of evaporation,
! here the limiting eliminates erroneous values
!
!
      FSW = MAX(MIN(1.,MAX((MIN(ETA_0,ETA_N(J,K))                       &
     &      *DT-ETA_DT(J,K))*1.E15,0.))                                 &
     &    ,MIN(1.,MAX((ETA_DT(J,K)-MAX(ETA_0,ETA_N(J,K))*DT)*1.E15,0.)))
      
      ETA_DT(J,K)=ETA_DT(J,K)*(1.-FSW) + FSW*(ETA_0+ETA_N(J,K))         &
     &  *DT*0.5
!      
! Average value of eta over the timestep is computed for diagnostics
! as is the estimated change in mass
!
      ETA_A(J,K) = ETA_DT(J,K)/DT
      TAU(J,K)   = EX2*ZZZ*ETA_DT(J,K)

    END subroutine GET_FORCING
!
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!
! ETA_INT is the integral of the difference between the vapor mixing
! ratio and the saturation vapor mixing ratio over the timestep.  It
! is calculated in a manner consistent with the dynamical and micro-
! physical forcings via the analytical/numerical method.
!
! X ......generic variable for designating difference between the
!         vapor mixing ratio and the saturation vapor mixing ratio
! X0 .....value of X at the end of the previous timestep
! XD .....value of X after dynamics
! C ......function related to dynamic forcing for an excess vapor
!         equation of the form dX/dt = D + CX
! DT .... length of the timestep
!
      FUNCTION ETA_INT(X0,XD,C,DT)

      IMPLICIT NONE

      REAL :: ETA_INT

      REAL,INTENT(IN):: X0,XD,DT,C
      REAL :: CI,D

      D=(XD-X0)/DT
      CI=1./C
      ETA_INT=CI*((D*CI+X0)*(EXP(C*DT)-1.) - D*DT)

      END function eta_int
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!
! ETA_NEW is the estimated value of ETA at the end of the timestep 
! given the forcings.  It is calculated in a manner consistent with 
! the dynamical and micro-physical forcings via the 
! analytical/numerical method.
!
! X ......generic variable for designating difference between the
!         vapor mixing ratio and the saturation vapor mixing ratio
! X0 .....value of X at the end of the previous timestep
! XD .....value of X after dynamics
! C ......function related to dynamic forcing for an excess vapor
!         equation of the form dX/dt = D + CX
! DT .... length of the timestep
!
      FUNCTION ETA_NEW(X0,XD,C,DT)

      IMPLICIT NONE

      REAL :: ETA_NEW

      REAL, INTENT(IN):: X0,XD,C,DT
      REAL :: CI,D

      D=(XD-X0)/DT
      CI=1./C
      ETA_NEW=(D*CI+X0)*EXP(C*DT) - D*CI

      END function ETA_NEW
!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!

end module
