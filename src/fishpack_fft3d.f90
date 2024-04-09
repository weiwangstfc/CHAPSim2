module fishpack_fft3d
  USE precision_mod
  include 'fftpack.f'
  implicit none

  INTEGER(4) :: LPEROD, MPEROD, NPEROD
  INTEGER(4) :: LP, MP, NP
  INTEGER(4) :: L, M, NG, NL,ML
  INTEGER(4) :: LDIMF, MDIMF
  INTEGER(4) :: WSZ
  REAL(WP)    :: C1, C2
  REAL(WP)    :: SCALX, SCALY
  REAL(WP), ALLOCATABLE    :: A(:), B(:), C(:), D(:), BB(:)
  REAL(WP), ALLOCATABLE    :: XRT(:),YRT(:), WX(:), WY(:)
  REAL(WP), ALLOCATABLE    :: FR(:,:,:), FK(:,:,:), T(:)
  REAL(WP), ALLOCATABLE    :: F_io   (:,:,:)
!----------------------------------------------------------------------------------------------------------
! processures
!----------------------------------------------------------------------------------------------------------
  public :: FISHPACK_POIS3D_INIT
  public :: FISHPACK_POIS3D_SIMPLE

  private :: TRID0
  private :: FFTPACK_XZ
  private :: FFTPACK_ROOT

contains
!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FISHPACK_POIS3D_INIT(dm)
    use parameters_constant_mod
    use udf_type_mod
    IMPLICIT NONE
    type(t_domain), intent(in) :: dm
  
    INTEGER(4) :: K
  
    IF( dm%ibcx(1, 1) == IBC_PERIODIC  .and. dm%ibcx(2, 1) == IBC_PERIODIC  )  LPEROD = 0
    IF( dm%ibcx(1, 1) == IBC_DIRICHLET .and. dm%ibcx(2, 1) == IBC_DIRICHLET )  LPEROD = 1
    IF( dm%ibcx(1, 1) == IBC_DIRICHLET .and. dm%ibcx(2, 1) == IBC_NEUMANN   )  LPEROD = 2
    IF( dm%ibcx(1, 1) == IBC_NEUMANN   .and. dm%ibcx(2, 1) == IBC_NEUMANN   )  LPEROD = 3
    IF( dm%ibcx(1, 1) == IBC_NEUMANN   .and. dm%ibcx(2, 1) == IBC_DIRICHLET )  LPEROD = 4
  
    IF( dm%ibcz(1, 1) == IBC_PERIODIC  .and. dm%ibcz(2, 1) == IBC_PERIODIC  )  LPEROD = 0
    IF( dm%ibcz(1, 1) == IBC_DIRICHLET .and. dm%ibcz(2, 1) == IBC_DIRICHLET )  LPEROD = 1
    IF( dm%ibcz(1, 1) == IBC_DIRICHLET .and. dm%ibcz(2, 1) == IBC_NEUMANN   )  LPEROD = 2
    IF( dm%ibcz(1, 1) == IBC_NEUMANN   .and. dm%ibcz(2, 1) == IBC_NEUMANN   )  LPEROD = 3
    IF( dm%ibcz(1, 1) == IBC_NEUMANN   .and. dm%ibcz(2, 1) == IBC_DIRICHLET )  LPEROD = 4
  
    NPEROD = 1
  
    LP = LPEROD + 1
    MP = MPEROD + 1
    NP = NPEROD + 1
  
    L = dm%nc(1) ! L: dim-x, NCL1_io
    M = dm%nc(3) ! M: dim-z, NCL3
    NG= dm%nc(2) ! N: dim-y, NCL2
    NL= dm%      !N2DO(MYID)
    ML= dm%      !N3DO(MYID)
  
    IF( (L .LE. 3) .or. (M .LE. 3) .or. (NG .LE. 3) ) &
    CALL ERRHDL('Dimensions in poisson solver should be larger than 3!',myid)
  
    C1 = dm%h2r(1) ! DXQI
    C2 = dm%h2r(3) ! DZQI
  
  ! y-dim global variable
    ALLOCATE (A (NG) ) ; A = 0.0_WP
    ALLOCATE (B (NG) ) ; B = 0.0_WP
    ALLOCATE (C (NG) ) ; C = 0.0_WP
    ALLOCATE (D (NG) ) ; D = 0.0_WP
    ALLOCATE (BB(NG) ) ; BB= 0.0_WP
  ! x-dim global variable
    ALLOCATE (XRT(L))  ; XRT = 0.0_WP
  ! z-dim global variable
    ALLOCATE (YRT(M))  ; YRT = 0.0_WP
    ALLOCATE (WX(WSZ)) ; WX = 0.0_WP
    ALLOCATE (WY(WSZ)) ; WY = 0.0_WP
    ALLOCATE (FR(L, M,  NL) ) ; FR = 0.0_WP
    ALLOCATE (FK(L, ML, NG) ) ; FK = 0.0_WP
    ALLOCATE (T(MAX0(L, M, NG)) ) ; T = 0.0_WP
    
    ALLOCATE ( F_io   (NCL1_io,NCL2,N3DO(0) )     )       ;  F_io= 0.0_WP

  
    DO K=1,NG
        A(K) = AMPH(K)
        B(K) = ACPH(K)
        C(K) = APPH(K)
    END DO
  
    CALL FFTPACK_ROOT
  
    RETURN
  
  END SUBROUTINE
!==========================================================================================================
!==========================================================================================================
SUBROUTINE FFTPACK_ROOT ! nothing changed.
    IMPLICIT NONE
    REAL(WP)    :: DUM
    REAL(WP)    :: DY, DJ
    REAL(WP),external :: PIMACH
    REAL(WP)    :: DX, DI
  
    INTEGER(4) :: I,J
    INTEGER(4) :: LRDEL
    INTEGER(4) :: MRDEL
    INTEGER(4) :: LR, MR
  
    PI = PIMACH(DUM)
    LR = L
    MR = M
!C
!C     GENERATE TRANSFORM ROOTS FOR X DIRECTION
!C
      LRDEL = ((LP-1)*(LP-3)*(LP-5))/3
      SCALX = DBLE(LR+LRDEL)
      DX = PI/(2.0_wp*SCALX)
      GO TO (108,103,101,102,101),LP
101   DI = 0.50_wp
      SCALX = 2.0_wp*SCALX
      GO TO 104
102   DI = 1.00_wp
      GO TO 104
103   DI = 0.00_wp
104   DO 105 I=1,LR
          XRT(I) = -4.0_wp*C1*(SIN((DBLE(I)-DI)*DX))**2
105   CONTINUE
      SCALX = 2.0_wp*SCALX
      GO TO (112,106,110,107,111),LP
106   CALL SINTI (LR,WX)
      GO TO 112
107   CALL COSTI (LR,WX)
      GO TO 112
108   XRT(1) = 0.0_wp
      XRT(LR) = -4.0_wp*C1
      DO 109 I=3,LR,2
          XRT(I-1) = -4.0_wp*C1*(SIN(DBLE((I-1))*DX))**2
          XRT(I) = XRT(I-1)
109   CONTINUE
      CALL RFFTI (LR,WX)
      GO TO 112
110   CALL SINQI (LR,WX)
      GO TO 112
111   CALL COSQI (LR,WX)
112   CONTINUE

!C
!C     GENERATE TRANSFORM ROOTS FOR Y DIRECTION (Z IN REAL)
!C  
    MRDEL = ((MP-1)*(MP-3)*(MP-5))/3
    SCALY = DBLE(MR+MRDEL)
    DY = PI/(2.0_wp*SCALY)
    GO TO (120,115,113,114,113),MP
113 DJ = 0.50_wp
    SCALY = 2.0_wp*SCALY
    GO TO 116
114 DJ = 1.00_wp
    GO TO 116
115 DJ = 0.00_wp
116 DO 117 J=1,MR
        YRT(J) = -4.0_wp*C2*(SIN((DBLE(J)-DJ)*DY))**2
117 CONTINUE
    SCALY = 2.0_wp*SCALY
    GO TO (124,118,122,119,123),MP
118 CALL SINTI (MR,WY)
    GO TO 124
119 CALL COSTI (MR,WY)
    GO TO 124
120 YRT(1) = 0.0_wp
    YRT(MR) = -4.0_wp*C2
    DO 121 J=3,MR,2
        YRT(J-1) = -4.0_wp*C2*(SIN(DBLE((J-1))*DY))**2
        YRT(J) = YRT(J-1)
121 CONTINUE
    CALL RFFTI (MR,WY)
    GO TO 124
122 CALL SINQI (MR,WY)
    GO TO 124
123 CALL COSQI (MR,WY)
124 CONTINUE

    RETURN  
    
  END SUBROUTINE
!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FISHPACK_POIS3D_SIMPLE(dm, dtmp, phi)
    IMPLICIT NONE
    INTEGER(4) :: IFWRD
    INTEGER(4) :: I, J, K, JJ
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(in) :: phi

!========RECONSTRUCT RHS TO FIT POIS3D.========
    DO I=1,L ! x-global
      DO K=1,M ! z-local
        DO J=1,NL ! y-local
          FR(I,K,J) = RHSLLPHI_io(I,J,K) ! third-dirction is local
        END DO
      END DO
    END DO     

!========FORWARD FFT IN X AND Z DIRECTION========
    IFWRD = 1
    CALL FFTPACK_XZ(IFWRD)

!========RESTORE RHSLLPHI========
    DO I=1,L
        DO K=1,M
            DO J=1,NL
                RHSLLPHI_io(I,J,K) = FR(I,K,J)
            END DO
        END DO
    END DO

!========TRANSPORT Y-DECOMP TO K-DECOMP========
    !CALL TRASP_Y2Z_RHSLLPHI_io
    CALL TRASP23_Y2Z(NCL1_IO, 1, N2DO(0), RHSLLPHI_io, F_IO)

!========CONSTRUCT DATA FOR TDMA========
    DO I=1,L
        DO K=1,N3DO(MYID)
            DO J=1,NG
                FK(I,K,J) = F_io(I,J,K)
            END DO
        END DO
    END DO 

!========TDMA IN Y DIRECTION FOR PART OF FR(:,PART,:)========
    DO I=1,L
      DO J=1,N3DO(MYID)
        JJ = KCL2G(J) 
        DO K=1,NG
          BB(K) = B(K)+XRT(I)/RCCI2(K)+YRT(JJ)
          T(K) = FK(I,J,K)
        END DO
    
        CALL TRID0
    
        DO K=1,NG
          FK(I,J,K) = T(K)
        END DO
      
      END DO
    END DO

! ========RE-CONSTRUCT DATA BACK========
    DO I=1,L
      DO K=1,N3DO(MYID)
        DO J=1,NG
          F_io(I,J,K) = FK(I,K,J) 
        END DO
      END DO
    END DO

! ========TRANSPORT Z-DECOMP TO Y-DECOMP========
    !CALL TRASP_Z2Y_RHSLLPHI_io  
    CALL TRASP23_Z2Y(NCL1_IO, 1, N2DO(0), RHSLLPHI_io, F_io)

! ========RESTORE RHSLLPHI========
    DO I=1,L
      DO K=1,M
        DO J=1,NL
          FR(I,K,J) = RHSLLPHI_io(I,J,K) 
        END DO
      END DO
    END DO 

!========BACKWARD FFT IN Z AND X DIRECTIONS========
    IFWRD = 2
    CALL FFTPACK_XZ(IFWRD)

!========SCALE THE CALCULATED VALUE========
    DO 167 I=1,L
      DO 166 J=1,M
        DO 165 K=1,NL
          FR(I,J,K) = FR(I,J,K)/(SCALX*SCALY)
165     CONTINUE
166   CONTINUE
167 CONTINUE

!========RE-STORE AND ASSIGN DATA TO DPH========
    DPH_IO = 0.0_WP
      DO I=1,L
        DO K=1,M
          DO J=1,NL
            DPH_io(I,J,K) = FR(I,K,J)
          END DO
        END DO
    END DO  
      
    RETURN
    
  END SUBROUTINE

!==========================================================================================================
!==========================================================================================================   
  SUBROUTINE FFTPACK_XZ(IFWRD)
!     FOR FFT TRANSFORM FROM SPACE TO WAVENUMBER, IFWRD =1, IS=  1
!     FOR FFT TRANSFORM FROM WAVENUMBER TO SPACE, IFWRD =2, IS= -1    
      IMPLICIT NONE
      INTEGER(4) :: IFWRD
      INTEGER(4) :: I,J,K
      INTEGER(4) :: LR, MR, NR
      
      LR = L
      MR = M
      NR = NL
         
      GO TO(125,142) IFWRD 
      
125 CONTINUE        
!
!     TRANSFORM X, in x-pencil
!      
    DO 141 J=1,MR ! z-global -> z-local 
      DO 140 K=1,NR ! y-local
        DO 126 I=1,LR ! x-global
            T(I) = FR(I,J,K)
126       CONTINUE
        GO TO (127,130,131,134,135),LP
127       GO TO (128,129),IFWRD
128       CALL RFFTF (LR,T,WX)
        GO TO 138
129       CALL RFFTB (LR,T,WX)
        GO TO 138
130       CALL SINT (LR,T,WX)
        GO TO 138
131       GO TO (132,133),IFWRD
132       CALL SINQF (LR,T,WX)
        GO TO 138
133       CALL SINQB (LR,T,WX)
        GO TO 138
134       CALL COST (LR,T,WX)
        GO TO 138
135       GO TO (136,137),IFWRD
136       CALL COSQF (LR,T,WX)
        GO TO 138
137       CALL COSQB (LR,T,WX)
138       CONTINUE
        DO 139 I=1,LR
            FR(I,J,K) = T(I)
139     CONTINUE
140   CONTINUE
141 CONTINUE
  GO TO (142,159),IFWRD
  
!C
!C     TRANSFORM Y
!C
142 CONTINUE
    DO 158 I=1,LR ! x-global -> x_local
      DO 157 K=1,NR !y-local
        DO 143 J=1,MR !z-global
            T(J) = FR(I,J,K)
143     CONTINUE
        GO TO (144,147,148,151,152),MP
144       GO TO (145,146),IFWRD
145       CALL RFFTF (MR,T,WY)
        GO TO 155
146       CALL RFFTB (MR,T,WY)
        GO TO 155
147       CALL SINT (MR,T,WY)
        GO TO 155
148       GO TO (149,150),IFWRD
149       CALL SINQF (MR,T,WY)
        GO TO 155
150       CALL SINQB (MR,T,WY)
        GO TO 155
151       CALL COST (MR,T,WY)
        GO TO 155
152       GO TO (153,154),IFWRD
153       CALL COSQF (MR,T,WY)
        GO TO 155
154       CALL COSQB (MR,T,WY)
155       CONTINUE
        DO 156 J=1,MR
            FR(I,J,K) = T(J)
156     CONTINUE
157   CONTINUE
158 CONTINUE
    GO TO (159,125),IFWRD
  
159 CONTINUE

  RETURN
  END SUBROUTINE 


!==========================================================================================================
!==========================================================================================================
  SUBROUTINE TRID0

    IMPLICIT NONE
    INTEGER(4) :: NR, MM1
    REAL(WP)    :: Z
    INTEGER(4) :: I,IP
    
    NR = NG
    MM1 = NR-1
    Z = 1.0_wp/BB(1)
    D(1) = C(1)*Z
    T(1) = T(1)*Z
    DO 101 I=2,MM1
        Z = 1.0_wp/(BB(I)-A(I)*D(I-1))
        D(I) = C(I)*Z
        T(I) = (T(I)-A(I)*T(I-1))*Z
101 CONTINUE
    Z = BB(NR)-A(NR)*D(MM1)
    IF (Z .NE. 0.0_wp) GO TO 102
    T(NR) = 0.0_wp
    GO TO 103
102 T(NR) = (T(NR)-A(NR)*T(MM1))/Z
103 CONTINUE
    DO 104 IP=1,MM1
        I = NR-IP
        T(I) = T(I)-D(I)*T(I+1)
104 CONTINUE
  RETURN
    
  END SUBROUTINE

end module