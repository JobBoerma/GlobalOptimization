! OBJECTIVE.f90 contains the simulation and objective 
!               (model-data distance) function
! --------------------------------------------------------------------
! --------------------------------------------------------------------

module OBJECTIVE
use, intrinsic :: iso_c_binding, only: c_int
    use myParams
    use genericParams
    Use nrtype
    implicit none
    PRIVATE
    PUBLIC objFun, dfovec, initial0, GetData


    interface

      function myFortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: myFortSleep
          integer (c_int), intent (in), VALUE :: seconds
      end function myFortSleep

   end interface

    real(8), parameter:: MISSING = -99.9
  REAL(DP) :: eps=epsilon(1.0_DP)
	REAL(SP), DIMENSION(:), POINTER :: fmin_fvecp

contains

FUNCTION objFun(theta)
! Sum of squares distance
! -------------------------------------------------------------------------
use genericParams
implicit none

    REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta
    REAL(DP),DIMENSION(p_nx)            :: pencons 
    REAL(DP)                            :: objFun0(1,1), objFun, penalty
    REAL(DP),DIMENSION(p_nmom)          :: Fout, FoutW
    REAL(DP),DIMENSION(p_nmom,1)        :: F,FW
    REAL(DP),DIMENSION(p_nmom,p_nmom)   :: W, Wsqr

    integer :: ii

    storeTheta=theta

    ! penalty
    pencons=10.0d5
    penalty=0.
    do ii=1,p_nx
        penalty = penalty + pencons(ii)*(min(0.D0,theta(ii)-p_bound(ii,1)))**2.D0
        penalty = penalty + pencons(ii)*(max(0.D0,theta(ii)-p_bound(ii,2)))**2.D0
    enddo

    call dfovec(p_nx,p_nmom,theta,Fout)

    ! Weighting matrix
    W=0.
    forall(ii=1:p_nmom) W(ii,ii)=1.
    forall(ii=1:NmomCS) W(ii,ii)=.8
    forall(ii=NmomCS+1:p_nmom) W(ii,ii)=.2

    F(:,1) = Fout
    FW = matmul(sqrt(W),F)
    FoutW = FW(:,1)

    objFun0=dot_product(FoutW,FoutW)
    objFun = objFun0(1,1)+penalty

END FUNCTION objFun

! -------------------------------------------------------------------------
SUBROUTINE GetData
use nrtype
use myParams
use genericParams
    integer:: i
    !allocate(Cx(p_nmom),gama(p_nmom))
    
    ! -----------------------------------------------------------------
    ! READ IN DATA - OUTPUT OF MATLAB CODE INPUT_FORTRAN.M
    !
    !       * Cx - data moments (detrended)
    !       * xt - normalized series of 1-year GDP growth
    !       * mt - 3-year mean earnings growth
    !       * gamma - scaling constant for moments ~ 0
    !
    OPEN(2, FILE = datadir // '/Cx.dat');        read(2,*) Cx_csA;      CLOSE(2);
    OPEN(2, FILE = datadir // '/V_ageprof.dat'); read(2,*) Cx_lc;       CLOSE(2);
    OPEN(2, FILE = datadir // '/xt.dat');        read(2,*) xt;          CLOSE(2);
    OPEN(2, FILE = datadir // '/mt1.dat');        read(2,*) mt1;          CLOSE(2);
    OPEN(2, FILE = datadir // '/gamma.dat');     read(2,*) gama_csA;    CLOSE(2);
    !
    ! -----------------------------------------------------------------

    Cx_cs = Cx_csA(1:p_nmom)
    gama_cs = gama_csA(1:p_nmom)
    gama_lc = 0.D0

    Cx = (/Cx_cs,Cx_lc/) 
    gama = (/gama_cs,gama_lc/)

    VAR25 = Cx_lc(1)
    ! id* series denote available years for the * change (* = sr,mr,lr)
    ! according to frequency in the PSID
    idsr=1
    idmr=1
    idlr=1
    if (country == 2) then
        i=20
        do while (i<=Tsim_sr)
            idsr(i) =  0
            i=i+2
        enddo
        i=18
        do while (i<=Tsim_mr)
            idmr(i) =  0
            i=i+2
        enddo
        i=16
        do while (i<=Tsim_lr)
            idlr(i) =  0
            i=i+2
        enddo
    endif

    ! -----------------------------------------------------------------
    ! GENERATE RANDOM VARIABLES
    !
    !   randva has dimension (nhhsim*5,Tlev_s,nsim)
    !   In each simulation nsim, a new matrix (N*5,T) of random vars 
    !   will be read. See dofvec below for more details.
    !
    call gen_rand
    !
    ! -----------------------------------------------------------------

END SUBROUTINE

! -------------------------------------------------------------------------
SUBROUTINE dfovec(np, nm, x, Fnout)
! Scaled distance (vector of nmoments size). For each moment i
!
!       Fnout(i)= [Csim(i)-Cx(i)] / [|Cx(i)|+gamma(i)]
!
!       where:
!       * Cx - vector of moments from data
!       * Csim - simulated vector of moments (average over nsim simulations)
!       * gamma - scaling factor for moments close to 0
!
! -------------------------------------------------------------------------
use nrtype
use myParams
use genericParams
IMPLICIT NONE

    INTEGER, INTENT(IN)     :: np, nm              
    REAL(DP), DIMENSION(np), INTENT(IN)  :: x
    REAL(DP), DIMENSION(nm,1),INTENT(OUT) :: Fnout
    
    ! Random variables: colums are periods, rows are
    ! 
    !       First nhhsim rows: uniform for mixture
    !       Next  nhhsim rows: normal for eta1
    !       Next  nhhsim rows: normal for eta2
    !       Next  nhhsim rows: normal for eta1
    !       Next  nhhsim rows: normal for epsilon 
    !
    ! 3rd dimension is one for each simulated economy
    !
    REAL(DP), DIMENSION(nhhsim*5,Tlev_s) :: randvar   ! Random variables, for a given simulation

    ! Simulation moments
    REAL(DP), DIMENSION(nm,nsim) ::C_THETA_M
    REAL(DP), DIMENSION(NmomCS_all+nage,nsim) ::C_THETA_M_ALL
    REAL(DP), DIMENSION(nm) ::C_THETA_BAR, C_x_abs
    REAL(DP), DIMENSION(NmomCS_all+nage) ::C_THETA_BAR_ALL
    logical:: chk

    !Variables for sleep call
    integer(c_int) :: mytime, dur

    INTEGER(i4b) :: i,i_m, i1, i2

    !Count number of function evaluations
     fe_counter=fe_counter+1

    ! -----------------------------------------------------------------
    ! SIMULATE NSIM TIMES
    !
    do i_m = 1,nsim
        randvar=randva(:,:,i_m)
        call f_moments(np,nm,x,randvar, xt, C_THETA_M(:,i_m),C_THETA_M_ALL(:,i_m) )
    enddo
    !
    ! DONE SIMULATING nsim TIMES
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! SCALED DISTANCE VECTOR FNOUT
    ! 
    C_THETA_BAR=sum(C_THETA_M,2)/real(nsim)
    C_THETA_BAR_ALL=sum(C_THETA_M_ALL,2)/real(nsim)
    C_x_abs= abs(Cx)
    FORALL (i=1:nm) Fnout(i,1)=(C_THETA_BAR(i)-Cx(i))/(C_x_abs(i)+gama(i))
    !                
    ! -- DONE CACULATING FNOUT
    ! -----------------------------------------------------------------

    mytime=3;
    dur=myFortSleep(mytime);

END SUBROUTINE dfovec

! -------------------------------------------------------------------------
SUBROUTINE f_moments(np, nm, x, randv, bcx, Cm, Cm_all)
! Simulates model and computes moments of income growth Cm
! -------------------------------------------------------------------------
use nrtype
use myParams
use genericParams
implicit none

    ! Arguments
    INTEGER, INTENT(IN)     :: np, nm
    REAL(DP), DIMENSION(np), INTENT(IN)  :: x
    REAL(DP), DIMENSION(Tlev_s), INTENT(IN)  :: bcx
    REAL(DP), DIMENSION(nhhsim*5,Tlev_s), INTENT(IN)  :: randv
    REAL(DP), DIMENSION(nm), INTENT(OUT)  :: Cm
    REAL(DP), DIMENSION(NmomCS_all+nage), INTENT(OUT)  :: Cm_all

    REAL(DP) :: Cm_lc(nage), Cm_cs(NmomCS)
    REAL(DP), DIMENSION(Tlev_s) :: gt

    ! Stochastic processes parameters
    real(dp), dimension(Tlev_s) :: mu_bar,mm1, mu_eta1, mu_eta2, mu_eta3, bc_x
    real(dp), dimension(nhhsim,Tlev_s) :: zerosNT, mu_eta1_m, mu_eta2_m, mu_eta3_m

    ! Income changes
    real(dp), dimension(nhhsim,Tsim_sr) :: INCCHANGE_SRtmp
    real(dp), dimension(nhhsim,Tsim_mr) :: INCCHANGE_MRtmp
    real(dp), dimension(nhhsim,Tsim_lr) :: INCCHANGE_LRtmp

    real(dp), dimension(nhhsim,Tdat_sr) :: INCCHANGE_SR 
    real(dp), dimension(nhhsim,Tdat_mr) :: INCCHANGE_MR 
    real(dp), dimension(nhhsim,Tdat_lr) :: INCCHANGE_LR 

    ! Moments
    real(dp), dimension(Tdat_sr) :: p10_SR, p50_SR, p90_SR, mea_SR
    real(dp), dimension(Tdat_mr) :: p10_MR, p50_MR, p90_MR, mea_MR
    real(dp), dimension(Tdat_lr) :: p10_LR, p50_LR, p90_LR, mea_LR
    real(dp) :: perc(3), pout(3)

    real(dp), dimension(NmomCS_all):: Cm_cs_all
    ! LC stuff
    real(dp) :: VARTRAN, VARPERM_1
    real(dp), dimension(nage) :: VARAGE_1, VARAGE_1_norm
    ! parameters
    real(dp) :: sig_eps, p1, p2, p3, sig_eta1, sig_eta2, sig_eta3, mu2, mu3, phi

    ! Random variables
    real(dp),dimension(nhhsim,Tlev_s) :: uniprob, &
                                       randEta1, randEta3, randEta2, &
                                       randEps, eps, eta1, eta2, eta3, eta
    integer(i4b), dimension(nhhsim,Tlev_s) :: mix1,mix2,mix3                                    
                                       
    integer :: i,ii,jj,i_t,i_a
    
    zerosNT=0.
    
    ! Parameters
    sig_eps = x(1)
    p1      = x(2)
    p2      = 1./2. - p1/2.
    p3      = p2
    mu2     = x(3)
    mu3     = x(4)
    sig_eta1= x(5)
    sig_eta2= x(6)
    sig_eta3= sig_eta2
    phi     = x(7)
    !m = x(8:8+Tlev_s-1)

    ! Aggregate risk series (bcx=-GDP log diff 1 year)
    bc_x = phi*bcx

    ! Read-in standard random series (first nhhsim lines are uniform, rest are
    ! std normal)
    i=0;   uniprob  = randv(nhhsim*i+1:nhhsim*(i+1),:); 
    i=i+1; randEta1 = randv(nhhsim*i+1:nhhsim*(i+1),:);
    i=i+1; randEta2 = randv(nhhsim*i+1:nhhsim*(i+1),:);
    i=i+1; randEta3 = randv(nhhsim*i+1:nhhsim*(i+1),:);
    i=i+1; randEps  = randv(nhhsim*i+1:nhhsim*(i+1),:);

    ! Create epsilon series
    eps = -(sig_eps**2.)/2. + sig_eps*randEps
    
    ! Create eta series
    mu_bar = -log(p1*exp(sig_eta1**2./2.)+ &
                  p2*exp(mu2-bc_x+(sig_eta2**2./2.))+ &
                  p3*exp(mu3-bc_x+(sig_eta3**2./2.))) 

    mm1=(/mt1,0.D0/)
    gt = mm1-mu_bar-p2*mu2-p3*mu3+(p2+p3)*bc_x
    mu_eta1 = mu_bar + gt
    mu_eta2 = mu_bar + mu2 - bc_x + gt
    mu_eta3 = mu_bar + mu3 - bc_x + gt

    ! mu_eta1_m is (nhhsim,Tlev_s)
    forall(i_t=1:Tlev_s) mu_eta1_m(:,i_t) = mu_eta1(i_t)
    forall(i_t=1:Tlev_s) mu_eta2_m(:,i_t) = mu_eta2(i_t)
    forall(i_t=1:Tlev_s) mu_eta3_m(:,i_t) = mu_eta3(i_t)

    eta1 = sig_eta1*randEta1 + mu_eta1_m
    eta2 = sig_eta2*randEta2 + mu_eta2_m
    eta3 = sig_eta3*randEta3 + mu_eta3_m

    mix1=0; mix2=0; mix3=0;
    WHERE(uniprob<=p1)
        mix1 = 1
    ELSEWHERE(uniprob<=(p1+p2))
        mix2 = 1
    ELSEWHERE
        mix3 = 1
    END WHERE
    
    eta = real(mix1)*eta1 + real(mix2)*eta2 + real(mix3)*eta3

    ! Income Changes
    !
    ! Dy[t,t+1] = Deps[t,t+1] + eta[t+1]
    !
    INCCHANGE_SRtmp = eps(:,1+L_sr:Tlev_s) - eps(:,1:Tlev_s-L_sr)
    do i = 1,L_sr
        INCCHANGE_SRtmp = INCCHANGE_SRtmp + eta(:,1+(i-1):Tlev_s-L_sr+(i-1))
    enddo
    INCCHANGE_MRtmp = eps(:,1+L_mr:Tlev_s) - eps(:,1:Tlev_s-L_mr)
    do i = 1,L_mr
        INCCHANGE_MRtmp = INCCHANGE_MRtmp + eta(:,1+(i-1):Tlev_s-L_mr+(i-1))
    enddo
    INCCHANGE_LRtmp = eps(:,1+L_lr:Tlev_s) - eps(:,1:Tlev_s-L_lr)
    do i = 1,L_lr
        INCCHANGE_LRtmp = INCCHANGE_LRtmp + eta(:,1+(i-1):Tlev_s-L_lr+(i-1))
    enddo

    ! get rid of changes not in the data for PSID
    ! id* series denote available years for the * change (* = sr,mr,lr)
    if (country == 2) then
        ii = 1; jj = 1;
        do while (ii<=Tsim_sr)
                if (idsr(ii)==0) ii=ii+1
                INCCHANGE_SR(:,jj) = INCCHANGE_SRtmp(:,ii)
                ii = ii+1; jj = jj+1;
        enddo
        ii = 1; jj = 1;
        do while (ii<=Tsim_mr)
                if (idmr(ii)==0) ii=ii+1
                INCCHANGE_MR(:,jj) = INCCHANGE_MRtmp(:,ii)
                ii = ii+1; jj = jj+1;
        enddo
        ii = 1; jj = 1;
        do while (ii<=Tsim_lr)
                if (idlr(ii)==0) ii=ii+1
                INCCHANGE_LR(:,jj) = INCCHANGE_LRtmp(:,ii)
                ii = ii+1; jj = jj+1;
        enddo
    else 
        INCCHANGE_SR = INCCHANGE_SRtmp(:,1:Tdat_sr)
        INCCHANGE_MR = INCCHANGE_MRtmp(:,1:Tdat_mr)
        INCCHANGE_LR = INCCHANGE_LRtmp(:,1:Tdat_lr)
    endif

    ! Moments
    perc=(/.1,.5,.9/)
    do i = 1,Tdat_sr
        mea_SR(i) = sum(INCCHANGE_SR(:,i))/size(INCCHANGE_SR(:,i))
        call percentile(INCCHANGE_SR(:,i),nhhsim,3,perc,pout)
        p10_SR(i)=pout(1)
        p50_SR(i)=pout(2)
        p90_SR(i)=pout(3)
    enddo
    do i = 1,Tdat_mr
        mea_MR(i) = sum(INCCHANGE_MR(:,i))/size(INCCHANGE_MR(:,i))
        call percentile(INCCHANGE_MR(:,i),nhhsim,3,perc,pout)
        p10_MR(i)=pout(1)
        p50_MR(i)=pout(2)
        p90_MR(i)=pout(3)
    enddo
    do i = 1,Tdat_lr
        mea_LR(i) = sum(INCCHANGE_LR(:,i))/size(INCCHANGE_LR(:,i))
        call percentile(INCCHANGE_LR(:,i),nhhsim,3,perc,pout)
        p10_LR(i)=pout(1)
        p50_LR(i)=pout(2)
        p90_LR(i)=pout(3)
    enddo

    !Life-cycle
    ! Version 1: time-invariant distribution of permanent shocks (imposing xt=0)
        ! variance of transitory and pemanent shock:
        VARTRAN = sig_eps
        call f_varperm(x,VARPERM_1)
        !VARPERM_1 = 1. 

        ! age profile of variance: age 25-60 
        VARAGE_1(1) = VARPERM_1
        do i_a = 2,nage
            VARAGE_1(i_a) = VARAGE_1(i_a-1)+VARPERM_1
        enddo
        VARAGE_1 = VARAGE_1 + VARTRAN
        
        !normalized version with VARAGE(1) = 0; corresponding to
        ! dummies
        VARAGE_1_norm = VARAGE_1 - VARTRAN - VARPERM_1
        Cm_lc=VAR25+VARAGE_1_norm
    
    ! Build matrix
    Cm_cs_all = (/p10_SR, p10_MR, p10_LR, &
               p50_SR, p50_MR, p50_LR, &
               p90_SR, p90_MR, p90_LR, &
               mea_SR, mea_MR, mea_LR/)

    Cm_cs = Cm_cs_all(1:p_nmom)

    Cm=(/Cm_cs,Cm_lc/)
    Cm_all=(/Cm_cs_all,Cm_lc/)

    open(unit=222,file=savedir//'/Cmout.dat',action='write')
	write(222,*) Cm_all
    close(222)
contains

        subroutine f_varperm(x,var)
        implicit none

        real(dp), dimension(p_nx), intent(in):: x
        real(dp), intent(out) :: var

        real(dp) :: MU1t,MU2t,MU3t,MUt
        real(dp) :: xsig_eps, xp1,xp2,xp3,xmu2,xmu3,xsig_eta1,xsig_eta2,xsig_eta3,xphi

    xsig_eps = x(1)
    xp1      = x(2)
    xp2      = 1./2. - p1/2.
    xp3      = p2
    xmu2     = x(3)
    xmu3     = x(4)
    xsig_eta1= x(5)
    xsig_eta2= x(6)
    xsig_eta3= sig_eta2
    xphi     = x(7)

        MU1t = -log( xp1*exp(xsig_eta1**2./real(2)) + xp2*exp(xmu2+xsig_eta2**2./real(2)) + xp3*exp(xmu3+xsig_eta3**2./real(2)) )
        MU2t = MU1t + xmu2
        MU3t = MU1t + xmu3

        MUt = xp1*MU1t + xp2*MU2t + xp3*MU3t
        var = xp1*(xsig_eta1**2. + MU1t**2.) + xp2*(xsig_eta2**2. + MU2t**2.) + xp3*(xsig_eta2**2. + MU3t**2.) - MUt**2.
end subroutine
      !-----------------------------------------------------------------------
      ! Calculate percentile
      !
      subroutine percentile(x, length, nl, per, oout)
      integer :: length,nl    ! elements # of x and per
      REAL(8),dimension(length):: x,xtos    ! input array, and its non-MISSING part
      REAL(8),dimension(nl)    :: per,oout  ! percentile level and results
      REAL(8):: bb,cc ! temporary variable
      integer :: nn,i,k,kp1  ! loop index and count #
      reaL(8) :: r

      ! check if percentile level is out of bound.
      if(maxval(per)>1..or.minval(per)<0.) stop 'ERROR in percentile: per out of range [0.,1.] !'

      ! choose non-MISSING input, but is it necessary? I think x is already chosen ......
      nn=0
      do i=1, length
        if(nomiss(x(i)))then
          nn=nn+1
          xtos(nn)=x(i)
        endif
      enddo
      if(nn.eq.0) then
        oout=MISSING
      else
        call sort(nn,xtos) ! ascending numerical order
        do i=1,nl
          r=per(i)*nn
          !bb=nn*per(i)+per(i)/3.+1./3.
          !cc=real(int(bb))
          !k=int(cc)
          !kp1=int(cc)+1
          k=floor(r+.5)
          kp1=k+1
          r=r-k;
          !if(k.ge.nn) then
           ! oout(i)=xtos(nn)
          !else
            !oout(i)=xtos(int(cc))+(bb-cc)*(xtos(int(cc)+1)-xtos(int(cc)))
            oout(i)=(.5-r)*xtos(k)+(.5+r)*xtos(kp1);
          !endif
        enddo
      endif

      return
      end subroutine percentile
!-----------------------------------------------------------------------
     !  Sort an array arr(1:n) into ascending numerical order using the Quicksort algorithm.
     !    n is element # of input array; arr is replace on output by its sorted rearrangement.
     !    Parameters: M is the size of subarrays sorted by straight insertion
     !    and NSTACK is the required auxiliary.
     !  from Numerical Recipes.
     !
      SUBROUTINE sort(n,arr)
      !use COMM !, only:DP
      INTEGER, PARAMETER :: M=7,NSTACK=50
      INTEGER  :: n, i,ir,j,jstack,k,l,istack(NSTACK)
      REAL(8) :: arr(n)
      REAL(8) :: a,temp

      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
         j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        !if(jstack.gt.NSTACK) then
         !   pause 'NSTACK too small in sort'
          !  endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE sort

      logical function ismiss(a)
      REAL(8) :: a
      REAL(8) :: rmiss
      rmiss=Missing+1.
      if(a.gt.rmiss) then
        ismiss=.FALSE.
      else
        ismiss=.TRUE.
      endif
      end function ismiss


      !  return TRUE if a is not MISSING
      !         FALSE if a is MISSING

      logical function nomiss(a)
      REAL(8) :: a
      REAL(8) :: rmiss
      rmiss=Missing+1.
      if(a.lt.rmiss) then
        nomiss=.FALSE.
      else
        nomiss=.TRUE.
      endif
      end function nomiss
!
! DONE WITH PERCENTILE FUNCTIONS
!-----------------------------------------------------------------------
END SUBROUTINE f_moments


!-----------------------------------------------------------------------
!FUNCTION funcv(g)
!! Scaled distance (vector of nmoments size). For each moment i
!!
!!       Mdis(i)= [Dmean(y(i+3)-y(i))-Mmean(y(i+3)-y(i))]
!!
!!       where:
!!       * Dmean
!!       * Mmean 
!!
!! -------------------------------------------------------------------------
!use nrtype
!use myParams
!use genericParams
!IMPLICIT NONE
!
!    REAL(SP), DIMENSION(Tlev_s), INTENT(IN)  :: g
!    REAL(DP), DIMENSION(Tsim_mr) ::funcv 
!
!    ! Stochastic processes parameters
!    real(dp), dimension(Tlev_s) :: mu_bar, mu_bar2, mu_eta1, mu_eta_t, mu_eta3, bc_x,mu_eta2
!    real(dp),dimension(Tsim_mr)::mu_eta_t3
!
!    ! parameters
!    real(dp) :: sig_eps, p1, p2, p3, sig_eta1, sig_eta2, sig_eta3, mu2, mu3, phi
!    INTEGER(i4b) :: i,i_m, i1, i2
!        
!    ! Parameters
!    sig_eps = storeTheta(1)
!    p1      = storeTheta(2)
!    p2      = 1./2. - p1/2.
!    p3      = p2
!    mu2     = storeTheta(3)
!    mu3     = storeTheta(4)
!    sig_eta1= storeTheta(5)
!    sig_eta2= storeTheta(6)
!    sig_eta3= sig_eta2
!    phi     = storeTheta(7)
!
!    ! Aggregate risk series (bcx=-GDP log diff 1 year)
!    bc_x = phi*xt
!
!    ! Create eta series
!    mu_bar = -log(p1*exp(sig_eta1**2./2.)+ &
!                  p2*exp(mu2-bc_x+(sig_eta2**2./2.))+ &
!                  p3*exp(mu3-bc_x+(sig_eta3**2./2.))) 
!    
!    mu_eta1 = mu_bar + g
!    mu_eta2 = mu_bar + mu2 - bc_x + g
!    mu_eta3 = mu_bar + mu3 - bc_x + g
!    mu_eta_t=p1*mu_eta1+p2*mu_eta2+p3*mu_eta3
!    mu_eta_t3=mu_eta_t(1:Tlev_s-3)+mu_eta_t(2:Tlev_s-2)+mu_eta_t(3:Tlev_s-1)
!    
!    forall (i=1:Tsim_mr) funcv(i)=mu_eta_t3(i)-mt(i)
!        funcv=dot_product(funcv,funcv)
!
!
!    END FUNCTION 

! Generate random sequences
!
subroutine gen_rand
use nrtype
use myParams
implicit none

real(dp), dimension(nhhsim*4,Tlev_s) :: normrand
real(dp), dimension(nhhsim,Tlev_s) :: unirand
!real(dp), dimension(nhhsim*5,Tsim,nsim) :: randva
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
integer :: i, j, k, nseed, seed0, m
!common randva

seed0=123456
call random_seed(size = nseed)
allocate(seed(nseed))
forall(m=1:nseed) seed(m) = (m-1)*100 + seed0
call random_seed(put = seed)
deallocate(seed)

do k=1,nsim        
    call random_number(unirand)
    randva(1:nhhsim,:,k)=unirand
    do i=1,nhhsim*4
            do j=1,Tlev_s
                    normrand(i,j)=random_normal()
            enddo
    enddo
    randva(nhhsim+1:nhhsim*5,:,k)=normrand

    open(unit=k+1000,file=savedir//'/randvar_sim'//achar(48+k)//'.dat',action='write')
        do i=1,nhhsim*5
                write(k+1000,*) randva(i,:,k)
        enddo
    close(k+1000)
enddo

contains

FUNCTION random_normal() RESULT(fn_val)
use nrtype
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL(dp) :: fn_val

!     Local variables
REAL(8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
REAL(8)      :: zero = 0.0, half = 0.5


!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

end subroutine 

! not sure about this one, Arun had it, and I don't want to mess up the algorithm
subroutine initial0
implicit none
call GetData
end subroutine

end MODULE objective

