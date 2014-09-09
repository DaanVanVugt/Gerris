! Sequential code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The main driver subroutine dagmg below provides a sequential implementation
! of the method presented in [1], where the used algorithms are described
! in detail.
! See the accompanying userguide for more details on how to use the software,
! and the README file for copyright notice and installation instructions.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!     Report GANMN 08-02
!
!-----------------------------------------------------------------------
      MODULE dagmg_mem
      SAVE
!---------------
! Parameters for subroutine dagmg that may be tuned by expert users
!
! INTEGER
!
!  maxlev is the maximal number of levels
!         (should be large enough - much larger than needed is armless)
!  real_len is the length of 1 REAL(kind(0.0d0)) in byte
!        (used only to display information on memory usage)
INTEGER, PARAMETER :: maxlev=40, real_len=8
!  nsmooth is the number of smoothing steps (pre + post);
!          each cycle includes (nsmooth+1)/2 pre-smoothing step(s)
!          and nsmooth/2 post-smoothing step(s)
!  nstep is the maximum number of coarsening steps
!        nstep<0 means that coarsening is stopped only according to
!        the matrix size, see parameter maxcoarsesize
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level)
INTEGER , PARAMETER :: nsmooth=2, nstep=-1, nlvcyc=0
!
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize,  it is factorized exactly and
!                 coarsening is stopped.
INTEGER , PARAMETER :: maxcoarsesize=400
!
! REAL
!
!  resi is the threshold t for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  checkdd is the threshold for strongly diagonally dominant rows,
!        see the first line of "Initialization" in Algorithm 2.1 in [1]
!  trswc is the Strong/Weak coupling threshold beta in Algorithms 2.1 & 2.2
!        in [1]
!  trspos is a threshold: if a row has a positive offidiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is kept
!         out of the aggregation
!  scalcg is the scaling factor applied to the coarse grid matrices
!
REAL (kind(0.0d0)) , PARAMETER :: resi=0.25
REAL (kind(0.0d0)) , PARAMETER :: checkdd=5.0, trswc=0.25, trspos=0.45, scalcg=1
!----------------
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!     Report GANMN 08-02
!---------------------------------------------
! Internal variables (do not modify)
!
TYPE InnerData
 REAL (kind(0.0d0)), DIMENSION(:), POINTER :: a
 INTEGER, DIMENSION(:), POINTER :: ja
 INTEGER, DIMENSION(:), POINTER :: ia
 REAL (kind(0.0d0)), DIMENSION(:), POINTER :: p
 INTEGER, DIMENSION(:), POINTER :: idiag
 INTEGER, DIMENSION(:), POINTER :: ind
 INTEGER, DIMENSION(:), POINTER :: inloc
 INTEGER, DIMENSION(:), POINTER :: ilstout
 INTEGER, DIMENSION(:), POINTER :: lstout
 INTEGER, DIMENSION(:), POINTER :: ilstin
 INTEGER, DIMENSION(:), POINTER :: lstin
END TYPE InnerData
!
        TYPE (InnerData) :: dt(maxlev)
        INTEGER :: nn(maxlev),kstat(3,maxlev)=0,innermax(maxlev)
        INTEGER :: nlev,nwrkcum,iout,nrst
        REAL (kind(0.0d0)) :: memi=0.0,memax=0.0,memr=0.0,mritr,rlenilen
        LOGICAL :: spd,wfo,wff,allzero
        INTEGER, PARAMETER :: IRANK=-9999
        END MODULE dagmg_mem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE dagmg(n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
          REAL (kind(0.0d0)) :: a(*),f(n),x(n)
          REAL (kind(0.0d0)) :: tol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The dimension of the matrix.
!
!  A       (input/output) REAL (kind(0.0d0)). Numerical values of the matrix.
!  IA      (input/output) INTEGER. Pointers for every row.
!  JA      (input/output) INTEGER. Column indices.
!
!              dagmg ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1). According what is written
!              above, dagmg assumes that some of these JA(K) (for
!              IA(I)<= K < IA(I+1) ) is equal to I with corresponding
!              A(K) carrying the value of the diagonal element,
!              which must be positive.
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order (The
!              matrix is mathematically the same, but stored in
!              different way).
!
!  F       (input/output) REAL (kind(0.0d0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB==0, 2, 3, 10 or 12
!
!  X       (input/output) REAL (kind(0.0d0)).
!          On input and if IJOB==10 or IJOB==12, initial guess
!             (for other values of IJOB, the default is used: the zero vector).
!          On output, the computed solution.
!
! IJOB     (input) INTEGER. Tells dagmg what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only (based on previous setup), no initial guess
!         12: solves only (based on previous setup), initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner on the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!   !!! IJOB==2,3,12 require that one has previously called dagmg with IJOB==1
!   !!!    It is mandatory to keep n, A, JA and IA unchanged between a call
!   !!!    to dagmg with IJOB==1 and all subsequent calls with IJOB== 2, 3 or 12
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!         !!  If NREST==1, Flexible CG is used instead of GCR
!             Should be used if and only if the matrix is
!             symmetric and positive definite.
!         !!  Significant only if IJOB == 0, 2, 3, 10 or 12.
!             If IJOB == 3, determines only the type of inner iterations
!             (FCG if NREST == 1, GCR otherwise).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB== 0, 2, 10 or 12.
!
!  TOL     (input) REAL (kind(0.0d0)).
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB== 0, 2, 10 or 12.
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Local variables
!
          INTEGER, SAVE :: nza
          REAL (kind(0.0d0)), SAVE :: cputm=0.0d0,eltm=0.0d0,flop=0.0d0,memh
          LOGICAL, SAVE :: preprocessed=.FALSE.,solve=.FALSE.
          INTEGER :: i,init
          REAL (kind(0.0d0)) :: cputmp,eltmp,resid
          INTEGER :: listrank,ifirstlistrank
!
          wfo=.TRUE.
          iout=iprint
          IF (iprint <= 0) THEN
             iout=6
             wfo=.FALSE.
          ELSE IF (iprint == 5) THEN
             iout=6
          END IF
!
          wff=wfo.AND.(IRANK<=0)
!
          IF (MOD(ijob,10) >= 2 .AND. .NOT.preprocessed) THEN
             WRITE (iout,1001) IRANK,ijob
             STOP
          END IF
!
          IF (ijob < 0 .AND. solve) GOTO 450
          IF (ijob < 0) GOTO 500
          IF (MOD(ijob,10) >= 2) GOTO 300
          IF (preprocessed) THEN
             CALL dagmg_relmem
             IF (.NOT.allzero) CALL dagmg_LAPACK(nn(nlev),a,ja,ia,f,-2,flop)
             preprocessed=.FALSE.
             solve=.FALSE.
             eltm=0.0d0
             cputm=0.0d0
             flop=0.0d0
             kstat=0
          END IF
          CALL dagmg_mestime(-1,0.0d0,0.0d0)
!
!       Initial settings
!
          nza=ia(n+1)-ia(1)
          IF (HUGE(n) > 1.0e10) THEN
             rlenilen=dble(8)/dble(real_len)
          ELSE
             rlenilen=dble(4)/dble(real_len)
          END IF
!
!
!-----  PREPROCESSING
!
          nlev=0
          IF (wfo) THEN
             WRITE(iout,900) IRANK
          END IF
          CALL dagmg_setup(1,n,a,ja,ia,listrank,ifirstlistrank)
          preprocessed=.TRUE.
          memh=memr+memi*rlenilen
!
          CALL dagmg_mestime(1,cputmp,eltmp)
          IF(wfo)THEN
               WRITE(iout,960) IRANK,memax/nza,real_len,&
                    &memax*real_len/(2**20)
               IF(MOD(ijob,10) == 1) THEN
                  WRITE(iout,961) IRANK,memh/nza,real_len,&
                       &memh*real_len/(2**20)
               END IF
!CPU_TIME: next line may be uncommented if implemented
!          (see subroutine mestime)
!              WRITE(iout,996) IRANK,cputmp
               WRITE(iout,997) IRANK,eltmp
               WRITE(iout,'()')
          END IF
!
          IF (MOD(ijob,10) == 1) RETURN
!
!-----  SOLUTION PROCESS
!
300       CONTINUE
          CALL dagmg_mestime(-2,0.0d0,0.0d0)
          resid=tol
          nrst=nrest
          init=MAX(IJOB/10,0)
          IF (nrst <= 0) nrst=10 ! Default restart paramater for GMRESR
          spd=.FALSE.  ! nrst=1 => system assumed SPD, use FCG instead of GCR
          IF (nrst == 1) spd=.TRUE.
   !
   !      DIRECT SOLVE IF ONLY ONE LEVEL
          IF (nlev == 1) THEN
   !
             x(1:n)=f(1:n)
             CALL dagmg_LAPACK(n,a,ja,ia,x,2,flop)
             mritr=0
   !
   !      END OF DIRECT SOLVE
   !
   !-     SINGLE APPLICATION OF MG PRECONDITIONER
          ELSE IF (MOD(ijob,10) >= 3) THEN
             IF (wfo) THEN
                WRITE(iout,901) IRANK
             END IF
             CALL dagmg_applyprec(n,f,x,a,ja,ia,flop)
             CALL dagmg_mestime(2,cputmp,eltmp)
             cputm=cputm+cputmp
             eltm=eltm+eltmp
             solve=.TRUE.
             RETURN
   !
   !-     END OF MG PREC
   !
   !-     GCR SOLUTION
          ELSE IF (nrst > 1) THEN
   !
             CALL dagmg_GCR(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !-     END OF GCR SOLVE
   !
   !-     FLEXIBLE CG SOLUTION
          ELSE
   !
             CALL dagmg_FlexCG(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !-     END OF FCG SOLVE
   !
          END IF
!
          CALL dagmg_mestime(2,cputm,eltm)
          solve=.FALSE.
          GOTO 455
450       CONTINUE
          IF (wfo) THEN
             WRITE(iout,'()')
             WRITE(iout,990) IRANK
             WRITE(iout,'()')
          END IF
455       CONTINUE
          IF (wfo) THEN
             IF (wff .AND. iter.NE.0) THEN
                DO i=2,nlev-1
                   WRITE(iout,955) i,kstat(3,i-1),kstat(3,i),         &
                        dble(kstat(3,i))/dble(kstat(3,i-1)),kstat(2,i)
                END DO
             END IF
             WRITE(iout,'()')
             WRITE(iout,954) IRANK,flop/dble(11*n+2*nza)
             WRITE(iout,962) IRANK,(memh+mritr)/nza,real_len, &
                  (memh+mritr)*real_len/(2**20)
!CPU_TIME: next line may be uncommented if implemented
!          (see subroutine mestime)
!            WRITE(iout,998) IRANK,cputm
             WRITE(iout,999) IRANK,eltm
             WRITE(iout,'()')
          END IF
!
          IF (MOD(ijob,10) > 0) RETURN
!
!----  RELEASE MEMORY
!
500       CONTINUE
          CALL dagmg_relmem
          IF (.NOT.allzero) CALL dagmg_LAPACK(n,a,ja,ia,f,-2,flop)
          preprocessed=.FALSE.
          solve=.FALSE.
          eltm=0.0d0
          cputm=0.0d0
          flop=0.0d0
          kstat=0
          IF (wfo) THEN
             WRITE (iout,902) IRANK
          END IF
!
          CONTINUE
!
          RETURN
900       FORMAT(i3,'*ENTERING AGMG **********************************',&
               '***************************')
901       FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902       FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
               '***************************')
954       FORMAT(i3,'*','  Equiv. number of CG iter.:',f9.2)
955       FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
                '   mean=',f7.2,'    max=',i3)
960       FORMAT(  i3,'*','         memory used (peak):',f9.2,        &
               ' real(',i1,') words per nnz (',f8.2,' Mb)')
961       FORMAT(  i3,'*','     memory still allocated:',f9.2,        &
               ' real(',i1,') words per nnz (',f8.2,' Mb)')
962       FORMAT(  i3,'*','   memory used for solution:',f9.2,        &
               ' real(',i1,') words per nnz (',f8.2,' Mb)')
990       FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996       FORMAT(i3,'*','           Setup time (CPU):   ',1PE10.2,     &
               ' seconds')
997       FORMAT(i3,'*','       Setup time (Elapsed):   ',1PE10.2,     &
               ' seconds')
998       FORMAT(i3,'*','        Solution time (CPU):   ',1PE10.2,     &
               ' seconds')
999       FORMAT(i3,'*','    Solution time (Elapsed):   ',1PE10.2,     &
               ' seconds')
1001      FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
               ' is not allowed')
        END SUBROUTINE dagmg
!-----------------------------------------------------------------------
        RECURSIVE SUBROUTINE dagmg_setup(l,n,a,ja,ia,listrank,ifl)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: l,n,ja(*),ia(n+1)
          REAL (kind(0.0d0)) :: a(*)
          INTEGER :: ifl,listrank(ifl:*)
          INTEGER :: nnext,ierr,nzan,nmax,i,j,k,llsto,lw,lwn,lw0,kw,nwprev
          INTEGER, POINTER, DIMENSION(:) :: iw,iw0,ja2,lcg,iwp
          REAL (kind(0.0d0)), ALLOCATABLE, DIMENSION(:)  :: a2
          INTEGER, SAVE    :: nlc(2),nlcp(2),nlc1(2),icum
          REAL (kind(0.0d0)), SAVE :: ngl(2),nglp(2),ngltot(2),nlctot(2),ngl1(2)
          LOGICAL, SAVE    :: slowcoarse=.FALSE.
          REAL (kind(0.0d0))       :: ops,ff,xsi,eta,checkddl,dum(2)
          CHARACTER(len=13) :: prtint
          REAL (kind(0.0d0)) :: fff
          INTEGER , parameter :: IONE=1
!
          nn(l)=n
          nlc(1)=n
          nlc(2)=ia(n+1)-ia(1)
!
          ngl=nlc
          IF ( l > 1 ) THEN
             nlctot=nlctot+nlc
             ngltot=ngltot+ngl
             IF (wfo) THEN
                IF (n > 9.9e10) THEN
                   WRITE(prtint(1:12),'(1pe12.5)') dble(n)
                ELSE
                   WRITE(prtint(1:12),'(i12)') n
                END IF
                WRITE(iout,920) prtint(1:12),dble(nlcp(1))/dble(n)
                IF (nlc(2) > 9.9e10) THEN
                   WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
                ELSE
                   WRITE(prtint(1:12),'(i12)') nlc(2)
                END IF
                WRITE(iout,921) prtint(1:12),dble(nlc(2))/dble(n),        &
                     dble(nlcp(2))/dble(nlc(2))
             END IF
             xsi=(6*1.0d0)/(10*1.0d0)
             eta=0.5d0
             ff=(ngl1(2)/ngl(2))*(xsi**(l-1))/dble(icum)
             IF (ff < 2.0d0-eta) THEN
                innermax(l)=1
             ELSE
                innermax(l)=2
             END IF
             icum=icum*innermax(l)
          ELSE
             nlctot=nlc
             ngltot=ngl
             ngl1=ngl
             nlc1=nlc
             icum=1
             innermax(1)=1
             slowcoarse=.FALSE.
             allzero=.FALSE.
             nglp(1)=0.0
             IF (wfo) THEN
                IF (n > 9.9e10) THEN
                   WRITE(prtint(1:12),'(1pe12.5)') dble(n)
                ELSE
                   WRITE(prtint(1:12),'(i12)') n
                END IF
                WRITE(iout,'()')
                WRITE(iout,918) prtint(1:12)
                IF (nlc(2) > 9.9e10) THEN
                   WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
                ELSE
                   WRITE(prtint(1:12),'(i12)') nlc(2)
                END IF
                WRITE(iout,919) prtint(1:12),dble(nlc(2))/dble(n)
                WRITE(iout,'()')
             END IF
          END IF
          IF( l == nstep+1  .OR. l == maxlev .OR.  ngl(1) <= maxcoarsesize &
               .OR. (slowcoarse .AND. 3*nglp(1) < 5*ngl(1))                &
               .OR. nglp(1) ==ngl(1) ) THEN
             nlev=l
             innermax(l)=0
          END IF
          slowcoarse=.FALSE.
          IF ( l > 1 .AND. 3*nglp(1) < 5*ngl(1) ) slowcoarse=.TRUE.
          nlcp=nlc
          nglp=ngl
          IF (n == 0) THEN
             nnext=0
             allzero=.TRUE.
             ! n=0 and nlev /= l possible with the parallel version
             IF (nlev /= l) THEN
                IF (wfo) THEN
                   WRITE(iout,914) IRANK,l+1
                END IF
                CALL dagmg_setup(l+1,nnext,a,ja,ia,listrank,ifl)
                RETURN
             END IF
             IF (allzero) GOTO 500 ! skip factorization on last level
             GOTO 300
          END IF
          lw=0
          lw0=0
          IF (l == nlev) GOTO 200
          lw=3*n+1
          lw0=n
          ALLOCATE(dt(l)%p(n),dt(l)%ind(n),dt(l)%idiag(n+1)           &
               ,lcg(4*n),a2(nlc(2)),ja2(nlc(2)),iw(lw),iw0(lw0))
          memi=memi+6*n+1+nlc(2)+lw+lw0
          memr=memr+n+nlc(2)
          memax=MAX(memax,memr+memi*rlenilen)
          iwp => iw
          kw=n+1
          IF (l>1) THEN
             checkddl=-1.0d0
          ELSE
             IF (wfo) THEN
                WRITE (iout,901) IRANK
             END IF
             IF (wff) THEN
                WRITE (iout,902) checkdd,trswc
                WRITE (iout,903) trspos
             END IF
             checkddl=checkdd
          END IF
          CALL dagmg_aggl4(n,nnext,a,ja,ia,lcg,dt(l)%ind,dt(l)%idiag,iw0, &
               iwp,iw(kw),dt(l)%p,trswc,a2,ja2,iw(n+kw),checkddl,trspos)
200       CONTINUE
          IF (lw > 0) DEALLOCATE(iw)
          memi=memi-lw
          lw=0
          IF (l == nlev) GOTO 300
          CALL dagmg_setsgs(n,a,ja,ia,dt(l)%p,dt(l)%idiag,a2,ja2,dt(l)%inloc)
          ALLOCATE(dt(l+1)%ia(nnext+1))
          memi=memi+nnext+1
          memax=MAX(memax,memr+memi*rlenilen)
          j=nnext
          IF (j > lw0) THEN
             IF (lw0 > 0) DEALLOCATE(iw0)
             ALLOCATE(iw0(j))
             memi=memi-lw0+j
             memax=MAX(memax,memr+memi*rlenilen)
             lw0=j
          END IF
          CALL dagmg_setacg(nnext,lcg,a,ja,ia,a2,ja2,dt(l+1)%ia,nzan,     &
               dt(l)%ind,iw0)
          DEALLOCATE(iw0,lcg)
          ALLOCATE(dt(l+1)%a(nzan),dt(l+1)%ja(nzan))
          memi=memi-lw0-4*n+nzan
          memr=memr+nzan
          memax=MAX(memax,memr+memi*rlenilen)
          dt(l+1)%a(1:nzan)=a2(1:nzan)
          dt(l+1)%ja(1:nzan)=ja2(1:nzan)
          DEALLOCATE(a2,ja2)
          memi=memi-nlc(2)
          memr=memr-nlc(2)
          IF (scalcg < 0.0d0) THEN
             ff=-scalcg*dble(nnext)/dble(n)
          ELSE
             ff=scalcg
          END IF
          IF (ff /= 1.0d0) THEN
             IF (wfo) THEN
                WRITE(iout,'()')
                WRITE(iout,915) IRANK,l+1,ff
             END IF
             CALL DSCAL(nzan,ff,dt(l+1)%a,IONE)
          ELSE
             IF (wfo) THEN
                WRITE(iout,'()')
                WRITE(iout,914) IRANK,l+1
             END IF
          END IF
          !
250       CONTINUE
          CALL dagmg_setup(l+1,nnext,dt(l+1)%a,dt(l+1)%ja,dt(l+1)%ia,   &
               listrank,nnext+1)
          GOTO 500
          !
300       CONTINUE
          !
          CALL dagmg_LAPACK(n,a,ja,ia,fff,1,ops)
          IF (wfo) THEN
             WRITE(iout,911) IRANK,ops/dble(11*nn(1)+2*nlc1(2))
             WRITE(iout,'()')
          END IF
          !
500       CONTINUE
          IF (l==1) THEN
             IF (wfo) THEN
                WRITE(iout,'()')
                WRITE(iout,954) nlctot(1)/dble(nlc1(1))
                WRITE(iout,955) nlctot(2)/dble(nlc1(2))
             END IF
             DO i=1,nlvcyc
                k=nlev-i
                IF (k > 2) innermax(k)=1
             END DO
             nwrkcum=0
             nwprev=0
             DO i=2,nlev-1
                IF (innermax(i) > 1)  THEN
                   nwrkcum=nwrkcum+MAX(nwprev,6*nn(i))
                   nwprev=nn(i)
                ELSE
                   nwrkcum=nwrkcum+MAX(nwprev,4*nn(i))
                   nwprev=0
                END IF
             END DO
             nwrkcum=nwrkcum+MAX(nwprev,2*nn(nlev))
          END IF
          !
          RETURN
901       FORMAT(i3,'*SETUP: Coarsening by double pairwise aggregation')
902       FORMAT('****       Strong diag. dom. trs:',f5.1,                 &
               ' ; Strong/Weak coupling trs:',f5.2)
903       FORMAT('****',18x,'Threshold for rows with large pos. offdiag.:',f5.2)
911       FORMAT(i3,'*','        Exact factorization:',f12.3,            &
               ' equiv. CG iterations')
914       FORMAT(i3,'*','                      Level:',I12)
915       FORMAT(i3,'*','                      Level:',I12,            &
               ' ;  applied scaling factor:',f6.3)
918       FORMAT('****','        Number of unknowns:', A12)
919       FORMAT('****','                 Nonzeros :', A12,                &
               ' (per row:',f7.2,')')
920       FORMAT('****','        Number of variables:',A12,              &
               '          (reduction ratio:',f5.2,')')
921       FORMAT('****','                   Nonzeros:',A12,              &
               ' (per row:',f4.1,    &
               '; red. ratio:',f5.2,')')
954       FORMAT('****',' relative complexity (grid):',f9.2)
955       FORMAT('****',' relative complexity (nnzs):',f9.2)
        END SUBROUTINE dagmg_setup
        !------------------------------------------------------------------
        SUBROUTINE dagmg_relmem
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: l
          DO l=1,nlev
             IF ( l>1 ) THEN
                IF ( nn(l-1)>0 ) THEN
                   DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%ia)
                END IF
             END IF
             IF (l<nlev .AND. nn(l)>0 )                                  &
                  DEALLOCATE(dt(l)%p,dt(l)%ind,dt(l)%idiag)
          END DO
          memi=0
          memr=0
          memax=0
          RETURN
        END SUBROUTINE dagmg_relmem
        !------------------------------------------------------------------
        SUBROUTINE dagmg_LAPACK(N,a,ja,ia,f,ijb,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: N,ia(n+1),ja(*),ijb
          REAL (kind(0.0d0)) :: a(*),f(n)
          REAL (kind(0.0d0)) :: flop
          !
          REAL (kind(0.0d0)), ALLOCATABLE, SAVE :: ac(:,:)
          INTEGER , ALLOCATABLE, SAVE :: ipiv(:)
          INTEGER, SAVE :: iflop
          INTEGER :: i,kk,inloc
          INTEGER  :: ierr
          INTEGER , parameter :: IONE=1
          !
          ierr=0
          IF (ijb == -2) THEN
             !
             DEALLOCATE (ac,ipiv)
             !
          ELSE IF (ijb == 1) THEN
             !
             ALLOCATE (ac(n,n),ipiv(n))
             memi=memi+n
             memr=memr+n*n
             memax=MAX(memax,memr+memi*rlenilen)
             ac=0.0d0
             DO i=1,n
                DO kk=ia(i),ia(i+1)-1
                   ac(i,ja(kk))=a(kk)
                END DO
             END DO
             CALL DGETRF(N,N,ac,N,ipiv,ierr)
             IF (ierr /= 0) THEN
                WRITE(iout, *) ' FATAL ERROR in GETRF: ierror=',ierr
                STOP
             END IF
             iflop=(2*n*n-n)
             flop=(2*1.0d0)/(3*1.0d0)*(dble(n)**3)
             !
          ELSE IF (ijb == 2) THEN
             !
             CALL DGETRS('N',N,IONE,ac,N,ipiv,f,N,ierr)
             IF (ierr /= 0) THEN
                WRITE(iout, *) ' FATAL ERROR in GETRS: ierror=',ierr
                STOP
             END IF
             flop=flop+dble(iflop)
             !
          END IF
          !
          RETURN
        END SUBROUTINE dagmg_LAPACK
        !-----------------------------------------------------------------------
        RECURSIVE SUBROUTINE dagmg_CGcorr(N,X,R,l,flop,w)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       ::  N, l, idum
          REAL (kind(0.0d0)) ::  flop
          REAL (kind(0.0d0)) ::  dum
          REAL (kind(0.0d0)) ::  X(N), R(N), w(*)
          INTEGER       ::  nnext, l1, RN, XN
          l1=l+1
          nnext=nn(l1)
          IF (nnext > 0) THEN
             XN=1
             RN=XN+nnext
             CALL dagmg_restaggl(N,nnext,R,w(RN),dt(l)%ind,flop)
             IF (l1 == nlev) THEN
                w(XN:XN+nnext-1)=w(RN:RN+nnext-1)
                CALL dagmg_LAPACK(nnext,dt(l1)%a,dt(l1)%ja,dt(l1)%ia,  &
                     w(XN),2,flop)
             ELSE IF ( innermax(l1) <= 1 ) THEN
                CALL dagmg_prec(nnext,w(XN),w(RN),dt(l1)%a,dt(l1)%ja,   &
                     dt(l1)%ia   &
                     ,l1,flop,w(RN+nnext))
                kstat(2,l+1)=MAX(kstat(2,l+1),1)
                kstat(3,l+1)=kstat(3,l+1)+1
                kstat(1,l+1)=kstat(1,l+1)+1
             ELSE IF (spd) THEN
                CALL dagmg_FlexCG_inner(nnext,w(XN),w(RN),l+1,flop)
             ELSE
                CALL dagmg_GCR_inner(nnext,w(XN),w(RN),l+1,flop)
             END IF
             CALL dagmg_prolaggl(N,nnext,X,w(XN),dt(l)%ind)
          ELSE
             X(1:N)=0.0d0
          END IF
          RETURN
        END SUBROUTINE dagmg_CGcorr
!-----------------------------------------------------------------------
        RECURSIVE SUBROUTINE dagmg_prec(N,X,B,a,ja,ia,l,flop,R)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: N, l
          REAL (kind(0.0d0)) ::  flop
          REAL (kind(0.0d0)) ::  B(N), X(N), R(N,*)
          INTEGER :: ja(*), ia(N+1)
          REAL (kind(0.0d0)) :: a(*)
          INTEGER :: is
          IF (nsmooth == 2) THEN
             CALL dagmg_sgsolve1(N,B,R(1,2),R,X,a,ja,ia,dt(l)%p,            &
                  dt(l)%idiag,flop, dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,&
                  dt(l)%lstin,dt(l)%ilstin)
             CALL dagmg_CGcorr(N,X,R,l,flop,R(1,3))
             CALL dagmg_sgsolve2(N,X,R(1,2),R,a,ja,ia,dt(l)%p,              &
                  dt(l)%idiag,flop,dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout, &
                  dt(l)%lstin,dt(l)%ilstin)
          ELSE
             X(1:N)=B(1:N)
             CALL dagmg_sgsolve(N, X, a, ja, ia, dt(l)%p, dt(l)%idiag, flop,&
                  dt(l)%inloc)
             DO is=2,(nsmooth+1)/2
                CALL dagmg_rescalc(N, X, R(1,2), B, a, ja, ia, flop,        &
                     dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,   &
                     dt(l)%ilstin)
                CALL dagmg_sgsolve(N,R(1,2),a, ja, ia, dt(l)%p, dt(l)%idiag,&
                     flop, dt(l)%inloc)
                X(1:N)=X(1:N)+R(1:N,2)
                flop=flop+dble(N)
             END DO
             CALL dagmg_rescalc(N, X, R, B, a, ja, ia, flop,                &
                  dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,      &
                  dt(l)%ilstin)
             CALL dagmg_CGcorr(N,R(1,2),R,l,flop,R(1,3))
             X(1:N)=X(1:N)+R(1:N,2)
             flop=flop+dble(N)
             DO is=1,nsmooth/2
                CALL dagmg_rescalc(N, X, R(1,2), B, a, ja, ia, flop,        &
                     dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,   &
                     dt(l)%ilstin)
                CALL dagmg_sgsolve(N,R(1,2),a, ja, ia, dt(l)%p, dt(l)%idiag,&
                     flop,dt(l)%inloc)
                X(1:N)=X(1:N)+R(1:N,2)
                flop=flop+dble(N)
             END DO
          END IF
          !
          RETURN
        END SUBROUTINE dagmg_prec
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_applyprec( N,f,X,a,ja,ia,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       :: N
          REAL (kind(0.0d0)) ::  flop
          REAL (kind(0.0d0)) ::  f(N), X(N)
          INTEGER       :: ja(*), ia(N+1)
          REAL (kind(0.0d0)) :: a(*)
          REAL (kind(0.0d0)), ALLOCATABLE :: S(:)
          !
          mritr=nwrkcum+2*N
          ALLOCATE (S(nwrkcum+2*N))
          CALL dagmg_prec(n,X,f,a,ja,ia,1,flop,S)
          kstat(3,1)=kstat(3,1)+1
          DEALLOCATE (S)
          !
          RETURN
        END SUBROUTINE dagmg_applyprec
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_FlexCG(N,f,X,ITER,RESID,a,ja,ia,init,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       :: N, ITER, init
          REAL (kind(0.0d0)) :: flop
          REAL (kind(0.0d0)) :: f(N), X(N)
          INTEGER       :: ja(*), ia(N+1)
          REAL (kind(0.0d0)) :: a(*)
          INTEGER       :: MAXIT, ierr
          REAL (kind(0.0d0)) :: TOL, BNORM, RESID, dum0
          REAL (kind(0.0d0)) :: ALPHA, BET0, RHO, dum(6)
          REAL (kind(0.0d0)), ALLOCATABLE :: S(:)
          REAL(kind(0.0d0)), external :: DDOT
          REAL(kind(0.0d0)), external :: DNRM2
          INTEGER , parameter :: IONE=1
          !
          mritr=nwrkcum+5*N
          ALLOCATE (S(n+1:nwrkcum+6*N))
          !
          flop=0.0d0
          kstat=0
          IF (wfo) THEN
             WRITE(iout,940) IRANK
          END IF
          IF (wff) THEN
             WRITE(iout,946) nsmooth
          END IF
          !
          TOL = RESID
          MAXIT = ITER
          RESID = 1.0d0
          ITER = 0
          dum(6) = DNRM2(N, f, IONE)**2
          flop=flop+dble(2*N)
          !
          !    compute initial residual
          !
          IF (init==1) THEN
             CALL dagmg_rescalc(N,x,f,f,a,ja,ia,flop,dt(1)%inloc,      &
                  dt(1)%lstout,dt(1)%ilstout,dt(1)%lstin,dt(1)%ilstin)
                  dum(5) = DNRM2(N, f, IONE)**2
             dum(2:3)=dum(5:6)
             BNORM=SQRT(dum(3))
             RESID=SQRT(dum(2))
             RESID=RESID/BNORM
             IF(wff.AND. (MAXIT <= 0 .OR. RESID <= TOL)) THEN
                WRITE(iout,900) 0, resid*bnorm, resid
             END IF
          END IF
          !
          !
          DO WHILE ( ITER < MAXIT .AND. RESID > TOL )
             ITER = ITER + 1
             !
             !        Preconditioner Solve.
             !
             CALL dagmg_prec(n,S(1+3*N),f,a,ja,ia,1,flop,S(1+4*N))
             !
             !        Compute direction vector.
             !
             IF ( ITER > 1 ) THEN
                dum(1) = - DDOT( N, S(1+N), IONE, S(1+3*N), IONE )
                BET0=dum(1)
                BET0=BET0/RHO
                CALL DAXPY( N, BET0, S(1+2*N), IONE, S(1+3*N), IONE )
                flop=flop+dble(4*N)
             ENDIF
             CALL DCOPY( N, S(1+3*N), IONE, S(1+2*N), IONE )
             CALL dagmg_matv(N,S(1+2*N),S(1+N),a,ja,ia,flop,dt(1)%inloc,   &
                  dt(1)%lstout,dt(1)%ilstout,dt(1)%lstin,dt(1)%ilstin)
             dum(4) =  DDOT( N, S(1+2*N), IONE, S(1+N), IONE )
             dum(5) =  DDOT(N,S(1+2*N),IONE,f,IONE)
             flop=flop+dble(4*N)
             IF (ITER==1) THEN
                dum(1:3)=dum(4:6)
                IF (init == 0) THEN
                   BNORM=SQRT(dum(3))
                END IF
                IF(wff) THEN
                  WRITE(iout,900) 0, resid*bnorm, resid
                END IF
             ELSE
                dum(1:2)=dum(4:5)
             END IF
             RHO=dum(1)
             ALPHA=dum(2)/RHO
             !
             IF (ITER == 1 .AND. init == 0) THEN
                CALL DCOPY(N,S(1+2*N),IONE,X,IONE)
                CALL DSCAL(N,ALPHA,X,IONE)
                flop=flop+dble(N)
             ELSE
                CALL DAXPY( N, ALPHA, S(1+2*N), IONE, X, IONE )
                flop=flop+dble(2*N)
             END IF
             CALL DAXPY( N, -ALPHA, S(1+N), IONE, f, IONE )
             dum0 = DNRM2(N,f,IONE)**2
             RESID=dum0
             RESID=SQRT(RESID)
             RESID=RESID/BNORM
             IF (wff) THEN
                WRITE(iout,900) iter, resid*bnorm, resid
             END IF
             flop=flop+dble(4*N)
          END DO
          !
          IF (wff) THEN
             WRITE(iout,952) iter
             WRITE(iout,'()')
          END IF
          !
          kstat(3,1)=ABS(iter)
          !
          DEALLOCATE(S)
          RETURN
900       FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
               '        Relat. res.=',e9.2)
940       FORMAT(i3,                                                     &
               '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
946       FORMAT(  '****          (with',i2,' SGS smoothing step per cycle)')
952       FORMAT('****  - Convergence reached in',I5,' iterations -')
          !
        END SUBROUTINE dagmg_FlexCG
!-----------------------------------------------------------------------
        RECURSIVE SUBROUTINE dagmg_FlexCG_inner( N,X,R,l,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       :: N, ITER, l, ierr
          REAL (kind(0.0d0)) ::  RESID, flop, BNORM, dum0
          REAL (kind(0.0d0)) ::  X(N), R(N,*)
          REAL (kind(0.0d0)) :: alpha1,alpha2,bet0,rho1,rho2,gamm0,dum(6)
          REAL(kind(0.0d0)), external :: DDOT
          REAL(kind(0.0d0)), external :: DNRM2
          INTEGER , parameter :: IONE=1
          !
          ! AT MOST 2 ITERATIONS
          !
          dum(6)=DNRM2(N, R, IONE)**2
          ITER = 1
          !
          !        Preconditioner Solve.
          !
          CALL dagmg_prec( N,X,R,dt(l)%a,dt(l)%ja,dt(l)%ia,l,flop,R(1,2))
          !
          CALL dagmg_matv(N,X,R(1,2),dt(l)%a,dt(l)%ja,dt(l)%ia,flop,      &
               dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,dt(l)%ilstin)
          !
          dum(4) =  DDOT(N,X,IONE,R(1,2),IONE)
          dum(5) =  DDOT(N,X,IONE,R,IONE)
          dum(1:3)=dum(4:6)
          BNORM=SQRT(dum(3))
          rho1=dum(1)
          alpha1=dum(2)
          !
          bet0=alpha1/rho1
          CALL DAXPY( N, -bet0, R(1,2), IONE, R, IONE )
          dum0=DNRM2(N,R,IONE)**2
          RESID=dum0
          RESID=SQRT(RESID)/BNORM
          IF (RESID <= resi) THEN
             CALL DSCAL(N,bet0,X,IONE)
             kstat(2,l)=MAX(kstat(2,l),iter)
             kstat(3,l)=kstat(3,l)+iter
             flop=flop+dble(11*N)
             RETURN
          END IF
          !
          ITER = 2
          !
          !         Preconditioner Solve.
          !
          CALL dagmg_prec(N,R(1,3),R,dt(l)%a,dt(l)%ja,dt(l)%ia,l,flop,R(1,4))
          !
          CALL dagmg_matv(N,R(1,3),R(1,4),dt(l)%a,dt(l)%ja,dt(l)%ia, flop,&
               dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,       &
               dt(l)%ilstin)
          !
          dum(4) = DDOT(N,R(1,2),IONE,R(1,3),IONE)
          dum(5) = DDOT(N,R(1,3),IONE,R,IONE)
          dum(6) = DDOT(N,R(1,3),IONE,R(1,4),IONE)
          dum(1:3)=dum(4:6)
          gamm0=dum(1)
          alpha2=dum(2)
          rho2 =dum(3)
          !
          rho2=rho2-gamm0*gamm0/rho1
          !
          bet0=(alpha1-alpha2*gamm0/rho2)/rho1
          CALL DSCAL( N, bet0, X, IONE )
          bet0=alpha2/rho2
          CALL DAXPY( N, bet0, R(1,3), IONE, X, IONE )
          RESID=RESID*RESID         ! crude estimation for statistic only
          flop=flop+dble(19*N)
          !
          kstat(2,l)=MAX(kstat(2,l),iter)
          kstat(3,l)=kstat(3,l)+iter
          IF ( RESID > resi )  kstat(1,l)=kstat(1,l)+1
          !
          RETURN
          !
        END SUBROUTINE dagmg_FlexCG_inner
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_GCR(N,f,X,ITER,RESID,a,ja,ia,init,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       :: N, ITER, init
          REAL (kind(0.0d0)) ::  flop
          REAL (kind(0.0d0)) ::  f(N), X(N)
          INTEGER       :: ja(*), ia(N+1)
          REAL (kind(0.0d0)) :: a(*)
          INTEGER   :: MAXIT,i,itm,kc,ku,kfree,iv,iR,ifree,irst,ierr
          REAL (kind(0.0d0)) ::  ALPHA, BET0
          REAL (kind(0.0d0)) ::  RESID, RESID2, BNORM2,  RHO, TRS
          REAL (kind(0.0d0)) ::  TOL, TOL2BNORM2, TOLT, dum0
          REAL (kind(0.0d0)) :: dum(6)
          REAL (kind(0.0d0)), ALLOCATABLE :: S(:)
          REAL(kind(0.0d0)), external :: DDOT
          REAL(kind(0.0d0)), external :: DNRM2
          INTEGER , parameter :: IONE=1
          INTEGER  :: itm1, m, info
          !
          mritr=nwrkcum+N*(2*nrst+2)+((nrst+1)*nrst)/2+nrst
          ALLOCATE (S(n+1:nwrkcum+N*(2*nrst+3)+((nrst+1)*nrst)/2+nrst))
          !
          flop=0.0d0
          kstat=0
          IF (wfo) THEN
             WRITE(iout,938) IRANK,nrst
          END IF
          IF (wff) THEN
             WRITE(iout,946) nsmooth
          END IF
          !
          m=MIN(nrst,ITER,N)
          kc=1
          ku=kc+m
          kfree=ku+m
          iR=0
          iv=iR+(m*(m+1))/2
          ifree=iv+m
          TRS=EPSILON(1.0d0)
          TRS=SQRT(SQRT(TRS))
          TOL = RESID
          TOL2BNORM2 = TOL
          MAXIT = ITER
          RESID2= 1.0d0
          ITER = 0
          dum(6) = DNRM2(N, f, IONE)**2
          flop=flop+dble(2*N)
          !
          !    compute initial residual
          !
          IF (init==1) THEN
             CALL dagmg_rescalc(N,x,f,f,a,ja,ia,flop,dt(1)%inloc,      &
                  dt(1)%lstout,dt(1)%ilstout,dt(1)%lstin,dt(1)%ilstin)
                  dum(5) = DNRM2(N, f, IONE)**2
             dum(2:3)=dum(5:6)
             BNORM2=dum(3)
             RESID2=dum(2)
             TOL2BNORM2 = TOL*TOL*BNORM2
             IF (wff.AND. (MAXIT <= 0 .OR. RESID2 <= TOL2BNORM2)) THEN
                WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
             END IF
          END IF
          !
          itm  = -1
          irst = 0
          DO WHILE ( ITER < MAXIT .AND. RESID2 > TOL2BNORM2 )
             itm  = itm  + 1
             ITER = ITER + 1
             !
             !      Restarting
             IF (itm == m) THEN
                CALL DTPTRS('U','N','U',m, IONE,S(1+iR+kfree*N),          &
                     S(1+iv+kfree*N),m,info)
                IF (irst == 0 .AND. init == 0) THEN
                   CALL DGEMV('N',N,m,1.0d0,S(1+ku*N), &
                        N,S(1+iv+kfree*N),IONE,0.0d0,&
                        X,IONE)
                   flop=flop+dble(2*m*N+m*(m+1))
                ELSE
                   CALL DGEMV('N',N,m,1.0d0,S(1+ku*N), &
                        N,S(1+iv+kfree*N),IONE,0.0d0,&
                        S(1+kc*N),IONE)
                   DO i=1,N
                      X(i)=X(i)+S(1+i-1+kc*N)
                   END DO
                   flop=flop+dble((2*m+1)*N+m*(m+1))
                END IF
                itm=0
                irst=irst+1
             END IF
             !
             !        Preconditioner Solve & MATVEC
             !
             CALL dagmg_prec(N,S(1+(ku+itm)*N),f,a,ja,ia,1,flop,        &
                  S(1+ifree+kfree*N))
             CALL dagmg_matv(N, S(1+(ku+itm)*N), S(1+(kc+itm)*N),       &
                  a, ja, ia, flop,dt(1)%inloc,dt(1)%lstout,            &
                  dt(1)%ilstout,dt(1)%lstin,dt(1)%ilstin)
             !
             !        Gram-Schmidt
             !
             IF (itm > 0) THEN
                DO i=0,itm-1
                   dum(1)=DDOT(N,S(1+(kc+i)*N),IONE,S(1+(kc+itm)*N),IONE)
                   bet0=dum(1)
                   bet0=bet0/S(1+iR+i+(i*(i+1))/2+kfree*N)
                   S(1+iR+i+(itm*(itm+1))/2+kfree*N)=bet0
                   CALL DAXPY(N,-bet0,S(1+(kc+i)*N),IONE,      &
                        S(1+(kc+itm)*N),IONE)
                   flop=flop+dble(4*N)
                END DO
             END IF
             !
             !           no normalisation: record norm instead
             !
             dum(4)=DNRM2(N,S(1+(kc+itm)*N),IONE)**2
             dum(5)=DDOT( N, S(1+(kc+itm)*N), IONE, f, IONE )
             IF (ITER == 1) THEN
                dum(1:3)=dum(4:6)
                IF (init == 0) THEN
                    BNORM2=dum(3)
                    RESID2=BNORM2
                    TOL2BNORM2=TOL*TOL*BNORM2
                END IF
                IF (wff) THEN
                   WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
                END IF
                TOLT = MAX(TOL2BNORM2,TRS*RESID2)
             ELSE
                dum(1:2)=dum(4:5)
             END IF
             rho=dum(1)
             alpha=dum(2)
             S(1+iR+itm+(itm*(itm+1))/2+kfree*N)=rho
             bet0=alpha/rho
             S(1+iv+itm+kfree*N)=bet0
             !
             CALL DAXPY( N, -bet0, S(1+(kc+itm)*N), IONE, f, IONE )
             flop=flop+dble(6*N)
             !
             RESID2 = RESID2 - alpha*alpha/rho
             IF (RESID2 <= TOLT) THEN
                dum0 = DNRM2(N,f,IONE)**2
                RESID2=dum0
                flop=flop+dble(2*N)
                TOLT = MAX(TOL2BNORM2,TRS*RESID2)
             END IF
             IF (wff)THEN
                WRITE(iout,900) iter,irst,SQRT(resid2),SQRT(resid2/bnorm2)
             END IF
             !
          END DO
          !
          IF (itm >= 0) THEN
             itm1=itm+1
             CALL DTPTRS('U','N','U',itm1, IONE,S(1+iR+kfree*N),     &
                  S(1+iv+kfree*N),m,info)
             IF (irst == 0 .AND. init == 0) THEN
                CALL DGEMV('N',N,itm1,1.0d0,S(1+ku*N), &
                     N,S(1+iv+kfree*N),IONE,0.0d0,   &
                     X,IONE)
                flop=flop+dble(2*(itm+1)*N+(itm+1)*(itm+2))
             ELSE
                CALL DGEMV('N',N,itm1,1.0d0,S(1+ku*N), &
                     N,S(1+iv+kfree*N),IONE,0.0d0,   &
                     S(1+kc*N),IONE)
                DO i=1,N
                   X(i)=X(i)+S(1+i-1+kc*N)
                END DO
                flop=flop+dble((2*(itm+1)+1)*N+(itm+1)*(itm+2))
             END IF
          END IF
          !
          RESID=SQRT(RESID2/BNORM2)
          IF (wff) THEN
             WRITE(iout,952) iter
             WRITE(iout,'()')
          END IF
          !
          kstat(3,1)=ABS(iter)
          !
          DEALLOCATE (S)
          !
          RETURN
900       FORMAT('****  Iter=',i5,' (',i2,' rest.)        Resid=',e9.2,    &
               '        Relat. res.=',e9.2)
938       FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
946       FORMAT(  '****          (with',i2,' SGS smoothing step per cycle)')
952       FORMAT('****  - Convergence reached in',I5,' iterations -')
          !
        END SUBROUTINE dagmg_GCR
!-----------------------------------------------------------------------
        RECURSIVE SUBROUTINE dagmg_GCR_inner( N,X,R,l,flop)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER       :: N, ITER, l, ierr
          REAL (kind(0.0d0)) ::  RESID,BNORM,flop,rho1,rho2
          REAL (kind(0.0d0)) ::  X(N), R(N,*)
          REAL (kind(0.0d0)) :: alpha1,alpha2,bet0,gamm0,dum(6)
          REAL(kind(0.0d0)), external :: DDOT
          REAL(kind(0.0d0)), external :: DNRM2
          INTEGER , parameter :: IONE=1
          !
          ! AT MOST 2 ITERATIONS
          !
          dum(6)=DNRM2(N, R, IONE)**2
          ITER = 1
          !
          !        Preconditioner Solve & MATVEC
          !
          CALL dagmg_prec(N,X,R,dt(l)%a,dt(l)%ja,dt(l)%ia,l,flop,R(1,2))
          CALL dagmg_matv(N, X, R(1,2),dt(l)%a,dt(l)%ja,dt(l)%ia,flop,    &
               dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,dt(l)%ilstin)
          !
          dum(4)=DNRM2(N,R(1,2),IONE)**2
          dum(5)=DDOT(N, R(1,2), IONE, R, IONE )
          dum(1:3)=dum(4:6)
          BNORM=dum(3)
          rho1=dum(1)
          alpha1=dum(2)
          !
          bet0=alpha1/rho1
          RESID=BNORM-alpha1*bet0
          !
          IF (RESID <= resi*resi*BNORM) THEN
             CALL DSCAL( N, bet0, X, IONE )
             flop=flop+dble(7*N)
             kstat(2,l)=MAX(kstat(2,l),iter)
             kstat(3,l)=kstat(3,l)+iter
             RESID=SQRT(MAX(0.0d0,RESID/BNORM))
             RETURN
          END IF
          !
          CALL DAXPY( N, -bet0, R(1,2), IONE, R, IONE )
          !
          ITER = 2
          !
          !           Preconditioner Solve & MATVEC
          !
          CALL dagmg_prec(N,R(1,3),R,dt(l)%a,dt(l)%ja,dt(l)%ia,l,flop,R(1,4))
          CALL dagmg_matv(N,R(1,3),R(1,4),dt(l)%a,dt(l)%ja,dt(l)%ia,flop, &
               dt(l)%inloc,dt(l)%lstout,dt(l)%ilstout,dt(l)%lstin,       &
               dt(l)%ilstin)
          !
          dum(4) = DDOT(N,R(1,2),IONE,R(1,4),IONE)
          dum(5) = DDOT(N,R(1,4),IONE,R,IONE)
          dum(6) = DNRM2(N,R(1,4),IONE)**2
          dum(1:3)=dum(4:6)
          gamm0=dum(1)
          alpha2=dum(2)
          rho2 =dum(3)
          !
          rho2=rho2-gamm0*gamm0/rho1
          bet0=(alpha1-alpha2*gamm0/rho2)/rho1
          CALL DSCAL( N, bet0, X, IONE )
          bet0=alpha2/rho2
          CALL DAXPY( N, bet0, R(1,3), IONE, X, IONE )
          RESID=RESID-alpha2*bet0
          RESID=SQRT(MAX(0.0d0,RESID/BNORM))
          !
          flop=flop+dble(17*N)
          !
          kstat(2,l)=MAX(kstat(2,l),iter)
          kstat(3,l)=kstat(3,l)+iter
          IF ( RESID > resi )  kstat(1,l)=kstat(1,l)+1
          !
          RETURN
          !
        END SUBROUTINE dagmg_GCR_inner
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_aggl2 (n, nc, a, ja, ia, lcg, ind, lpair, deg  &
             , next, prev, first, odmax, trs, checkdd, trspos, maxdg)
          IMPLICIT NONE
          INTEGER :: n, nc, ja(*), ia(n+1), lcg(4*n), ind(n), lpair(n), deg(n)
          INTEGER :: next(n), prev(n), first(0:n-1), maxdg
          REAL(kind(0.0d0)) :: a(*), odmax(n), vald
          REAL(kind(0.0d0)) :: checkdd, trspos
          REAL(kind(0.0d0)) :: trs, val, tent, odm, valp
          INTEGER :: mindg, i, j, jj, k, ipair, dg, isel, nmark, nm1, nm2, kd
          nmark = 0
          nm1 = - n - 1
          nm2 = - n - 2
          maxdg = 0
          nc = 0
          IF (checkdd.GT.0.0d0) THEN
             DO i = 1, n
                j = ia (i)
                jj = ia (i + 1) - 1
                dg = jj - j
                val = 0.0d0
                odm = 0.0d0
                valp= 0.0d0
                DO k = j, jj
                   IF (ja (k) .EQ.i) THEN
                      vald = a(k)
                   ELSE
                      odm = MAX(odm,ABS(a(k)))
                      valp= MAX(valp,a(k))
                      val = val + ABS(a(k))
                   ENDIF
                ENDDO
                IF (dg.EQ.0 .OR. ABS(vald).GT.checkdd*val) THEN
                   ind (i) = nm1
                   nmark = nmark + 1
                   odmax (i) = 0.0d0
                   deg (i) = nm1
                ELSE IF (valp.GT.trspos*vald) THEN
                   ind(i) = nm2
                   deg (i) = 0
                   odmax (i) = - trs * odm
                ELSE
                   ind (i) = 0
                   deg (i) = 0
                   odmax (i) = - trs * odm
                ENDIF
             ENDDO
          ELSE
             DO i = 1, n
                j = ia (i)
                jj = ia (i + 1) - 1
                dg = jj - j
                odm = 0.0d0
                valp= 0.0d0
                DO k = j, jj
                   IF (ja (k) .EQ.i) THEN
                      vald = a (k)
                   ELSE
                      odm = MAX(odm,ABS(a(k)))
                      valp= MAX(valp,a(k))
                   ENDIF
                ENDDO
                IF (dg.EQ.0) THEN
                   ind (i) = nm1
                   nmark = nmark + 1
                   odmax (i) = 0.0d0
                   deg (i) = nm1
                ELSE IF (valp.GT.trspos*vald) THEN
                   ind(i) = nm2
                   deg (i) = 0
                ELSE
                   ind (i) = 0
                   deg (i) = 0
                   odmax (i) = - trs * odm
                ENDIF
             ENDDO
          ENDIF
          DO i = 1, n
             DO k = ia (i), ia (i + 1) - 1
                IF ( (ind (i).EQ.0 .OR. ind(i).EQ.nm2)            &
                   .AND. a(k).LT.odmax(i)             &
                   .AND. ja (k).NE.i  ) THEN
                   dg = deg (ja (k) ) + 1
                   deg (ja (k) ) = dg
                   maxdg = MAX(maxdg,dg)
                ENDIF
             ENDDO
          ENDDO
          DO i = 0, maxdg
             first (i) = 0
          ENDDO
          DO i = n, 1, - 1
             IF (ind(i).EQ.0 .OR. ind(i).EQ.nm2) THEN
                dg = deg(i)
                IF (first(dg) .GT. 0) prev ( first(dg) ) = i
                next (i) = first (dg)
                prev (i) = 0
                first (dg) = i
             ENDIF
          ENDDO
          DO WHILE (nmark.LT.n)
             mindg = - 1
             ipair = 0
             DO WHILE (ipair.EQ.0)
                mindg = mindg + 1
                !           if (mindg .gt. maxdg) stop 'isel'
                IF (first (mindg) .GT.0) THEN
                   ipair = first (mindg)
                   first (mindg) = next (ipair)
                   IF (next (ipair) .GT.0) prev (next (ipair) ) = 0
                ENDIF
             ENDDO
             nc = nc + 1
             ind (ipair) = nc
             isel = 0
             IF (ind(ipair) .EQ. nm2) GOTO 20
             val = 0.0d0
             DO i = ia (ipair), ia (ipair + 1) - 1
                j = ja (i)
                   IF (ind(j).EQ.0 .AND. a(i).LT.odmax(ipair)) THEN
                      tent = a(i)
                      IF (tent.LT.1.0001*val) THEN
                         isel = j
                         val = tent
                      ENDIF
                   ENDIF
             ENDDO
 20          CONTINUE
             IF (isel.EQ.0) THEN
                lcg (nc) = ipair
                lpair (nc) = 0
                nmark = nmark + 1
             ELSE
                ind (ipair) = - nc
                ind (isel) = nc
                lcg (nc) = isel
                nmark = nmark + 2
                lpair (nc) = ipair
                dg = deg(isel)
                IF (prev (isel) .GT.0) THEN
                   next (prev (isel) ) = next (isel)
                ELSE
                   first (dg) = next (isel)
                ENDIF
                IF (next (isel) .GT.0) prev (next (isel) ) = prev (isel)
                DO i = ia (isel), ia (isel + 1) - 1
                   j = ja (i)
                      IF ( (ind (j).EQ.0 .OR. ind(j).EQ.nm2)            &
                           .AND. a(i).LT.odmax(isel) ) THEN
                         dg = deg (j)
                         dg = dg - 1
                         deg (j) = dg
!                         IF (dg .LT. maxdeg) THEN
                            IF (prev (j) .GT.0) THEN
                              next (prev (j) ) = next (j)
                            ELSE
                               first (dg+1) = next (j)
                            ENDIF
                            IF (next (j) .GT.0) prev (next (j) ) = prev (j)
                            IF (first (dg) .GT.0) prev (first (dg) ) = j
                            next (j) = first (dg)
                            prev (j) = 0
                            first (dg) = j
!                         ENDIF
                      END IF
                ENDDO
             ENDIF
             DO i = ia (ipair), ia (ipair + 1) - 1
                j = ja (i)
                   IF ( (ind (j).EQ.0 .OR. ind(j).EQ.nm2)            &
                        .AND. a(i).LT.odmax(ipair) ) THEN
                      dg = deg (j)
                      dg = dg - 1
                      deg (j) = dg
!                      IF (dg .LT. maxdeg) THEN
                         IF (prev (j) .GT.0) THEN
                            next (prev (j) ) = next (j)
                         ELSE
                            first (dg+1) = next (j)
                         ENDIF
                         IF (next (j) .GT.0) prev (next (j) ) = prev (j)
                         IF (first (dg) .GT.0) prev (first (dg) ) = j
                         next (j) = first (dg)
                         prev (j) = 0
                         first (dg) = j
!                      ENDIF
                   ENDIF
             ENDDO
          ENDDO
          RETURN
        END SUBROUTINE dagmg_aggl2
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_aggl4 (n, nc, a, ja, ia, lcg, ind, lcg1,      &
             lpair1, lcg2, lpair2, odmax, trs, a2, ja2, ia2, checkdd, trspos)
          IMPLICIT NONE
          INTEGER :: n, nc, ja(*), ia(n+1), lcg(4*n), ind(n)
          INTEGER :: lpair1(n), lcg1 (n), lcg2(n), lpair2(n)
          INTEGER :: ia2(n+1), ja2(*)
          REAL(kind(0.0d0)) :: a(*), odmax (n), a2(*), vald, wextmx, wextmn
          REAL(kind(0.0d0)) :: trs, checkdd, trspos, wextmax, wextmin
          INTEGER :: nc1, i, jj, jc, jcol, kb, jpos, nz, maxdg
          INTEGER, POINTER, DIMENSION(:) :: ifirst
          CALL dagmg_aggl2 (n, nc1, a, ja, ia, lcg1, ind, lpair1, lcg2,  &
               lpair2, ia2, lcg, odmax, trs, checkdd, trspos, maxdg)
          nz = 0
          ia2 (1) = 1
          DO i = 1, nc1
             lcg (i) = 0
          ENDDO
          DO i = 1, nc1
             jc=nc1+1
             jj = lcg1 (i)
             DO kb = ia (jj), ia (jj+1) - 1
                jcol = ja (kb)
                   jcol = ABS (ind (jcol) )
!                   IF (jcol.NE.i.AND.jcol.LE.nc1) THEN
                   IF (jcol.LE.nc1) THEN
                      jpos = lcg (jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz + 1
                         ja2 (nz) = jcol
                         lcg (jcol) = nz
                         a2 (nz) = a (kb)
                      ELSE
                         a2 (jpos) = a2 (jpos) + a (kb)
                      ENDIF
                   ENDIF
             ENDDO
             jj = lpair1 (i)
             IF (jj.GT.0) THEN
                DO kb = ia (jj), ia (jj + 1) - 1
                   jcol = ja (kb)
                      jcol = ABS (ind (jcol) )
!                      IF (jcol.NE.i.AND.jcol.LE.nc1) THEN
                      IF (jcol.LE.nc1) THEN
                         jpos = lcg (jcol)
                         IF (jpos.EQ.0) THEN
                            nz = nz + 1
                            ja2 (nz) = jcol
                            lcg (jcol) = nz
                            a2 (nz) = a (kb)
                         ELSE
                            a2 (jpos) = a2 (jpos) + a (kb)
                         ENDIF
                      ENDIF
                ENDDO
             ENDIF
             DO kb = ia2 (i), nz
                lcg (ja2 (kb) ) = 0
             ENDDO
             ia2 (i + 1) = nz + 1
          ENDDO
          IF (maxdg+1 .LE. ia(n+1)-nz-1) THEN
             CALL dagmg_aggl2 (nc1, nc, a2, ja2, ia2, lcg2, lcg, lpair2,&
                  lcg(n+1), lcg(2*n+1), lcg(3*n+1), ja2(nz+1),         &
                  odmax, trs, -1.0d0, trspos, maxdg)
          ELSE
             ALLOCATE(ifirst(maxdg+1))
             CALL dagmg_aggl2 (nc1, nc, a2, ja2, ia2, lcg2, lcg, lpair2,&
               lcg(n+1), lcg(2*n+1), lcg(3*n+1), ifirst,               &
               odmax, trs, -1.0d0, trspos, maxdg)
             DEALLOCATE(ifirst)
         END IF
          DO i = 1, n
             jc = ind (i)
             jcol = ABS (jc)
             IF (jcol.LE.nc1) THEN
                IF (jc.GT.0) THEN
                   ind (i) = lcg (jcol)
                ELSE
                   ind (i) = - ABS (lcg (jcol) )
                ENDIF
             ENDIF
          ENDDO
          CALL dagmg_setlcg (nc, lcg, lcg1, lpair1, lcg2, lpair2)
          RETURN
        END SUBROUTINE dagmg_aggl4
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_setlcg (nc, lcg, lcg1, lpair1, lcg2, lpair2)
          IMPLICIT NONE
          INTEGER nc, lcg(4,nc), i, j1, j2
          INTEGER lpair1(*), lcg1(*), lcg2(*), lpair2 (*)
          DO i = 1, nc
             j1 = lcg2 (i)
             j2 = lpair2 (i)
             lcg (1, i) = lcg1 (j1)
             lcg (2, i) = lpair1 (j1)
             IF (j2.GT.0) THEN
                lcg (3, i) = lcg1 (j2)
                lcg (4, i) = lpair1 (j2)
             ELSE
                lcg (3, i) = 0
                lcg (4, i) = 0
             ENDIF
          ENDDO
          RETURN
        END SUBROUTINE dagmg_setlcg
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_setsgs (n, a, ja, ia, p, idiag, w, iw, inloc)
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), idiag(n), iw(*), inloc(n)
          REAL(kind(0.0d0)) :: a(*), p(n), w(*), t
          INTEGER :: i, j, k, ipos, nzu
          DO i = 1, n
             !       find diag, save upper part, concatenate lower part
             ipos = ia (i)
             nzu = 0
             DO k = ia (i), ia (i + 1) - 1
                j = ja (k)
                IF (j.GT.i) THEN
                   nzu = nzu + 1
                   w (nzu) = a (k)
                   iw (nzu) = j
                ELSEIF (j.LT.i) THEN
                   a (ipos) = a (k)
                   ja (ipos) = j
                   ipos = ipos + 1
                ELSE
                   p (i) = a (k)
                ENDIF
             ENDDO
             !
             !       copy back diagonal entry
             idiag (i) = ipos
             a (ipos) = p (i)
             ja (ipos) = i
             !
             !       copy back upper part
             DO k = 1, nzu
                ipos = ipos + 1
                a (ipos) = w (k)
                ja (ipos) = iw (k)
             ENDDO
             !
             !      save inverse
             t=p(i)
             p (i) = 1.0d0 / t
          ENDDO
          RETURN
        END SUBROUTINE dagmg_setsgs
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_matv (n, x, y, a, ja, ia, flop,                &
             inloc, lstout, ilstout, lstin, ilstin )
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier
          INTEGER :: inloc(n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
          REAL(kind(0.0d0)) :: x(n), y(n), a(*), t
          REAL(kind(0.0d0)) :: flop
          DO i = 1, n
             k1 = ia (i)
             t = a (k1) * x (ja (k1) )
             k2 = ia (i+1)
             DO kk = k1 + 1, k2 - 1
                t = t + a (kk) * x (ja (kk) )
             ENDDO
             y (i) = t
          ENDDO
          flop = flop + dble(2 * (ia (n + 1) - 1) - n)
          RETURN
        END SUBROUTINE dagmg_matv
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_sgsolve1(n, b, t1, r, v, a, ja, ia, p, idiag, flop,&
             inloc, lstout, ilstout, lstin, ilstin )
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), idiag(n), kk, j, k2, i, ier
          INTEGER :: inloc (n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
          REAL(kind(0.0d0)) :: b(n), t1(n), r(n), v(n)
          REAL(kind(0.0d0)) :: a(*), p(n), t
          REAL(kind(0.0d0)) :: flop
          !a
          t1 (1) = p (1) * b (1)
          v (1) = b (1)
          DO kk = 2, n
             t = b (kk)
             DO j = ia (kk), idiag (kk) - 1
                t = t - a (j) * t1 (ja (j) )
             ENDDO
             v (kk) = t
             t1 (kk) = p (kk) * t
          ENDDO
          !b
          DO kk = n - 1, 1, - 1
             t = 0.0d0
             k2 = ia (kk+1)
             DO j = idiag (kk) + 1, k2 - 1
                t = t - a (j) * t1 (ja (j) )
             ENDDO
             t1 (kk) = t1 (kk) + p (kk) * t
          ENDDO
          !c
          !...
          r (1) = b (1) - v (1)
          DO kk = 2, n
             t = b (kk) - v (kk)
             DO j = ia (kk), idiag (kk) - 1
                t = t - a (j) * t1 (ja (j) )
             ENDDO
             r (kk) = t
          ENDDO
          !...
          flop = flop + dble(3 * (ia (n + 1) - 1) )
          RETURN
        END SUBROUTINE dagmg_sgsolve1
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_sgsolve2(n, x, t1, t3, a, ja, ia, p, idiag, flop, &
             inloc, lstout, ilstout, lstin, ilstin )
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), idiag(n), kk, j, i, k2, ier
          INTEGER :: inloc (n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
          REAL(kind(0.0d0)) :: x(n), t1(n), t3(n), a(*), p(n), t
          REAL(kind(0.0d0)) :: flop
          ! input : x=t2 ; converted to t1+t2
          ! output: x=t1+t1+t3
          !a
          x (n) = x (n) + t1 (n)
          t3 (n) = 0.0d0
          DO kk = n - 1, 1, - 1
             x (kk) = x (kk) + t1 (kk)
             t = 0.0d0
             DO j = idiag (kk) + 1, ia (kk + 1) - 1
                t = t + a (j) * x (ja (j) )
             ENDDO
             t3 (kk) = t
          ENDDO
          !b
          t3 (1) = p (1) * t3 (1)
          DO kk = 2, n
             t = t3 (kk)
             DO j = ia (kk), idiag (kk) - 1
                t = t - a (j) * t3 (ja (j) )
             ENDDO
             t3 (kk) = p (kk) * t
          ENDDO
          !c
          DO kk = 1, n
             t3 (kk) = x (kk) + t3 (kk)
          ENDDO
          !d
          DO kk = n - 1, 1, - 1
             t = 0.0d0
             k2 = ia (kk+1)
             DO j = idiag (kk) + 1, k2 - 1
                t = t - a (j) * t3 (ja (j) )
             ENDDO
             t3 (kk) = t3 (kk) + p (kk) * t
          ENDDO
          !e
          DO kk = 1, n
             x (kk) = x (kk) + t1 (kk) - t3 (kk)
          ENDDO
          flop = flop + dble(3 * (ia (n + 1) - 1 + n) )
          RETURN
        END SUBROUTINE dagmg_sgsolve2
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_sgsolve (n, x, a, ja, ia, p, idiag, flop, inloc)
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), idiag(n)
          INTEGER ::  kk, j, k2, inloc(n)
          REAL(kind(0.0d0)) :: x(n), a(*), p(n), t
          REAL(kind(0.0d0)) :: flop
          x (1) = p (1) * x (1)
          DO kk = 2, n
             t = x (kk)
             DO j = ia (kk), idiag (kk) - 1
                t = t - a (j) * x (ja (j) )
             ENDDO
             x (kk) = p (kk) * t
          ENDDO
          DO kk = n - 1, 1, - 1
             t = 0.0d0
             k2 = ia (kk+1)
             DO j = idiag (kk) + 1, k2 - 1
                t = t - a (j) * x (ja (j) )
             ENDDO
             x (kk) = x (kk) + p (kk) * t
          ENDDO
          flop = flop + dble(2 * (ia (n + 1) - 1) )
          RETURN
        END SUBROUTINE dagmg_sgsolve
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_rescalc (n, x, y, b, a, ja, ia, flop,              &
             inloc, lstout, ilstout, lstin, ilstin )
          USE dagmg_mem
          IMPLICIT NONE
          INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier
          INTEGER :: inloc(n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
          REAL(kind(0.0d0)) :: x(n), y(n), b(n), a (*), t
          REAL(kind(0.0d0)) :: flop
          DO i = 1, n
             k1 = ia (i)
             t = b (i) - a (k1) * x (ja (k1) )
             k2 = ia (i+1)
             DO kk = k1 + 1, k2 - 1
                t = t - a (kk) * x (ja (kk) )
             ENDDO
             y (i) = t
          ENDDO
          flop = flop + dble(2 * (ia (n + 1) - 1) )
          RETURN
        END SUBROUTINE dagmg_rescalc
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_prolaggl (n, nc, V, B, ind)
          IMPLICIT NONE
          INTEGER :: n, nc, ind (n), k, i
          REAL(kind(0.0d0)) :: V (n), B (nc)
          DO i = 1, n
             k = ABS (ind (i) )
             IF (k.LE.nc) THEN
                V (i) = B (k)
             ELSE
                V (i) = 0.0d0
             ENDIF
          ENDDO
          RETURN
        END SUBROUTINE dagmg_prolaggl
!-----------------------------------------------------------------------
        SUBROUTINE dagmg_restaggl (n, nc, V, B, ind, flop)
          IMPLICIT NONE
          INTEGER :: n, nc, ind (n), k, i
          REAL(kind(0.0d0)) :: V (n), B (nc)
          REAL(kind(0.0d0)) :: flop
          B(1:nc)=0.0d0
          DO i = 1, n
             k = ABS (ind (i) )
             IF (k.LE.nc) B (k) = B (k) + V (i)
          ENDDO
          flop = flop + dble(n)
          RETURN
        END SUBROUTINE dagmg_restaggl
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE dagmg_setacg (nc, lcg, a, ja, ia, ac, jac, iac, nzac,  &
           ind, iw)
        IMPLICIT NONE
        INTEGER :: nc, lcg(4,nc), ja(*), ia(*), jac(*)
        INTEGER :: iac(nc+1), nzac
        INTEGER :: ind(*), iw(*)
        REAL(kind(0.0d0)) :: a(*), ac(*)
        INTEGER :: i, kk, jj, jc, kb, jcol, jpos
        DO i = 1, nc
           iw (i) = 0
        ENDDO
        nzac = 0
        iac (1) = 1
        DO 10 i = 1, nc
           DO 8 kk = 1, 4
              jj = lcg (kk, i)
              IF (jj.NE.0) THEN
                 DO kb = ia (jj), ia (jj + 1) - 1
                    jc = ja (kb)
                    jcol = ABS (ind (jc) )
                    IF (jcol.LE.nc) THEN
                       jpos = iw (jcol)
                       IF (jpos.EQ.0) THEN
                          nzac = nzac + 1
                          jac (nzac) = jcol
                          iw (jcol) = nzac
                          ac (nzac) = a (kb)
                       ELSE
                          ac (jpos) = ac (jpos) + a (kb)
                       ENDIF
                    ENDIF
                 ENDDO
              END IF
8          END DO
           DO kb = iac (i), nzac
              iw (jac (kb) ) = 0
           ENDDO
           iac (i + 1) = nzac + 1
10      END DO
        RETURN
      END SUBROUTINE dagmg_setacg
!*Routine for bacward compatibility************************************
      SUBROUTINE agmg(n,a,ja,ia,f,ijob,iprint,nrest,iter,tol)
          IMPLICIT NONE
          INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
          REAL (kind(0.0d0)) :: a(*),f(n),x(n)
          REAL (kind(0.0d0)) :: tol
          CALL dagmg(n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol)
          CALL DCOPY(n,x,1,f,1)
          RETURN
      END SUBROUTINE agmg
!*Timings***************************************************************
      SUBROUTINE dagmg_mestime (id,cputm,eltm)
        IMPLICIT NONE
        INTEGER, SAVE :: cpt_init(10)=-1,cpt_fin,cpt_max,freq,cpt
        REAL, SAVE :: t1(10), t2
        REAL(kind(0.0d0)) :: cputm,eltm
        INTEGER :: id
        IF (id>0) THEN
           !Next line may be uncommented if FORTRAN 95 function
           !CPU_TIME is implemented
           !   CALL CPU_TIME(t2)
           CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
           !
           cpt = cpt_fin - cpt_init(id)
           IF (cpt_fin < cpt_init(id)) cpt = cpt + cpt_max
           eltm = dble(cpt) / freq
           cputm = dble(t2 - t1(id))
           !
        ELSE
           !
           CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
           !Next line may be uncommented if FORTRAN 95 function
           !CPU_TIME is implemented
           !   CALL CPU_TIME(t1(-id))
           !
        END IF
        RETURN
      END SUBROUTINE dagmg_mestime
!***********************************************************************
