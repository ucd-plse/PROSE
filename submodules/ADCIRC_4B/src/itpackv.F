!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 41.09                                    *
!                                                                             *
!  mod history                                                                *
!  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
!                                          output buffer flush (m.cobb)       *
!                                          3D fixes (k.dresback)              *
!                                          drop MNPROC in fort.15 (t.campbell)*
!                                          various bug fixes in RBCs          *
!                                          ZSURFBUOY/BCPG calc                *
!  v40.02m002 - 12/22 - jjw/vjp - Vic suggested to avoid compiler conflicts   *
!****************************************************************************** 
! 
      MODULE ITPACKV
!
!  vjp  9/19/99
!
!-----------------------------------------------------------------------------
!     version -  itpackv 2d (january 1990)
!
!     code written by - david kincaid, roger grimes, john respess
!                       center for numerical analysis
!                       university of texas
!                       austin, tx  78712
!                       (512) 471-1242
!-----------------------------------------------------------------------------
!

#ifdef CMPI
      USE MESSENGER                !   MPI interface for padcirc
!!#else
!!      USE SIZES
#endif

#ifdef GPTL
      USE gptl
#endif
      implicit none

      logical :: foo


!--------------------data declarations end here-------------------------------c


      CONTAINS

      subroutine jcg (n,ndim,maxnz,jcoef,coef,rhs,u,iwksp,nw,wksp,
     *                iparm,rparm,ier)
      implicit none
      integer n,ndim,maxnz,nw,ier,ib1,ib2,ib3,ib4,ib5,nb,n3,i,nbo,
     *  itmax1,loop,idgts
      integer j, nandex, lnans, gnans, iimax !jgf46.00 add
      integer jcoef(ndim,maxnz),iwksp(3*n),iparm(12)
      real(kind=8) timi1,timj1,temp,time1,time2,timi2,timj2,digit1,
     *   digit2,tol
      real(kind=8) coef(ndim,maxnz),rhs(n),u(n),wksp(nw),rparm(12)
      real(kind=8), parameter :: c512_0 = 512.0, c0_0 = 0.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!     itpackv 2d main routine  jcg  (jacobi conjugate gradient)
!     each of the main routines --
!           jcg, jsi, sor, ssorcg, ssorsi, rscg, rssi
!     can be used independently of the others
!
! ... function --
!
!          jcg drives the jacobi conjugate gradient algorithm.
!
! ... parameter list --
!
!          n      input integer.  dimension of the matrix.
!          ndim   row dimension of jcoef and coef arrays in calling
!                   routine
!          maxnz  maximum number of nonzeros per row
!          jcoef  integer array for sparse matrix representation.
!          coef   array for sparse matrix representation.
!                 jcoef and coef use the ellpack data structure.
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution. on output, it contains
!                 the latest estimate to the solution.
!          iwksp  integer vector workspace of length 3*n
!          nw     input integer.  length of available wksp.  on output,
!                 iparm(8) is amount used.
!          wksp   vector used for working space.  jacobi conjugate
!                 gradient needs this to be in length at least
!                 4*n + 4*itmax.  here itmax = iparm(1) is the
!                 maximum allowable number of iterations.
!          iparm  integer vector of length 12.  allows user to specify
!                 some integer parameters which affect the method.
!          rparm  vector of length 12. allows user to specify some
!                 parameters which affect the method.
!          ier    output integer.  error flag.
!
! ... jcg module references --
!
!         from itpackv    chgcon, determ, dfault, echall,
!                         eigvns, eigvss, eqrt1s, iterm ,
!                         itjcg , parcon, permat, peror,
!                         pervec, pjac  , pmult , prbndx, pstop ,
!                         sbelm , scal  , unscal,
!                         vout  , zbrent
!          system         abs, log10, amax0, amax1, mod, sqrt
!
! ... local itpackv references --
!
!          echall, itjcg , permat,
!          peror, pervec, pjac  , prbndx, sbelm , scal  , unscal
!
!     version -  itpackv 2d (january 1990)
!
!     code written by - david kincaid, roger grimes, john respess
!                       center for numerical analysis
!                       university of texas
!                       austin, tx  78712
!                       (512) 471-1242
!
!     for additional details on the
!          (a) routine    see toms article 1982
!          (b) algorithm  see cna report 150
!
!     based on theory by - david young, david kincaid, lou hageman
!
!     reference the book - applied iterative methods
!                          l. hageman, d. young
!                          academic press, 1981
!
!     **************************************************
!     *               important note                   *
!     *                                                *
!     *      when installing itpackv routines on a     *
!     *  different computer, reset some of the values  *
!     *  in  subroutne dfault.   most important are    *
!     *                                                *
!     *   srelpr      machine relative precision       *
!     *   rparm(1)    stopping criterion               *
!     *                                                *
!     *   also change system-dependent routine         *
!     *                                                *
!     **************************************************
!
! ... variables in common block - itcom1
!
!     in     - iteration number
!     is     - iteration number when parameters last changed
!     isym   - symmetric/nonsymmetric case switch
!     itmax  - maximum number of iterations allowed
!     level  - level of output control switch
!     nout   - output unit number
!
! ... variables in common block - itcom2
!
!     adapt  - fully adaptive procedure switch
!     betadt - switch for adaptive determination of beta
!     caseii - adaptive procedure case switch
!     halt   - stopping test switch
!     partad - partially adaptive procedure switch
!
! ... variables in common block - itcom3
!
!     bdelnm - two norm of b times delta-super-n
!     betab  - estimate for the spectral radius of lu matrix
!     cme    - estimate of largest eigenvalue
!     delnnm - inner product of pseudo-residual at iteration n
!     delsnm - inner product of pseudo-residual at iteration s
!     ff     - adaptive procedure damping factor
!     gamma  - acceleration parameter
!     omega  - overrelaxation parameter for sor and ssor
!     qa     - pseudo-residual ratio
!     qt     - virtual spectral radius
!     rho    - acceleration parameter
!     rrr    - adaptive parameter
!     sige   - parameter sigma-sub-e
!     sme    - estimate of smallest eigenvalue
!     specr  - spectral radius estimate for ssor
!     srelpr - machine relative precision
!     stptst - stopping parameter
!     udnm   - two norm of u
!     zeta   - stopping criterion
!
! ... initialize common blocks
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      ier = 0
      level = iparm(2)
      nout = iparm(4)
      if (level.ge.1) write (nout,10)
   10 format (///1x,'i t p a c k      j c g      ')
      if (iparm(1).le.0) go to 370
!     if (iparm(11).eq.0) timj1 = timer(0.0)
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,rparm,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!      temp = 500.0*srelpr  !jgf46.00 commented out
      temp = c512_0*srelpr !jgf46.00 added
      if (zeta.ge.temp) go to 30
      if (level.ge.1) write (nout,20) zeta,srelpr,temp
   20 format (/1x,'warning rparm(1) =',e10.3,' (zeta)'/1x,'srelpr =',e10.3/1x,'zeta reset to ',e10.3)
      zeta = temp
   30 continue
      time1 = rparm(9)
      time2 = rparm(10)
      digit1 = rparm(11)
      digit2 = rparm(12)
!
! ... verify n
!
      if (n.gt.0) go to 50
      ier = 11
      if (level.ge.0) write (nout,40) n
   40 format (/1x,'error: called from itpackv routine jcg invalid matrix dimension, n =',i8)
      go to 370
!
! ... scale linear system, u, and rhs by the square root of the
! ... diagonal elements.
!
   50 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call scal (n,ndim,maxnz,jcoef,coef,rhs,u,wksp,ier)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif 
      if (ier.eq.0) go to 70
      if (level.ge.0) write (nout,60) ier
   60 format (/1x,'error'//1x,'called from itpackv routine jcg error detected ier = ',i5)
      go to 370
!
! ... remove rows and columns if requested
!
   70 continue
      if (iparm(10).eq.0) go to 80
      tol = rparm(8)
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call sbelm (n,ndim,maxnz,jcoef,coef,rhs,wksp,tol)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
! ... initialize wksp base addresses.
!
   80 ib1 = 1
      ib2 = ib1+n
      ib3 = ib2+n
      ib4 = ib3+n
      ib5 = ib4+n
!
! ... permute to  red-black system if requested
!
      nb = iparm(9)
      if (nb.ge.0) go to 110
      if (nb.le.-2) go to 170
      n3 = n*3
      do 90 i = 1,n3
         iwksp(i) = 0
   90 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call prbndx (n,ndim,maxnz,jcoef,iwksp,iwksp(ib2),nb,level,nout,
     *   ier)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      if (ier.eq.0) go to 110
      if (level.ge.0) write (nout,100) ier,nb
  100 format (/1x,'error: called from itpackv routine jcg ier = ',i5,' iparm(9) = ',i5,' (nb)')
      go to 350
  110 if (nb.ge.0.and.nb.le.n) go to 130
      ier = 14
      if (level.ge.0) write (nout,120) ier,nb
  120 format (/1x,'error: called from itpackv routine jcg. ier = ',i5,' iparm(9) =',i5,' (nb)')
      go to 350
  130 if (nb.ne.0.and.nb.ne.n) go to 150
      nbo = nb
      nb = n/2
      if (level.ge.2) write (nout,140) nbo,nb
  140 format (/10x,' nb = ',i5,' implies matrix is diagonal nb reset to ',i5)
!
! ... permute matrix and rhs
!
  150 if (level.ge.2) write (nout,160) nb
  160 format (/10x,'order of black subsystem = ',i5,' (nb)')
      if (iparm(9).ge.0) go to 170
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call permat (n,ndim,maxnz,jcoef,coef,iwksp,wksp,iwksp(ib3:ib3+n))
      call pervec (n,iwksp,rhs,wksp)
      call pervec (n,iwksp,u,wksp)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
! ... check for sufficient workspace.
!
  170 iparm(8) = 4*n+4*itmax
      if (nw.ge.iparm(8)) go to 190
      ier = 12
      if (level.ge.0) write (nout,180) nw,iparm(8)
  180 format (/1x,'error: called from itpackv routine jcg not enough workspace at ',i10/1x, ' iparm(8) =',i10)
      go to 330
!
  190 continue
      if (level.le.2) go to 220
      write (nout,200)
  200 format (///1x,'in the following, rho and gamma are acceleration parameters')
      if (adapt) write (nout,210)
  210 format (1x,'cme is the estimate of the largest eigenvalue of the jacobi matrix')
  220 continue
!     if (iparm(11).eq.0) timi1 = timer(0.0)
!
! ... compute initial pseudo-residual
!
      do 230 i = 1,nw
         wksp(i) = c0_0 !jgf46.00 added d0
  230 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call scopy (n,rhs,1,wksp(ib2:),1)
      call pjac (n,ndim,maxnz,jcoef,coef,u,wksp(ib2:))
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      do 240 i = 1,n
         wksp(n+i) = wksp(n+i)-u(i)
  240 continue
!
! ... iteration sequence
!
      itmax1 = itmax+1
      do 260 loop = 1,itmax1
         in = loop-1
         if (mod(in,2).eq.1) go to 250
!
! ... code for the even iterations.
!
!     u           = u(in)             wksp(ib2) = del(in)
!     wksp(ib1)   = u(in-1)           wksp(ib3) = del(in-1)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
         call itjcg (n,ndim,maxnz,jcoef,coef,u,wksp(ib1:),wksp(ib2:),
     *      wksp(ib3:),wksp(ib4:),wksp(ib5:))
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
         if (halt) go to 290
         go to 260
!
! ... code for the odd iterations.
!
!     u           = u(in-1)           wksp(ib2) = del(in-1)
!     wksp(ib1)   = u(in)             wksp(ib3) = del(in)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
  250    call itjcg (n,ndim,maxnz,jcoef,coef,wksp(ib1:),u,wksp(ib3:),
     *      wksp(ib2:),wksp(ib4:),wksp(ib5:))
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
         if (halt) go to 290
  260 continue
!
!..... itmax has been reached
!
      if (iparm(11).ne.0) go to 270
!     timi2 = timer(0.0)
!     time1 = timi2-timi1
  270 ier = 13
      if (level.ge.1) write (nout,280) itmax
  280 format (/1x,'warning: in itpackv routine jcg failure to converge in',i5,' iterations')
      if (iparm(3).eq.0) rparm(1) = stptst
      go to 320
!
! ... method has converged
!
  290 if (iparm(11).ne.0) go to 300
!     timi2 = timer(0.0)
!     time1 = timi2-timi1
  300 if (level.ge.1) write (nout,310) in
  310 format (/1x,'jcg  has converged in ',i5,' iterations')
!
! ... put solution into u if not already there.
!
  320 continue
      if (mod(in,2).eq.1) then
#ifdef GPTL
            gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
            call scopy (n,wksp,1,u,1)
#ifdef GPTL
            gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      end if
!
! ... un-permute matrix,rhs, and solution
!
  330 if (iparm(9).ne.-1) go to 340
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call permat (n,ndim,maxnz,jcoef,coef,iwksp(ib2:),wksp(ib4:),
     *   iwksp(ib3:))
      call pervec (n,iwksp(ib2:),rhs,wksp(ib4:))
      call pervec (n,iwksp(ib2:),u,wksp(ib4:))
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      if (ier.eq.12) go to 350
!
! ... optional error analysis
!
  340 idgts = iparm(12)
      if (idgts.lt.0) go to 350
      if (iparm(2).le.0) idgts = 0
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call peror (n,ndim,maxnz,jcoef,coef,rhs,u,wksp,digit1,digit2,
     *   idgts)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
! ... unscale the matrix, solution, and rhs vectors.
!
  350 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      call unscal (n,ndim,maxnz,jcoef,coef,rhs,u,wksp)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
!
! ... set return parameters in iparm and rparm
!
      iparm(8) = iparm(8)-4*(itmax-in)
!      if (iparm(11).ne.0) go to 360 !jgf46.00 commented out
!     timj2 = timer(0.0)
!     time2 = timj2-timj1
!      time2 = 0.0                   !jgf46.00 commented out
  360 if (iparm(3).ne.0) go to 370
      iparm(1) = in
      iparm(9) = nb
      rparm(2) = cme
      rparm(3) = sme
!      rparm(9) = time1              !jgf46.00 commented out
!      rparm(10) = time2             !jgf46.00 commented out
      rparm(11) = digit1
      rparm(12) = digit2
!
  370 continue
      if (level.ge.3) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
            call echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,
     *   rparm,2)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::jcg", gptl_handle)
#endif
      end if
      if (ier.eq.0.and.level.ge.1) write (nout,380)
  380 format (/1x,'execution successful')
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::jcg", gptl_handle)
#endif
      return
      end subroutine


      subroutine itjcg (n,ndim,maxnz,jcoef,coef,u,u1,d,d1,dtwd,tri)
      use global,only:nperseg
      implicit none
!
      integer n,ndim,maxnz,i,j
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz)
      real(kind=8) u(n),u1(n),d(n),d1(n),dtwd(n),tri(3*n)
      real(kind=8) gamold,rhoold,rhotmp,dnrm,con,dtnrm,c1,c2,c3,c4
      real(kind=8) dumy1(1),dumy2(1),del3nrms(3),unorm
      real(kind=8) vic1, vic2 !jgf46.00 added
      logical q1
      real(kind=8), parameter :: c0_0 = 0.0, c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... itjcg performs one iteration of the jacobi conjugate gradient
!     algorithm.  it is called by jcg.
!
! ... parameter list --
!
!          n      input integer.  dimension of the matrix.
!          ndim   row dimension of jcoef and coef arrays in calling
!                   routine
!          maxnz  maximum number of nonzeros per row
!          jcoef  integer sparse matrix representation
!          coef   sparse matrix representation
!          u      input vector.  contains the value of the
!                 solution vector at the end of in iterations.
!          u1     input/output vector.  on input, it contains
!                 the value of the solution at the end of the in-1
!                 iteration.  on output, it will contain the newest
!                 estimate for the solution vector.
!          d      input vector.  contains the pseudo-residual
!                 vector after in iterations.
!          d1     input/output vector.  on input, d1 contains
!                 the pseudo-residual vector after in-1 iterations.  on
!                 output, it will contain the newest pseudo-residual
!                 vector.
!          dtwd   array.  used in the computations of the
!                 acceleration parameter gamma and the new pseudo-
!                 residual.
!          tri    array.  stores the tridiagonal matrix associated
!                 with the eigenvalues of the conjugate gradient
!                 polynomial.
!
! ... local itpackv references --
!          chgcon, iterm , parcon, pjac  , pstop 
!
!     description of variables in common blocks in routine jcg
!
! ... compute new estimate for cme if adapt = .true.
!
      save
      if (adapt) then
            call chgcon (itmax,tri,gamold,rhoold,1)
      end if
!
!   ...UPDATE PSEUDO-RESIDUAL VECTOR "d" OF SIZE "n"
!   ...BEFORE PERFORMING ONE JACOBI ITERATION
!

#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif

#ifdef CMPI
      dumy1(1) = c0_0 !jgf46.00 added
      dumy2(1) = c0_0 !jgf46.00 added
      if(nperseg.gt.0)then
          CALL UPDATER_W_PERBC(d,dumy1,dumy2,1) ;
      else
          CALL UPDATER(d,dumy1,dumy2,1)
      endif
#endif

      do 10 i = 1,n
         dtwd(i) = c0_0
   10 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      call pjac (n,ndim,maxnz,jcoef,coef,d,dtwd)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif

!       
! ... test for stopping
!
       if (q1 .OR.  ((in .gt.5) .AND. 
     *        (mod(in,5).ne. 0))) then

#ifdef CMPI
        call ps2dots(n,d,dtwd,del3nrms)
        delnnm = del3nrms(1)
        dtnrm  = del3nrms(2)
        unorm  = del3nrms(3)
#else
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      delnnm = sdot(n,d,1,d,1)
      dtnrm = sdot(n,d,1,dtwd,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      unorm = c1_0 !jgf46.00 added
#endif

        else

#ifdef CMPI
        call ps3dots(n,d,dtwd,u,del3nrms)
        delnnm = del3nrms(1)
        dtnrm  = del3nrms(2)
        unorm  = del3nrms(3)
#else
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      delnnm = sdot(n,d,1,d,1)
      dtnrm = sdot(n,d,1,dtwd,1)
      unorm = sdot(n,u,1,u,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif
#endif
       end if

      dnrm = delnnm
      con = cme
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      call pstop_nrms (n,unorm,dnrm,con,1,q1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif
      if (halt) go to 50
 


       ! ... compute rho and gamma - acceleration parameters
 
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
   20 call parcon (dtnrm,c1,c2,c3,c4,gamold,rhoold,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::itjcg", gptl_handle)
#endif
!
! ... compute u(in+1) and d(in+1)
!
   30 do 40 i = 1,n
         u1(i) = c1*d(i)+c2*u(i)+c3*u1(i)
         d1(i) = c1*dtwd(i)+c4*d(i)+c3*d1(i)
   40 continue

!
! ... output intermediate information
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::itjcg", gptl_handle)
#endif
   50 call iterm (n,ndim,maxnz,coef,u,dtwd,1)

      return
      end subroutine


      subroutine peror (n,ndim,maxnz,jcoef,coef,rhs,u,work,
     *                   digit1,digit2,idgts)
      implicit none
      integer i,n,ndim,maxnz,idgts
      integer jcoef(ndim,maxnz)
      real(kind=8)  digit1,digit2,bnrm,rnrm,temp
      real(kind=8) coef(ndim,maxnz),rhs(n),u(n),work(n)
      real(kind=8), parameter :: c0_0 = 0.0, c2_0 = 2.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!     peror computes the residual, r = rhs - a*u.  the user
!     also has the option of printing the residual and/or the
!     unknown vector depending on idgts.
!
! ... parameter list --
!
!          n      dimension of matrix
!          ndim   row dimension of jcoef and coef in calling routine
!          maxnz  maximum number of nonzeros per row
!          jcoef  integer array of sparse matrix representation
!          coef   array of sparse matrix representation
!          rhs    right hand side of matrix problem
!          u      latest estimate of solution
!          work   workspace vector of length 2*n
!          digit1 output - measure of accuracy of stopping test
!          digit2 output - measure of accuracy of solution
!          idgts   parameter controlling level of output
!                    if idgts < 1 or idgts > 4, then no output.
!                            = 1, then number of digits is printed, pro-
!                                 vided level .ge. 1
!                            = 2, then solution vector is printed, pro-
!                                 vided level .ge. 1
!                            = 3, then residual vector is printed, pro-
!                                 vided level .ge. 1
!                            = 4, then both vectors are printed, pro-
!                                 vided level .ge. 1
!
! ... local itpackv references --
!
!          pmult , vout
!
! ... specifications for arguments
!
!
!     description of variables in common block in main routine
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
      digit1 = c0_0
      digit2 = c0_0
      if (n.le.0) go to 70
!
      digit1 = -log10(abs(srelpr))
      if (stptst.gt.c0_0) digit1 = -log10(abs(stptst))
      do 10 i = 1,n
         work(i) = rhs(i)/coef(i,1)
   10 continue

#ifdef CMPI
      bnrm = psdot(n,work,work)
#else
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
      bnrm = sdot(n,work,1,work,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
#endif

      if (bnrm.eq.c0_0) go to 30
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
      call pmult (n,ndim,maxnz,jcoef,coef,u,work)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
      do 20 i = 1,n
         work(i) = (rhs(i)-work(i))/coef(i,1)
   20 continue

#ifdef CMPI
      rnrm = psdot(n,work,work)
#else
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
      rnrm = sdot(n,work,1,work,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
#endif


      temp = rnrm/bnrm
      if (temp.eq.c0_0) go to 30
      digit2 = -log10(abs(temp))/c2_0
      go to 40
!
   30 digit2 = -log10(abs(srelpr))
!
   40 if ((idgts.lt.1).or.(level.le.0)) go to 70
      write (nout,50) digit1,digit2
   50 format (/10x,'approx. no. of digits in stopping tes', f5.1, f5.1,'(digit2)')
!
      if (idgts.le.1.or.idgts.gt.4) go to 70
      if (idgts.ge.3) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
            call vout (n,work,1,nout)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
      end if
      do 60 i = 1,n
         work(i) = u(i)*coef(i,1)
   60 continue
      if (idgts.ne.3) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
            call vout (n,work,2,nout)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::peror", gptl_handle)
#endif
      end if
!
   70 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::peror", gptl_handle)
#endif
      return
      end subroutine

      subroutine pstop (n,u,dnrm,ccon,iflag,q1)
      implicit none
      integer n,iflag
      logical q1
      real(kind=8) u(n),dnrm,con,ccon,uold,tr,tl
      real(kind=8), parameter :: c1_0 = 1.0, c1000_0 = 1000.0, c0_0 = 0.0, c2_0 = 2.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!     pstop performs a test to see if the iterative
!     method has converged to a solution inside the error
!     tolerance, zeta.
!
! ... parameter list --
!
!          n      order of system
!          u      present solution estimate
!          dnrm   inner product of pseudo-residuals at preceding
!                    iteration
!          con    stopping test parameter (= ccon)
!          iflag  stopping test integer flag
!                    iflag = 0,  sor iteration zero
!                    iflag = 1,  non-rs method
!                    iflag = 2,  rs method
!          q1     stopping test logical flag
!
!
!     description of variables in common block in main routine
!

#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::pstop", gptl_handle)
#endif
      con = ccon
      halt = .false.
!
!     special procedure for zeroth iteration
!
      if (in.ge.1) go to 10
      q1 = .false.
      udnm = c1_0
      stptst = c1000_0
      if (iflag.le.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pstop", gptl_handle)
#endif
      return
      end if
!
! ... test if udnm needs to be recomputed
!
   10 continue
      if (q1) go to 20
      if ((in.gt.5).and.(mod(in,5).ne.0)) go to 20
      uold = udnm

#ifdef CMPI
      udnm = psdot(n,u,u)
#else
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pstop", gptl_handle)
#endif
      udnm = sdot(n,u,1,u,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::pstop", gptl_handle)
#endif
#endif

      if (udnm.eq.c0_0) udnm = c1_0
      if ((in.gt.5).and.(abs(udnm-uold).le.udnm*zeta)) q1 = .true.
!
! ... compute stopping test
!
   20 tr = sqrt(udnm)
      tl = c1_0
      if (con.eq.c1_0) go to 40
      if (iflag.eq.2) go to 30
      tl = sqrt(dnrm)
      tr = tr*(c1_0-con)
      go to 40
   30 tl = sqrt(c2_0*dnrm)
      tr = tr*(c1_0-con*con)
   40 stptst = tl/tr
      if (tl.ge.tr*zeta) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pstop", gptl_handle)
#endif
      return
      end if
      halt = .true.
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pstop", gptl_handle)
#endif
      return
      end subroutine


      subroutine pstop_nrms (n,unrm,dnrm,ccon,iflag,q1)
      implicit none
      integer n,iflag
      logical q1
      real(kind=8) unrm,dnrm,con,ccon,uold,tr,tl
      real(kind=8), parameter :: c1_0 = 1.0, c1000_0 = 1000.0, c0_0 = 0.0, c2_0 = 2.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!     pstop performs a test to see if the iterative
!     method has converged to a solution inside the error
!     tolerance, zeta.
!
! ... parameter list --
!
!          n      order of system
!          u      present solution estimate
!          dnrm   inner product of pseudo-residuals at preceding
!                    iteration
!          con    stopping test parameter (= ccon)
!          iflag  stopping test integer flag
!                    iflag = 0,  sor iteration zero
!                    iflag = 1,  non-rs method
!                    iflag = 2,  rs method
!          q1     stopping test logical flag
!
!
!     description of variables in common block in main routine
!
#ifdef GPTL
      gptl_ret = gptlstart_handle( "::ITPACKV::pstop_nrms", gptl_handle)
#endif
      con = ccon
      halt = .false.
!
!     special procedure for zeroth iteration
!
      if (in.ge.1) go to 10
      q1 = .false.
      udnm = c1_0
      stptst = c1000_0
      if (iflag.le.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle( "::ITPACKV::pstop_nrms", gptl_handle)
#endif
      return
      end if
!
! ... test if udnm needs to be recomputed
!
   10 continue
      if (q1) go to 20
      if ((in.gt.5).and.(mod(in,5).ne.0)) go to 20
      uold = udnm

      udnm = unrm

      if (udnm.eq.c0_0) udnm = c1_0 !jgf46.00 added d0 to 0.0
      if ((in.gt.5).and.(abs(udnm-uold).le.udnm*zeta)) q1 = .true.
!
! ... compute stopping test
!
   20 tr = sqrt(udnm)
      tl = c1_0
      if (con.eq.c1_0) go to 40
      if (iflag.eq.2) go to 30
      tl = sqrt(dnrm)
      tr = tr*(c1_0-con)
      go to 40
   30 tl = sqrt(c2_0*dnrm)
      tr = tr*(c1_0-con*con)
   40 stptst = tl/tr
      if (tl.ge.tr*zeta) then
#ifdef GPTL
      gptl_ret = gptlstop_handle( "::ITPACKV::pstop_nrms", gptl_handle)
#endif
      return
      end if
      halt = .true.
!
#ifdef GPTL
      gptl_ret = gptlstop_handle( "::ITPACKV::pstop_nrms", gptl_handle)
#endif
      return
      end subroutine


      subroutine chgcon (ldt,tri,gamold,rhoold,ibmth)
      implicit none
      integer ldt,ibmth,ip,ier
      real(kind=8) tri(ldt)
      real(kind=8) gamold,rhoold,cmold,start,end1
      real(kind=8), parameter :: c0_0 = 0.0, c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... chgcon computes the new estimate for the largest eigenvalue
!     for conjugate gradient acceleration.
!
! ... parameter list --
!
!          ldt    leading dimension of tri
!          tri    tridiagonal matrix associated with the eigenvalues
!                    of the conjugate gradient polynomial
!          gamold
!            and
!          rhoold previous values of acceleration parameters
!          ibmth  indicator of basic method being accelerated by cg
!                      ibmth = 1,  jacobi
!                            = 2,  reduced system
!                            = 3,  ssor
!
! ... local itpackv references --
!
!          eigvns, eigvss
      save
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      go to (10,20,30), ibmth
!
! ... jacobi conjugate gradient
!
   10 start = cme
      ip = in
      go to 40
!
! ... reduced system cg
!
   20 start = cme**2
      ip = in
      go to 40
!
! ... ssor cg
!
   30 if (adapt) start = spr
      if (.not.adapt) start = specr
      ip = in-is
!
! ... define the matrix
!
   40 if (ip.ge.2) go to 60
      if (ip.eq.1) go to 50
!
! ... ip = 0
!
      end1 = c0_0   !jgf46.00 added d0
      cmold = c0_0  !jgf46.00 added d0
      go to 110
!
! ... ip = 1
!
   50 end1 = c1_0-c1_0/gamma
      tri(1) = end1
      tri(1+ldt) = c0_0
      go to 110
!
! ... ip > 1
!
   60 if (abs(start-cmold).le.zeta*start) go to 120
      cmold = start
!
! ... compute the largest eigenvalue
!
      tri(ip) = c1_0-c1_0/gamma
      tri(ip+ldt) = (c1_0-rho)/(rho*rhoold*gamma*gamold)
      if (isym.ne.0) go to 80
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      end1 = eigvss(ip,tri,start,zeta,itmax,ier)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      if (ier.eq.0) go to 100
      if (level.ge.2) write (nout,70) ier
   70 format (/10x,'difficulty in computation of maximum eigenvalue ier =',i5)
      go to 100
   80 continue
      !end1 = eigvns(ldt,ip,tri,tri(1,3),tri(1,4),ier)
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      end1 = eigvns(ldt,ip,tri,tri(1+ldt*2:ldt*4),tri(1+ldt*3:ldt*4),ier)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      if (ier.eq.0) go to 100
      if (level.ge.2) write (nout,90) ier
   90 format (/10x,'difficulty in computation of maximum eigenvalue ier =',i5)
  100 continue
      if (ier.ne.0) go to 130
!
! ... set spectral radius for the various methods
!
  110 if (ibmth.eq.1) cme = end1
      if (ibmth.eq.2) cme = sqrt(abs(end1))
      if (ibmth.eq.3.and.adapt) spr = end1
      if (ibmth.eq.3.and..not.adapt) specr = end1
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      return
!
! ... relative change in cme is less than zeta.  therefore stop
!     changing.
!
  120 adapt = .false.
      partad = .false.
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      return
!
! ... estimate for cme > one.  therefore need to stop adaptive
!     procedure and keep old value of cme.
!
  130 adapt = .false.
      partad = .false.
      if (level.ge.2) write (nout,140) in,start
  140 format (/10x,'estimate of maximum eigenvalue',i5/10x, 'final estimate =',e15.7/)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::chgcon", gptl_handle)
#endif
      return
      end  subroutine

      real(kind=8) function determ (ldt,n,tri,wk1,wk2,xlmda)
      implicit none
      integer ldt,n,i,l
      real(kind=8) xlmda
      real(kind=8) tri(3*n),wk1(n),wk2(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::determ", gptl_handle)
#endif
!
!     determ computes the determinant of a symmetric
!     tridiagonal matrix given by tri. det(tri - xlmda*i) = 0
!
! ... parameter list --
!
!          ldt    leading dimension of array tri
!          n      order of tridiagonal system
!          tri    symmetric tridiagonal matrix of order n
!          wk1,   workspace vectors of length n
!           wk2
!          xlmda  argument for characteristic equation
!
!
      do 10 i = 1,n
         wk1(i) = tri(i)-xlmda
   10 continue
      wk2(n) = wk1(n)
      wk2(n-1) = wk1(n-1)*wk2(n)+tri(n+ldt)
      if (n.eq.2) go to 30
!
! ... beginning of loop
!
      do 20 l = n-2,1,-1
         wk2(l) = wk1(l)*wk2(l+1)+tri(l+1+ldt)*wk2(l+2)
   20 continue
!
!     wk2(1) = solrn (n,wk1(-1),-1,tri(0,2),-1,wk2,-1)
!
! ... determinant computed
!
   30 determ = wk2(1)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::determ", gptl_handle)
#endif
      return
      end  function


      subroutine dfault (iparm,rparm)
      implicit none
      integer iparm(12)
      real(kind=8) rparm(12), temp
      real(kind=8), parameter :: c1_0 = 1.0, c0_5 = 0.5, c2_0 = 2.0, c5e_6 = 5e-6, c0_0 = 0.0, c0_75 = 0.75, c0_25 = 0.25, c128_0 = 128.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... dfault sets the default values of iparm and rparm.
!
! ... parameter list --
!
!          iparm
!           and
!          rparm  arrays specifying options and tolerances
!
! ... specifications for arguments
!
!
!     description of variables in common blocks in main routine
!
!     srelpr  - computer precision (approx.)
!     if installer of package does not know srelpr value,
!     an approximate value can be determined from a simple
!     fortran program such as
!
!     srelpr = 1.0e0
!   2 srelpr = 0.5e0*srelpr
!     temp = srelpr + 1.0e0
!     if (temp .gt. 1.0e0)  go to 2
!     srelpr = 2.0e0*srelpr
!     write (6,3) srelpr
!   3 format (5x,e15.8)
!     stop
!     end
!
!     some values are-
!
!     srelpr = 7.1e-15   for cray x-mp, y-mp  (approx.) 2**-47
!          = 1.49e-8   for dec 10  (approx.) 2**-26
!          = 1.192e-7  for vax 11/780 (approx) 2**-23
!          = 1.192e-7  for Sun Spark Station 2 (approx) 2**-23
!          = 1.192e-7  for IBM RISC 6000 2 (approx) 2**-23
!          = 4.768e-7  for ibm 370/158
!          = 0.5960E-7 for Lahey fortran on ALR (486 PC) rl 12/92
!          = 0.5960E-7 for Vax 6000
!             *** should be changed for other machines ***
!
!     to facilitate convergence, rparm(1) should be set to
!          500.*srelpr or larger
!
!vjp--Determine macheps 
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::dfault", gptl_handle)
#endif
      srelpr = c1_0
    2 srelpr = c0_5*srelpr
      temp = srelpr + c1_0
      if (temp .gt. c1_0)  go to 2
      srelpr = c2_0*srelpr
!
      iparm(1) = 100
      iparm(2) = -1 !jgf46.00 changed to -1, was 0
      iparm(3) = 0
      iparm(4) = 6
      iparm(5) = 0
      iparm(6) = 1
      iparm(7) = 1
      iparm(8) = 0
      iparm(9) = -2
      iparm(10) = 0
      iparm(11) = 0
      iparm(12) = 3 !jgf46.00 changed to -3, was 0
!
      rparm(1) = c5e_6 !jgf46.00 was 512.e0*srelpr
      rparm(2) = c0_0
      rparm(3) = c0_0
      rparm(4) = c0_75
      rparm(5) = c1_0
      rparm(6) = c0_0
      rparm(7) = c0_25
      rparm(8) = c128_0*srelpr !jgf46.00 was 100.0e0*srelpr
      rparm(9) =  c0_0
      rparm(10) = c0_0
      rparm(11) = c0_0
      rparm(12) = c0_0
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::dfault", gptl_handle)
#endif
      return
      end  subroutine

      subroutine echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,rparm,icall)
      implicit none
      integer n,i,j,ndim,maxnz,icall
      integer jcoef(ndim,maxnz),iparm(12)
      real(kind=8) coef(ndim,maxnz),rhs(n),rparm(12)
      real(kind=8), parameter :: c0_0 = 0.0, c1_0 = 1.0, c_1_0 = -1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... echall initializes the itpackv common blocks from the
! ... information contained in iparm and rparm. echall also prints the
! ... values of all the parameters in iparm and rparm.
!
! ... parameter list --
!
!          iparm
!           and
!          rparm  arrays of parameters specifying options and
!                    tolerances
!          icall  indicator of which parameters are being printed
!                    icall = 1,  initial parameters
!                    icall = 2,  final parameters
!
! ... specifications for arguments
!
!     description of variables in common blocks in main routine
!
      if (icall.ne.1) go to 120
!
! ... initialize itpackv common
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::echall", gptl_handle)
#endif
      zeta = rparm(1)
      cme = rparm(2)
      sme = rparm(3)
      ff = rparm(4)
      omega = rparm(5)
      specr = rparm(6)
      betab = rparm(7)
      itmax = iparm(1)
      level = iparm(2)
      isym = iparm(5)
!
      adapt = .false.
      partad = .false.
      betadt = .false.
      if (iparm(6).eq.1.or.iparm(6).eq.3) adapt = .true.
      if (iparm(6).eq.1) betadt = .true.
      if (iparm(6).eq.2) partad = .true.
!
      caseii = .false.
      if (iparm(7).eq.2) caseii = .true.
      if (caseii) sme = -cme
      if (.not.caseii.and.sme.eq.c0_0) sme = c_1_0 !jgf46.00 added d0
      spr = sme
!
! ... set rest of common variables to zero
!
      in = 0
      is = 0
      halt = .false.
      bdelnm = c0_0 !jgf46.00 added d0 to all these
      delnnm = c0_0
      delsnm = c0_0
      gamma = c0_0
      qa = c0_0
      qt = c0_0
      rho = c0_0
      rrr = c0_0
      sige = c0_0
      stptst = c0_0
      udnm = c0_0
!
      if (level.le.4) go to 100
!
!     this section of echall causes printing of the linear system and
!     the iterative parameters
!
      write (nout,10)
   10 format (///5x,'the linear system is as follows')
      write (nout,20)
   20 format (/2x,'jcoef array')
      do 30 i = 1,n
         write (nout,40) (jcoef(i,j),j=1,maxnz)
   30 continue
   40 format (1x,8(1x,i8))
      write (nout,50)
   50 format (/2x,'coef array')
      do 60 i = 1,n
         write (nout,70) (coef(i,j),j=1,maxnz)
   60 continue
   70 format (1x,5(2x,g14.6))
      write (nout,80)
   80 format (/2x,'rhs array')
      write (nout,90) (rhs(i),i=1,n)
   90 format (1x,5g16.6)
  100 if (level.le.2) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::echall", gptl_handle)
#endif
      return
      end if
      write (nout,110)
  110 format (///5x,'initial iterative parameters'/)
      go to 140
  120 write (nout,130)
  130 format (///5x,'final iterative parameters'/)
  140 write (nout,150) iparm(1),level,iparm(3),nout,isym,iparm(6)
  150 format (10x,'iparm(1)  =',i15,4x,i15,4x,i15,4x,i15,4x,i15,4x,i15,4x)
      write (nout,160) iparm(7),iparm(8),iparm(9),iparm(10),iparm(11),
     *   iparm(12)
  160 format (10x,'iparm(7)  =',i15,4x,i15,4x,i15,4x,i15,4x,i15,4x,i15,4x)
      write (nout,170) zeta,cme,sme,ff,omega,specr
  170 format (10x,'rparm(1)  =',e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x)
      write (nout,180) betab,rparm(8),rparm(9),rparm(10),rparm(11),
     *   rparm(12)
  180 format (10x,'rparm(7)  =',e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x,e15.8,4x)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::echall", gptl_handle)
#endif
      return
      end  subroutine


      real(kind=8) function eigvns (ldt,n,tri,d,e2,ier)
      implicit none
      integer ldt,n,i,ier
      real(kind=8) tri(3*n),d(n),e2(n)
      real(kind=8), parameter :: c0_0 = 0.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::eigvns", gptl_handle)
#endif
!
! ... eigvns computes the largest eigenvalue of a symmetric
!     tridiagonal matrix for conjugate gradient acceleration.
!
! ... parameter list --
!
!          ldt    leading dimension of tri
!          n      order of tridiagonal system
!          tri    symmetric tridiagonal matrix of order n
!          d      array for eqrt1s (negative diagonal elements)
!          e2     array for eqrt1s (super diagonal elements)
!          ier    error flag -- on return, ier=0 indicates that
!                    the largest eigenvalue of tri was found.
!
! ... local itpackv references --
!
!          eqrt1s
!
! ... specifications for arguments
!
!
      eigvns = c0_0
!
      d(1) = -tri(1)
      do 10 i = 2,n
         d(i) = -tri(i)
         e2(i) = abs(tri(i+ldt))
   10 continue
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::eigvns", gptl_handle)
#endif
      call eqrt1s (d,e2,n,1,0,ier)
      eigvns = -d(1)
!
      return
      end function 


      real(kind=8) function eigvss (n,tri,start,zeta,itmax,ier)
      implicit none
      integer n,itmax,itmp,ier,maxfn,nsig
      real(kind=8) tri(3*n),start,eps,a,b,zeta
      real(kind=8), parameter :: c1_0 = 1.0, c0_0 = 0.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::eigvss", gptl_handle)
#endif
!
! ... eigvss computes the largest eigenvalue of a symmetric
!     tridiagonal matrix for conjugate gradient acceleration.
!     modified imsl routine zbrent used.
!
! ... parameter list --
!
!          n      order of tridiagonal system
!          tri    symmetric tridiagonal matrix of order n
!          start  initial lower bound of interval containing root
!          zeta   stopping criteria
!          ier    error flag -- on return, ier = 0 indicates that
!                    the largest eigenvalue of tri was found.
!
! ... local itpackv references --
!
!          zbrent
!
! ... specifications for arguments
!
!
      eigvss = c0_0
!      itmp = int(-log10(abs(zeta))) !jgf46.00 commented out
!vjp 2/5/06 added kind parameter for ibm platform
#ifdef IBM
      itmp = int(-log10(abs(zeta)),KIND(c0_0)) !jgf46.00 added
#else
      itmp = int(-log10(abs(zeta))) !jgf46.16 added
#endif
      nsig = max0(itmp,4)
      maxfn = max0(itmax,50)
      eps = c0_0
      a = start
      b = c1_0
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::eigvss", gptl_handle)
#endif
      call zbrent (n,tri,eps,nsig,a,b,maxfn,ier)
      eigvss = b
!
      return
      end  function


      subroutine eqrt1s (d,e2,n,m,isw,ier)
      implicit none
      integer n,m,i,j,k,isw,ier,ii,k1
      real(kind=8) d(n),e2(n),err,s,tot,p,q,qp,r,delta,f,ep,dlam
      real(kind=8), parameter :: c0_0 = 0.0, c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!   modified imsl routine name   - eqrt1s
!
!-----------------------------------------------------------------------
!
!   computer            - cdc/single
!
!   latest revision     - june 1, 1980
!
!   purpose             - smallest or largest m eigenvalues of a
!                           symmetric tridiagonal matrix
!
!   usage               - call eqrt1s (d,e2,n,m,isw,ier)
!
!   arguments    d      - input vector of length n containing
!                           the diagonal elements of the matrix.  the
!                           computed eigenvalues replace the first m
!                           components of the vector d in non-
!                           decreasing sequence, while the remaining
!                           components are lost.
!                e2     - input vector of length n containing
!                           the squares of the off-diagonal elements
!                           of the matrix.  input e2 is destroyed.
!                n      - input scalar containing the order of the
!                           matrix.
!                m      - input scalar containing the number of
!                           smallest eigenvalues desired (m is
!                           less than or equal to n).
!                isw    - input scalar meaning as follows -
!                           isw=1 means that the matrix is known to be
!                             positive definite.
!                           isw=0 means that the matrix is not known
!                             to be positive definite.
!                ier    - error parameter. (output)
!                           warning error
!                             ier = 601 indicates that successive
!                               iterates to the k-th eigenvalue were not
!                               monotone increasing. the value k is
!                               stored in e2(1).
!                           terminal error
!                             ier = 602 indicates that isw=1 but matrix
!                               is not positive definite
!
!   precision/hardware  - single and double/h32
!                       - single/h36,h48,h60
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   remarks      as written, the routine computes the m smallest
!                eigenvalues. to compute the m largest eigenvalues,
!                reverse the sign of each element of d before and
!                after calling the routine. in this case, isw must
!                equal zero.
!
!   copyright           - 1980 by imsl, inc. all rights reserved.
!
!   warranty            - imsl warrants only that imsl testing has been
!                           applied to this code. no other warranty,
!                           expressed or implied, is applicable.
!
!-----------------------------------------------------------------------
!
!
!                                  srelpr = machine precision
!                                  first executable statement
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::eqrt1s", gptl_handle)
#endif
      ier = 0
      dlam = c0_0
      err = c0_0
      s = c0_0
!
!                                  look for small sub-diagonal entries
!                                  define initial shift from lower
!                                  gerschgorin bound.
!
      tot = d(1)
      q = c0_0
      j = 0
      do 30 i = 1,n
         p = q
         if (i.eq.1) go to 10
         if (p.gt.srelpr*(abs(d(i))+abs(d(i-1)))) go to 20
   10    e2(i) = c0_0
!
!                                  count if e2(i) has underflowed
!
   20    if (e2(i).eq.c0_0) j = j+1
         q = c0_0
         if (i.ne.n) q = sqrt(abs(e2(i+1)))
         tot = min(d(i)-p-q,tot)
   30 continue
      if (isw.eq.1.and.tot.lt.c0_0) go to 50
      do 40 i = 1,n
         d(i) = d(i)-tot
   40 continue
      go to 60
   50 tot = c0_0
   60 do 190 k = 1,m
!
!                                  next qr transformation
!
   70    tot = tot+s
         delta = d(n)-s
         i = n
         f = abs(srelpr*tot)
         if (dlam.lt.f) dlam = f
         if (delta.gt.dlam) go to 90
         if (delta.ge.(-dlam)) go to 160
         ier = 602
         if (level.ge.1) write (nout,80)
   80    format (/1x,'warning')
         go to 200
!
!                                  replace small sub-diagonal squares
!                                  by zero to reduce the incidence of
!                                  underflows
!
   90    if (k.eq.n) go to 110
         k1 = k+1
         do 100 j = k1,n
            if (e2(j).le.(srelpr*(d(j)+d(j-1)))**2) e2(j) = 0.0
  100    continue
  110    f = e2(n)/delta
         qp = delta+f
         p = c1_0
         if (k.eq.n) go to 140
         k1 = n-k
         do 130 ii = 1,k1
            i = n-ii
            q = d(i)-s-f
            r = q/qp
            p = p*r+c1_0
            ep = f*r
            d(i+1) = qp+ep
            delta = q-ep
            if (delta.gt.dlam) go to 120
            if (delta.ge.(-dlam)) go to 160
            ier = 602
            if (level.ge.1) write (nout,80)
            go to 200
  120       f = e2(i)/q
            qp = delta+f
            e2(i+1) = qp*ep
  130    continue
  140    d(k) = qp
         s = qp/p
         if (tot+s.gt.tot) go to 70
         ier = 601
         e2(1) = k
         if (level.ge.1) write (nout,150) k
  150    format (/1x,'warning: not monotonic',i10/1x)
!
!                                  set error -- irregular end
!                                  deflate minimum diagonal element
!
         s = c0_0
         i = ismin(n-k+1,d(k:),1)
         delta = min(qp,d(i))
!
!                                  convergence
!
  160    if (i.lt.n) e2(i+1) = e2(i)*f/qp
         if (i.eq.k) go to 180
         do 170 j = i-1,k,-1
            d(j+1) = d(j)-s
            e2(j+1) = e2(j)
  170    continue
  180    d(k) = tot
         err = err+abs(delta)
         e2(k) = err
  190 continue
      if (ier.eq.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::eqrt1s", gptl_handle)
#endif
          return
      end if
  200 continue
  210 return
      end  subroutine

      subroutine iterm (n,ndim,maxnz,coef,u,wk,imthd)
      implicit none
      integer n,i,ip,imthd,ndim,maxnz
      real(kind=8) coef(ndim,maxnz),u(n),wk(n),qtff
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!     iterm produces the iteration summary line at the end
!     of each iteration. if level .ge. 4, the latest approximation
!     to the solution will be printed.
!
! ... parameter list --
!
!          n      order of system or, for reduced system
!                    routines, order of black subsystem
!          coef   iteration matrix
!          u      solution estimate
!          wk     work array of length n
!          imthd  indicator of method
!                    imthd = 1,  jcg
!                    imthd = 2,  jsi
!                    imthd = 3,  sor
!                    imthd = 4,  ssorcg
!                    imthd = 5,  ssorsi
!                    imthd = 6,  rscg
!                    imthd = 7,  rssi
!
! ... specifications for arguments
!
!
! ... print various parameters after each iteration
!
      if (level.lt.2) return
      go to (10,100,140,170,50,10,100), imthd
   10 if (in.gt.0) go to 30
!
! ... print header for jcg and rscg
!
      write (nout,20)
   20 format (' iterations... '//)
!
! ... print summary line
!
   30 write (nout,40) in,stptst,cme,rho,gamma
   40 format (3x,i5,3x,5e13.5)
      if (level.ge.4) go to 200
!
      return
!
   50 if (in.gt.0) go to 70
!
! ... print header for ssor-si
!
      write (nout,60)
   60 format (' iterations... '//)
!
! ... print summary line
!
   70 ip = in-is
      if (imthd.eq.7) ip = 2*ip
      if (ip.lt.3) go to 80
      qtff = qt**ff
      write (nout,40) in,stptst,qa,qtff,rho,gamma
      if (level.ge.4) go to 200
      return
!
   80 write (nout,90) in,stptst,rho,gamma
   90 format (3x,i5,3x,e13.5,26x,2e13.5)
      if (level.ge.4) go to 200
      return
!
  100 if (in.gt.0) go to 120
!
! ... print header for j-si and rs-si
!
      write (nout,110)
  110 format (' iterations... '//)
!
! ... print summary line
!
  120 ip = in-is
      if (imthd.eq.7) ip = 2*ip
      if (ip.lt.3) go to 130
      qtff = qt**ff
      write (nout,40) in,stptst,qa,qtff,rho
      if (level.ge.4) go to 200
      return
!
  130 write (nout,90) in,stptst,rho
      if (level.ge.4) go to 200
      return
!
! ... print various parameters after each iteration for sor.
!
  140 if (in.gt.0) go to 160
!
! ... print header for sor
!
      write (nout,150)
  150 format (' iterations... '//)
!
! ... print summary line for sor
!
  160 continue
      write (nout,40) in,stptst,cme,omega,specr
      if (level.ge.4) go to 200
!
      return
!
! ... print various parameters after each iteration for ssor-cg.
!
  170 if (in.gt.0) go to 190
!
! ... print header for ssor-cg
!
      write (nout,180)
  180 format (' iterations... '//)
!
! ... print summary line for ssor-cg
!
  190 continue
      write (nout,40) in,stptst,specr,spr,rho,gamma
      if (level.ge.4) go to 200
      return
!
  200 if (imthd.gt.5) go to 220
      write (nout,210) in
  210 format (/1x,2x,'estimate of solution at iteration ',i5)
      go to 240
  220 write (nout,230) in
  230 format (/1x,2x,'estimate of solution at black points ',
     *   'at iteration ',i5)
  240 do 250 i = 1,n
         wk(i) = u(i)*coef(i,1)
  250 continue
      write (nout,260) (wk(i),i=1,n)
  260 format (1x,5g16.7)
      write (nout,270)
  270 format (//)
!
      return
      end  subroutine

      subroutine parcon (dtnrm,c1,c2,c3,c4,gamold,rhotmp,ibmth)
      implicit none
      integer ip,ibmth
      real(kind=8) dtnrm,c1,c2,c3,c4,gamold,rhotmp,rhoold
      real(kind=8), parameter :: c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... parcon computes acceleration parameters for conjugate gradient
!     acceleration methods.
!
! ... parameter list --
!
!          dtnrm  inner product of residuals
!          c1     output -- rho*gamma
!          c2     output -- rho
!          c3     output -- 1-rho
!          c4     output -- rho*(1-gamma)
!          gamold output -- value of gamma at preceding iteration
!          rhotmp last estimate for value of rho
!          ibmth  indicator of basic method being accelerated by cg
!                      ibmth = 1,   jacobi
!                            = 2,   reduced system
!                            = 3,   ssor
!
!     description of variables in common blocks in main routine
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::parcon", gptl_handle)
#endif
      ip = in-is
!
! ... set rhoold and gamold
!
      rhoold = rho
      gamold = gamma
!
! ... compute gamma (in+1) 
!
! ... for jacobi or reduced system cg
!
!      if (ibmth.le.2) gamma = 1.0e0/(1.0e0-dtnrm/delnnm) !jgf46.00 comm. out
      gamma = c1_0/(c1_0-dtnrm/delnnm) !jgf46.00 added
!
! ... for ssor cg
!
!      if (ibmth.eq.3) gamma = delnnm/dtnrm !jgf46.00 commented out
!
! ... compute rho (in+1)
!
      rho = c1_0
      if (ip.eq.0) go to 20
!      if (isym.eq.0) go to 10                 !jgf46.00 commented out
!      rho = 1.0e0/(1.0e0-gamma*rhotmp/delsnm) !jgf46.00 commented out
!      go to 20                                !jgf46.00 commented out
   10 rho = c1_0/(c1_0-gamma*delnnm/(gamold*delsnm*rhoold))
!
! ... compute constants c1, c2, c3, and c4
!
   20 delsnm = delnnm
      rhotmp = rhoold
      c1 = rho*gamma
      c2 = rho
      c3 = c1_0-rho
      c4 = rho*(c1_0-gamma)
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::parcon", gptl_handle)
#endif
      return
      end  subroutine



      subroutine permat (n,ndim,maxnz,jcoef,coef,p,work,iwork)
      implicit none
!
! ... permat takes the sparse matrix representation
!     of the matrix stored in the arrays jcoef and coef and
!     permutes both rows and columns, overwriting the previous
!     structure.
!
! ... parameter list --
!
!          n         order of system
!          ndim      row dimension of arrays jcoef and coef in
!                       the calling routine
!          maxnz     maximum number of nonzero entries per row
!          jcoef     integer array for data
!          coef      array for data structure coefficients
!          p         permutation vector
!          work      workspace of length n
!          iwork     integer workspace of length n
!
! ... it is assumed that the i-th entry of the permutation vector
!     p indicates the row the i-th row gets mapped into.  (i.e.
!     if ( p(i) = j ) row i gets mapped into row j)
!
!     *** note ***  this routine is to be called after routine scal.
!
! ... specifications for arguments
!
      integer ndim, maxnz, n, i, j
      integer jcoef(ndim,maxnz),p(ndim*maxnz),iwork(n)
      real(kind=8) coef(ndim,maxnz),work(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
      if (n.le.0) return
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::permat", gptl_handle)
#endif
      do 50 j = 1,maxnz
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::permat", gptl_handle)
#endif
         call scopy (n,coef(1:,j),1,work,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::permat", gptl_handle)
#endif
         do 10 i = 1,n
            iwork(i) = jcoef(i,j)
   10    continue
         do 20 i = 1,n
            coef(p(i),j) = work(i)
   20    continue
         do 30 i = 1,n
            jcoef(p(i),j) = iwork(i)
   30    continue
         do 40 i = 1,n
            jcoef(i,j) = p(jcoef(i,j))
   40    continue
   50 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::permat", gptl_handle)
#endif
      return
      end  subroutine


      subroutine pervec (n,p,v,work)
      implicit none
      integer i,n
      integer p(n)
      real(kind=8) v(n),work(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! ... pervec permutes a vector as dictated by the permutation
! ... vector p.  if p(i) = j, then v(j) gets v(i).
!
! ... parameter list --
!
!          n       length of vectors p, v, and work
!          p       integer permutation vector
!          v       vector to be permuted
!          work    workspace vector of length n
!
!
      call scopy (n,v,1,work,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::pervec", gptl_handle)
#endif
      do 10 i = 1,n
         v(p(i)) = work(i)
   10 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pervec", gptl_handle)
#endif
      return
      end  subroutine


      subroutine pjac (n,ndim,maxnz,jcoef,coef,u,rhs)
      implicit none
!      integer n,ndim,maxnz,maxm1 !jgf46.00 commented out
      integer i,j,n,ndim,maxnz    !jgf46.00 added
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz),u(n),rhs(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::pjac", gptl_handle)
#endif
!
! ... pjac performs one jacobi iteration.
!
! ... parameter list --
!
!         n       dimension of matrix
!         ndim    row dimension of jcoef and coef arrays in calling
!                   routine
!         maxnz   maximum number of nonzeros per row
!         jcoef   integer data structure for coefficient columns
!         coef    data structure for array coefficients
!         u       estimate of solution of a matrix problem
!         rhs     on input -- contains the right hand side of the
!                             matrix problem
!                 on output -- contains b*u + rhs  where b = i - a
!                              and a has been scaled to have a unit
!                              diagonal
!
!
!      maxm1 = maxnz-1  !jgf46.00 commented out
!      call ymasx2 (ndim,n,maxm1,coef(1,2),jcoef(1,2),rhs,u)!jgf46.00 comm.out
!     jgf46.00 Begin add.
      do j = 2,maxnz ! jgf50.03: st3 from maxnz-1, Bug fix, 06.02.2010 
      do i = 1,n
         rhs(i) = rhs(i) - coef(i,j)*u(jcoef(i,j))
      enddo
      enddo
!     jgf46.00 End add.
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pjac", gptl_handle)
#endif
      return
      end  subroutine


      subroutine pmult (n,ndim,maxnz,jcoef,coef,b,c)
      implicit none
!      integer n,ndim,maxnz,maxm1 !jgf46.00 commented out
      integer i,j,n,ndim,maxnz    !jgf46.00 added
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz),b(n),c(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! ... pmult computes c = a*b, a matrix-vector product.  matrix
!     a is assumed to be stored in the coef, jcoef ellpack
!     data structure and all entries in the column array jcoef
!     are assumed to be between 1 and n, inclusive.
!     a is assumed to have a unit diagonal.
!
! ... parameter list --
!
!          n        dimension of matrix
!          ndim     row dimension of coef and jcoef in calling routine
!          maxnz    maximum number of nonzeros per row
!          jcoef    integer array for coefficient columns
!          coef     array for coefficients
!          b        multiplying vector of length n
!          c        product vector of length n
!
!
      call scopy (n,b,1,c,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::pmult", gptl_handle)
#endif
!      maxm1 = maxnz-1 !jgf46.00 commented out
!      call ypasx2 (ndim,n,maxm1,coef(1,2),jcoef(1,2),c,b) !jgf46.00 comm.out
!     jgf46.00 Begin add.
      do j = 2,maxnz !jgf50.03: st3 from maxnz-1, Bug fix, 06.02.2010
      do i = 1,n
         c(i) = c(i) - coef(i,j)*b(jcoef(i,j))
      enddo
      enddo
!     jgf46.00 End add.
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::pmult", gptl_handle)
#endif
      return
      end  subroutine

      subroutine prbndx (n,ndim,maxnz,jcoef,p,ip,nblack,level,nout,ier)
      implicit none
      integer i,n,ndim,maxnz,nblack,level,nout,ier,ibgn,next,last,
     *  j,k,nxttyp,jcol,l,nred
      integer jcoef(ndim,maxnz),p(n),ip(n)
      integer first,old,young,curtyp,type
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::prbndx", gptl_handle)
#endif
!
!**************************************************************
!
!     prbndx computes the red-black permutation
!     vectors p ( and its inverse ip ) if possible.
!
!     the algorithm is to mark the first node as red (arbitrary).
!     all of its adjacent nodes are marked black and placed in
!     a stack.  the remainder of the code pulls the first node
!     off the top of the stack and tries to type its adjacent nodes.
!     the typing of the adjacent point is a five way case statement
!     which is well commented below (see do loop 100).
!
!     the array p is used both to keep track of the color of a node
!     (red node is positive, black is negative) but also the father
!     node that caused the color marking of that point.  since
!     complete information on the adjacency structure is hard to come
!     by this forms a link to enable the color change of a partial
!     tree when a recoverable color conflict occurs.
!
!     the array ip is used as a stack to point to the set of nodes
!     left to be typed that are known to be adjacent to the current
!     father node.
!
!     *** note ***  this routine is to be called after routine scal.
!
!*********************************************************************
!
!     input parameters --
!
!        n      number of nodes.  (integer, scalar)
!
!        ndim   row dimension of jcoef in calling routine.
!
!        maxnz  maximum number of nonzeros per row
!
!        jcoef  array of column indices.  it is assumed
!               that for every row where only one element is
!               stored that element corresponds to the diagonal
!               entry.  the diagonal must be the first entry stored.
!                 (integer, arrays)
!
!        level  switch for printing
!
!        nout   output tape number
!
!     output parameters --
!
!        nblack number of black nodes.  number of red nodes is
!               n - nblack.  (integer, scalar)
!
!        p, ip  permutation and inverse permutation vectors.
!               (integer, arrays each of length n)
!
!        ier    error flag. (integer, scalar)
!
!               ier = 0, normal return.  indexing performed
!                        successfully
!               ier = 201, red-black indexing not possible.
!
!********************************************************************
!
!
      ier = 0
      if (n.le.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::prbndx", gptl_handle)
#endif
            return
      endif
      do 10 i = 1,n
         p(i) = 0
         ip(i) = 0
   10 continue
!
! ... handle the first set of points until some adjacent points
! ... are found
!
      first = 1
!
   20 p(first) = first
      if (maxnz.gt.1) go to 40
!
! ... search for next entry that has not been marked
!
      if (first.eq.n) go to 120
      ibgn = first+1
      do 30 i = ibgn,n
         if (p(i).ne.0) go to 30
         first = i
         go to 20
   30 continue
      go to 120
!
! ... first set of adjacent points found
!
   40 next = 1
      last = 1
      ip(1) = first
!
! ... loop over labeled points indicated in the stack stored in
! ... the array ip
!
   50 k = ip(next)
      curtyp = p(k)
      nxttyp = -curtyp
      do 100 j = 2,maxnz
         jcol = jcoef(k,j)
         if (jcol.eq.k) go to 100
         type = p(jcol)
!
!==================================================================
!
!     the following is a five way case statement dealing with the
!     labeling of the adjacent node.
!
! ... case i.  if the adjacent node has already been labeled with
!              label equal to nxttyp, then skip to the next adjacent
!              node.
!
         if (type.eq.nxttyp) go to 100
!
! ... case ii.  if the adjacent node has not been labeled yet label
!               it with nxttyp and enter it in the stack
!
         if (type.ne.0) go to 60
         last = last+1
         ip(last) = jcol
         p(jcol) = nxttyp
         go to 100
!
! ... case iii.  if the adjacent node has already been labeled with
!                opposite color and the same father seed, then there
!                is an irrecoverable color conflict.
!
   60    if (type.eq.curtyp) go to 140
!
! ... case iv.  if the adjacent node has the right color and a different
!               father node, then change all nodes of the youngest fathe
!               node to point to the oldest father seed and retain the
!               same colors.
!
         if (type*nxttyp.lt.1) go to 80
         old = min0(iabs(type),iabs(nxttyp))
         young = max0(iabs(type),iabs(nxttyp))
         do 70 l = young,n
            if (iabs(p(l)).eq.young) p(l) = isign(old,p(l))
   70    continue
         curtyp = p(k)
         nxttyp = -curtyp
         go to 100
!
! ... case v.  if the adjacent node has the wrong color and a different
!              father node, then change all nodes of the youngest father
!              node to point to the oldest father node along with
!              changing their colors.  since until this time the
!              youngest father node tree has been independent no other
!              color conflicts will arise from this change.
!
   80    old = min0(iabs(type),iabs(nxttyp))
         young = max0(iabs(type),iabs(nxttyp))
         do 90 l = young,n
            if (iabs(p(l)).eq.young) p(l) = isign(old,-p(l))
   90    continue
         curtyp = p(k)
         nxttyp = -curtyp
!
! ... end of case statement
!
!==================================================================
!
  100 continue
!
! ... advance to next node in the stack
!
      next = next+1
      if (next.le.last) go to 50
!
! ... all nodes in the stack have been removed
!
! ... check for nodes not labeled.  if any are found
! ... start the labeling process again at the first
! ... node found that is not labeled.
!
      ibgn = first+1
      do 110 i = ibgn,n
         if (p(i).ne.0) go to 110
         first = i
         go to 20
  110 continue
!
!===================================================================
!
! ... all nodes are now typed either red or black
!
! ... generate permutation vectors
!
  120 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::prbndx", gptl_handle)
#endif
      call whenige (n,p,1,0,ip,nred)
      call whenilt (n,p,1,0,ip(nred+1),nblack)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::prbndx", gptl_handle)
#endif
      do 130 i = 1,n
         p(ip(i)) = i
  130 continue
!
! ... successful red-black ordering completed
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::prbndx", gptl_handle)
#endif
      return
!
! ...... type conflict
!
  140 ier = 201
      if (level.ge.0) write (nout,150)
  150 format ('error: red-black indexing not possible'//)
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::prbndx", gptl_handle)
#endif
      return
      end  subroutine

      subroutine sbelm (n,ndim,maxnz,jcoef,coef,rhs,work,tol)
      implicit none
      integer i,j,jcol,n,ndim,maxnz
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz),rhs(n),work(n),tol
      real(kind=8), parameter :: c0_0 = 0.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::sbelm", gptl_handle)
#endif
!
! ... sbelm is designed to remove rows of the matrix for which
! ... all off-diagonal elements are very small (less than tol).
! ... this is to take care of matrices arising from finite
! ... element discretizations of partial differential equations
! ... with dirichlet boundary conditions.  any such rows and
! ... corresponding columns are then eliminated (set to the
! ... identity after correcting the rhs).
! ... *** note ***  this routine is to be called after routine scal.
!
! ... parameter list --
!
!         n       dimension of matrix
!         ndim    row dimension of arrays jcoef and coef in the
!                    calling program
!         maxnz   maximum number of nonzero entries per row
!         jcoef   integer array of matrix representation
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         work    work array of length n
!         tol     tolerance factor
!
!
      if (n.le.0.or.maxnz.lt.2) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::sbelm", gptl_handle)
#endif
            return
      end if
!
! ... find maximum off-diagonal elements in absolute value.
!
      do 10 i = 1,n
         work(i) = c0_0
   10 continue
      do 30 j = 2,maxnz
         do 20 i = 1,n
            work(i) = max(work(i),abs(coef(i,j)))
   20    continue
   30 continue
!
! ... eliminate desired rows and columns.
!
      do 60 j = 2,maxnz
         do 50 i = 1,n
            if (work(i).lt.tol) go to 40
            jcol = jcoef(i,j)
            if (work(jcol).ge.tol) go to 50
            rhs(i) = rhs(i)-coef(i,j)*rhs(jcol)
   40       coef(i,j) = c0_0
            jcoef(i,j) = i
   50    continue
   60 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::sbelm", gptl_handle)
#endif
      return
      end  subroutine

      subroutine scal (n,ndim,maxnz,jcoef,coef,rhs,u,work,ier)
      implicit none
      integer i,j,n,ndim,maxnz,nsgncg,ier
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz),rhs(n),u(n),work(n),save
      real(kind=8), parameter :: c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
! ... scal scales original matrix to a unit diagonal matrix.  rhs
! ... and u vectors are scaled accordingly.  the data
! ... structure is adjusted to have diagonal entries in
! ... column 1.  zero entries in jcoef array are changed to
! ... positive integers between 1 and n.
!
! ... parameter list --
!
!         n       dimension of matrix
!         ndim    row dimension of arrays jcoef and coef in the
!                    calling program
!         maxnz   maximum number of nonzero entries per row
!         jcoef   integer array of matrix representation
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         work    work array of length n
!         ier     error flag -- on return, nonzero values mean
!                    401 -- zero diagonal element
!                    402 -- nonexistent diagonal element
!
!
!     description of variables in common block in main routine
!
! ... check for positive diagonal entries for each row.
! ... put diagonal entries in column 1.  replace zeros in
! ... row i of jcoef with i.
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::scal", gptl_handle)
#endif
      ier = 0
      nsgncg = 0
      if (n.le.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scal", gptl_handle)
#endif
            return
      endif
      do 110 i = 1,n
         if (jcoef(i,1).eq.i) go to 50
         if (maxnz.lt.2) go to 20
         do 10 j = 2,maxnz
            if (jcoef(i,j).eq.i) go to 40
   10    continue
!
! ... fatal error -- no diagonal entry for row i.
!
   20    ier = 402
         if (level.ge.0) write (nout,30) i
   30    format (//1x,'error: no diagonal entry in row',i10)
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scal", gptl_handle)
#endif
         return
!
! ... shift row i so that diagonal element is in column 1.
!
   40    save = coef(i,j)
         coef(i,j) = coef(i,1)
         jcoef(i,j) = jcoef(i,1)
         coef(i,1) = save
         jcoef(i,1) = i
!
! ... check sign of diagonal entry.  if negative, change signs of
! ... all row coefficients and corresponding rhs element.
!
   50    if (coef(i,1)) 60 , 90 , 110
   60    do 70 j = 1,maxnz
            coef(i,j) = -coef(i,j)
   70    continue
         rhs(i) = -rhs(i)
         nsgncg = nsgncg+1
         if (level.ge.5) write (nout,80) i
   80    format (//1x,'note: equation ',i10,' has been negated')
         go to 110
!
! ... fatal error -- zero diagonal element for row i.
!
   90    ier = 401
         if (level.ge.0) write (nout,100) i
  100    format (//1x,'error: diagonal entry in row ',i10,' is zero')
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scal", gptl_handle)
#endif
         return
  110 continue
!
! ... change zero elements of jcoef array.
!
      if (maxnz.lt.2) go to 140
      do 130 j = 2,maxnz
         do 120 i = 1,n
            if (jcoef(i,j).le.0) jcoef(i,j) = i
  120    continue
  130 continue
!
! ... scale rhs and u arrays.  store reciprocal square roots
! ... of diagonal entries in column 1 of coef.
!
  140 do 150 i = 1,n
         work(i) = sqrt(coef(i,1))
  150 continue
      do 160 i = 1,n
         u(i) = u(i)*work(i)
  160 continue
      do 170 i = 1,n
         work(i) = c1_0/work(i)
  170 continue
      call scopy (n,work,1,coef(:,1),1)
      do 180 i = 1,n
         rhs(i) = rhs(i)*work(i)
  180 continue
!
! ... scale matrix.
!
      if (maxnz.lt.2) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scal", gptl_handle)
#endif
      return
      end if
      do 200 j = 2,maxnz
         do 190 i = 1,n
            coef(i,j) = coef(i,j)*work(i)*work(jcoef(i,j))
  190    continue
  200 continue
!
! ... adjust isym if the  0 .lt. nsgncg .lt. n
!
      if (nsgncg.gt.0.and.nsgncg.lt.n) isym = 1
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scal", gptl_handle)
#endif
      return
      end  subroutine


      subroutine unscal (n,ndim,maxnz,jcoef,coef,rhs,u,work)
      implicit none
      integer i,j,n,ndim,maxnz
      integer jcoef(ndim,maxnz)
      real(kind=8) coef(ndim,maxnz),rhs(n),u(n),work(n)
      real(kind=8), parameter :: c1_0 = 1.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! ... unscal reverses the scaling done in routine scal.
!
! ... parameter list --
!
!         n       dimension of matrix
!         ndim    row dimension of arrays jcoef and coef in the
!                    calling program
!         maxnz   maximum number of nonzero entries per row
!         jcoef   integer array of matrix representation
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         work    work array of length n
!
!
! ... unscale u and rhs arrays.
!
      call scopy (n,coef(:,1),1,work,1)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::unscal", gptl_handle)
#endif
      do 10 i = 1,n
         u(i) = u(i)*work(i)
   10 continue
      do 20 i = 1,n
         work(i) = c1_0/work(i)
   20 continue
      do 30 i = 1,n
         rhs(i) = rhs(i)*work(i)
   30 continue
!
! ... unscale matrix.
!
      if (maxnz.lt.2) go to 80
      do 50 j = 2,maxnz
         do 40 i = 1,n
            coef(i,j) = coef(i,j)*work(i)*work(jcoef(i,j))
   40    continue
   50 continue
!
! ... put original zeros back in icoef array and restore unscaled
! ... diagonal entries to column one.
!
      do 70 j = 2,maxnz
         do 60 i = 1,n
            if (jcoef(i,j).eq.i) jcoef(i,j) = 0
   60    continue
   70 continue
   80 do 90 i = 1,n
         coef(i,1) = work(i)**2
   90 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::unscal", gptl_handle)
#endif
      return
      end  subroutine


      subroutine vout (n,v,iswt,nout)
      implicit none
!
!     vout effects printing of residual and solution
!     vectors - called from peror
!
! ... parameter list --
!
!          v      vector of length n
!          iswt   labelling information
!          nout   output device number
!
! ... specifications for arguments
!
      integer n, iswt, nout, k, kupper, i, j, jm1
      real(kind=8)  v(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::vout", gptl_handle)
#endif
!
!        if (n .le. 0) return
!
      kupper = min0(n,4)
      if (iswt.eq.1) write (nout,10)
   10 format (//5x,'residual vector')
      if (iswt.eq.2) write (nout,20)
   20 format (//5x,'solution vector')
      write (nout,30) (i,i=1,kupper)
   30 format (10x,4i15)
      write (nout,40)
   40 format (10x,65('-')/)
!
      do 60 j = 1,n,4
         kupper = min0(j+3,n)
         jm1 = j-1
         write (nout,50) jm1,(v(k),k=j,kupper)
   50    format (4x,i5,'+  ',4e15.5)
   60 continue
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::vout", gptl_handle)
#endif
      return
      end  subroutine

      subroutine zbrent (n,tri,eps,nsig,a,b,maxfn,ier)
      USE SIZES
      implicit none
      integer n,ier,nsig,maxfn,ib3,ib4,ic
      real(kind=8) a,b,c,d,e,eps,p,q,r,s,t,fa,fb,fc,tol,rm,rone,temp
      real(kind=8) tri(3*n)
      real(kind=8), parameter :: c10_0 = 10.0, c0_0 = 0.0, c0_1 = 0.1, c0_5 = 0.5, c1_0 = 1.0, c3_0 = 3.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
! *** begin -- itpackv common
!
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
!
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
!
      real(kind=8) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
!
! *** end   -- itpackv common
!
!   modified imsl routine name   - zbrent
!
!-----------------------------------------------------------------------
!
!   computer            - cdc/single
!
!   latest revision     - january 1, 1978
!
!   purpose             - zero of a function which changes sign in a
!                           given interval (brent algorithm)
!
!   usage               - call zbrent (f,eps,nsig,a,b,maxfn,ier)
!
!   arguments    tri    - a tridiagonal matrix of order n
!                eps    - first convergence criterion (input).  a root,
!                           b, is accepted if abs(f(b)) is less than or
!                           equal to eps.  eps may be set to zero.
!                nsig   - second convergence criterion (input).  a root,
!                           b, is accepted if the current approximation
!                           agrees with the true solution to nsig
!                           significant digits.
!                a,b    - on input, the user must supply two points, a
!                           and b, such that f(a) and f(b) are opposite
!                           in sign.
!                           on output, both a and b are altered.  b
!                           will contain the best approximation to the
!                           root of f. see remark 1.
!                maxfn  - on input, maxfn should contain an upper bound
!                           on the number of function evaluations
!                           required for convergence.  on output, maxfn
!                           will contain the actual number of function
!                           evaluations used.
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 501 indicates the algorithm failed to
!                             converge in maxfn evaluations.
!                           ier = 502 indicates f(a) and f(b) have the
!                             same sign.
!
!   precision/hardware  - single and double/h32
!                       - single/h36,h48,h60
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   remarks  1.  let f(x) be the characteristic function of the matrix
!                tri evaluated at x. function determ evaluates f(x).
!                on exit from zbrent, when ier=0, a and b satisfy the
!                following,
!                f(a)*f(b) .le. 0,
!                abs(f(b)) .le. abs(f(a)), and
!                either abs(f(b)) .le. eps or
!                abs(a-b) .le. max(abs(b),0.1)*10.0**(-nsig).
!                the presence of 0.1 in this error criterion causes
!                leading zeroes to the right of the decimal point to be
!                counted as significant digits. scaling may be required
!                in order to accurately determine a zero of small
!                magnitude.
!            2.  zbrent is guaranteed to reach convergence within
!                k = (log((b-a)/d)+1.0)**2 function evaluations where
!                  d=min(over x in (a,b) of
!                    max(abs(x),0.1)*10.0**(-nsig)).
!                this is an upper bound on the number of evaluations.
!                rarely does the actual number of evaluations used by
!                zbrent exceed sqrt(k). d can be computed as follows,
!                  p = amin1(abs(a),abs(b))
!                  p = amax1 (0.1,p)
!                  if ((a-0.1)*(b-0.1).lt.0.0) p = 0.1
!                  d = p*10.0**(-nsig)
!
!   copyright           - 1977 by imsl, inc. all rights reserved.
!
!   warranty            - imsl warrants only that imsl testing has been
!                           applied to this code. no other warranty,
!                           expressed or implied, is applicable.
!
!-----------------------------------------------------------------------
!
!
!     description of variables in common block in main routine
!
! ... local itpackv references --
!
!          determ
!
!                                  first executable statement
!
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      ier = 0
      ib3 = 2*itmax+1
      ib4 = 3*itmax+1
      t = c10_0**(-nsig) !jgf46.00 added d0
      ic = 2
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      fa = determ(itmax,n,tri,tri(ib3:ib3+n),tri(ib4:ib4+n),a)
      fb = determ(itmax,n,tri,tri(ib3:ib3+n),tri(ib4:ib4+n),b)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      s = b
!
!                                  test for same sign
!
!      if (fa*fb.gt.0.0) go to 110 jgf46.00 commented out
!     jgf46.00 Begin add.
! cvjp 2/10/06  replaced because produced NaNs
!     if (fa*fb.gt.0.0e0) go to 110
            
      if ( (fa > c0_0 .and. fb > c0_0) .or.
     &     (fa < c0_0 .and. fb < c0_0) ) then
        go to 110
      endif
!     jgf46.00 End add.

   10 c = a
      fc = fa
      d = b-c
      e = d
   20 if (abs(fc).ge.abs(fb)) go to 30
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   30 continue

#ifdef REAL8
      tol = t*max(abs(b),c0_1)
#else
      tol = t*max(abs(b),c0_1)
#endif

      rm = (c-b)*c0_5
!
!                                  test for first convergence criteria
!
      if (abs(fb).le.eps) go to 80
!
!                                  test for second convergence criteria
!
      if (abs(c-b).le.tol) go to 80
!
!                                  check evaluation counter
!
      if (ic.ge.maxfn) go to 90
!
!                                  is bisection forced
!
      if (abs(e).lt.tol) go to 60
      if (abs(fa).le.abs(fb)) go to 60
      s = fb/fa
      if (a.ne.c) go to 40
!
!                                  linear interpolation
!
      p = (c-b)*s
      q = c1_0-s
      go to 50
!
!                                  inverse quadratic interpolation
!
   40 q = fa/fc
      r = fb/fc
      rone = r-c1_0
      p = s*((c-b)*q*(q-r)-(b-a)*rone)
      q = (q-c1_0)*rone*(s-c1_0)
   50 if (p.gt.c0_0) q = -q
      if (p.lt.c0_0) p = -p
      s = e
      e = d
!
!                                  if abs(p/q).ge.75*abs(c-b) then
!                                     force bisection
!
      if (p+p.ge.c3_0*rm*q) go to 60
!
!                                  if abs(p/q).ge..5*abs(s) then force
!                                     bisection. s = the value of p/q
!                                     on the step before the last one
!
      if (p+p.ge.abs(s*q)) go to 60
      d = p/q
      go to 70
!
!                                  bisection
!
   60 e = rm
      d = e
!
!                                  increment b
!
   70 a = b
      fa = fb
      temp = d
!jjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
      if (abs(temp).le.c0_5*tol) temp = sign(0.5e0*tol,rm)
#else
      if (abs(temp).le.c0_5*tol) temp = sign(0.5e0*tol,rm)
#endif
      b = b+temp
      s = b
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      fb = determ(itmax,n,tri,tri(ib3:ib3+n),tri(ib4:ib4+n),s)
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      ic = ic+1
!      if (fb*fc.le.0.0) go to 20 !jgf46.00 commented out
!     jgf46.00 Begin add. 
!vjp 2/10/06 replace because produced NaNs
!     if (fb*fc.le.0.0e0) go to 20

      if (fb == c0_0 .or. fc == c0_0) goto 20

      if ( (fb > c0_0 .and. fc < c0_0) .or.
     &     (fb < c0_0 .and. fc > c0_0) ) then
        go to 20
      endif
!     jgf46.00 End add.

      go to 10
!
!                                  convergence of b
!
   80 a = c
      maxfn = ic
      go to 130
!
!                                  maxfn evaluations
!
   90 ier = 501
      a = c
      maxfn = ic
      if (level.ge.1) write (nout,100) maxfn
  100 format (/1x,'warning:  algorithm failed to converge in',i6,' iterations ')
      go to 130
!
!                                  terminal error - f(a) and f(b) have
!                                  the same sign
!
  110 ier = 502
      maxfn = ic
      if (level.ge.1) write (nout,120)
  120 format (/1x,'warning: in itpackv routine zbrent f(a) and f(b) have same sign   ')
  130 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::zbrent", gptl_handle)
#endif
      return
      end  subroutine

!     jgf46.00 deleted the following subroutine
!      subroutine ypasx2 (ndim,n,m,a,ja,y,x)

!     jgf46.00 deleted the following subroutine
!      subroutine ymasx2 (ndim,n,m,a,ja,y,x)


      subroutine whenige (n,p,inc,itarg,ip,npt)
      implicit none
      integer n,inc,itarg,npt,i
      integer p(n), ip(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::whenige", gptl_handle)
#endif
!
      npt = 0
      do 10 i = 1,n
         if (p(i) .lt. itarg) go to 10
         npt = npt + 1
         ip(npt) = i
 10   continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::whenige", gptl_handle)
#endif
      return
      end subroutine


      subroutine whenilt (n,p,inc,itarg,ip,npt)
      implicit none
      integer n,inc,itarg,npt,i
      integer p(n), ip(n)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::whenilt", gptl_handle)
#endif
!
      npt = 0
      do 10 i = 1,n
         if (p(i) .ge. itarg) go to 10
         npt = npt + 1
         ip(npt) = i
 10   continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::whenilt", gptl_handle)
#endif
      return
      end subroutine


!
!*****************************************************************
!
!  Modified BLAS for ITPACK2D
!
!*****************************************************************
!

      real(kind=8) function sdot(n,sx,incx,sy,incy)
      implicit none
      integer i,n,incx,incy
      real(kind=8) sx(n),sy(n)
      REAL(kind=8) ddot
      real(kind=8), parameter :: c0_0 = 0.0
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::sdot", gptl_handle)
#endif
!
      ddot = c0_0
      sdot = ddot
!
      if (n.gt.0) then  
        do i = 1,n
           ddot = ddot + sx(i)*sy(i)
        enddo 
        sdot = ddot
      endif
!
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::sdot", gptl_handle)
#endif
      return
      end function

      integer function ismin (n,sx,incx)
      implicit none
      integer n,incx,ns,i,ii
      real(kind=8) sx(n*incx),smin,xval
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
!
!     find smallest index of minimum value of single precision sx.
!
      ismin = 0
      if (n.le.0) return
      ismin = 1
      if (n.le.1) return
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::ismin", gptl_handle)
#endif
      if (incx.eq.1) go to 30
!
!        code for increments not equal to 1.
!
      smin = sx(1)
      ns = n*incx
      ii = 1
      do 20 i = 1,ns,incx
         xval = sx(i)
         if (xval.ge.smin) go to 10
         ismin = ii
         smin = xval
   10    ii = ii+1
   20 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::ismin", gptl_handle)
#endif
      return
!
!        code for increments equal to 1.
!
   30 smin = sx(1)
      do 40 i = 2,n
         xval = sx(i)
         if (xval.ge.smin) go to 40
         ismin = i
         smin = xval
   40 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::ismin", gptl_handle)
#endif
      return
      end function



      subroutine  scopy(n,sx,incx,sy,incy)
      implicit none
      integer n,incx,incy,ix,iy,i,m,mp1,ns
      real(kind=8) sx(n*incx),sy(n*incy)
#ifdef GPTL
      integer :: gptl_ret, gptl_handle = 0
#endif
#ifdef GPTL
      gptl_ret = gptlstart_handle("::ITPACKV::scopy", gptl_handle)
#endif
!
!     copy single precision sx to single precision sy.
!
      if(n.le.0) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scopy", gptl_handle)
#endif
            return
      end if
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scopy", gptl_handle)
#endif
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop so remaining vector length is a multiple of 7.
!
   20 m = n - (n/7)*7
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) then
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scopy", gptl_handle)
#endif
      return
      end if
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scopy", gptl_handle)
#endif
      return
!
!        code for equal, positive, nonunit increments.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sx(i)
   70     continue
#ifdef GPTL
      gptl_ret = gptlstop_handle("::ITPACKV::scopy", gptl_handle)
#endif
      return
      end subroutine


      END MODULE ITPACKV
