      program CW263
c
c   Stream function wave theory code:
c   * automatic selection of order
c   * uniform current
c   * returns proportion of limiting height
c   j.r.chaplin@soton.ac.uk
c
c   Subroutine cw260 solves the wave.  Kinematics are then available
c   through subtroutine kmts as illustrated below in the main program.
c   Stores results in CW263.PSI in a format compatible with that used
c      by CW6.FOR (24/3/99)
c
      character ans*1
c
   10 write(*,'(a)') ' Water depth (m)   = '
      read*,d
      write(*,'(a)') ' Period      (s)   = '
      read*,t
      write(*,'(a)') ' Wave height (m)   = '
      read*,h
      write(*,'(a)') ' Current     (m/s) = '
      read*,u
      nverb= 1
c
      call cw260(d,t,h,u,nverb,n,el)
c
   20 write(*,'(/a)')
     :  ' (H)orizontal, (V)ertical, (S)urface, (N)ew wave, (Q)uit : '
      read(*,'(a)') ans
      goto (20,21,22,23,10,29) (index('HhVvSsNnQq',ans)+3)/2
c
c   Horizontal
c
   21 write(*,'(a)') ' y (m) = '
      read*,yy
      npt= 21
      write(*,3)
      do 31 i=1,npt
      xx= el*(i-1)/float(npt-1)
      tt= 0.0
      call kmts(n,xx,yy,tt,uu,vv,ut,vt,du,dv,etah)
      ans= ' '
      if (yy.gt.etah) ans='*'
   31 write(*,'(f9.3,8f8.3,1x,a)') xx,yy,uu,vv,ut,vt,du,dv,etah,ans
      goto 20
c
c   Vertical
c
   22 write(*,'(a)') ' x/L = '
      read*,xl
      xx= xl*el
      yy= 0.0
      tt= 0.0
      call kmts(n,xx,yy,tt,uu,vv,ut,vt,du,dv,etah)
      npt= 21
      write(*,3)
      do 32 i=1,npt
      yy= (d+etah)*(npt-i)/float(npt-1)-d
      call kmts(n,xx,yy,tt,uu,vv,ut,vt,du,dv,etah)
   32 write(*,'(f9.3,8f8.3,1x,1a)') xx,yy,uu,vv,ut,vt,du,dv,etah
      goto 20
c
c   Free surface
c
   23 npt= 21
      write(*,3)
      do 33 i=1,npt
      xx= el*(i-1)/float(npt-1)
      yy= h
      tt= 0.0
      call kmts(n,xx,yy,tt,uu,vv,ut,vt,du,dv,etah)
   33 write(*,'(f9.3,8f8.3)') xx,etah,uu,vv,ut,vt,du,dv,etah
      goto 20
c
   29 stop
    3 format('      x       y       u       v      ut   ',
     :  '   vt      du      dv      eta')
      end
c
c
c
      subroutine cw260(zd,zt,zh,zu,nverb,nfun,zel)
c
c  Input:  zd=depth; zt=period; zh=height; u=current; nverb=verbosity
c  Output: nfun=order; zel=wavelength
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      real zd,zt,zh,zu,zel
      double precision k
      character itl*79,datim*22
      common /one/ d,t,h,u,k
     :       /two/ eta(nmax),c(nmax),amp(0:nmax)
c
      pi= 4*atan(1.0)
c
      d= zd
      t= zt
      hw= zh
      u= zu
c
c  Get required solution order
c
      call wavecel(t,d,u,tr,cel)
      call limit(hw,d,tr,rat,1)
      dl0= d/(9.81*t*t/(2*pi))
      a= 0.86/sqrt(dl0)
      b= 7+2.2*log(dl0)
      cc= 2.7-3*log(dl0)
      nve= nint((a+b*rat+cc*rat**2)/2)*2
      nve= nve+2
      hb= hw/rat
c
      k= 2*pi/(cel*t)
c
c  Start with height=<Hb/2
c    First step up orders; then heights ...
c
      if (nverb.ne.0) write(*,'(/2a)')
     :   '      d       T       H       U    order  iter  ',
     :   'rms error  code      L'
      nfun = 6
      h= min(hw,0.5*hb)
      do 10 i=0,nmax
   10 amp(i)= 0.0
      amp(1)= h/2
   12 continue
      call cw261(nfun,iter,fsumsq,ifail)
      if (ifail.eq.0) then
         el= 2*pi/k
         if (nverb.ne.0) write(*,'(a,4f8.3,i5,i7,1p1e12.3,0p,i5,f10.3)')
     :      ' ',d,t,h,u,nfun,iter,fsumsq,ifail,el
      else
         if (nverb.ne.0) write(*,'(a,4f8.3,i5,i7,1p1e12.3,0p,i5)')
     :      ' ',d,t,h,u,nfun,iter,fsumsq,ifail
         stop
      endif
      if (nfun.lt.nve) then
         nfun= nfun+2
         goto 12
      endif
c
      if (hw.gt.0.5*hb) then
         fac= 1.1
   11    hm= h
         if (h*fac.gt.hw) then
            h= hw
            ilast= 1
         else
            h= h*fac
            ilast= 0
         endif
         do 14 i=1,nfun-1
   14    amp(i)= (h/hm)*amp(i)
         call cw261(nfun,iter,fsumsq,ifail)
         if (ifail.eq.0) then
            el= 2*pi/k
            if (nverb.ne.0)
     :         write(*,'(a,4f8.3,i5,i7,1p1e12.3,0p,i5,f10.3)')
     :         ' ',d,t,h,u,nfun,iter,fsumsq,ifail,el
         else
            if (nverb.ne.0) write(*,'(a,4f8.3,i5,i7,1p1e12.3,0p,i5)')
     :         ' ',d,t,h,u,nfun,iter,fsumsq,ifail
            stop
         endif
         if (ilast.ne.1) then
            fac= 0.995*fac
            goto 11
         endif
      endif
      zel= 2*pi/k
c
c$$$      open(12,file='cw263.psi')
c$$$      itl= 'Solved by CW263'
c$$$      datim= 'Space for date & time:'
c$$$      ver= 5.01
c$$$      el0= 9.81*zt**2/(2*pi)
c$$$      hl0= zh/el0
c$$$      dl0= zd/el0
c$$$      wl0= zel/el0
c$$$      write(12,101) itl,datim,ver,hl0,dl0,wl0,zh,zd,zt,nfun-1,zu
c$$$      write(12,102) (eta(i),c(i+1),i=1,nfun-1),eta(nfun)
c$$$      write(12,107) (amp(i)/zh,i=0,nfun-1)
c$$$      close(12)
c
  101 format(a/a,f10.2/3f16.10/3f16.10,i5,f16.10)
  102 format(1p2e25.16)
  107 format(1p1e25.16)
c
      return
      end
c
c
c
      subroutine cw261(nfun,iter,fsumsq,ifail)
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      double precision k
      dimension x(nmax),f(nmax),etas(2*nmax),bb(0:nmax)
      common /one/ d,t,h,u,k
     :       /two/ eta(nmax),c(nmax),amp(0:nmax)
c
      pi= 4*atan(1.0)
c
      do 10 i=1,nfun-2
      th= (i-1)*pi/(nfun-1)
      x(i)= 0.0
      do 10 j=1,nfun-1
   10 x(i)= x(i)+cos(j*th)*amp(j)
      x(nfun-1)= k
      jverb= 0
      call gaf(nfun,nfun-1,x,f,fsumsq,jverb,iter,ifail)
      if (ifail.eq.1) return
      k= x(nfun-1)
      do 11 i=1,nfun
   11 etas(i)= eta(i)
      do 12 i=1,nfun-1
   12 etas(nfun+i)= eta(nfun-i)
      call four(etas,2*nfun-2,amp,bb,nfun-1)
      amp(nfun)= 0.0
      bb(nfun)= 0.0
      return
      end
c
c
c
      subroutine lsfun(nfun,x,f,jac,sq,ifail)
c
c  Computes function errors and Jacobian
c  x(1:nfun-1) = unknowns,  f(1:nfun) = functions
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      double precision k,jac,ks
      dimension x(nmax),f(nmax),jac(nmax,nmax)
      dimension dex(nmax,nmax),th(nmax),s(nmax),
     :	chcs(nmax,nmax),chsn(nmax,nmax),shcs(nmax,nmax),shsn(nmax,nmax),
     :	fu(nmax,nmax),df(nmax,nmax),a(nmax),da(nmax,nmax),
     :  dax(nmax,nmax),dcx(nmax,nmax),dfx(nmax,nmax),
     :	us(nmax),vs(nmax),dux(nmax,nmax),dvx(nmax,nmax)
      common /one/ d,t,h,u,k
     :       /two/ eta(nmax),c(nmax),amp(0:nmax)
c
      pi= 4*atan(1.0d0)
      g= 9.81d0
      ifail= 0
      om= 2*pi/t
c
c  Get surface elevations
c
      do 40 i=1,nfun-2
   40 eta(i)= x(i)
      eta(nfun)= eta(1)-h
      sm= 0.0
      do 10 i=2,nfun-2
   10 sm= sm+eta(i)
      eta(nfun-1)= -sm-(eta(1)+eta(nfun))/2
      do 11 i=2,nfun
      if (eta(i).gt.eta(i-1)+5*h/nfun) then
	 ifail= 1
c        write(*,'(7f11.3)') (eta(m),m=1,nfun-1)
	 return
      endif
   11 continue
      k= x(nfun-1)
c
c  Get dex(i,j) = d(eta(i))/d(eta(j)), 1 <= j <= nfun-2
c
      do 18 i=1,nfun
      dex(i,nfun-1)= 0.0d0
      do 18 j=1,nfun-2
      if (i.lt.nfun-1) then
	 if (i.eq.j) then
	    dex(i,j)= 1.0d0
	 else
	    dex(i,j)= 0.0d0
	 endif
      else if (i.eq.nfun-1) then
	 dex(i,j)= -1.0d0
      else if (i.eq.nfun) then
	 if (j.eq.1) then
	    dex(i,j)= 1.0d0
	 else
	    dex(i,j)= 0.0d0
	 endif
      endif
   18 continue
c
c  Set some useful functions
c
      do 23 i=1,nfun
      th(i)= pi*(i-1)/float(nfun-1)
      s(i)= d+eta(i)
      ks= k*s(i)
      do 23 n=1,nfun
      ch= cosh(n*ks)
      sh= sinh(n*ks)
      cs= cos(n*th(i))
      sn= sin(n*th(i))
      chcs(i,n)= ch*cs
      shsn(i,n)= sh*sn
      chsn(i,n)= ch*sn
   23 shcs(i,n)= sh*cs
c
c  Get normalised stream function coefficients a(j)
c
      do 12 i=1,nfun
      fu(i,1)= k*s(i)
      df(i,1)= k
      do 12 j=2,nfun
      n= j-1
      fu(i,j)= shcs(i,n)
   12 df(i,j)= chcs(i,n)*n*k
      call trans2(nfun,fu,df,a,da)
c
c  Get dax(i,j) = d(a(i))/d(eta(j)), 1 <=j <= nfun-2
c
      do 41 i=1,nfun
      do 41 j=1,nfun-2
      dij= 0.0
      do 42 n=1,nfun
   42 dij= dij+da(i,n)*dex(n,j)
   41 dax(i,j)= dij
c
c   Get dax(i,nfun-1) = d(a(i))/dk
c
      do 19 i=1,nfun
      sm= 0.0
      do 13 m=1,nfun
   13 sm= sm+da(i,m)*s(m)/k
   19 dax(i,nfun-1)= sm
c
c  De-normalise to get correct reverse mean flow
c
      a1k2= a(1)*k**2
      r= (u*k-om)/a1k2
      do 24 i=1,nfun
      c(i)= a(i)*r
      do 24 j=1,nfun-1
      drx= -(u*k-om)*dax(1,j)*(k/a1k2)**2
      if (j.eq.nfun-1) then
         drx= drx+u/a1k2-2*a(1)*k*(u*k-om)/a1k2**2
      endif
   24 dcx(i,j)= r*dax(i,j)+a(i)*drx
c
c  Get surface velocities and functions
c
      do 14 i=1,nfun
      su= k*c(1)
      sv= 0.0d0
      do 15 n=1,nfun-1
      su= su+n*k*c(n+1)*chcs(i,n)
   15 sv= sv+n*k*c(n+1)*shsn(i,n)
      us(i)= su
   14 vs(i)= sv
      sm= 0.0d0
      do 16 i=1,nfun
      f(i)= (us(i)**2+vs(i)**2)/(2*g)+eta(i)
   16 sm= sm+f(i)
      sm= sm/float(nfun)
      sq= 0.0d0
      do 17 i=1,nfun
      f(i)= f(i)-sm
   17 sq= sq+f(i)**2
      sq= sqrt(sq/nfun)/h
c
c  Get d(u(i))/d(eta(j)) and d(u(i))/dk
c
      do 20 i=1,nfun
      do 21 j=1,nfun-2
      su= k*dcx(1,j)
      sv= 0.0d0
      do 22 n=1,nfun-1
      su= su+n*k*dcx(n+1,j)*chcs(i,n)+(n*k)**2*c(n+1)*shcs(i,n)*dex(i,j)
   22 sv= sv+n*k*dcx(n+1,j)*shsn(i,n)+(n*k)**2*c(n+1)*chsn(i,n)*dex(i,j)
      dux(i,j)= su
   21 dvx(i,j)= sv
      j= nfun-1
      su= c(1)+k*dcx(1,j)
      sv= 0.0d0
      do 25 n=1,nfun-1
      su= su+(c(n+1)+k*dcx(n+1,j))*n*chcs(i,n)+
     :	    n**2*k*s(i)*c(n+1)*shcs(i,n)
   25 sv= sv+(c(n+1)+k*dcx(n+1,j))*n*shsn(i,n)+
     :	    n**2*k*s(i)*c(n+1)*chsn(i,n)
      dux(i,j)= su
   20 dvx(i,j)= sv
c
c  Get derivatives of functions
c
      do 44 i=1,nfun
      do 44 j=1,nfun-1
   44 dfx(i,j)= (us(i)*dux(i,j)+vs(i)*dvx(i,j))/g+dex(i,j)
c
      do 45 j=1,nfun-1
      sm= 0.0
      do 46 i=1,nfun
   46 sm= sm+dfx(i,j)
      sm= sm/nfun
      do 45 i=1,nfun
   45 jac(i,j)= dfx(i,j)-sm
c
      return
      end
c
c
c
      subroutine trans2(n,f,df,a,da)
c
c  Gramm-Schmidt orthogonalisation with derivatives
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      dimension f(nmax,nmax),df(nmax,nmax),a(nmax),da(nmax,nmax)
      dimension c(nmax,nmax),b(nmax),cr(nmax,nmax),g(nmax,nmax),
     :   dly(nmax),dc(nmax,nmax),br(nmax,nmax),db(nmax)
      do 10 i=1,n
      do 10 j=1,n
   10 c(i,j)= 0.0d0
c
      st= 0.0d0
      do 54 i=1,n
   54 st= st+f(i,1)**2
      st= dsqrt(st)
c
      do 16 i=1,n
   16 g(i,1)= f(i,1)/st
      c(1,1)= st
      do 18 k=2,n
      km= k-1
      do 29 j=1,km
      st= 0.0d0
      do 20 i=1,n
   20 st= st+g(i,j)*f(i,k)
   29 c(j,k)= st
      c(k,k)= 1.0d0
      do 21 i=1,n
      st= 0.0d0
      do 22 j=1,km
   22 st= st+c(j,k)*g(i,j)
   21 g(i,k)= f(i,k)-st
      st= 0.0d0
      do 50 i=1,n
   50 st= st+g(i,k)**2
      st= dsqrt(st)
      do 52 i=1,n
   52 g(i,k)= g(i,k)/st
      c(k,k)= st
      do 200 j=1,km
      st= 0.0d0
      do 202 i=1,n
  202 st= st+g(i,j)*g(i,k)
  200 db(j)= st
      do 204 i=1,n
      st= 0.0d0
      do 206 j=1,km
  206 st= st+db(j)*g(i,j)
  204 g(i,k)=g(i,k)-st
   18 continue
      do 24 j=1,n
      st= 0.0d0
      do 26 i=1,n
   26 st= st+g(i,j)
   24 b(j)= st
      st= 0.0d0
      do 208 i=1,n
      sb= 0.0d0
      do 210 j=1,n
  210 sb= sb+b(j)*g(i,j)
  208 st= st+(sb-1.0d0)**2
      call trinv(n,c,cr)
      do 28 i=1,n
      st= 0.0d0
      do 30 j=1,n
   30 st= st+cr(i,j)*b(j)
   28 a(i)= st
      st= 0.0d0
      do 72 i=1,n
      sb= 0.0d0
      do 74 j=1,n
   74 sb= sb+a(j)*f(i,j)
   72 st= st+(sb-1.0d0)**2
      st= dsqrt(st/float(n))
      do 80 k=1,n
      do 81 i=1,n
      do 81 j=1,n
   81 dc(i,j)= 0.0d0
      dly(1)= 2.0d0*f(k,1)*df(k,1)
      dc(1,1)= 0.5d0*dly(1)/c(1,1)
      do 82 j=2,n
   82 dc(1,j)= (f(k,j)*df(k,1)+df(k,j)*f(k,1))/c(1,1)
     :	    -0.5d0*dly(1)*c(1,j)/c(1,1)**2
      do 84 l=2,n
      lm= l-1
      st= 0.0d0
      do 86 m=1,lm
   86 st= st+c(m,l)*dc(m,l)
      dly(l)= 2.0d0*(f(k,l)*df(k,l)-st)
      dc(l,l)= 0.5d0*dly(l)/c(l,l)
      if (l.eq.n) goto 84
      lp= l+1
      do 85 m=lp,n
      st= 0.0d0
      do 90 j=1,lm
   90 st= st+c(j,m)*dc(j,l)+dc(j,m)*c(j,l)
   85 dc(l,m)= (f(k,m)*df(k,l)+df(k,m)*f(k,l))/c(l,l)
     :	    -st/c(l,l)-0.5d0*dly(l)*c(l,m)/c(l,l)**2
   84 continue
      do 92 l=1,n
      do 92 m=1,n
      st= 0.0d0
      do 96 j=1,n
   96 st= st+g(l,j)*dc(j,m)
   92 br(l,m)= -st
      do 94 m=1,n
   94 br(k,m)= br(k,m)+df(k,m)
      do 32 j=1,n
      sb= 0.0d0
      do 34 i=1,n
      st= 0.0d0
      do 36 l=1,n
   36 st= st+br(i,l)*cr(l,j)
   34 sb= sb+st
   32 db(j)= sb
      do 38 i=1,n
      st= db(i)
      do 40 j=1,n
   40 st= st-dc(i,j)*a(j)
   38 dly(i)= st
      do 42 i=1,n
      st= 0.0d0
      do 44 j=1,n
   44 st= st+cr(i,j)*dly(j)
   42 da(i,k)= st
   80 continue
      return
      end
c
c
c
      subroutine trinv(n,c,w)
c
c  3-diagonal system solver
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      dimension c(nmax,nmax),w(nmax,nmax)
      do 10 i=1,n
      do 10 j=1,n
   10 w(i,j)= 0.0d0
      do 12 jj=1,n
      j= n+1-jj
      w(j,j)= 1.0d0
      if (j.eq.1) goto 12
      k= j-1
      do 13 ii=1,k
      i= k-ii+1
      sum= 0.0d0
      do 14 l=1,ii
      m= n+2-jj-l
   14 sum= sum+c(i,m)*w(m,j)/c(m,m)
   13 w(i,j)= -sum
   12 continue
      do 16 i=1,n
      sum= c(i,i)
      do 16 j=1,n
   16 w(i,j)= w(i,j)/sum
      return
      end
c
c
c
      subroutine wavecel(ta,d,u,tr,c)
c
c  Linear theory C by series approximation for waves on a current u
c    d  = still water depth
c    ta = absolute period (fixed reference frame)
c    tr = relative period (reference frame moving with the current)
c    c  = celerity (reference frame moving with the current)
c
      implicit double precision (a-h,o-z)
      pi= 4*atan(1.0d0)
      g= 9.81
      sigma= 2*pi/ta
      y= sigma*sigma*d/g
      a= 1.0/(1.0+y*(0.6667+y*(0.3556+y*(0.1608+y*(0.06321+
     :	y*(0.02744+y*0.01))))))
      c= dsqrt((d*g)/(y+a))
      if (abs(u).lt.1.0d-6) then
	 tr= ta
	 return
      else
	 el= c*ta
	 iter= 0
   10	 tr= el/(el/ta-u)
         elp= (g*tr**2/(2*pi))*tanh(2*pi*d/el)
	 del= elp-el
	 el= el+del*0.4
	 if (abs(del/el).gt.1.0d-6) then
	    iter= iter+1
	    if (iter.eq.100) then
               write(*,'(a)') ' WAVECEL error'
	       stop
	    endif
	    goto 10
	 endif
	 tr= el/(el/ta-u)
	 c= el/tr
	 return
      endif
      end
c
c
c
      subroutine slpds(a,b,n,cc)
c
c  Linear system solution
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      dimension a(nmax,nmax),b(nmax),cc(nmax)
      n1= n-1
      do 100 k=1,n1
      c= a(k,k)
      k1= k+1
      if (abs(c).lt.1.0d-10) then
         write(*,'(a,i5)') ' Matrix error 1: ',k
         stop
      endif
      do 4 j=k1,n
    4 a(k,j)= a(k,j)/c
      b(k)= b(k)/c
      do 10 i=k1,n
      c= a(i,k)
      do 5 j=k1,n
    5 a(i,j)= a(i,j)-c*a(k,j)
   10 b(i)= b(i)-c*b(k)
  100 continue
      if (abs(a(n,n)).lt.1.0d-10) then
         write(*,'(a,i5)') ' Matrix error 2: ',n
         stop
      endif
      b(n)= b(n)/a(n,n)
      do 200 l=1,n1
      k= n-l
      k1= k+1
      do 200 j=k1,n
  200 b(k)= b(k)-a(k,j)*b(j)
      do 300 i=1,n
  300 cc(i)= b(i)
      return
      end
c
c
c
      subroutine gaf(nf,nv,xc,fvecc,fsumsq,nverb,iter,ifail)
c
c   Non-linear system error minimisation
c   ifail =  0: OK
c            1: looks hopeless
c            2: poor convergence
c
      implicit double precision (a-h,o-z)
      parameter (nmax=25)
      dimension xc(nmax),fvecc(nmax)
      dimension fjacc(nmax,nmax),
     :   aa(nmax,nmax),bb(nmax),cc(nmax),xcm(nmax)
c
      iter = 0
      icalls= 0
      orf= 2.0/nf
      fsumsm= 100.0
c
   15 do 16 i=1,nv
   16 xcm(i)= xc(i)
   19 call lsfun(nf,xc,fvecc,fjacc,fsumsq,ifl)
      icalls= icalls+1
      if ((fsumsq.gt.fsumsm.and.iter.gt.1).or.ifl.ne.0) then
         if (orf.lt.0.05) then
            ifail= 1
            return
         endif
         orf= 0.8*orf
         iter= max(iter-1,0)
         do 18 i=1,nv
   18    xc(i)= xcm(i)
         fsumsm= 100.0
         goto 19
      endif
      call monit(nf,fvecc,icalls,nverb)
c
      do 10 i=1,nv
      do 10 j=1,nv
      ra= 0.0d0
      do 11 l=1,nf
   11 ra= ra+fjacc(l,i)*fjacc(l,j)
   10 aa(i,j)= ra
c
      do 12 i=1,nv
      ra= 0.0d0
      do 13 l=1,nf
   13 ra= ra+fvecc(l)*fjacc(l,i)
   12 bb(i)= -ra
c
      call slpds(aa,bb,nv,cc)
c
      dxmax= 0.0d0
      do 14 i=1,nv
      dxmax= max(dxmax,abs(cc(i)))
   14 xc(i)= xc(i)+orf*cc(i)
c
      iter= iter+1
      fsumsm= fsumsq
      orf= min(1.0d0,orf*1.1)
c
      if (iter.ge.50.and.fsumsq.lt.1.0d-4) then
         ifail= 2
         return
      endif
      if (iter.ge.50) then
         ifail= 1
         return
      endif
      if (fsumsq.gt.1.0d-6) goto 15
      ifail= 0
      return
      end
c
c
c
      subroutine monit(nfun,f,icalls,nverb)
      implicit double precision (a-h,o-z)
      double precision k
      parameter (nmax=25)
      dimension f(nmax)
      common /one/ d,t,h,u,k
     :       /two/ eta(nmax),c(nmax),amp(0:nmax)
c
c  Outputs monitoring data
c
      if (nverb.eq.0) return
      sm= 0.0d0
      do 10 i=1,nfun
   10 sm=sm+f(i)**2
      sq= sqrt(sm/nfun)/h
      write(*,'(/i11,1p1e11.3)') icalls,sq
      write(*,'(1p7e11.3)') (eta(i),i=1,nfun)
      return
      end
c
c
c
      subroutine limit(h,d,t,rat,nverb)
      double precision h,d,t,rat
c
c  Estimates H/(limiting height for d and t) = rat
c
      dimension dl0(18),hl0(18)
      data dl0/2,0.578,0.440,0.356,0.293,0.243,0.201,0.166,0.1359,
     :         0.1100,0.0876,0.0686,0.0524,0.0390,0.0277,0.01879,
     :         0.01168,0.00638/
      data hl0/0.1682,0.1665,0.1613,0.1531,0.1423,0.1298,0.1159,
     :         0.1017,0.0873,0.0735,0.0605,0.0487,0.0380,0.0289,
     :         0.0208,0.01440,0.00911,0.00501/
      pi= 4*atan(1.0)
      el0= 9.81*t**2/(2*pi)
      ha= h/el0
      da= d/el0
      if (da.gt.dl0(1)) then
         rat= ha/hl0(1)
      else if (da.lt.dl0(18)) then
         rat= ha/(0.8*da)
      else
         do 10 i=2,18
         if (dl0(i).lt.da) goto 11
   10    continue
   11    x1= log(dl0(i))
         x2= log(dl0(i-1))
         y1= log(hl0(i))
         y2= log(hl0(i-1))
         r= (log(da)-x1)/(x2-x1)
         hb= exp(y1+r*(y2-y1))
         rat= ha/hb
      endif
c
c      if (nverb.ne.0.or.rat.gt.1.0) then
c         write(*,'(a,f5.3)') ' H/Hb              = ',rat
c      endif
      if (rat.gt.1.0) stop
c
      return
      end
c
c
c
      subroutine four(f,n,a,b,nb)
c
c  Fourier analysis
c
      implicit double precision (a-h,o-z)
      dimension a(0:nb),b(0:nb),f(n)
      pi= 4*atan(1.0d0)
      rn= 2.0d0/n
      t= 2*pi/n
      c= cos(t)
      s= sin(t)
      vk= 0.0d0
      vl= -1.0d0
      do 10 k=0,nb 
      t= c*vk  
      ck= t-vl 
      vl= vk   
      vk= ck+t 
      t= ck+ck 
      ul= 0.0d0
      um= f(n) 
      do 12 mm=3,n 
      m= n+2-mm
      u0= ul   
      ul= um   
   12 um= t*ul-u0+f(m) 
      a(k)= (ck*um-ul+f(1))*rn 
   10 b(k)= s*vl*um*rn 
      a(0)= a(0)*0.5d0
      if (2*nb.ne.n) return
      a(nb)= a(nb)*0.5d0
      b(nb)= 0.0d0
      return   
      end
c
c
c
      subroutine kmts(nfun,xx,yy,tt,uu,vv,ut,vt,du,dv,etah)
c
c   Computes
c        horizontal and vertical velocity components (u,v)
c        horizontal and vertical local acceleration components (ut,vt)
c        horizontal and vertical total acceleration components (du,dv)
c        water surface elevation (etah)
c   at t=tt; x=xx; y=yy.
c   If yy>eta kinematics are returned at the free surface
c
      implicit double precision (a-h,o-z)
      real xx,yy,tt,uu,vv,ut,vt,du,dv,etah
      parameter (nmax=25)
      double precision k,ks
      common /one/ d,t,h,u,k
     :       /two/ eta(nmax),c(nmax),amp(0:nmax)
      pi=4*atan(1.0d0)
      om= 2*pi/t
      theta= k*xx-om*tt
c
      etah= 0.0
      do 11 i=1,nfun-1
   11 etah= etah+cos(i*theta)*amp(i)
c
      ks= k*(d+min(yy,etah))
      s1= 0.0
      s2= 0.0
      s3= 0.0
      s4= 0.0
c
      do 10 i=1,nfun-1
      ip= i+1
      ch= cosh(i*ks)
      sh= sinh(i*ks)
      cs= cos(i*theta)
      sn= sin(i*theta)
      s1= s1+i*ch*cs*c(ip)
      s2= s2+i*sh*sn*c(ip)
      s3= s3+i*i*ch*sn*c(ip)
   10 s4= s4+i*i*sh*cs*c(ip)
      uu=  u+k*s1
      vv=  k*s2
      ut=  k*om*s3
      vt= -k*om*s4
      ux= -k*k*s3
      vx=  k*k*s4
      uy=  vx
      vy= -ux
      du= ut+uu*ux+vv*uy
      dv= vt+uu*vx+vv*vy
c
      return
      end

