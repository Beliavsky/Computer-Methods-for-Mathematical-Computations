      subroutine decomp(ndim,n,a,cond,ipvt,work)
c
      integer ndim,n
      double precision a(ndim,n),cond,work(n)
      integer ipvt(n)
c
c     decomposes a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     use solve to compute solutions to linear systems.
c
c     input..
c
c        ndim = declared row dimension of the array containing  a.
c
c        n = order of the matrix.
c
c        a = matrix to be triangularized.
c
c     output..
c
c        a  contains an upper triangular matrix  u  and a permuted
c          version of a lower triangular matrix  i-l  so that
c          (permutation matrix)*a = l*u .
c
c        cond = an estimate of the condition of  a .
c           for the linear system  a*x = b, changes in  a  and  b
c           may cause changes  cond  times as large in  x .
c           if  cond+1.0 .eq. cond , a is singular to working
c           precision.  cond  is set to  1.0d+32  if exact
c           singularity is detected.
c
c        ipvt = the pivot vector.
c           ipvt(k) = the index of the k-th pivot row
c           ipvt(n) = (-1)**(number of interchanges)
c
c     work space..  the vector  work  must be declared and included
c                   in the call.  its input contents are ignored.
c                   its output contents are usually unimportant.
c
c     the determinant of a can be obtained on output by
c        det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n).
c
      double precision ek, t, anorm, ynorm, znorm
      integer nm1, i, j, k, kp1, kb, km1, m
      double precision dabs, dsign
c
      ipvt(n) = 1
      if (n .eq. 1) go to 80
      nm1 = n - 1
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         t = 0.0d0
         do 5 i = 1, n
            t = t + dabs(a(i,j))
    5    continue
         if (t .gt. anorm) anorm = t
   10 continue
c
c     gaussian elimination with partial pivoting
c
      do 35 k = 1,nm1
         kp1= k+1
c
c        find pivot
c
         m = k
         do 15 i = kp1,n
            if (dabs(a(i,k)) .gt. dabs(a(m,k))) m = i
   15    continue
         ipvt(k) = m
         if (m .ne. k) ipvt(n) = -ipvt(n)
         t = a(m,k)
         a(m,k) = a(k,k)
         a(k,k) = t
c
c        skip step if pivot is zero
c
         if (t .eq. 0.0d0) go to 35
c
c        compute multipliers
c
         do 20 i = kp1,n
             a(i,k) = -a(i,k)/t
   20    continue
c
c        interchange and eliminate by columns
c
         do 30 j = kp1,n
             t = a(m,j)
             a(m,j) = a(k,j)
             a(k,j) = t
             if (t .eq. 0.0d0) go to 30
             do 25 i = kp1,n
                a(i,j) = a(i,j) + a(i,k)*t
   25        continue
   30    continue
   35 continue
c
c     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
c     estimate obtained by one step of inverse iteration for the
c     small singular vector.  this involves solving two systems
c     of equations, (a-transpose)*y = e  and  a*z = y  where  e
c     is a vector of +1 or -1 chosen to cause growth in y.
c     estimate = (1-norm of z)/(1-norm of y)
c
c     solve (a-transpose)*y = e
c
      do 50 k = 1, n
         t = 0.0d0
         if (k .eq. 1) go to 45
         km1 = k-1
         do 40 i = 1, km1
            t = t + a(i,k)*work(i)
   40    continue
   45    ek = 1.0d0
         if (t .lt. 0.0d0) ek = -1.0d0
         if (a(k,k) .eq. 0.0d0) go to 90
         work(k) = -(ek + t)/a(k,k)
   50 continue
      do 60 kb = 1, nm1
         k = n - kb
         t = 0.0d0
         kp1 = k+1
         do 55 i = kp1, n
            t = t + a(i,k)*work(i)
   55    continue
         work(k) = t + work(k)
         m = ipvt(k)
         if (m .eq. k) go to 60
         t = work(m)
         work(m) = work(k)
         work(k) = t
   60 continue
c
      ynorm = 0.0d0
      do 65 i = 1, n
         ynorm = ynorm + dabs(work(i))
   65 continue
c
c     solve a*z = y
c
      call solve(ndim, n, a, work, ipvt)
c
      znorm = 0.0d0
      do 70 i = 1, n
         znorm = znorm + dabs(work(i))
   70 continue
c
c     estimate condition
c
      cond = anorm*znorm/ynorm
      if (cond .lt. 1.0d0) cond = 1.0d0
      return
c
c     1-by-1
c
   80 cond = 1.0d0
      if (a(1,1) .ne. 0.0d0) return
c
c     exact singularity
c
   90 cond = 1.0d+32
      return
      end


      subroutine solve(ndim, n, a, b, ipvt)
c
      integer ndim, n, ipvt(n)
      double precision a(ndim,n),b(n)
c
c   solution of linear system, a*x = b .
c   do not use if decomp has detected singularity.
c
c   input..
c
c     ndim = declared row dimension of array containing a .
c
c     n = order of matrix.
c
c     a = triangularized matrix obtained from decomp .
c
c     b = right hand side vector.
c
c     ipvt = pivot vector obtained from decomp .
c
c   output..
c
c     b = solution vector, x .
c
      integer kb, km1, nm1, kp1, i, k, m
      double precision t
c
c     forward elimination
c
      if (n .eq. 1) go to 50
      nm1 = n-1
      do 20 k = 1, nm1
         kp1 = k+1
         m = ipvt(k)
         t = b(m)
         b(m) = b(k)
         b(k) = t
         do 10 i = kp1, n
             b(i) = b(i) + a(i,k)*t
   10    continue
   20 continue
c
c     back substitution
c
      do 40 kb = 1,nm1
         km1 = n-kb
         k = km1+1
         b(k) = b(k)/a(k,k)
         t = -b(k)
         do 30 i = 1, km1
             b(i) = b(i) + a(i,k)*t
   30    continue
   40 continue
   50 b(1) = b(1)/a(1,1)
      return
      end

      subroutine spline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end

      double precision function seval(n, u, x, y, b, c, d)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end

      subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
c
      double precision fun, a, b, abserr, relerr, result, errest, flag
      integer nofun
c
c   estimate the integral of fun(x) from a to b
c   to a user provided tolerance.
c   an automatic adaptive routine based on
c   the 8-panel newton-cotes rule.
c
c   input ..
c
c   fun     the name of the integrand function subprogram fun(x).
c   a       the lower limit of integration.
c   b       the upper limit of integration.(b may be less than a.)
c   relerr  a relative error tolerance. (should be non-negative)
c   abserr  an absolute error tolerance. (should be non-negative)
c
c   output ..
c
c   result  an approximation to the integral hopefully satisfying the
c           least stringent of the two error tolerances.
c   errest  an estimate of the magnitude of the actual error.
c   nofun   the number of function values used in calculation of result.
c   flag    a reliability indicator.  if flag is zero, then result
c           probably satisfies the error tolerance.  if flag is
c           xxx.yyy , then  xxx = the number of intervals which have
c           not converged and  0.yyy = the fraction of the interval
c           left to do when the limit on  nofun  was approached.
c
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),f(16),x(16),fsave(8,30),xsave(8,30)
      double precision dabs,dmax1
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
c
c   ***   stage 1 ***   general initialization
c   set constants.
c
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
c
c   trouble when nofun reaches nofin
c
      w0 =   3956.0d0 / 14175.0d0
      w1 =  23552.0d0 / 14175.0d0
      w2 =  -3712.0d0 / 14175.0d0
      w3 =  41984.0d0 / 14175.0d0
      w4 = -18160.0d0 / 14175.0d0
c
c   initialize running sums to zero.
c
      flag = 0.0d0
      result = 0.0d0
      cor11  = 0.0d0
      errest = 0.0d0
      area   = 0.0d0
      nofun = 0
      if (a .eq. b) return
c
c   ***   stage 2 ***   initialization for first interval
c
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev  = 0.0d0
      f0 = fun(x0)
      stone = (b - a) / 16.0d0
      x(8)  =  (x0  + x(16)) / 2.0d0
      x(4)  =  (x0  + x(8))  / 2.0d0
      x(12) =  (x(8)  + x(16)) / 2.0d0
      x(2)  =  (x0  + x(4))  / 2.0d0
      x(6)  =  (x(4)  + x(8))  / 2.0d0
      x(10) =  (x(8)  + x(12)) / 2.0d0
      x(14) =  (x(12) + x(16)) / 2.0d0
      do 25 j = 2, 16, 2
         f(j) = fun(x(j))
   25 continue
      nofun = 9
c
c   ***   stage 3 ***   central calculation
c   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
c   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
c
   30 x(1) = (x0 + x(2)) / 2.0d0
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.0d0
         f(j) = fun(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0d0
      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6))
     1  + w3*(f(3)+f(5))  +  w4*f(4)) * step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))
     1  + w3*(f(11)+f(13)) + w4*f(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
c
c   ***   stage 4 *** interval convergence test
c
      esterr = dabs(qdiff) / 1023.0d0
      tolerr = dmax1(abserr,relerr*dabs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
c
c   ***   stage 5   ***   no convergence
c   locate next interval.
c
   50 nim = 2*nim
      lev = lev+1
c
c   store right hand elements for future use.
c
      do 52 i = 1, 8
         fsave(i,lev) = f(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
c
c   assemble left hand elements for immediate use.
c
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         f(2*j+18) = f(j+9)
         x(2*j+18) = x(j+9)
   55 continue
      go to 30
c
c   ***   stage 6   ***   trouble section
c   number of function values is about to exceed limit.
c
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
c
c   current level is levmax.
c
   62 flag = flag + 1.0d0
c
c   ***   stage 7   ***   interval converged
c   add contributions into running sums.
c
   70 result = result + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0d0
c
c   locate next interval.
c
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
c
c   assemble elements required for the next interval.
c
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8
         f(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      go to 30
c
c   ***   stage 8   ***   finalize and return
c
   80 result = result + cor11
c
c   make sure errest not less than roundoff level.
c
      if (errest .eq. 0.0d0) return
   82 temp = dabs(result) + errest
      if (temp .ne. dabs(result)) return
      errest = 2.0d0*errest
      go to 82
      end

      subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c     written by h.a.watts and l.f.shampine
c                   sandia laboratories
c                  albuquerque,new mexico
c
c    rkf45 is primarily designed to solve non-stiff and mildly stiff
c    differential equations when derivative evaluations are inexpensive.
c    rkf45 should generally not be used when the user is demanding
c    high accuracy.
c
c abstract
c
c    subroutine  rkf45  integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
c              where the y(i) are given at t .
c    typically the subroutine is used to integrate from t to tout but it
c    can be used as a one-step integrator to advance the solution a
c    single step in the direction of tout.  on return the parameters in
c    the call list are set for continuing the integration. the user has
c    only to call rkf45 again (and perhaps define a new value for tout).
c    actually, rkf45 is an interfacing routine which calls subroutine
c    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
c    computes an approximate solution over one step.
c
c    rkf45  uses the runge-kutta-fehlberg (4,5)  method described
c    in the reference
c    e.fehlberg , low-order classical runge-kutta formulas with stepsize
c                 control , nasa tr r-315
c
c    the performance of rkf45 is illustrated in the reference
c    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
c                 differential equations-the state of the art ,
c                 sandia laboratories report sand75-0182 ,
c                 to appear in siam review.
c
c
c    the parameters represent-
c      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated
c      y(*) -- solution vector at t
c      t -- independent variable
c      tout -- output point at which solution is desired
c      relerr,abserr -- relative and absolute error tolerances for local
c            error test. at each step the code requires that
c                 abs(local error) .le. relerr*abs(y) + abserr
c            for each component of the local error and solution vectors
c      iflag -- indicator for status of integration
c      work(*) -- array to hold information internal to rkf45 which is
c            necessary for subsequent calls. must be dimensioned
c            at least  3+6*neqn
c      iwork(*) -- integer array used to hold information internal to
c            rkf45 which is necessary for subsequent calls. must be
c            dimensioned at least  5
c
c
c  first call to rkf45
c
c    the user must provide storage in his calling program for the arrays
c    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
c    declare f in an external statement, supply subroutine f(t,y,yp) and
c    initialize the following parameters-
c
c      neqn -- number of equations to be integrated.  (neqn .ge. 1)
c      y(*) -- vector of initial conditions
c      t -- starting point of integration , must be a variable
c      tout -- output point at which solution is desired.
c            t=tout is allowed on the first call only, in which case
c            rkf45 returns with iflag=2 if continuation is possible.
c      relerr,abserr -- relative and absolute local error tolerances
c            which must be non-negative. relerr must be a variable while
c            abserr may be a constant. the code should normally not be
c            used with relative error control smaller than about 1.e-8 .
c            to avoid limiting precision difficulties the code requires
c            relerr to be larger than an internally computed relative
c            error parameter which is machine dependent. in particular,
c            pure absolute error is not permitted. if a smaller than
c            allowable value of relerr is attempted, rkf45 increases
c            relerr appropriately and returns control to the user before
c            continuing the integration.
c      iflag -- +1,-1  indicator to initialize the code for each new
c            problem. normal input is +1. the user should set iflag=-1
c            only when one-step integrator control is essential. in this
c            case, rkf45 attempts to advance the solution a single step
c            in the direction of tout each time it is called. since this
c            mode of operation results in extra computing overhead, it
c            should be avoided unless needed.
c
c
c  output from rkf45
c
c      y(*) -- solution at t
c      t -- last point reached in integration.
c      iflag = 2 -- integration reached tout. indicates successful retur
c                   and is the normal mode for continuing integration.
c            =-2 -- a single successful step in the direction of tout
c                   has been taken. normal mode for continuing
c                   integration one step at a time.
c            = 3 -- integration was not completed because relative error
c                   tolerance was too small. relerr has been increased
c                   appropriately for continuing.
c            = 4 -- integration was not completed because more than
c                   3000 derivative evaluations were needed. this
c                   is approximately 500 steps.
c            = 5 -- integration was not completed because solution
c                   vanished making a pure relative error test
c                   impossible. must use non-zero abserr to continue.
c                   using the one-step integration mode for one step
c                   is a good way to proceed.
c            = 6 -- integration was not completed because requested
c                   accuracy could not be achieved using smallest
c                   allowable stepsize. user must increase the error
c                   tolerance before continued integration can be
c                   attempted.
c            = 7 -- it is likely that rkf45 is inefficient for solving
c                   this problem. too much output is restricting the
c                   natural stepsize choice. use the one-step integrator
c                   mode.
c            = 8 -- invalid input parameters
c                   this indicator occurs if any of the following is
c                   satisfied -   neqn .le. 0
c                                 t=tout  and  iflag .ne. +1 or -1
c                                 relerr or abserr .lt. 0.
c                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8
c      work(*),iwork(*) -- information which is usually of no interest
c                   to the user but necessary for subsequent calls.
c                   work(1),...,work(neqn) contain the first derivatives
c                   of the solution vector y at t. work(neqn+1) contains
c                   the stepsize h to be attempted on the next step.
c                   iwork(1) contains the derivative evaluation counter.
c
c
c  subsequent calls to rkf45
c
c    subroutine rkf45 returns with all information needed to continue
c    the integration. if the integration reached tout, the user need onl
c    define a new tout and call rkf45 again. in the one-step integrator
c    mode (iflag=-2) the user must keep in mind that each step taken is
c    in the direction of the current tout. upon reaching tout (indicated
c    by changing iflag to 2),the user must then define a new tout and
c    reset iflag to -2 to continue in the one-step integrator mode.
c
c    if the integration was not completed but the user still wants to
c    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
c    the relerr parameter has been adjusted appropriately for continuing
c    the integration. in the case of iflag=4 the function counter will
c    be reset to 0 and another 3000 function evaluations are allowed.
c
c    however,in the case iflag=5, the user must first alter the error
c    criterion to use a positive value of abserr before integration can
c    proceed. if he does not,execution is terminated.
c
c    also,in the case iflag=6, it is necessary for the user to reset
c    iflag to 2 (or -2 when the one-step integration mode is being used)
c    as well as increasing either abserr,relerr or both before the
c    integration can be continued. if this is not done, execution will
c    be terminated. the occurrence of iflag=6 indicates a trouble spot
c    (solution is changing rapidly,singularity may be present) and it
c    often is inadvisable to continue.
c
c    if iflag=7 is encountered, the user should use the one-step
c    integration mode with the stepsize determined by the code or
c    consider switching to the adams codes de/step,intrp. if the user
c    insists upon continuing the integration with rkf45, he must reset
c    iflag to 2 before calling rkf45 again. otherwise,execution will be
c    terminated.
c
c    if iflag=8 is obtained, integration can not be continued unless
c    the invalid input parameters are corrected.
c
c    it should be noted that the arrays work,iwork contain information
c    required for subsequent integration. accordingly, work and iwork
c    should not be altered.
c
c
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr,abserr,work(1)
c     if compiler checks subscripts, change work(1) to work(3+6*neqn)
c
      external f
c
      integer k1,k2,k3,k4,k5,k6,k1m
c
c
c     compute indices for the splitting of the work array
c
      k1m=neqn+1
      k1=k1m+1
      k2=k1+neqn
      k3=k2+neqn
      k4=k3+neqn
      k5=k4+neqn
      k6=k5+neqn
c
c     this interfacing routine merely relieves the user of a long
c     calling list via the splitting apart of two working storage
c     arrays. if this is not compatible with the users compiler,
c     he must use rkfs directly.
c
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),work(k1m),
     1          work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),
     2          work(k6+1),iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
c
      return
      end
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,
     1                f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c
c     rkfs integrates a system of first order ordinary differential
c     equations as described in the comments for rkf45 .
c     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
c     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
c     internally by the code and appear in the call list to eliminate
c     local retention of variables between calls. accordingly, they
c     should not be altered. items of possible interest are
c         yp - derivative of solution vector at t
c         h  - an appropriate stepsize to be used for the next step
c         nfe- counter on the number of derivative function evaluations
c
c
      logical hfaild,output
c
      integer  neqn,iflag,nfe,kop,init,jflag,kflag
      double precision  y(neqn),t,tout,relerr,abserr,h,yp(neqn),
     1  f1(neqn),f2(neqn),f3(neqn),f4(neqn),f5(neqn),savre,
     2  savae
c
      external f
c
      double precision  a,ae,dt,ee,eeoet,esttol,et,hmin,remin,rer,s,
     1  scale,tol,toln,u26,epsp1,eps,ypk
c
      integer  k,maxnfe,mflag
c
      double precision  dabs,dmax1,dmin1,dsign
c
c     remin is the minimum acceptable value of relerr.  attempts
c     to obtain higher accuracy with this subroutine are usually
c     very expensive and often unsuccessful.
c
      data remin/1.d-12/
c
c
c     the expense is controlled by restricting the number
c     of function evaluations to be approximately maxnfe.
c     as set, this corresponds to about 500 steps.
c
      data maxnfe/3000/
c
c
c     check input parameters
c
c
      if (neqn .lt. 1) go to 10
      if ((relerr .lt. 0.0d0)  .or.  (abserr .lt. 0.0d0)) go to 10
      mflag=iabs(iflag)
      if ((mflag .eq. 0) .or. (mflag .gt. 8)) go to 10
      if (mflag .ne. 1) go to 20
c
c     first call, compute machine epsilon
c
      eps = 1.0d0
    5 eps = eps/2.0d0
      epsp1 = eps + 1.0d0
      if (epsp1 .gt. 1.0d0) go to 5
      u26 = 26.0d0*eps
      go to 50
c
c     invalid input
   10 iflag=8
      return
c
c     check continuation possibilities
c
   20 if ((t .eq. tout) .and. (kflag .ne. 3)) go to 10
      if (mflag .ne. 2) go to 25
c
c     iflag = +2 or -2
      if ((kflag .eq. 3) .or. (init .eq. 0)) go to 45
      if (kflag .eq. 4) go to 40
      if ((kflag .eq. 5)  .and.  (abserr .eq. 0.0d0)) go to 30
      if ((kflag .eq. 6)  .and.  (relerr .le. savre)  .and.
     1    (abserr .le. savae)) go to 30
      go to 50
c
c     iflag = 3,4,5,6,7 or 8
   25 if (iflag .eq. 3) go to 45
      if (iflag .eq. 4) go to 40
      if ((iflag .eq. 5) .and. (abserr .gt. 0.0d0)) go to 45
c
c     integration cannot be continued since user did not respond to
c     the instructions pertaining to iflag=5,6,7 or 8
   30 stop
c
c     reset function evaluation counter
   40 nfe=0
      if (mflag .eq. 2) go to 50
c
c     reset flag value from previous call
   45 iflag=jflag
      if (kflag .eq. 3) mflag=iabs(iflag)
c
c     save input iflag and set continuation flag value for subsequent
c     input checking
   50 jflag=iflag
      kflag=0
c
c     save relerr and abserr for checking input on subsequent calls
      savre=relerr
      savae=abserr
c
c     restrict relative error tolerance to be at least as large as
c     2*eps+remin to avoid limiting precision difficulties arising
c     from impossible accuracy requests
c
      rer=2.0d0*eps+remin
      if (relerr .ge. rer) go to 55
c
c     relative error tolerance too small
      relerr=rer
      iflag=3
      kflag=3
      return
c
   55 dt=tout-t
c
      if (mflag .eq. 1) go to 60
      if (init .eq. 0) go to 65
      go to 80
c
c     initialization --
c                       set initialization completion indicator,init
c                       set indicator for too many output points,kop
c                       evaluate initial derivatives
c                       set counter for function evaluations,nfe
c                       estimate starting stepsize
c
   60 init=0
      kop=0
c
      a=t
      call f(a,y,yp)
      nfe=1
      if (t .ne. tout) go to 65
      iflag=2
      return
c
c
   65 init=1
      h=dabs(dt)
      toln=0.
      do 70 k=1,neqn
        tol=relerr*dabs(y(k))+abserr
        if (tol .le. 0.) go to 70
        toln=tol
        ypk=dabs(yp(k))
        if (ypk*h**5 .gt. tol) h=(tol/ypk)**0.2d0
   70 continue
      if (toln .le. 0.0d0) h=0.0d0
      h=dmax1(h,u26*dmax1(dabs(t),dabs(dt)))
      jflag=isign(2,iflag)
c
c
c     set stepsize for integration in the direction from t to tout
c
   80 h=dsign(h,dt)
c
c     test to see if rkf45 is being severely impacted by too many
c     output points
c
      if (dabs(h) .ge. 2.0d0*dabs(dt)) kop=kop+1
      if (kop .ne. 100) go to 85
c
c     unnecessary frequency of output
      kop=0
      iflag=7
      return
c
   85 if (dabs(dt) .gt. u26*dabs(t)) go to 95
c
c     if too close to output point,extrapolate and return
c
      do 90 k=1,neqn
   90   y(k)=y(k)+dt*yp(k)
      a=tout
      call f(a,y,yp)
      nfe=nfe+1
      go to 300
c
c
c     initialize output point indicator
c
   95 output= .false.
c
c     to avoid premature underflow in the error tolerance function,
c     scale the error tolerances
c
      scale=2.0d0/relerr
      ae=scale*abserr
c
c
c     step by step integration
c
  100 hfaild= .false.
c
c     set smallest allowable stepsize
c
      hmin=u26*dabs(t)
c
c     adjust stepsize if necessary to hit the output point.
c     look ahead two steps to avoid drastic changes in the stepsize and
c     thus lessen the impact of output points on the code.
c
      dt=tout-t
      if (dabs(dt) .ge. 2.0d0*dabs(h)) go to 200
      if (dabs(dt) .gt. dabs(h)) go to 150
c
c     the next successful step will complete the integration to the
c     output point
c
      output= .true.
      h=dt
      go to 200
c
  150 h=0.5d0*dt
c
c
c
c     core integrator for taking a single step
c
c     the tolerances have been scaled to avoid premature underflow in
c     computing the error tolerance function et.
c     to avoid problems with zero crossings,relative error is measured
c     using the average of the magnitudes of the solution at the
c     beginning and end of a step.
c     the error estimate formula has been grouped to control loss of
c     significance.
c     to distinguish the various arguments, h is not permitted
c     to become smaller than 26 units of roundoff in t.
c     practical limits on the change in the stepsize are enforced to
c     smooth the stepsize selection process and to avoid excessive
c     chattering on problems having discontinuities.
c     to prevent unnecessary failures, the code uses 9/10 the stepsize
c     it estimates will succeed.
c     after a step failure, the stepsize is not allowed to increase for
c     the next attempted step. this makes the code more efficient on
c     problems having discontinuities and more effective in general
c     since local extrapolation is being used and extra caution seems
c     warranted.
c
c
c     test number of derivative function evaluations.
c     if okay,try to advance the integration from t to t+h
c
  200 if (nfe .le. maxnfe) go to 220
c
c     too much work
      iflag=4
      kflag=4
      return
c
c     advance an approximate solution over one step of length h
c
  220 call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
c
c     compute and test allowable tolerances versus local error estimates
c     and remove scaling of tolerances. note that relative error is
c     measured with respect to the average of the magnitudes of the
c     solution at the beginning and end of the step.
c
      eeoet=0.0d0
      do 250 k=1,neqn
        et=dabs(y(k))+dabs(f1(k))+ae
        if (et .gt. 0.0d0) go to 240
c
c       inappropriate error tolerance
        iflag=5
        return
c
  240   ee=dabs((-2090.0d0*yp(k)+(21970.0d0*f3(k)-15048.0d0*f4(k)))+
     1                        (22528.0d0*f2(k)-27360.0d0*f5(k)))
  250   eeoet=dmax1(eeoet,ee/et)
c
      esttol=dabs(h)*eeoet*scale/752400.0d0
c
      if (esttol .le. 1.0d0) go to 260
c
c
c     unsuccessful step
c                       reduce the stepsize , try again
c                       the decrease is limited to a factor of 1/10
c
      hfaild= .true.
      output= .false.
      s=0.1d0
      if (esttol .lt. 59049.0d0) s=0.9d0/esttol**0.2d0
      h=s*h
      if (dabs(h) .gt. hmin) go to 200
c
c     requested error unattainable at smallest allowable stepsize
      iflag=6
      kflag=6
      return
c
c
c     successful step
c                        store solution at t+h
c                        and evaluate derivatives there
c
  260 t=t+h
      do 270 k=1,neqn
  270   y(k)=f1(k)
      a=t
      call f(a,y,yp)
      nfe=nfe+1
c
c
c                       choose next stepsize
c                       the increase is limited to a factor of 5
c                       if step failure has just occurred, next
c                          stepsize is not allowed to increase
c
      s=5.0d0
      if (esttol .gt. 1.889568d-4) s=0.9d0/esttol**0.2d0
      if (hfaild) s=dmin1(s,1.0d0)
      h=dsign(dmax1(s*dabs(h),hmin),h)
c
c     end of core integrator
c
c
c     should we take another step
c
      if (output) go to 300
      if (iflag .gt. 0) go to 100
c
c
c     integration successfully completed
c
c     one-step mode
      iflag=-2
      return
c
c     interval mode
  300 t=tout
      iflag=2
      return
c
      end
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c    fehl integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt=f(t,y(1),---,y(neqn))
c    where the initial values y(i) and the initial derivatives
c    yp(i) are specified at the starting point t. fehl advances
c    the solution over the fixed step h and returns
c    the fifth order (sixth order accurate locally) solution
c    approximation at t+h in array s(i).
c    f1,---,f5 are arrays of dimension neqn which are needed
c    for internal storage.
c    the formulas have been grouped to control loss of significance.
c    fehl should be called with an h not smaller than 13 units of
c    roundoff in t so that the various independent arguments can be
c    distinguished.
c
c
      integer  neqn
      double precision  y(neqn),t,h,yp(neqn),f1(neqn),f2(neqn),
     1  f3(neqn),f4(neqn),f5(neqn),s(neqn)
c
      double precision  ch
      integer  k
c
      ch=h/4.0d0
      do 221 k=1,neqn
  221   f5(k)=y(k)+ch*yp(k)
      call f(t+ch,f5,f1)
c
      ch=3.0d0*h/32.0d0
      do 222 k=1,neqn
  222   f5(k)=y(k)+ch*(yp(k)+3.0d0*f1(k))
      call f(t+3.0d0*h/8.0d0,f5,f2)
c
      ch=h/2197.0d0
      do 223 k=1,neqn
  223   f5(k)=y(k)+ch*(1932.0d0*yp(k)+(7296.0d0*f2(k)-7200.0d0*f1(k)))
      call f(t+12.0d0*h/13.0d0,f5,f3)
c
      ch=h/4104.0d0
      do 224 k=1,neqn
  224   f5(k)=y(k)+ch*((8341.0d0*yp(k)-845.0d0*f3(k))+
     1                            (29440.0d0*f2(k)-32832.0d0*f1(k)))
      call f(t+h,f5,f4)
c
      ch=h/20520.0d0
      do 225 k=1,neqn
  225   f1(k)=y(k)+ch*((-6080.0d0*yp(k)+(9295.0d0*f3(k)-
     1         5643.0d0*f4(k)))+(41040.0d0*f1(k)-28352.0d0*f2(k)))
      call f(t+h/2.0d0,f1,f5)
c
c     compute approximate solution at t+h
c
      ch=h/7618050.0d0
      do 230 k=1,neqn
  230   s(k)=y(k)+ch*((902880.0d0*yp(k)+(3855735.0d0*f3(k)-
     1        1371249.0d0*f4(k)))+(3953664.0d0*f2(k)+
     2        277020.0d0*f5(k)))
c
      return
      end

      double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end

      double precision function fmin(ax,bx,f,tol)
      double precision ax,bx,f,tol
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result ( .ge. 0.0d0)
c
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a minimum
c
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
c  then fmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
      double precision  fu,fv,fw,fx,x
      double precision  dabs,dsqrt,dsign
c
c  c is the squared inverse of the golden ratio
c
      c = 0.5d0*(3. - dsqrt(5.0d0))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
      eps = 1.0d00
   10 eps = eps/2.0d00
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d00) go to 10
      eps = dsqrt(eps)
c
c  initialization
c
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0d0
      fx = f(x)
      fv = fx
      fw = fx
c
c  main loop starts here
c
   20 xm = 0.5d0*(a + b)
      tol1 = eps*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
c
c  check stopping criterion
c
      if (dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) go to 90
c
c is golden-section necessary
c
      if (dabs(e) .le. tol1) go to 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d00*(q - r)
      if (q .gt. 0.0d0) p = -p
      q =  dabs(q)
      r = e
      e = d
c
c  is parabola acceptable
c
   30 if (dabs(p) .ge. dabs(0.5d0*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
c
c  a parabolic interpolation step
c
      d = p/q
      u = x + d
c
c  f must not be evaluated too close to ax or bx
c
      if ((u - a) .lt. tol2) d = dsign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = dsign(tol1, xm - x)
      go to 50
c
c  a golden-section step
c
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e
c
c  f must not be evaluated too close to x
c
   50 if (dabs(d) .ge. tol1) u = x + d
      if (dabs(d) .lt. tol1) u = x + dsign(tol1, d)
      fu = f(u)
c
c  update  a, b, v, w, and x
c
      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
c
c  end of main loop
c
   90 fmin = x
      return
      end

      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,scale,anorm
      double precision dsqrt,dmax1,dabs,dsign
      logical matu,matv
c
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input.
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output.
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     this is a modified version of a routine from the eispack
c     collection by the nats project
c
c     modified to eliminate machep
c
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      anorm = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    anorm = dmax1(anorm,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
c     .......... for k=n step -1 until 1 do -- ..........
  510 do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            if (dabs(rv1(l)) + anorm .eq. anorm) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            if (dabs(w(l1)) + anorm .eq. anorm) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            if (dabs(f) + anorm .eq. anorm) go to 565
            g = w(i)
            h = dsqrt(f*f+g*g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0d0 * h * y)
         g = dsqrt(f*f+1.0d0)
         f = ((x - z) * (x + z) + h * (y / (f + dsign(g,f)) - h)) / x
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = dsqrt(f*f+h*h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = dsqrt(f*f+h*h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end

      real function urand(iy)
      integer  iy
c
c      urand is a uniform random number generator based  on  theory  and
c  suggestions  given  in  d.e. knuth (1969),  vol  2.   the integer  iy
c  should be initialized to an arbitrary integer prior to the first call
c  to urand.  the calling program should  not  alter  the  value  of  iy
c  between  subsequent calls to urand.  values of urand will be returned
c  in the interval (0,1).
c
      integer  ia,ic,itwo,m2,m,mic
      double precision  halfm
      real  s
      double precision  datan,dsqrt
      data m2/0/,itwo/2/
      if (m2 .ne. 0) go to 20
c
c  if first entry, compute machine integer word length
c
      m = 1
   10 m2 = m
      m = itwo*m2
      if (m .gt. m2) go to 10
      halfm = m2
c
c  compute multiplier and increment for linear congruential method
c
      ia = 8*idint(halfm*datan(1.d0)/8.d0) + 5
      ic = 2*idint(halfm*(0.5d0-dsqrt(3.d0)/6.d0)) + 1
      mic = (m2 - ic) + m2
c
c  s is the scale factor for converting to floating point
c
      s = 0.5/halfm
c
c  compute next random number
c
   20 iy = iy*ia
c
c  the following statement is for computers which do not allow
c  integer overflow on addition
c
      if (iy .gt. mic) iy = (iy - m2) - m2
c
      iy = iy + ic
c
c  the following statement is for computers where the
c  word length for addition is greater than for multiplication
c
      if (iy/2 .gt. m2) iy = (iy - m2) - m2
c
c  the following statement is for computers where integer
c  overflow affects the sign bit
c
      if (iy .lt. 0) iy = (iy + m2) + m2
      urand = float(iy)*s
      return
      end
