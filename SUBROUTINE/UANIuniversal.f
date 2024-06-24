! ----------------------------------------------------------------------  !
!! Universal UANISOHYPER_INV for constitutive neural networks            !!
! ----------------------------------------------------------------------  !
!  When using, please cite: 
!  "A universal material subroutine for soft matter system" 
!  M. Peirlinck, J.A. Hurtado, M.K. Rausch, A. Buganza Tepole, E. Kuhl
! ----------------------------------------------------------------------  !
!
      subroutine uanisohyper_inv (aInv, ua, zeta, nFibers, nInv,
     *                            ui1, ui2, ui3, temp, noel,
     *                            cmname, incmpFlag, ihybFlag,
     *                            numStatev, statev,
     *                            numFieldv, fieldv, fieldvInc,
     *                            numProps, props)
c
c     UANISOHYPER_INV for constitutive neural networks
c		- including zeroth invariant manipulation layer
c		- including anisotropic terms
c		- including mixed invariants
c		- including compressibility
c      
      include 'aba_param.inc'
      include 'aba_tcs_param.inc'
c
      character *80 cmname
      dimension aInv(nInv), ua(2), zeta(nFibers*(nFibers-1)/2)
      dimension ui1(nInv), ui2(nInv*(nInv+1)/2)
      dimension ui3(nInv*(nInv+1)/2), statev(numStatev)
      dimension fieldv(numFieldv), fieldvInc(numFieldv)
      dimension props(numProps)
c
      dimension aInv0(15)
c
      parameter ( zero = 0.d0, one = 1.d0, three = 3.d0 )
c
c     for table collection, parameter tables, property tables
      character*80 cTableColl(n_tcsC_TC_size)
      dimension jTableColl(n_tcsI_TC_size)
      character*80 ptName
      parameter (maxParams = 180)
      character*80 cParams(maxParams)
      dimension iParamsDataType(maxParams), iParams(maxParams)
      dimension rParams(maxParams)
c
c     Local arrays to keep track of coupling terms due to derived invariants
c
      parameter (maxNumDerInv = 99) ! Maximum number of possible derived invariants
      parameter (maxNumDerInvUsed = 10) ! Maximum used in a given material definition
      dimension IndxDerived(maxNumDerInv,2)
      dimension a_Derived(maxNumDerInvUsed,15) 
      dimension j_Derived(maxNumDerInvUsed,15) 
c
c     Process Network coefficients
c
      ptName = ''
c
      jErrorTC = 0
      jError = 0
      numParams = 0
      numRows = 0
      numNodes = 0
      numDerived = 0
      call queryTableCollection(jTableColl, cTableColl, jErrorTC)
      if (jErrorTC.eq.0) then
         ptName = 'UNIVERSAL_TAB'
         call queryParameterTable(ptName, numParams, numRows, jError)
         if (jError.eq.0) then 
            numNodes = numRows
         end if
c     Check for derived invariants
         ptName = 'MIXED_INV'
         call queryParameterTable(ptName, numParams, numRows, jError)
         if (jError.eq.0) then 
            numDerived = numRows
            do jDer = 1, numDerived
               call getParameterTableRow(ptName, jDer, numParams, 
     *              iParamsDataType, iParams, rParams, cParams, jError)
               kDerInv = iParams(1)
               IndxDerived(kDerInv,1) = jDer
c     -- set auxiliary arrays
               nDeriv = 0
               do kInv = 1, 15
                  aik = rParams(1+kInv)
                  if ( aik .ne.zero ) then
                     nDeriv = nDeriv + 1
                     a_Derived(jDer,nDeriv) = aik
                     j_Derived(jDer,nDeriv) = kInv
                  end if 
               end do
               IndxDerived(kDerInv,2) = nDeriv
            end do
         end if
      end if
c Initialize strain energy function and array of derivatives
      ua(1) = zero
      do kInv = 1, nInv
         indx2 = indx(kInv,kInv)
         ui1(kInv)=zero
         ui2(indx2)=zero
      end do
c    
c Set array of invariants in reference configuration, aInv0
c
      aInv0(1) = three
      aInv0(2) = three
      do kInv = 3, nInv
         aInv0(kInv) = one
      end do
      if (nFibers.gt.1) then
         aInv0(6) = zeta(1)
         aInv0(7) = zeta(1)
         if (nFibers.gt.2) then
            aInv0(10) = zeta(2)
            aInv0(11) = zeta(2)
            aInv0(12) = zeta(3)
            aInv0(13) = zeta(3)
         end if
      end if
c
c Add contribution from each active neuron
c
      ptName = 'UNIVERSAL_TAB'
      do jRow = 1, numNodes
         call getParameterTableRow(ptName, jRow, numParams, 
     *        iParamsDataType, iParams, rParams, cParams, jError)
         kInv = iParams(1)
         kf0  = iParams(2)
         kf1  = iParams(3)
         kf2  = iParams(4)
         w0   = rParams(5)
         w1   = rParams(6)
         w2   = rParams(7)
c     Set shifted invariant 
         if ( kInv .gt. 100 ) then
            kDerInv = mod(kInv,100)
            jDer = IndxDerived(kDerInv,1) 
            nDeriv = IndxDerived(kDerInv,2)
            xInv = zero
            do j = 1, nDeriv 
               aij = a_Derived(jDer,j)
               jInv = j_Derived(jDer,j)
               xInv = xInv + aij*(aInv(jInv) - aInv0(jInv))
            end do
         else 
            xInv = aInv(kInv) - aInv0(kInv)
         end if            
c     Add contribution from this neuron
         lUdevUpdate = 1
         if (kInv.eq.3) then
            lUdevUpdate = 0
         end if
c
         UA1 = zero 
         UA2 = zero 
         duDi = zero 
         d2uDiDi = zero
         UI3 = zero 
c
         call uCANN(xInv, kf0, w0, kf1, w1, kf2, w2,
     *        lUdevUpdate, UA1, UA2, duDi, d2uDiDi )
c
         ua(1) = ua(1) + UA1
         ua(2) = ua(2) + UA2
         if ( kInv .gt. 100 ) then
            kDerInv = mod(kInv,100)
            jDer = IndxDerived(kDerInv,1) 
            nDeriv = IndxDerived(kDerInv,2)
            do j1 = 1, nDeriv 
               aij1 = a_Derived(jDer,j1)
               j1Inv = j_Derived(jDer,j1)
               ui1(j1Inv) = ui1(j1Inv) + duDi * aij1
               do j2 = j1, nDeriv 
                  aij2 = a_Derived(jDer,j2)
                  j2Inv = j_Derived(jDer,j2)
                  indx12 = indx(j1Inv,j2Inv) 
                  ui2(indx12) = ui2(indx12) + d2uDiDi * aij1 * aij2
               end do
            end do
         else 
            ui1(kInv) = ui1(kInv) + duDi
            ui2(indx(kInv,kInv)) = ui2(indx(kInv,kInv)) + d2uDiDi
         end if
c
      end do
c
      return
      end
c
c
      subroutine uCANN(xInv, kf0, w0, kf1, w1, kf2, w2,
     *        lUdevUpdate, UA1, UA2, duDi, d2uDiDi)
c
      include 'aba_param.inc'
c
c     Process zeroth layer
      call uCANN_h0(kf0, xInv, f0, df0, ddf0 )
c     Process first layer
      call uCANN_h1(kf1, w0, f0, f1, df1, ddf1)
c     Process second layer
      call uCANN_h2(kf2, w1, f1, f2, df2, ddf2)

      if (lUdevUpdate.eq.1) then
         UA2 = UA2+w2*f2
      end if 
      UA1 = UA1 + w2*f2
      duDi = duDi+w2*df2*df1*df0
      d2uDiDi = d2uDiDi+w2*((ddf2*df1**2+df2*ddf1)*df0**2+df2*df1*ddf0)
      return
      end
c
c
      subroutine uCANN_h0(kf, x, f, df, ddf)
c
      include 'aba_param.inc'
c
      parameter ( zero = 0.d0, one = 1.d0, half = 0.5d0 )
c
c     Branch based on function type
      if ( kf .eq. 1 ) then
c     f(x) = x 
         f = x
         df = one
         ddf = zero
      else if ( kf .eq. 2 ) then
c     f(x) = <x> = (|x|+x)/2
         f = max(x,zero)
         df = half*(sign(one,x)+one)
         ddf = zero
      else if ( kf .eq. 3 ) then
c     f(x) = |x| 
         f = abs(x)
         df = sign(one,x) 
         ddf = zero
      end if
      
      return
      end
c
c
      subroutine uCANN_h1(kf, w, x, f, df, ddf)
c
      include 'aba_param.inc'
c
      parameter ( zero = 0.d0, two = 2.d0 )
c
c     Branch based on function type
      if ( kf.eq. 1 ) then
c     f(x) = w*x 
         f = w*x
         df = w
         ddf = zero
      else if ( kf.eq. 2 ) then
c     f(x) = (w*x)^2 
         f = w**2 * x**2
         df = w**2 * two*x 
         ddf = w**2 * two
      else if ( kf.ge. 3 ) then
c     f(x) = (w*x)^kf
         f = w**kf * x**kf
         df = kf * w**kf * x**(kf-1) 
         ddf = kf * (kf-1) * w**kf * x**(kf-2)
      end if 
c     
      return
      end
c
c
      subroutine uCANN_h2(kf, w, x, f, df, ddf)
c
      include 'aba_param.inc'
c
      parameter ( zero = 0.d0, one = 1.d0 )
c
c     Branch based on function type
      if ( kf .eq. 1 ) then
c     f(x) = w*x 
         f = w*x
         df = w
         ddf = zero
      else if ( kf.eq. 2 ) then
c     f(x) = exp(w*x)-1 
         f = exp(w*x)-one
         df = w * exp(w*x)
         ddf = w * df
      else if ( kf.eq. 3 ) then
c     f(x) = -ln(1-w*x)
         f = -log(one-w*x)
         df = w / (one-w*x)
         ddf = df**2 
      end if 
c      
      return
      end
c
c
c Maps index from Square to Triangular storage of symmetric matrix 
c
      integer function indx( i, j )
c
      include 'aba_param.inc'
c
      ii = min(i,j)
      jj = max(i,j)
c
      indx = ii + jj*(jj-1)/2
c
      return
      end