program kTemp
!
! Copyright 2019 Stephen Connolly
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! Programs tested using Microsoft Visual Studio 2010 and Intel Parallel Studio XE 2013.
!
!   This program computes the numerical error between implementations
!    of a hyperelastic constitutive model in terms of Cauchy-Green 
!    invariants and also equivalently defined in terms of principal stretches.
!
!   Defined in the spatial configuration in terms of the Kirchhoff stress
!    and the Oldroyd rate of the Cauchy stress.
!
!   The invariant implementation is computed in the body of the program
!    while the principal stretch implementation is defined in the 
!    subroutine "Lambda".
!
!   User-input: deformation gradient, material parameters, derivatives
!    (in terms of both the strain invariants and principal stretches),
!    and other constitutive model specific variables.
!
!   The do loop related to "k3" may be used if the error is to be computed
!    for a range of perturbed deformations gradients.
!
!   This subroutine uses the Mooney-Rivlin constitutive model.
!
!    Use breakpoints or print statements for outputs
!
!    General arrays are defined
!
      implicit none
      real(8) J,PI,W,U,SSE,I1,I2,Ibar1,Ibar2,trb2,&
          dWdI1,dWdI2,dWdI11,dWdI22,dWdI12,kc,kd,ka1,ka2,&
          kb1,kb2,kb3,kb4,kde1,kde2,kde3,kde4,kde5,kde6,kde7,kde8
      real(8), dimension(3,3) :: defg,F,FT,b,b2,Q,R1,R2,R3
      real(8), dimension(6) :: bv,b2v,Iden,kiso,VSC,CST
      real(8), dimension(6,6) :: bvdy,bdyb2,b2dyb,bdyId,Iddyb,b2dyb2,&
          Iddyb2,b2dyId,Idy,Id4,bv4,a_vol,a_iso,ketens
      real(8) la_strs(6),la_etens(6,6),test,error(18),error1,error2
      integer knshr,k1,k2,k3
!
      real(8) zero,one,two,three,four,five,six,third,nine,half,ninth
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,&
          five=5d0,six=6d0,nine=9d0,third=(1d0/3d0),half=(1d0/2d0),&
          ninth=(one/nine))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Define the number of material parameters and their symbolic notation
!
      real(8) kprops(2),C10,C01
!
!   Numerically assign material constants
!
      kprops(1) = 1d0
      kprops(2) = 0.5d0
!
!    Assign material coefficients
!
      C10 = kprops(1)
      C01 = kprops(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Do loop to perturb the deformation gradient
!
      do k3=1,18
          test = 1d0*10**(0d0-(k3*1d0))
!
!    Numerically assign values for testing program
!
!    Fully define deformation gradient defg(3,3)
!
      defg(1,1) = 0.25d0
      defg(1,2) = 0d0
      defg(1,3) = 0d0
      defg(2,1) = 0d0
      defg(2,2) = 2d0
      defg(2,3) = 0d0
      defg(3,1) = 0d0
      defg(3,2) = 0d0
      defg(3,3) = 2d0
!
!    Define defg(3,3) as a homogeneous deformation
!
!    Stretch magnitude
!
!      defg(1,1) = 0.5d0
!
!      defg(1,2) = 0d0
!      defg(1,3) = 0d0
!      defg(2,1) = 0d0
!      defg(2,3) = 0d0
!      defg(3,1) = 0d0
!      defg(3,2) = 0d0
!    Uniaxial
!      defg(2,2) = (one/(defg(1,1)**half)) !+ 1d0*10**(0d0-(k3*1d0)) !perturbation
!      defg(3,3) = (one/(defg(1,1)**half)) !- 1d0*10**(0d0-(k3*1d0)) !perturbation
!    Pure Shear
!      defg(2,2) = one
!      defg(3,3) = (one/defg(1,1))
!    Equibiaxial
!      defg(2,2) = defg(1,1)
!      defg(3,3) = (one/((defg(1,1))**two))
!
!    Numerically assign knshr(integer)
!
      knshr = 3
!
      PI=4.D0*ATAN(1.D0)
!
      R1(1,1) = cos(PI/4d0)
      R1(2,1) = sin(PI/4d0)
      R1(2,2) = cos(PI/4d0)
      R1(1,2) = -sin(PI/4d0)
      R1(3,1) = 0d0
      R1(3,2) = 0d0
      R1(1,3) = 0d0
      R1(2,3) = 0d0
      R1(3,3) = 1d0
!
      R2(1,1) = cos(PI/3d0)
      R2(1,3) = sin(PI/3d0)
      R2(3,3) = cos(PI/3d0)
      R2(3,1) = -sin(PI/3d0)
      R2(1,2) = 0d0
      R2(2,1) = 0d0
      R2(2,3) = 0d0
      R2(3,2) = 0d0
      R2(2,2) = 1d0
!
      R3(2,2) = cos(PI/6d0)
      R3(3,2) = sin(PI/6d0)
      R3(3,3) = cos(PI/6d0)
      R3(2,3) = -sin(PI/6d0)
      R3(1,2) = 0d0
      R3(2,1) = 0d0
      R3(1,3) = 0d0
      R3(3,1) = 0d0
      R3(1,1) = 1d0
!
      Q=matmul(R1,matmul(R2,R3))
!
      F = matmul(Q,defg)
!
!    Calculate the volume ratio J, determinant of F
!
      call bdet(F, knshr, J)
!
!    Calculate the transpose of the deformation gradient
!
      FT = transpose(F)
!
!    Calculate left Cauchy-Green tensor: b(3,3)=F*(F^T)
!    
      b = matmul(F, FT)
!
!    Calculate b squared: b2(3,3)=bb
!    
      b2 = matmul(b, b)
!
!    Isochoric b in Voigt notation, bv(6)
!
      bv(1) = b(1,1)*J**(-two/three)
      bv(2) = b(2,2)*J**(-two/three)
      bv(3) = b(3,3)*J**(-two/three)
      bv(4) = b(1,2)*J**(-two/three)
      bv(5) = b(1,3)*J**(-two/three)
      bv(6) = b(2,3)*J**(-two/three)
!
!    Isochoric b2 in Voigt notation, b2v(6)
!
      b2v(1) = b2(1,1)*J**(-four/three)
      b2v(2) = b2(2,2)*J**(-four/three)
      b2v(3) = b2(3,3)*J**(-four/three)
      b2v(4) = b2(1,2)*J**(-four/three)
      b2v(5) = b2(1,3)*J**(-four/three)
      b2v(6) = b2(2,3)*J**(-four/three)
!
!    Calculate the isochoric first invariant of b, Ibar1
!
      Ibar1    = bv(1) + bv(2) + bv(3)
!
!    Calculate the trace of b2, trb2
!
      trb2 = b2v(1) + b2v(2) + b2v(3)
!
!    Calculate the isochoric second invariant of b, Ibar2
!
      Ibar2 = half*((Ibar1**two)-trb2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Define uncoupled strain energy density function, its derivatives
!     and any additional variables
!
!    Note: additional variables must be declared if "implicit none" is
!           is used
!
!     Define the isochoric energy, W
!
      W = C10*(Ibar1-three) + C01*(Ibar2-three)
!
!     Define the Volumetric energy, U
!
      U = zero
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      SSE = W + U
!
!    Calculate isochoric stress coefficient, dWdI1=dW/dIbar1
!
      dWdI1 = C10
!
!    Calculate isochoric stress coefficient, dWdI2=dW/dIbar2
!
      dWdI2 = C01
!
!    Calculate isochoric elasticity coefficient, dWdI11=d^2W/dIbar1^2
!
      dWdI11 = zero
!
!    Calculate isochoric elasticity coefficient, dWdI22=d^2W/dIbar2^2
!
      dWdI22 = zero
!
!    Calculate isochoric elasticity coefficient, dWdI12=d^2W/dIbar1dIbar2
!
      dWdI12 = zero
!
!    Calculate volumetric stress coefficient, kc=dU/dJ
!
      kc = zero
!
!    Calculate volumetric elasticity coefficient, kd=d^2U/dJ^2
!
      kd = zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!                    CALCULATE THE KIRCHHOFF STRESS
!***********************************************************************
!
!    Define the second-order identitiy tensor in Voigt notation, Iden(6)
!
      do k1=1,6
       if (k1.lt.4) then
        Iden(k1) = 1d0
       else
        Iden(k1) = 0d0
       end if
      end do
!
!    Calculate the isochoric stress coefficient, ka=2dW/dIbar1
!
      ka1 = (two*dWdI1 + two*Ibar1*dWdI2)
      ka2 = (two*dWdI2)
!
!    Calculate the isochoric 2nd P-K Stress stress, kiso(6)
!
      do k1=1,6
       kiso(k1) = ka1*((bv(k1))-(third*Ibar1*Iden(k1))) &
           - ka2*((b2v(k1))-(third*trb2*Iden(k1)))
      end do
!
!    Volumetric Kirchhoff Stress Magnitude: VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the Kirchhoff Stress Tensor, CST(6)
!
      do k1=1,6
       CST(k1) = kiso(k1) + VSC(k1)
      end do
!
!***********************************************************************
!                   CALCULATE THE ELASTICITY MODULI
!***********************************************************************
!
!     Calculate the dyadic products and symmetric dyadic products
!
      call m_dyad(bv, bv, bvdy)
      call m_dyad(bv, b2v, bdyb2)
      call m_dyad(b2v, bv, b2dyb)
      call m_dyad(bv, Iden, bdyId)
      call m_dyad(Iden, bv, Iddyb)
      call m_dyad(b2v, b2v, b2dyb2)
      call m_dyad(Iden, b2v, Iddyb2)
      call m_dyad(b2v, Iden, b2dyId)
      call m_dyad(Iden, Iden, Idy)
      call m_symdy(Iden, Iden, Id4)
      call m_symdy(bv, bv, bv4)
!
!    Calculate the constants kb1 to kb4
!
      kb1 = (four*(dWdI2+dWdI11)) + (four*(Ibar1**two)*dWdI22) + &
            (8d0*Ibar1*dWdI12)
      kb2 = four*((Ibar1*dWdI22)+(dWdI12))
      kb3 = four*dWdI22
      kb4 = four*dWdI2
!
!    Calculate the isochoric elasticity coefficients, kde1 to kde8
!
      kde1 = kb1
      kde2 = -kb2
      kde3 = (-third*kb1*Ibar1) + (third*kb2*trb2) - &
             ((two/three)*ka1)
      kde4 = kb3
      kde5 = (third*kb2*Ibar1) - (third*kb3*trb2) + (third*kb4) + &
             ((two/three)*ka2)
      kde6 = (1d0/9d0)*((kb1*Ibar1*Ibar1) - (two*kb2*Ibar1*trb2) + &
          (kb3*(trb2**two)) - (kb4*trb2) + (two*ka1*Ibar1) - &
          (two*ka2*trb2))
      kde7 = (2d0/3d0)*((ka1*Ibar1) - (ka2*trb2))
      kde8 = -kb4
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_iso(k1,k2) = ((kde1*(bvdy(k1,k2))) + &
              (kde2*(bdyb2(k1,k2)+b2dyb(k1,k2))) + &
              (kde3*(bdyId(k1,k2)+Iddyb(k1,k2))) + &
              (kde4*(b2dyb2(k1,k2))) + &
              (kde5*(Iddyb2(k1,k2)+b2dyId(k1,k2))) + &
              (kde6*(Idy(k1,k2))) + &
              (kde7*(Id4(k1,k2))) + &
              (kde8*(bv4(k1,k2))))
       end do
      end do
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = ((J*(kc+(J*kd))*Idy(k1,k2))) - &
                       (two*kc*J*Id4(k1,k2))
       end do
      end do
!
      ketens = a_iso + a_vol
!
      call lambda(F,kprops,la_strs,la_etens)
!
      error(k3) = zero
      error1 = zero
      error2 = zero
      do k1=1,6
       do k2=1,6
        error1 = error1 + ((la_etens(k1,k2) - KETENS(k1,k2))**2)
        error2 = error2 + ((KETENS(k1,k2))**2)
       end do
      end do
!
      error(k3) = (error1**half) / (error2**half)
!
      end do
!
end program kTemp
!***********************************************************************
    subroutine lambda(kdefg,props,strs,etens)
!
!   This subroutine computes the spatial Kirchhoff stress and elasticity 
!    tensor in terms of the Oldroyd rate of the Kirchhoff stress with the
!    constitutive model defined in terms of principal stretches.
!
      implicit none
      real(8) J,tol1,kc,kd,PI,W,U,SSE
      real(8), dimension(3) :: la2,la,la2bar,labar,dWdLa,B_a
      real(8), dimension(3,3) :: kdefg,F,FT,b,na,Id,d2WLa2,Ga,Q,R1,R2,R3
      real(8), dimension(6) :: bv,n11,n22,n33,Iden,kiso,VSC,CST,n12,&
        n13,n23
      real(8), dimension(6,6) :: Idy,Id4,a_vol,n11dy11,n22dy11,n33dy11,&
        n11dy22,n22dy22,n33dy22,n11dy33,n22dy33,n33dy33,n12dy12,&
        n13dy13,n23dy23,a_iso,ketens
      real(8) Ibar1,strs(6),etens(6,6)
      integer knshr,k1,k2
!
      real(8) zero,one,two,three,four,five,six,third,nine,half,ninth
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,&
          five=5d0,six=6d0,nine=9d0,third=(1d0/3d0),half=(1d0/2d0),&
          ninth=(one/nine))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Define the number of material parameters and their symbolic notation
!
      real(8) props(2),C10,C01
!
!    Assign material coefficients
!
      C10 = props(1)
      C01 = props(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Numerically assign knshr(integer)
!
      knshr = 3
!
      F = kdefg
!
!    Calculate the volume ratio J, determinant of F
!
      call bdet(F, knshr, J)
!
!    Calculate the transpose of the deformation gradient
!
      FT = transpose(F)
!
!    Calculate left Cauchy-Green tensor: b(3,3)=F*(F^T)
!    
      b = matmul(F, FT)
!
!    Store b in Voigt notation, bv(6)
!
      bv(1) = b(1,1)
      bv(2) = b(2,2)
      bv(3) = b(3,3)
      bv(4) = b(1,2)
      bv(5) = b(1,3)
      bv(6) = b(2,3)
!
!    Calculate eigenvalues la2(a) and eigenvectors n_a(3,3) of b
!     (Note: Jacobian algorithm destroys upper triangular components 12,13,23)
!
      call DSYEVJ3(b,na,la2)
!
!    Calculate the eigenvalues' square-root: la(3)
!
      la(1) = (la2(1)**half)
      la(2) = (la2(2)**half)
      la(3) = (la2(3)**half)
!
!    Calculate the isochoric eigenvalues: la2bar(3)
!
      la2bar(1) = la2(1)*(J**(-two/three))
      la2bar(2) = la2(2)*(J**(-two/three))
      la2bar(3) = la2(3)*(J**(-two/three))
!
!    Calculate their square-root: labar(3)
!
      labar(1) = (la2bar(1)**half)
      labar(2) = (la2bar(2)**half)
      labar(3) = (la2bar(3)**half)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Define uncoupled strain energy density function, its derivatives
!     and any additional variables
!
!    Note: additional variables must be declared if "implicit none" is
!           is used
!
!     Define the isochoric energy, W
!
      W = C10*(la2bar(1)+la2bar(2)+la2bar(3)-three) + &
          C01*((1/la2bar(1))+(1/la2bar(2))+(1/la2bar(3))-three)
!
!     Define the Volumetric energy, U
!
      U = zero
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      SSE = W + U
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch la(k1)
!
      do k1=1,3
      dWdLa(k1) = two*C10*labar(k1) - two*C01*(labar(k1)**-three)
      end do
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch la(k1)
!
      do k1=1,3
       do k2=1,3
        if (k1.eq.k2) then
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a = b
!
         d2WLa2(k1,k2) = two*C10 + 6d0*C01*(labar(k1)**-4d0)
!
        else
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a /= b
!
         d2WLa2(k1,k2) = zero
!
        end if
       end do
      end do
!
!    Calculate volumetric stress coefficient, kc=dU/dJ
!
      kc = zero
!
!    Calculate volumetric elasticity coefficient, kd=d^2U/dJ^2
!
      kd = zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!                    CALCULATE THE KIRCHHOFF STRESS
!***********************************************************************
!
!    Calculate the stress coefficients B_a(3)
!
      do k1=1,3
          B_a(k1) = dWdLa(k1)*labar(k1) - third*(dWdLa(1)*labar(1) + &
              dWdLa(2)*labar(2) + dWdLa(3)*labar(3))
      end do
!
!    Calculate the eigenvector dyadic products (eigenvalue bases)
!
!    n11
      n11(1) = na(1,1)*na(1,1)
      n11(2) = na(2,1)*na(2,1)
      n11(3) = na(3,1)*na(3,1)
      n11(4) = na(1,1)*na(2,1)
      n11(5) = na(1,1)*na(3,1)
      n11(6) = na(2,1)*na(3,1)
!    n22
      n22(1) = na(1,2)*na(1,2)
      n22(2) = na(2,2)*na(2,2)
      n22(3) = na(3,2)*na(3,2)
      n22(4) = na(1,2)*na(2,2)
      n22(5) = na(1,2)*na(3,2)
      n22(6) = na(2,2)*na(3,2)
!    n33
      n33(1) = na(1,3)*na(1,3)
      n33(2) = na(2,3)*na(2,3)
      n33(3) = na(3,3)*na(3,3)
      n33(4) = na(1,3)*na(2,3)
      n33(5) = na(1,3)*na(3,3)
      n33(6) = na(2,3)*na(3,3)
!
!    Define the second-order identitiy tensor, Id(3,3)
!
      do k1=1,3
       do k2=1,3
        if (k1.eq.k2) then
         Id(k1,k2) = 1d0
        else
         Id(k1,k2) = 0d0
        end if
       end do
      end do
!
!    Define the second-order identitiy tensor in Voigt notation, Iden(6)
!
      do k1=1,6
          if (k1.lt.4) then
              Iden(k1) = 1d0
          else
              Iden(k1) = 0d0
          end if
      end do
!
!    Calculate the isochoric Kirchhoff stress, kiso(6)
!
      do k1=1,6
          kiso(k1) = (B_a(1)*n11(k1))+(B_a(2)*n22(k1))+(B_a(3)*n33(k1))
      end do
!
!    Calculate the volumetric Kirchhoff Stress, VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the Kirchhoff Stress Tensor, CST(6)
!
      do k1=1,6
       CST(k1) = kiso(k1) + VSC(k1)
      end do
!
!***********************************************************************
!                   CALCULATE THE ELASTICITY MODULI
!***********************************************************************
!
!    Calculate the isochoric elasticity coefficient, Ga(3,3)
!
      do k1=1,3
       do k2=1,3
        Ga(k1,k2) = ((d2WLa2(k1,k2)*labar(k1)*labar(k2))+(dWdLa(k1) &
         *Id(k1,k2)*labar(k1))) + &
         (ninth*( &
         ((d2WLa2(1,1)*la2bar(1))+(dWdLa(1)*labar(1))) + &
         ((d2WLa2(2,2)*la2bar(2))+(dWdLa(2)*labar(2))) + &
         ((d2WLa2(3,3)*la2bar(3))+(dWdLa(3)*labar(3))) + &
         (two*( &
         ((d2WLa2(1,2))*labar(1)*labar(2)) + &
         ((d2WLa2(1,3))*labar(1)*labar(3)) + &
         ((d2WLa2(2,3))*labar(2)*labar(3)))))) - &
         (third*( &
         ((d2WLa2(k1,1)*labar(k1)*labar(1))+(dWdLa(1) &
         *Id(k1,1)*labar(1))) + &
         ((d2WLa2(k1,2)*labar(k1)*labar(2))+(dWdLa(2) &
         *Id(k1,2)*labar(2))) + &
         ((d2WLa2(k1,3)*labar(k1)*labar(3))+(dWdLa(3) &
         *Id(k1,3)*labar(3))) + &
         ((d2WLa2(1,k2)*labar(1)*labar(k2))+(dWdLa(k2) &
         *Id(1,k2)*labar(k2))) + &
         ((d2WLa2(2,k2)*labar(2)*labar(k2))+(dWdLa(k2) &
         *Id(2,k2)*labar(k2))) + &
         ((d2WLa2(3,k2)*labar(3)*labar(k2))+(dWdLa(k2) &
         *Id(3,k2)*labar(k2)))))
       end do
      end do
!
!    Calculate the "average" eigenvector dyads
!
!    na12:
      n12(1) = (na(1,1)*na(1,2))
      n12(2) = (na(2,1)*na(2,2))
      n12(3) = (na(3,1)*na(3,2))
      n12(4) = 0.5d0*((na(1,1)*na(2,2)) + (na(1,2)*na(2,1)))
      n12(5) = 0.5d0*((na(1,1)*na(3,2)) + (na(1,2)*na(3,1)))
      n12(6) = 0.5d0*((na(2,1)*na(3,2)) + (na(2,2)*na(3,1)))
!    na13:
      n13(1) = (na(1,1)*na(1,3))
      n13(2) = (na(2,1)*na(2,3))
      n13(3) = (na(3,1)*na(3,3))
      n13(4) = 0.5d0*((na(1,1)*na(2,3)) + (na(1,3)*na(2,1)))
      n13(5) = 0.5d0*((na(1,1)*na(3,3)) + (na(1,3)*na(3,1)))
      n13(6) = 0.5d0*((na(2,1)*na(3,3)) + (na(2,3)*na(3,1)))
!    na23:
      n23(1) = (na(1,2)*na(1,3))
      n23(2) = (na(2,2)*na(2,3))
      n23(3) = (na(3,2)*na(3,3))
      n23(4) = 0.5d0*((na(1,2)*na(2,3)) + (na(1,3)*na(2,2)))
      n23(5) = 0.5d0*((na(1,2)*na(3,3)) + (na(1,3)*na(3,2)))
      n23(6) = 0.5d0*((na(2,2)*na(3,3)) + (na(2,3)*na(3,2)))
!
!    Calculate the fourth-order eigenvector dyadic products
!
      call m_dyad(n11, n11, n11dy11)
      call m_dyad(n22, n11, n22dy11)
      call m_dyad(n33, n11, n33dy11)
      call m_dyad(n11, n22, n11dy22)
      call m_dyad(n22, n22, n22dy22)
      call m_dyad(n33, n22, n33dy22)
      call m_dyad(n11, n33, n11dy33)
      call m_dyad(n22, n33, n22dy33)
      call m_dyad(n33, n33, n33dy33)
!
      call m_dyad(n12, n12, n12dy12)
!
      call m_dyad(n13, n13, n13dy13)
!
      call m_dyad(n23, n23, n23dy23)
!
!    Calculate the dyad of Iden and Iden, Idy(6,6)
!
      call m_dyad(Iden, Iden, Idy)
!
!    Calculate the symmetric dyad of Iden and Iden, Id4(6,6)
!
      call m_symdy(Iden, Iden, Id4)
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = ((J*(kc+(J*kd))*Idy(k1,k2))) - &
                       (two*kc*J*Id4(k1,k2))
       end do
      end do
!
!    Compute the isochoric tangent modulus: a_iso(6,6)
!
!     First calculate terms independent of eigenvalue similarity
!
      do k1=1,6
       do k2=1,6
        a_iso(k1,k2) = &
            ((Ga(1,1)*n11dy11(k1,k2)) + (Ga(1,2)*n11dy22(k1,k2)) + &
             (Ga(1,3)*n11dy33(k1,k2)) + (Ga(2,1)*n22dy11(k1,k2)) + &
             (Ga(2,2)*n22dy22(k1,k2)) + (Ga(2,3)*n22dy33(k1,k2)) + &
             (Ga(3,1)*n33dy11(k1,k2)) + (Ga(3,2)*n33dy22(k1,k2)) + &
             (Ga(3,3)*n33dy33(k1,k2)))- &
             (two*((B_a(1)*n11dy11(k1,k2))+(B_a(2)*n22dy22(k1,k2)) + &
             (B_a(3)*n33dy33(k1,k2))))
       end do
      end do
!
!    Define eigenvalue similarity tolerance
!
      tol1 = 1d-6     !With tolerance value
      !tol1 = 1d-0    !With L’Hôpital’s rule
      !tol1 = 0d-50   !Without L’Hôpital’s rule
!
!     If La(1)=la(2) AND If La(2)=la(3)
!
      if ((abs(la(1)-la(2))).lt.tol1) then
       if ((abs(la(2)-la(3))).lt.tol1) then
           PRINT *, "1=2=3"
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ &
             ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + &
             ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + &
!
             ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + &
             ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + &
!
             ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + &
             ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If La(1)=la(2) AND If La(1)=la(3)
!
       else if ((abs(la(1)-la(3))).lt.tol1) then
           PRINT *, "1=2;1=3"
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ &
             ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + &
             ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + &
!
             ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + &
             ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + &
!
             ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + &
             ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If La(1)=la(2)
!
       else
           PRINT *, "1=2"
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + &
             ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + &
             ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + &
!
             ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*&
             (two*n13dy13(k1,k2))) + &
             ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*&
             (two*n13dy13(k1,k2))) + &
!
             ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*&
             (two*n23dy23(k1,k2))) + &
             ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*&
             (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If La(1)=la(3) AND If La(2)=la(3)
!
      else if((abs(la(1)-la(3))).lt.tol1) then
        if ((abs(la(2)-la(3))).lt.tol1) then
           PRINT *, "1=3,2=3"
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ &
             ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + &
             ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + &
!
             ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + &
             ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + &
!
             ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + &
             ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If La(1)=la(3)
!
        else
          PRINT *, "1=3"
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + &
             ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*&
             (two*n12dy12(k1,k2))) + &
             ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*&
             (two*n12dy12(k1,k2))) + &
!
             ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + &
             ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - &
             (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + &
!
             ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*&
             (two*n23dy23(k1,k2))) + &
             ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*&
             (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If La(2)=la(3)
!
      else if((abs(la(2)-la(3))).lt.tol1) then
          PRINT *, "2=3"
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + &
             ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*&
             (two*n12dy12(k1,k2))) + &
             ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*&
             (two*n12dy12(k1,k2))) + &
!
             ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*&
             (two*n13dy13(k1,k2))) + &
             ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*&
             (two*n13dy13(k1,k2))) + &
!
             ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - &
             (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + &
             ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - &
             (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
        end do
       end do
!
!    If eigenvalues are unique
!
      else
          PRINT *, "unique"
       do k1=1,6
        do k2=1,6
         a_iso(k1,k2) = a_iso(k1,k2) + &
             ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*&
             (two*n12dy12(k1,k2))) + &
             ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*&
             (two*n12dy12(k1,k2))) + &
!
             ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*&
             (two*n13dy13(k1,k2))) + &
             ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*&
             (two*n13dy13(k1,k2))) + &
!
             ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*&
             (two*n23dy23(k1,k2))) + &
             ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*&
             (two*n23dy23(k1,k2)))
!
        end do
       end do
      end if
!
      ketens = a_iso + a_vol
!
      strs  = CST
      etens = ketens
!
    end subroutine lambda
!***********************************************************************
      subroutine bdet(mat_A, nshr, DET)
!
!    This subroutine calculates the determinant [DET] of a 3x3 
!    matrix [mat_A].
!
      implicit none
      real(8) :: DET, mat_A(3,3)
      integer, intent(IN) :: nshr
!
      DET = mat_A(1,1)*mat_A(2,2)*mat_A(3,3) - &
           mat_A(1,2)*mat_A(2,1)*mat_A(3,3)
      if (nshr.eq.3) then
          DET = DET + mat_A(1,2)*mat_A(2,3)*mat_A(3,1)&
               + mat_A(1,3)*mat_A(2,1)*mat_A(3,2)&
               - mat_A(1,1)*mat_A(2,3)*mat_A(3,2)&
               - mat_A(1,3)*mat_A(2,2)*mat_A(3,1)
      end if
      return
    end subroutine bdet
!***********************************************************************
      subroutine m_dyad(mat_B, mat_C, dyad_D)
!
!    This subroutine calculates the dyadic product [dyad_D] of two
!    symmetric 2nd-order tensors [mat_B & mat_C] written in Voigt form.
!
      implicit none
      real(8) :: mat_B(6), mat_C(6), dyad_D(6,6)
      integer :: k1,k2
!
      do k1=1, 6
          do k2=1, 6
              dyad_D(k1,k2) = mat_B(k1)*mat_C(k2)
          end do
      end do
      return
      end subroutine m_dyad
!***********************************************************************
     subroutine m_symdy(mat_B, mat_C, symdy_D)
!
!    This subroutine calculates the symmetric dyadic product [symdy_D]
!    of two symmetric 2nd-order tensors [mat_B & mat_C] in Voigt 
!    notation for an isotropic material: D=(1/2)(Bik(x)Cjl + Bil(x)Cjk)
!
      implicit none
      real(8) mat_B(6),mat_C(6),mat_B2(6,6),mat_C2(6,6),symdy_D(6,6)
      integer k1,k2
!
      mat_B2(1,1) = mat_B(1)*mat_C(1)
      mat_B2(1,2) = mat_B(4)*mat_C(4)
      mat_B2(1,3) = mat_B(5)*mat_C(5)
      mat_B2(1,4) = mat_B(1)*mat_C(4)
      mat_B2(1,5) = mat_B(1)*mat_C(5)
      mat_B2(1,6) = mat_B(4)*mat_C(5)
      mat_B2(2,2) = mat_B(2)*mat_C(2)
      mat_B2(2,3) = mat_B(6)*mat_C(6)
      mat_B2(2,4) = mat_B(4)*mat_C(2)
      mat_B2(2,5) = mat_B(4)*mat_C(6)
      mat_B2(2,6) = mat_B(2)*mat_C(6)
      mat_B2(3,3) = mat_B(3)*mat_C(3)
      mat_B2(3,4) = mat_B(5)*mat_C(6)
      mat_B2(3,5) = mat_B(5)*mat_C(3)
      mat_B2(3,6) = mat_B(6)*mat_C(3)
      mat_B2(4,4) = mat_B(1)*mat_C(2)
      mat_B2(4,5) = mat_B(1)*mat_C(6)
      mat_B2(4,6) = mat_B(4)*mat_C(6)
      mat_B2(5,5) = mat_B(1)*mat_C(3)
      mat_B2(5,6) = mat_B(4)*mat_C(3)
      mat_B2(6,6) = mat_B(2)*mat_C(3)
      do k1=1, 6
        do k2=1, k1-1
          mat_B2(k1,k2) = mat_b2(k2,k1)
        end do
      end do
!
      mat_C2(1,1) = mat_B(1)*mat_C(1)
      mat_C2(1,2) = mat_B(4)*mat_C(4)
      mat_C2(1,3) = mat_B(5)*mat_C(5)
      mat_C2(1,4) = mat_B(4)*mat_C(1)
      mat_C2(1,5) = mat_B(5)*mat_C(1)
      mat_C2(1,6) = mat_B(5)*mat_C(4)
      mat_C2(2,2) = mat_B(2)*mat_C(2)
      mat_C2(2,3) = mat_B(6)*mat_C(6)
      mat_C2(2,4) = mat_B(2)*mat_C(4)
      mat_C2(2,5) = mat_B(6)*mat_C(4)
      mat_C2(2,6) = mat_B(6)*mat_C(2)
      mat_C2(3,3) = mat_B(3)*mat_C(3)
      mat_C2(3,4) = mat_B(6)*mat_C(5)
      mat_C2(3,5) = mat_B(3)*mat_C(5)
      mat_C2(3,6) = mat_B(3)*mat_C(6)
      mat_C2(4,4) = mat_B(4)*mat_C(4)
      mat_C2(4,5) = mat_B(5)*mat_C(4)
      mat_C2(4,6) = mat_B(5)*mat_C(2)
      mat_C2(5,5) = mat_B(5)*mat_C(5)
      mat_C2(5,6) = mat_B(5)*mat_C(6)
      mat_C2(6,6) = mat_B(6)*mat_C(6)
      do k1=1, 6
        do k2=1, k1-1
          mat_C2(k1,k2) = mat_C2(k2,k1)
        end do
      end do
!
      do k1=1, 6
          do k2=1, 6
              symdy_D(k1,k2) = (0.5d0*(mat_B2(k1,k2)+mat_C2(k1,k2)))
          end do
      end do
      return
    end subroutine m_symdy
!***********************************************************************
! ----------------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
! ----------------------------------------------------------------------------
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
! ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using the Jacobi algorithm.
! The upper triangular part of A is destroyed during the calculation,
! the diagonal elements are read but not destroyed, and the lower
! triangular elements are not referenced at all.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
!     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

!     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
!     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

!     Initialize Q to the identitity matrix
!     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

!     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

!     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
!     Main iteration loop
      DO 40 I = 1, 50
!       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

!       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) &
                          .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
!             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
!             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

!             Update eigenvectors
!             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
    END SUBROUTINE
! End of subroutine DSYEVJ3
!***********************************************************************