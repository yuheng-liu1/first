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
!   This program provides a template for computing the numerical error 
!    between implementations for a user-defined hyperelastic constitutive
!    model in terms of Cauchy-Green invariants and also equivalently 
!    defined in terms of principal stretches.
!
!   Defined in the material configuration in terms of the 2nd Piola-Kirchhoff
!    stress and the material elasticity tensor.
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
!    Use breakpoints or print statements for outputs
!
!    General arrays are defined
!
      implicit none
      real(8) J,PI,W,U,SSE,I1,I2,Ibar1,Ibar2,trC2,&
          dWdI1,dWdI2,dWdI11,dWdI22,dWdI12,kc,kd,ka1,ka2,&
          kb1,kb2,kb3,kb4,kde1,kde2,kde3,kde4,kde5,kde6,kde7,kde8
      real(8), dimension(3,3) :: defg,F,FT,C,C2,Cin,Q,R1,R2,R3
      real(8), dimension(6) :: Cv,Cinv,Iden,kiso,VSC,k2pk
      real(8), dimension(6,6) :: Iddy,CdyId,IddyC,CidyId,IddyCi,CdyC,&
          CidyC,CdyCi,Cvdy,Cv4,Idy4,a_vol,a_iso,ketens
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
      real(8) kprops(n),kn1,kn2,etc...
!
!   Numerically assign material constants
!
      kprops(n)  = "insert magnitude"d0
      kprops(n+1)  = "insert magnitude"d0
!
!    Assign material coefficients
!
      "symbol 1" = kprops(1)
      "symbol 2" = kprops(2)
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
!      defg(1,1) = 1.5d0
!      defg(1,2) = 0.2d0
!      defg(1,3) = 0.2d0
!      defg(2,1) = 0d0
!      defg(2,2) = 0.8d0
!      defg(2,3) = 0.2d0
!      defg(3,1) = 0d0
!      defg(3,2) = 0d0
!      defg(3,3) = 0.8d0
!
!    Define defg(3,3) as a homogeneous deformation
!
!    Stretch magnitude
!
      defg(1,1) = 2d0
!
!      defg(1,2) = 0d0
!      defg(1,3) = 0d0
!      defg(2,1) = 0d0
!      defg(2,3) = 0d0
!      defg(3,1) = 0d0
!      defg(3,2) = 0d0
!    Uniaxial
      defg(2,2) = (one/(defg(1,1)**half)) + 1d0*10**(0d0-(k3*1d0)) !perturbation
      defg(3,3) = (one/(defg(1,1)**half)) - 1d0*10**(0d0-(k3*1d0)) !perturbation
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
!      F = defg
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
      C = matmul(FT, F)
!
!    Calculate the inverse of C
!
      call kin(C, Cin)
!
!    Calculate C squared: C2(3,3)=CC
!    
      C2 = matmul(C, C)
!
!    Store C in Voigt notation, Cv(6)
!
      Cv(1) = C(1,1)
      Cv(2) = C(2,2)
      Cv(3) = C(3,3)
      Cv(4) = C(1,2)
      Cv(5) = C(1,3)
      Cv(6) = C(2,3)
!
!    Store Cinv in Voigt notation, Cv(6)
!
      Cinv(1) = Cin(1,1)
      Cinv(2) = Cin(2,2)
      Cinv(3) = Cin(3,3)
      Cinv(4) = Cin(1,2)
      Cinv(5) = Cin(1,3)
      Cinv(6) = Cin(2,3)
!
!    Calculate the first invariant of C, I1
!
      I1 = C(1,1) + C(2,2) + C(3,3)
!
!    Calculate the trace of C, trC2
!
      trC2 = C2(1,1) + C2(2,2) + C2(3,3)
!
!    Calculate the second invariant of C, I2
!
      I2 = half*((I1**two)-trC2)
!
!    Calculate the first & second invariants of C, I1 & I2
!
      Ibar1 = (J**(-two/three))*I1
      Ibar2 = (J**(-four/three))*I2
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
      W = "isochoric energy"
!
!     Defivne the Volumetric energy, U
!
      U = "volumetric energy"
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      SSE = W + U
!
!    Calculate isochoric stress coefficient, dWdI1=dW/dIbar1
!
      dWdI1 = "isochoric stress coefficient"
!
!    Calculate isochoric stress coefficient, dWdI2=dW/dIbar2
!
      dWdI2 = "isochoric stress coefficient"
!
!    Calculate isochoric elasticity coefficient, dWdI11=d^2W/dIbar1^2
!
      dWdI11 = "isochoric elasticity coefficient"
!
!    Calculate isochoric elasticity coefficient, dWdI22=d^2W/dIbar2^2
!
      dWdI22 = "isochoric elasticity coefficient"
!
!    Calculate isochoric elasticity coefficient, dWdI12=d^2W/dIbar1dIbar2
!
      dWdI12 = "isochoric elasticity coefficient"
!
!    Calculate volumetric stress coefficient, kc=dU/dJ
!
      kc = "volumetric stress coefficient"
!
!    Calculate volumetric elasticity coefficient, kd=d^2U/dJ^2
!
      kd = "volumetric elasticity coefficient"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!              CALCULATE THE 2ND PIOLA-KIRCHHOFF STRESS
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
       kiso(k1) = ka1*(J**(-two/three))*(Iden(k1)-(third*I1*Cinv(k1))) - &
            ka2*(J**(-four/three))*(Cv(k1)-(third*trC2*Cinv(k1)))
      end do
!
!    Calculate the volumetric Kirchoff Stress, VSC(6)
!
      do k1=1,6
       VSC(k1) = kc*J*Cinv(k1)
      end do
!
!    Calculate the total 2nd P-K Stress Tensor, CST(6)
!
      do k1=1,6
       k2pk(k1) = kiso(k1) + VSC(k1)
      end do
!
!***********************************************************************
!                   CALCULATE THE ELASTICITY MODULI
!***********************************************************************
!
!     Calculate the dyadic products and symmetric dyadic products
!
      call m_dyad(Iden, Iden, Iddy)
      call m_dyad(Cv, Iden, CdyId)
      call m_dyad(Iden, Cv, IddyC)
      call m_dyad(Cinv, Iden, CidyId)
      call m_dyad(Iden, Cinv, IddyCi)
      call m_dyad(Cv, Cv, CdyC)
      call m_dyad(Cinv, Cv, CidyC)
      call m_dyad(Cv, Cinv, CdyCi)
      call m_dyad(Cinv, Cinv, Cvdy)
      call m_symdy(Cinv, Cinv, Cv4)
      call m_symdy(Iden, Iden, Idy4)
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
      kde1 = kb1*J**(-four/three)
      kde2 = -kb2*J**(-two)
      kde3 = (-third*kb1*I1*J**(-four/three)) + &
          (third*kb2*trC2*J**(-two)) - ((two/three)*ka1*J**(-two/three))
      kde4 = kb3*J**(-8d0/3d0)
      kde5 = (third*kb2*I1*J**(-two)) - &
          (third*kb3*trC2*J**(-8d0/3d0)) + &
          (third*kb4*J**(-4d0/3d0)) + ((two/three)*ka2*J**(-four/three))
      kde6 = (1d0/9d0)*((kb1*I1*I1*J**(-4d0/3d0)) - &
          (two*kb2*I1*trC2*J**(-two)) + &
          (kb3*(trC2**two)*J**(-8d0/3d0)) - &
          (kb4*trC2*J**(-4d0/3d0)) + &
          (two*ka1*I1*J**(-2d0/3d0)) - &
          (two*ka2*trC2*J**(-4d0/3d0)))
      kde7 = (2d0/3d0)*((ka1*I1*J**(-2d0/3d0)) - &
          (ka2*trC2*J**(-4d0/3d0)))
      kde8 = -kb4*J**(-4d0/3d0)
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_iso(k1,k2) = ((kde1*(Iddy(k1,k2))) + &
              (kde2*(CdyId(k1,k2)+IddyC(k1,k2))) + &
              (kde3*(CidyId(k1,k2)+IddyCi(k1,k2))) + &
              (kde4*(CdyC(k1,k2))) + &
              (kde5*(CidyC(k1,k2)+CdyCi(k1,k2))) + &
              (kde6*(Cvdy(k1,k2))) + &
              (kde7*(Cv4(k1,k2))) + &
              (kde8*(Idy4(k1,k2))))
       end do
      end do
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = (J*(kc+(kd*J))*Cvdy(k1,k2)) - &
                       (two*kc*J*(Cv4(k1,k2)))
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
!   This subroutine computes the 2nd Piola-Kirchhoff stress and elasticity 
!    tensor defined in terms of principal stretches.
!
      implicit none
      real(8) J,tol1,kc,kd,PI,W,U,SSE
      real(8), dimension(3) :: la2,la,la2bar,labar,dWdLa,B_a
      real(8), dimension(3,3) :: kdefg,F,FT,C,Cinv,Na,Id,d2WLa2,Ga,&
        Q,R1,R2,R3
      real(8), dimension(6) :: Cv,N11,N22,N33,Iden,kiso,VSC,k2pk,N12,&
        N13,N23
      real(8), dimension(6,6) :: Cvdy,Cv4,a_vol,N11dy11,N22dy11,N33dy11,&
        N11dy22,N22dy22,N33dy22,N11dy33,N22dy33,N33dy33,N12dy12,&
        N13dy13,N23dy23,a_iso,ketens
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
      real(8) props(kn),kn1,kn2,etc...
!
!    Assign material coefficients
!
      "symbol 1" = kprops(1)
      "symbol 2" = kprops(2)
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
      C = matmul(FT, F)
!
!    Calculate the inverse of C
!
      call kin(C, Cinv)
!
!    Store Cinv in Voigt notation, Cv(6)
!
      Cv(1) = Cinv(1,1)
      Cv(2) = Cinv(2,2)
      Cv(3) = Cinv(3,3)
      Cv(4) = Cinv(1,2)
      Cv(5) = Cinv(1,3)
      Cv(6) = Cinv(2,3)
!
!    Calculate eigenvalues la2(a) and eigenvectors n_a(3,3) of b
!     (Note: Jacobian algorithm destroys upper triangular components 12,13,23)
!
      call DSYEVJ3(C,Na,la2)
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
!     Calculate the isochoric first invariant
!
      Ibar1 = La2bar(1) + La2bar(2) + La2bar(3)
!
!     Define the isochoric energy, W
!
      W = "isochoric energy"
!
!     Define the Volumetric energy, U
!
      U = "volumetric energy"
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      SSE = W + U
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch labar(k1)
!
      do k1=1,3
      dWdLa(k1) = "first derivative"
      end do
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch labar(k1)
!
      do k1=1,3
       do k2=1,3
        if (k1.eq.k2) then
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a = b
!
         d2WLa2(k1,k2) = "second derivative for a=b"
!
        else
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a /= b
!     (a not equal b)            
!
         d2WLa2(k1,k2) = "second derivative for a/=b"
!
        end if
       end do
      end do
!
!    Calculate volumetric stress coefficient, kc=dU/dJ
!
      kc = "volumetric stress coefficient"
!
!    Calculate volumetric elasticity coefficient, kd=d^2U/dJ^2
!
      kd = "volumetric elasticity coefficient"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!              CALCULATE THE 2ND PIOLA-KIRCHHOFF STRESS
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
!    N11
      N11(1) = Na(1,1)*Na(1,1)
      N11(2) = Na(2,1)*Na(2,1)
      N11(3) = Na(3,1)*Na(3,1)
      N11(4) = Na(1,1)*Na(2,1)
      N11(5) = Na(1,1)*Na(3,1)
      N11(6) = Na(2,1)*Na(3,1)
!    N22
      N22(1) = Na(1,2)*Na(1,2)
      N22(2) = Na(2,2)*Na(2,2)
      N22(3) = Na(3,2)*Na(3,2)
      N22(4) = Na(1,2)*Na(2,2)
      N22(5) = Na(1,2)*Na(3,2)
      N22(6) = Na(2,2)*Na(3,2)
!    N33
      N33(1) = Na(1,3)*Na(1,3)
      N33(2) = Na(2,3)*Na(2,3)
      N33(3) = Na(3,3)*Na(3,3)
      N33(4) = Na(1,3)*Na(2,3)
      N33(5) = Na(1,3)*Na(3,3)
      N33(6) = Na(2,3)*Na(3,3)
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
!    Calculate the isochoric 2nd P-K Stress stress, kiso(6)
!
      do k1=1,6
          kiso(k1) = ((one/la2(1))*B_a(1)*N11(k1)) + &
                     ((one/la2(2))*B_a(2)*N22(k1)) + &
                     ((one/la2(3))*B_a(3)*N33(k1))
      end do
!
!    Calculate the volumetric 2nd P-K Stress, VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Cv(k1)
      end do
!
!    Calculate the total 2nd P-K Stress Tensor, CST(6)
!
      do k1=1,6
       k2pk(k1) = kiso(k1) + VSC(k1)
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
!    Na12:
      N12(1) = (Na(1,1)*Na(1,2))
      N12(2) = (Na(2,1)*Na(2,2))
      N12(3) = (Na(3,1)*Na(3,2))
      N12(4) = 0.5d0*((Na(1,1)*Na(2,2)) + (Na(1,2)*Na(2,1)))
      N12(5) = 0.5d0*((Na(1,1)*Na(3,2)) + (Na(1,2)*Na(3,1)))
      N12(6) = 0.5d0*((Na(2,1)*Na(3,2)) + (Na(2,2)*Na(3,1)))
!    Na13:
      N13(1) = (Na(1,1)*Na(1,3))
      N13(2) = (Na(2,1)*Na(2,3))
      N13(3) = (Na(3,1)*Na(3,3))
      N13(4) = 0.5d0*((Na(1,1)*Na(2,3)) + (Na(1,3)*Na(2,1)))
      N13(5) = 0.5d0*((Na(1,1)*Na(3,3)) + (Na(1,3)*Na(3,1)))
      N13(6) = 0.5d0*((Na(2,1)*Na(3,3)) + (Na(2,3)*Na(3,1)))
!    Na23:
      N23(1) = (Na(1,2)*Na(1,3))
      N23(2) = (Na(2,2)*Na(2,3))
      N23(3) = (Na(3,2)*Na(3,3))
      N23(4) = 0.5d0*((Na(1,2)*Na(2,3)) + (Na(1,3)*Na(2,2)))
      N23(5) = 0.5d0*((Na(1,2)*Na(3,3)) + (Na(1,3)*Na(3,2)))
      N23(6) = 0.5d0*((Na(2,2)*Na(3,3)) + (Na(2,3)*Na(3,2)))
!
!    Calculate the fourth-order eigenvector dyadic products
!
      call m_dyad(N11, N11, N11dy11)
      call m_dyad(N22, N11, N22dy11)
      call m_dyad(N33, N11, N33dy11)
      call m_dyad(N11, N22, N11dy22)
      call m_dyad(N22, N22, N22dy22)
      call m_dyad(N33, N22, N33dy22)
      call m_dyad(N11, N33, N11dy33)
      call m_dyad(N22, N33, N22dy33)
      call m_dyad(N33, N33, N33dy33)
!
      call m_dyad(N12, N12, N12dy12)
!
      call m_dyad(N13, N13, N13dy13)
!
      call m_dyad(N23, N23, N23dy23)
!
!     Calculate the symmetric dyad of Cv and Cv
!
      call m_symdy(Cv, Cv, Cv4)
!
!     Calculate the dyad of Cv and Cv
!
      call m_dyad(Cv, Cv, Cvdy)
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = (J*(kc+(kd*J))*Cvdy(k1,k2)) - &
                        (two*kc*J*(Cv4(k1,k2)))
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
            (((one/la2(1))*(one/la2(1))*Ga(1,1)*N11dy11(k1,k2)) + &
             ((one/la2(1))*(one/la2(2))*Ga(1,2)*N11dy22(k1,k2)) + &
             ((one/la2(1))*(one/la2(3))*Ga(1,3)*N11dy33(k1,k2)) + &
             ((one/la2(2))*(one/la2(1))*Ga(2,1)*N22dy11(k1,k2)) + &
             ((one/la2(2))*(one/la2(2))*Ga(2,2)*N22dy22(k1,k2)) + &
             ((one/la2(2))*(one/la2(3))*Ga(2,3)*N22dy33(k1,k2)) + &
             ((one/la2(3))*(one/la2(1))*Ga(3,1)*N33dy11(k1,k2)) + &
             ((one/la2(3))*(one/la2(2))*Ga(3,2)*N33dy22(k1,k2)) + &
             ((one/la2(3))*(one/la2(3))*Ga(3,3)*N33dy33(k1,k2)))- &
             two*(((one/(la2(1)**two))*B_a(1)*N11dy11(k1,k2)) + &
             ((one/(la2(2)**two))*B_a(2)*N22dy22(k1,k2)) + &
             ((one/(la2(3)**two))*B_a(3)*N33dy33(k1,k2)))
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
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(1))*(one/la2(2)))*Ga(1,2))) * &
             two*(N12dy12(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(2))*(one/la2(1)))*Ga(2,1))) * &
             two*(N12dy12(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(1))*(one/la2(3)))*Ga(1,3))) * &
             two*(N13dy13(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(3))*(one/la2(1)))*Ga(3,1))) * &
             two*(N13dy13(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(2))*(one/la2(3)))*Ga(2,3))) * &
             two*(N23dy23(k1,k2))) + &
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(3))*(one/la2(2)))*Ga(3,2))) * &
             two*(N23dy23(k1,k2)))
          end do
         end do
!
!     If La2(1)=la2(2) AND If La2(1)=la2(3)
!
       else if ((abs(la(1)-la(3))).lt.tol1) then
           PRINT *, "1=2;1=3"
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ &
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(1))*(one/la2(2)))*Ga(1,2))) * &
             two*(N12dy12(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(2))*(one/la2(1)))*Ga(2,1))) * &
             two*(N12dy12(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(1))*(one/la2(3)))*Ga(1,3))) * &
             two*(N13dy13(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(3))*(one/la2(1)))*Ga(3,1))) * &
             two*(N13dy13(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(2))*(one/la2(3)))*Ga(2,3))) * &
             two*(N23dy23(k1,k2))) + &
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(3))*(one/la2(2)))*Ga(3,2))) * &
             two*(N23dy23(k1,k2)))
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
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(1))*(one/la2(2)))*Ga(1,2))) * &
             two*(N12dy12(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(2))*(one/la2(1)))*Ga(2,1))) * &
             two*(N12dy12(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(1)/la2(1)))/(la2(3)-la2(1)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(3)/la2(3)))/(la2(1)-la2(3)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(2)/la2(2)))/(la2(3)-la2(2)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2))) + &
             ((((B_a(2)/la2(2))-(B_a(3)/la2(3)))/(la2(2)-la2(3)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2)))
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
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(1))*(one/la2(2)))*Ga(1,2))) * &
             two*(N12dy12(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(2))*(one/la2(1)))*Ga(2,1))) * &
             two*(N12dy12(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(1))*(one/la2(3)))*Ga(1,3))) * &
             two*(N13dy13(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(3))*(one/la2(1)))*Ga(3,1))) * &
             two*(N13dy13(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(2))*(one/la2(3)))*Ga(2,3))) * &
             two*(N23dy23(k1,k2))) + &
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(3))*(one/la2(2)))*Ga(3,2))) * &
             two*(N23dy23(k1,k2)))
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
             ((((B_a(2)/la2(2))-(B_a(1)/la2(1)))/(la2(2)-la2(1)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(2)/la2(2)))/(la2(1)-la2(2)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(1))*(one/la2(3)))*Ga(1,3))) * &
             two*(N13dy13(k1,k2))) + &
             ((((one/(la2(1)**2))*((half*Ga(1,1))-B_a(1))) - &
             (half*((one/la2(3))*(one/la2(1)))*Ga(3,1))) * &
             two*(N13dy13(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(2)/la2(2)))/(la2(3)-la2(2)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2))) + &
             ((((B_a(2)/la2(2))-(B_a(3)/la2(3)))/(la2(2)-la2(3)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2)))
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
             ((((B_a(2)/la2(2))-(B_a(1)/la2(1)))/(la2(2)-la2(1)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(2)/la2(2)))/(la2(1)-la2(2)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(1)/la2(1)))/(la2(3)-la2(1)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(3)/la2(3)))/(la2(1)-la2(3)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
!
             ((((one/(la2(3)**2))*((half*Ga(3,3))-B_a(3))) - &
             (half*((one/la2(2))*(one/la2(3)))*Ga(2,3))) * &
             two*(N23dy23(k1,k2))) + &
             ((((one/(la2(2)**2))*((half*Ga(2,2))-B_a(2))) - &
             (half*((one/la2(3))*(one/la2(2)))*Ga(3,2))) * &
             two*(N23dy23(k1,k2)))
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
             ((((B_a(2)/la2(2))-(B_a(1)/la2(1)))/(la2(2)-la2(1)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(2)/la2(2)))/(la2(1)-la2(2)))*&
             (N12dy12(k1,k2)+N12dy12(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(1)/la2(1)))/(la2(3)-la2(1)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
             ((((B_a(1)/la2(1))-(B_a(3)/la2(3)))/(la2(1)-la2(3)))*&
             (N13dy13(k1,k2)+N13dy13(k1,k2))) + &
!
             ((((B_a(3)/la2(3))-(B_a(2)/la2(2)))/(la2(3)-la2(2)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2))) + &
             ((((B_a(2)/la2(2))-(B_a(3)/la2(3)))/(la2(2)-la2(3)))*&
             (N23dy23(k1,k2)+N23dy23(k1,k2)))
!
        end do
       end do
      end if
!
      ketens = a_iso + a_vol
!
      strs  = k2pk
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
      subroutine kin(mat_B, kinv_B)
!
!    This subroutine calculates the inverse [DET] of a 3x3 
!    matrix [mat_A].
!
      implicit none
      real(8) :: mat_B(3,3), cof(3,3), inv_B(3,3), kdet, kinv_B(3,3)
      integer kint
!
      kint = 3
!
!    Calculate the determinant of mat_B, kdet
!
      call bdet(mat_B, kint, kdet)
!
!    Calculate the cofactor matrix
!
      cof(1,1) = +(mat_B(2,2)*mat_B(3,3)-mat_B(2,3)*mat_B(3,2))
      cof(1,2) = -(mat_B(2,1)*mat_B(3,3)-mat_B(2,3)*mat_B(3,1))
      cof(1,3) = +(mat_B(2,1)*mat_B(3,2)-mat_B(2,2)*mat_B(3,1))
      cof(2,1) = -(mat_B(1,2)*mat_B(3,3)-mat_B(1,3)*mat_B(3,2))
      cof(2,2) = +(mat_B(1,1)*mat_B(3,3)-mat_B(1,3)*mat_B(3,1))
      cof(2,3) = -(mat_B(1,1)*mat_B(3,2)-mat_B(1,2)*mat_B(3,1))
      cof(3,1) = +(mat_B(1,2)*mat_B(2,3)-mat_B(1,3)*mat_B(2,2))
      cof(3,2) = -(mat_B(1,1)*mat_B(2,3)-mat_B(1,3)*mat_B(2,1))
      cof(3,3) = +(mat_B(1,1)*mat_B(2,2)-mat_B(1,2)*mat_B(2,1))
!
!    Transpose the cofactor matrix and divide by the determinant to give
!     the inverse of mat_B, kinv_B
!
      kinv_B = transpose(cof)/kdet
!
      return
      end subroutine kin
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