! **********************************************************************
!
! Copyright 2019 Stephen Connolly
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! Programs tested using Abaqus/Standard v2016 HF6, Microsoft Visual Studio 2010 
!      and Intel Parallel Studio XE 2013.
!
!    UMAT for extended-tube model to compute the Cauchy stress (stress), 
!     spatial elasticity modulus (DDSDDE) (in terms of the Jaumann-rate
!     of the Cauchy stress), and strain energy density (SSE) for a 
!     hyperelastic constitutive model in terms of isochoric principal
!     stretches using explicit computation of the eigenvectors.
!
!    User input required in "la_sub" and "kstress"
!
!    Stephen Connolly, December 2019.
!      Tested and implemented in Abaqus 2016 HF6
!      Cannot be used for plane stress
!
! **********************************************************************
!      
!     OPTIONAL ADDITIONAL VARIABLES TO BE DEFINED BY USER
!
! **********************************************************************
!
      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, 
     1  ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, 
     2  dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, 
     3  nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, 
     4  npt, layer, kspt, kstep, kinc)
!
      include 'aba_param.inc'
!
      character*80 cmname
!
      real(8) stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3, 3),
     3 dfgrd0(3, 3), dfgrd1(3, 3), time(2)
!
!    Define arrays for axisymmetric and 3D dummy tangents
!
      real(8) axi(6,6), k3d(4,4)
!
      if (nshr.eq.1) then
          call la_sub(dfgrd1,nshr,props,stress,axi,DDSDDE,SSE)
      else if (nshr.eq.3) then
          call la_sub(dfgrd1,nshr,props,stress,DDSDDE,k3d,SSE)
      end if
!
!     Call dummy for debugging purposes
!
!      if (NPT.EQ.1) then
!       if (NOEL.EQ.1) then
!         read(*,*) dummy
!       end if
!      end if
!
      return
      end subroutine umat
!
!***********************************************************************
      subroutine la_sub(kdefg,knshr,kPROPS,CST,kEtens,axtens,kSSE)
!
!    This program computes the stress (stress), tangent modulus (KETens)
!     and strain energy density (SSE) for a hyperelastic material model
!     defined in terms of principal stretches; using analytically 
!     derived expressions.
!
!    Inputs: 
!     deformation gradient:   kdefg(3,3)
!     number of shear terms:  knshr (integer)
!     material properties:    kPROPS(N)
!
!    Outputs:
!     Cauchy stress tensor:   CST(6)
!     3D Elasticity Tensor:   kEtens(6,6)
!     Axi Elasticity Tensor:  axtens(4,4)
!     Specific Elastic SE:    kSSE (real)
!
      implicit none
      real(8) J,tol1,kc,kd,la2_max,PI,W,U,kSSE
      real(8), dimension(3) :: la2,la2bar,la,labar,dWdLa,B_a
      real(8), dimension(3,3) :: kdefg,F,FT,b,na,Ga,Id,d2WLa2
      real(8), dimension(6) :: bv,n11,n22,n33,Iden,kiso,VSC,CST,n12,
     1   n13,n23
      real(8), dimension(6,6) :: Idy,Id4,a_vol,n11dy11,n22dy11,n33dy11,
     1   n11dy22,n22dy22,n33dy22,n11dy33,n22dy33,n33dy33,n12dy12,
     2   n13dy13,n23dy23,a_iso,CSTdyId,IddyCST,KETENS
      real(8), dimension(4,4) :: axtens
      integer knshr, k1, k2
!
!    Numerical Parameters
!
      real(8) zero, one, two, three, four, nine, half, third, ninth
      parameter(zero=0d0, one=1d0, two=2d0, three=3d0, four=4d0,
     1          nine=9d0, half=5d-1, third=1d0/3d0, ninth=1d0/9d0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Define the number of material parameters and their symbolic notation
!
      real(8) kprops(4),Gc,Ke,Be,De,I1,Ibar1
!
!    Assign material coefficients
!
      Gc = kprops(1)
      Ke = kprops(2)
      Be = kprops(3)
      De = kprops(4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Define tolerance for L’Hôpital’s rule
!
      tol1 = 1d-6
!
!    Calculate the determinant of deformation gradient: J(real)
!
      F = kdefg
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
!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Define uncoupled strain energy density function, its derivatives
!     and any additional variables
!
!    Note: additional variables must be declared if "implicit none" is
!           is used
!
!     Extended-tube model requires first invariant, I1
!
      I1 = bv(1) + bv(2) + bv(3)
      Ibar1 = I1*J**(-two/three)
!
!     Define the isochoric energy, W
!
      W = (((Gc/two)*((((one-(De**two))*(Ibar1-three))/(one-((de**two) 
     1     *(Ibar1-three))))+(LOG(one-((De**two)*(Ibar1-three)))))) + 
     2     ((two*Ke)/(Be**two))*(((labar(1)**(-Be))-one)+ 
     3     ((labar(2)**(-Be))-one)+((labar(3)**(-Be))-one)))
!
!     Define the Volumetric energy, U
!
      U = zero
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      kSSE = W + U
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch labar(k1)
!
      do k1=1,3
      dWdLa(k1) = ((labar(k1)*Gc*(one+((Ibar1-three)*(De**4))- 
     1       (two*(De**2)))) / ((((De**2)*(Ibar1-three))-one)**2)) 
     2       - ((two*Ke*(labar(k1)**(-Be-one)))/Be)
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
         d2WLa2(k1,k2) = ((Gc*(-one+((Ibar1-three)*((-two*la2bar(k1))
     1        +Ibar1-three)*De**6d0)+(((6d0*la2bar(k1))-(three*Ibar1)+
     2        nine)*De**4d0)+(((-four*la2bar(k1))+Ibar1-one)*De**2d0)))/
     3        ((((De**2d0)*(Ibar1-three))-one)**3d0)) + 
     4        ((two*Ke*(labar(k1)**(-Be-two))*(Be+one))/Be)
!
        else
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a /= b
!     (a not equal b)            
!
         d2WLa2(k1,k2) = -(two*Gc*labar(k1)*labar(k2)*(De**2)*(two+ 
     1       ((Ibar1-three)*(De**4))-(three*De**2)))/((((De**2)*(Ibar1 
     2        -three))-one)**3)
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
!                   CALCULATE THE CAUCHY STRESS
!***********************************************************************
!
!    Calculate the stress coefficients B_a(3)
!
      do k1=1,3
          B_a(k1) = dWdLa(k1)*labar(k1) - third*(dWdLa(1)*labar(1) +
     1         dWdLa(2)*labar(2) + dWdLa(3)*labar(3))
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
!    Calculate the volumetric Kirchhoff Stress, VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the Cauchy Stress Tensor, CST(6)
!
      do k1=1,6
       CST(k1) = (one/J)*(kiso(k1) + VSC(k1))
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
        Ga(k1,k2) = ((d2WLa2(k1,k2)*labar(k1)*labar(k2))+(dWdLa(k1) 
     &    *Id(k1,k2)*labar(k1))) + 
     &    (ninth*( 
     &    ((d2WLa2(1,1)*la2bar(1))+(dWdLa(1)*labar(1))) + 
     &    ((d2WLa2(2,2)*la2bar(2))+(dWdLa(2)*labar(2))) + 
     &    ((d2WLa2(3,3)*la2bar(3))+(dWdLa(3)*labar(3))) + 
     &    (two*( 
     &    ((d2WLa2(1,2))*labar(1)*labar(2)) + 
     &    ((d2WLa2(1,3))*labar(1)*labar(3)) + 
     &    ((d2WLa2(2,3))*labar(2)*labar(3)))))) - 
     &    (third*( 
     &    ((d2WLa2(k1,1)*labar(k1)*labar(1))+(dWdLa(1) 
     &    *Id(k1,1)*labar(1))) + 
     &    ((d2WLa2(k1,2)*labar(k1)*labar(2))+(dWdLa(2) 
     &    *Id(k1,2)*labar(2))) + 
     &    ((d2WLa2(k1,3)*labar(k1)*labar(3))+(dWdLa(3) 
     &    *Id(k1,3)*labar(3))) + 
     &    ((d2WLa2(1,k2)*labar(1)*labar(k2))+(dWdLa(k2) 
     &    *Id(1,k2)*labar(k2))) + 
     &    ((d2WLa2(2,k2)*labar(2)*labar(k2))+(dWdLa(k2) 
     &    *Id(2,k2)*labar(k2))) + 
     &    ((d2WLa2(3,k2)*labar(3)*labar(k2))+(dWdLa(k2) 
     &    *Id(3,k2)*labar(k2)))))
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
!    Compute the isochoric tangent modulus: a_iso(6,6)
!
!     First calculate terms independent of eigenvalue similarity
!
      do k1=1,6
       do k2=1,6
        a_iso(k1,k2) = 
     &       ((Ga(1,1)*n11dy11(k1,k2)) + (Ga(1,2)*n11dy22(k1,k2)) + 
     &        (Ga(1,3)*n11dy33(k1,k2)) + (Ga(2,1)*n22dy11(k1,k2)) + 
     &        (Ga(2,2)*n22dy22(k1,k2)) + (Ga(2,3)*n22dy33(k1,k2)) + 
     &        (Ga(3,1)*n33dy11(k1,k2)) + (Ga(3,2)*n33dy22(k1,k2)) + 
     &        (Ga(3,3)*n33dy33(k1,k2)))- 
     &        (two*((B_a(1)*n11dy11(k1,k2))+(B_a(2)*n22dy22(k1,k2)) + 
     &        (B_a(3)*n33dy33(k1,k2))))
       end do
      end do
!
!     If la(1)=la(2) AND If la(2)=la(3)
!
      if ((abs(la(1)-la(2))).lt.tol1) then
       if ((abs(la(2)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(2) AND If la(1)=la(3)
!
       else if ((abs(la(1)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(2)
!
       else
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If la(1)=la(3) AND If la(2)=la(3)
!
      else if((abs(la(1)-la(3))).lt.tol1) then
        if ((abs(la(2)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(3)
!
        else
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If la(2)=la(3)
!
      else if((abs(la(2)-la(3))).lt.tol1) then
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
        end do
       end do
!
!    If eigenvalues are unique
!
      else
       do k1=1,6
        do k2=1,6
         a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
!
        end do
       end do
      end if
!
!    Compute the volumetric tangent modulus: a_vol(6,6)
!    
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = ((J*(kc+(J*kd))*Idy(k1,k2))) -
     1                  (two*kc*J*Id4(k1,k2))
       end do
      end do
!
!    Calculate symmetric dyadic products (geometric tangent contributions)
!
      call m_symdy(Iden,CST,IddyCST)
      call m_symdy(CST,Iden,CSTdyId)
!
!     Calculate the spatial tangent modulus, kEtens(6,6)
!
      do k1=1,6
       do k2=1,6
        kEtens(k1,k2) = (1/J)*(a_iso(k1,k2) + a_vol(k1,k2)) +
     1            (IddyCST(k1,k2)+CSTdyId(k1,k2))
                  
       end do
      end do
!
!     If axisymmetric or plane strain, use axtens(4,4)
!
      do k1=1,4
       do k2=1,4
        axtens(k1,k2) = kEtens(k1,k2)
       end do
      end do
!
      return
      end subroutine la_sub
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
      DET = mat_A(1,1)*mat_A(2,2)*mat_A(3,3) -
     1      mat_A(1,2)*mat_A(2,1)*mat_A(3,3)
      if (nshr.eq.3) then
          DET = DET + mat_A(1,2)*mat_A(2,3)*mat_A(3,1)
     1          + mat_A(1,3)*mat_A(2,1)*mat_A(3,2)
     2          - mat_A(1,1)*mat_A(2,3)*mat_A(3,2)
     3          - mat_A(1,3)*mat_A(2,2)*mat_A(3,1)
      end if
      return
      end subroutine bdet
!**********************************************************************
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
!**********************************************************************
      subroutine m_symdy(mat_B, mat_C, symdy_D)
!
!    This subroutine calculates the symmetric dyadic product [symdy_D]
!    of two symmetric 2nd-order tensors [mat_B & mat_C] in Voigt 
!    notation for an isotropic material.
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
!**********************************************************************
* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
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

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
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
              
*             Apply Jacobi transformation
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

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
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
* End of subroutine DSYEVJ3
* ----------------------------------------------------------------------------