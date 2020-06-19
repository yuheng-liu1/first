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
!    UMAT computes the Cauchy stress (stress), spatial elasticity
!     modulus (DDSDDE) (in terms of the Jaumann-rate of the Cauchy stress),
!     and strain energy density (SSE) for a Mooney-Rivlin hyperelastic constitutive 
!     model in terms of isochoric Cauchy-Green invariants using explicit 
!     computation of the eigenvectors.
!
!    User input required in "I1I2" and "kstress"
!
!    Stephen Connolly, 2019.
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
!    Skip calculation if "dummy" step
!
      if(dtime.eq.zero) return
!
      if (nshr.eq.1) then
          call I1I2(dfgrd1,nshr,props,stress,axi,DDSDDE,SSE)
      else if (nshr.eq.3) then
          call I1I2(dfgrd1,nshr,props,stress,DDSDDE,k3d,SSE)
      end if
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
      subroutine I1I2(kdefg,knshr,kPROPS,kST,kEtens,axtens,kSSE)
!
!   This program computes the stress (stress), tangent modulus (KETens)
!    and strain energy density (SSE) for a hyperelastic material model
!    defined in terms of principal invariants; using analytically 
!    derived expressions.
!
!    Inputs: 
!     deformation gradient:   kdefg(3,3)
!     number of shear terms:  k1nshr (integer)
!     material properties:    kPROPS(N)
!
!    Outputs:
!     Cauchy stress tensor:   kST(6)
!     3D Elasticity Tensor:   kEtens(6,6)
!     Axi Elasticity Tensor:  axtens(4,4)
!     Specific Elastic SE:    k1SSE (real)
!
      implicit none
      real(8) J,PI,W,U,kSSE,I1,I2,Ibar1,Ibar2,trb2,
     1    dWdI1,dWdI2,dWdI11,dWdI22,dWdI12,kc,kd,ka1,ka2,
     2    kb1,kb2,kb3,kb4,kde1,kde2,kde3,kde4,kde5,kde6,kde7,kde8
      real(8), dimension(3,3) :: kdefg,F,FT,b,b2,Q,R1,R2,R3
      real(8), dimension(6) :: bv,b2v,Iden,kiso,VSC,kST
      real(8), dimension(6,6) :: bvdy,bdyb2,b2dyb,bdyId,Iddyb,b2dyb2,
     1    Iddyb2,b2dyId,Idy,Id4,bv4,a_vol,a_iso,CSTdyId,IddyCST,ketens
      real(8), dimension(4,4) :: axtens
      integer knshr,k1,k2
!
      real(8) zero,one,two,three,four,five,six,third,nine,half,ninth
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,
     1    five=5d0,six=6d0,nine=9d0,third=(1d0/3d0),half=(1d0/2d0),
     2    ninth=(one/nine))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  1: USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Define the number of material parameters and their symbolic notation
!
      real(8) kprops(2),C10,C01
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
!    Store defg in F(3,3)
!
      F = kdefg
!
!    Calculate the determinant of deformation gradient: J
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
!    Store the isochoric Finger tensor in Voigt notation, bv(6)
!
      bv(1) = b(1,1)*J**(-two/three)
      bv(2) = b(2,2)*J**(-two/three)
      bv(3) = b(3,3)*J**(-two/three)
      bv(4) = b(1,2)*J**(-two/three)
      bv(5) = b(1,3)*J**(-two/three)
      bv(6) = b(2,3)*J**(-two/three)
!
!    Store the isochoric Finger tensor in Voigt notation, bv(6)
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
!!!!!!!!!!!!!!!!!!!!!!!!!  2: USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!                    CALCULATE THE CAUCHY STRESS
!***********************************************************************
!
!    Define the Krocker delta (Identity matrix), Iden(3,3)
!
      do k1=1,3
       Iden(k1) = one
      end do
      do k1=4,6
       Iden(k1) = zero
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
       kiso(k1) = ka1*((bv(k1))-(third*Ibar1*Iden(k1)))
     1      - ka2*((b2v(k1))-(third*trb2*Iden(k1)))
      end do
!
!    Volumetric Cauchy Stress Magnitude: VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the Cauchy Stress Tensor, CST(6)
!
      do k1=1,6
       kST(k1) = (1/J)*kiso(k1) + VSC(k1)
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
      kb1 = (four*(dWdI2+dWdI11)) + (four*(Ibar1**two)*dWdI22) + 
     1       (8d0*Ibar1*dWdI12)
      kb2 = four*((Ibar1*dWdI22)+(dWdI12))
      kb3 = four*dWdI22
      kb4 = four*dWdI2
!
!    Calculate the isochoric elasticity coefficients, kde1 to kde8
!
      kde1 = kb1
      kde2 = -kb2
      kde3 = (-third*kb1*Ibar1) + (third*kb2*trb2) - 
     1        ((two/three)*ka1)
      kde4 = kb3
      kde5 = (third*kb2*Ibar1) - (third*kb3*trb2) + (third*kb4) + 
     1        ((two/three)*ka2)
      kde6 = (1d0/9d0)*((kb1*Ibar1*Ibar1) - (two*kb2*Ibar1*trb2) + 
     1     (kb3*(trb2**two)) - (kb4*trb2) + (two*ka1*Ibar1) - 
     2     (two*ka2*trb2))
      kde7 = (2d0/3d0)*((ka1*Ibar1) - (ka2*trb2))
      kde8 = -kb4
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_iso(k1,k2) = ((kde1*(bvdy(k1,k2))) + 
     1         (kde2*(bdyb2(k1,k2)+b2dyb(k1,k2))) + 
     2         (kde3*(bdyId(k1,k2)+Iddyb(k1,k2))) + 
     3         (kde4*(b2dyb2(k1,k2))) + 
     4         (kde5*(Iddyb2(k1,k2)+b2dyId(k1,k2))) + 
     5         (kde6*(Idy(k1,k2))) + 
     6         (kde7*(Id4(k1,k2))) + 
     7         (kde8*(bv4(k1,k2))))
       end do
      end do
!
!    Calculate the volumetric contribution a_vol
!
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = ((J*(kc+(J*kd))*Idy(k1,k2))) -
     1                  (two*kc*J*Id4(k1,k2))
       end do
      end do
!
!    Convert to Jaumann-rate of Cauchy Stress
!
      call m_symdy(kST,Iden,CSTdyId)
      call m_symdy(Iden,kST,IddyCST)
!
      do k1=1,6
       do k2=1,6
        KETENS(k1,k2) = (one/J)*(a_iso(k1,k2) + a_vol(k1,k2)) +
     1                   CSTdyId(k1,k2) + IddyCST(k1,k2)
       end do
      end do
!
!     If axisymmetric, use axtens(4,4)
!
      do k1=1,4
          do k2=1,4
              axtens(k1,k2) = kEtens(k1,k2)
          end do
      end do
!
      return
      end subroutine I1I2
!
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
!
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
!
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