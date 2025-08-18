      module GlobalVar
       implicit none
         integer, parameter :: Nel = 17459, Nip = 8  ! NEED TO BE CHANGED FOR EACH SIMULATIONS
         double precision, target :: GPold(Nel,Nip,4),GPnew(Nel,Nip,4),
     & Volint(Nel,Nip,1)
         double precision incem,Iteration_U
      end module        
      
CC#########################Octree module################################            
      module OctreeAlg
      implicit none
!Define point type
      type point
      real(kind = 8), pointer :: tx, ty, tz, tv, tb
      end type point
      

!Define octree point type
      type octree_node
      type(point), dimension(:), pointer :: points
      type(octree_node), dimension(:), pointer :: children
      real(kind = 8) :: x_min, x_max, y_min, y_max, z_min, z_max
      logical leaf
      end type octree_node
      
      
      type(point) :: center
      real(kind = 8) ::  grad, sumwt, radius, cLt, mw0, ca2
      type(octree_node), pointer :: root
  
!Recursively construct an octree struture
      Contains
      
      recursive subroutine build_octree(node)
      type(octree_node) :: node
      type(point), dimension(:), pointer ::  points,points_1, points_2, 
     &        points_3, points_4,points_5, points_6,points_7, points_8
      integer :: i, j, k, num_np(8)
      real(kind = 8) :: x_mid, y_mid, z_mid
      real(kind = 8) :: dis, weight,  evpw
      
      x_mid = (node%x_min + node%x_max) / 2.0d0
      y_mid = (node%y_min + node%y_max) / 2.0d0  
      z_mid = (node%z_min + node%z_max) / 2.0d0            
                
      if (size(node%points)>=0 .and. size(node%points)<50) then
          
          node%leaf = .true.
    
      else if (size(node%points)>=50) then

          node%leaf = .false.
     
       allocate(node%children(8)) 
       
       call allocate_points_in_range(node%points, x_mid, y_mid, z_mid, 
     &      points_1, points_2, points_3, points_4,points_5,points_6,
     &      points_7,points_8, num_np)

       deallocate(node%points)
       
       do i = 1, 8
          allocate(node%children(i)%points(num_np(i))) 
       end do

       node%children(1)%points = points_1
       node%children(2)%points = points_2
       node%children(3)%points = points_3
       node%children(4)%points = points_4
       node%children(5)%points = points_5
       node%children(6)%points = points_6
       node%children(7)%points = points_7
       node%children(8)%points = points_8
       
       node%children(1)%x_min = node%x_min
       node%children(1)%x_max = x_mid
       node%children(1)%y_min = node%y_min
       node%children(1)%y_max = y_mid
       node%children(1)%z_min = node%z_min
       node%children(1)%z_max = z_mid
       
       node%children(2)%x_min = node%x_min
       node%children(2)%x_max = x_mid
       node%children(2)%y_min = y_mid
       node%children(2)%y_max = node%y_max
       node%children(2)%z_min = node%z_min
       node%children(2)%z_max = z_mid
       
       node%children(3)%x_min = x_mid
       node%children(3)%x_max = node%x_max
       node%children(3)%y_min = node%y_min
       node%children(3)%y_max = y_mid
       node%children(3)%z_min = node%z_min
       node%children(3)%z_max = z_mid
       
       node%children(4)%x_min = x_mid
       node%children(4)%x_max = node%x_max
       node%children(4)%y_min = y_mid
       node%children(4)%y_max = node%y_max
       node%children(4)%z_min = node%z_min
       node%children(4)%z_max = z_mid
       
       node%children(5)%x_min = node%x_min
       node%children(5)%x_max = x_mid
       node%children(5)%y_min = node%y_min
       node%children(5)%y_max = y_mid
       node%children(5)%z_min = z_mid
       node%children(5)%z_max = node%z_max
       
       node%children(6)%x_min = node%x_min
       node%children(6)%x_max = x_mid
       node%children(6)%y_min = y_mid
       node%children(6)%y_max = node%y_max
       node%children(6)%z_min = z_mid
       node%children(6)%z_max = node%z_max
       
       node%children(7)%x_min = x_mid
       node%children(7)%x_max = node%x_max
       node%children(7)%y_min = node%y_min
       node%children(7)%y_max = y_mid
       node%children(7)%z_min = z_mid
       node%children(7)%z_max = node%z_max
       
       node%children(8)%x_min = x_mid
       node%children(8)%x_max = node%x_max
       node%children(8)%y_min = y_mid
       node%children(8)%y_max = node%y_max
       node%children(8)%z_min = z_mid
       node%children(8)%z_max = node%z_max
       
       do i = 1, 8
           call build_octree(node%children(i))
       end do
     
      end if  
      
     
      end subroutine build_octree
     
     
     
      recursive subroutine search_octree(node, x_min, x_max, y_min, 
     &                                   y_max, z_min, z_max)
      type(octree_node), intent(in) :: node
      real(kind = 8), intent(in) :: x_min,x_max,y_min,y_max,z_min,z_max
      real(kind = 8) ::  dis, weight
      integer :: i  
      
!The source code will be publicly accessible via GitHub upon article publication.      

      end subroutine search_octree
     
      subroutine allocate_points_in_range(points, x_mid, y_mid, z_mid, 
     &      points_1, points_2, points_3, points_4,points_5,points_6,
     &      points_7,points_8, num_np)
      type(point),dimension(:), pointer :: points, points_1, points_2, 
     &     points_3, points_4, points_5, points_6, points_7, points_8
      real(kind = 8) :: x_mid, y_mid, z_mid
      integer :: num_points, i, j, num_np(8)

!The source code will be publicly accessible via GitHub upon article publication.
      
      end 

      
      function distance(x, y, z) result(dis)
      real(kind = 8), intent(in) :: x, y, z
      real(kind = 8) :: dis
      
      dis = sqrt ((center%tx - x)**2.0d0 + (center%ty - y)**2.0d0+ 
     &           (center%tz - z)**2.0d0)

      end function distance

      end module OctreeAlg      
CC#########################Octree module################################     
      

      module ModelParam
       implicit none     
       double precision EK,miu,fphi,pi,c0,Epsilon,betaS,Psi,As
       double precision Ftol,Newton,rho,zeta,Maxit,cLt
       integer index_s,index_nl 
       
       contains 
          subroutine Initialize(para,npara)
        !===================================
            double precision para(npara)
            integer npara

! model parameters and algorithm parameters
           EK = para(1)
           miu = para(2)
           c0 = para(3)
           fphi = para(4)
           Psi = para(5)
           As = para(6)
           Epsilon = para(7)
           betaS = para(8)
           pi = para(9) 
           Maxit = para(10)
           rho = para(11)
           zeta = para(12)
           cLt = para(13)
           index_s = para(14)
           index_nl = para(15)
           
          return
          end subroutine Initialize      
      
          end module  

      SUBROUTINE USDFLD (FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      Use GlobalVar
      INCLUDE 'ABA_PARAM.INC'
C
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      CALL GETVRM('IVOL',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     $MATLAYO,LACCFLA)
      FIELD(1) = ARRAY(1)
      Volint(noel,npt,1) = FIELD(1)

C
      RETURN
      end    


      SUBROUTINE UMAT (STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,
     1 DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,
     3 COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     4 KSPT,KSTEP,KINC)
c      INCLUDE 'ABA_PARAM.INC'
!       
      use OctreeAlg
      use GlobalVar
       implicit none
       
       CHARACTER*80 CMNAME

      Integer NTENS, NDI, NSHR, NSTATV, NPROPS, NOEL, NPT,
     1 LAYER, KSPT, KSTEP, KINC, INITTENSION,DTIME
      DOUBLE PRECISION  stress(ntens),statev(nstatv),
     1  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),jstep(4),
     2  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     3  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), DFGRD1(3,3),
     4  SSE,SPD,SCD,RPL,DRPLDT,TEMP,DTEMP,PNEWDT,CELENT

       double precision EK,miu,fphi,pi,c0,Epsilon,betaS,Psi,As,
     $ rho,zeta,Maxit,alpha
       double precision ftol,BK,GK,dphi,gi,normg,ft,
     $ damage,damage_i,dp,epn,epnw,p,q,theta,f_ep
       
       double precision Ivol(6,6),Isym(6,6),delta(6,1),De(6,6),
     $ CTO(6,6),g(8,1),xnew(8,1),xold(8,1),Jacob(8,8),IJacob(8,8),
     $ s(6,1),para(24),dstra(6,1),st(6,1),dx(8,1),A21(1,6),f_s(6,1)
         
      double precision mw,ca1,index_s,index_nl,Amw,Aca1
      integer K1,K2,K3,Iteration,Newton,npara
      
      GPnew(noel,npt,1:3) = Coords(1:3)
      if (time(2).LT.1.D-10) then
          statev = 0.d0
          statev(1:6) = stress(1:6)
          statev(7) = 1.d-8
          GPnew(noel,npt,4) = statev(7)
          GPold=GPnew 
      endif

! determine whether a new incremental step has been initiated
      if ((Noel.Eq.1).and.(Npt.eq.1)) then
          if  (incem.LT.kinc) then  
               GPold = GPnew
               incem = kinc
               Iteration_U = 0.d0
          else
               Iteration_U = Iteration_U + 1.d0
          endif     
      endif
      
      
! Define algorithm parameters
      Epsilon = 0.1d0
      betaS = 0.999d0      
      pi = 3.1415926535897934d0
      Maxit = 20.d0   
      rho = 1.D-4
      zeta = 0.1d0      
      Newton = 20
      Iteration = 1
      npara = 24
      Ftol = 1.D-6
! Read model parameters
      EK = props(1)   
      miu = props(2)
      c0 = props(3)
      fphi = props(4)*pi/180.d0
      Psi = props(5)*pi/180.d0
      As = props(6)
      cLt = props(7)    
      index_s = props(8)
      index_nl = props(9)
!Store material parameters in the para() array
      para(1:6) = props(1:6)
      para(4) = fphi    
      para(5) = Psi    
      para(7) = Epsilon
      para(8) = betaS
      para(9) = pi
      para(10) = Maxit
      para(11) = rho
      para(12) = zeta
      para(13) = cLt 
      para(14) = index_s
      para(15) = index_nl  
! Calculate the elastic stiffness matrix and Read the effective stress at step n
      De = 0.d0
      BK = EK/3.d0/(1.d0-2.d0*miu)
      GK = EK/2.d0/(1.d0+miu)
      call TenMat1(Ivol,Isym,delta)
      De = (3.0d0*BK - 2.0d0*GK)*Ivol+2.0d0*GK*Isym
      stress(1:6) = statev(1:6) 
      s = 0.0d0
      dstra = 0.0d0
      call DimIn(stress,dstran,s,dstra,NTENS)    
      epn = statev(7)
! calculate elastic trial stress and perform loading and unloading judgments
      st = s + matmul(De,dstra)
      xold(1:6,1)=s(1:6,1)
      xold(7,1) = epn
      xold(8,1) = 0.d0
      xnew = xold
      xnew(1:6,1) = st(1:6,1)
      ft = 0.d0
      call yieldfun(xnew,para,npara,ft,f_ep,f_s)
!      xnew = xold
      call cgfun(xnew,xold,para,npara,De,dstra,9,gi,g)
      normg = norm2(g)  
C##***********************************Nonlocal computation*************************************##C
!      call nlfunTra(mw,ca1,ntens,noel,npt,para,npara,time,
!     & nstatv,statev)

      call nlfunOct(mw,ca1,statev,nstatv,COORDS,time)

! ca1 represents the contribution of neighbourhood integration points to nonlocal computation
      epnw = (1.d0-mw)*statev(7) + ca1
      epnw = max(epnw,1.D-8)
      call dam_fun(damage,dp,epnw,para,npara)
!      damage = statev(11)
      if (ft.LE.Ftol) then 
        CTO = De*(1.d0-damage)
        normg = 0.d0
        statev(1:6) = st(1:6,1)
        xnew(1:6,1) = st(1:6,1)*(1.d0-damage)
      else    
      do while ((normg.GT.Ftol).and.(Iteration .LT. Newton))
         call Jacobfun(xnew,xold,para,npara,De,dstra,Jacob)
         call inverse(Jacob,IJacob,8)
         dx = -matmul(IJacob,g)
         call LSM(xnew,xold,dx,para,npara,De,dstra,alpha,g)
         xnew = xnew + alpha*dx 
         call cgfun(xnew,xold,para,npara,De,dstra,9,gi,g)
         normg = norm2(g)
         Iteration = Iteration +1
 !        pause
      enddo
      
      call Jacobfun(xnew,xold,para,npara,De,dstra,Jacob)
      call inverse(Jacob,IJacob,8)
      CTO = 0.0d0
      CTO(1:6,1:6) = IJacob(1:6,1:6)
      CTO = matmul(CTO,De) 
      A21(1,1:6) = IJacob(7,1:6)
      s(1:6,1) = xnew(1:6,1)
      
! Calculate nonlocal damage internal variable and damage variable
        epnw = (1.d0-mw)*xnew(7,1) + ca1
        epnw = max(epnw,1.D-8)
        call dam_fun(damage,dp,epnw,para,npara)
! Calculate the consistent tangent stiffness matrix and nominal stress, and store the effective stress        
        dp = (1.d0-mw)*dp 
        CTO = (1.d0-damage)*CTO - dp*matmul(s,A21)
        statev(1:6) = xnew(1:6,1)
        xnew(1:6,1) = xnew(1:6,1)*(1.d0-damage)     
      end if 

      call dam_fun(damage_i,dp,xnew(7,1),para,npara)
      statev(7) = max(xnew(7,1),1.D-8)
      statev(8) = Iteration
      statev(9) = normg         
      statev(10) = damage_i
      statev(11) = damage
! Record the damage/plastic internal variable of the integration point 
      GPnew(Noel,Npt,4) = statev(7)      
      call DimOut(DDSDDE,stress,CTO,xnew,NTENS)   
      
      if ((Iteration .GE. Newton).or.isnan(normg).or.isnan(norm2(CTO))
     2.or.isnan(norm2(xnew)).or.(Iteration_U .GE. 50.d0)) then
         PNEWDT = 0.25d0
! Give an unreasonable stress and stiffness
         CTO = -1.d60
         stress = 1.d60
         GPnew  = GPold
         Iteration_U = 0.d0
!         call XIT
      end if        

      Return
       end 

 !The source code will be publicly accessible via GitHub upon article publication.      