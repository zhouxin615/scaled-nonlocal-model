!************************************************************************************
!* An unconstrained stress updating algorithm with the line search method for       *
!* elastoplastic soil models                                                        *
!* Copyright (C) Zhou, X., Lu, D.C., Su, C.C., Gao, Z.W., Du, X.L.                  *
!* All Rights Reserved.                                                             *
!* Unauthorized copying of this file, via any medium is strictly prohibited         *
!* Copyright (C) 2021 Zhou,X. (zhouxin615@126.com) Lu,D.C. (dechun@bjut.edu.cn)     *
!* The code can be used for academic research on the premise that the paper         *
!* "Zhou, X., Lu, D.C., Zhao, J.D., Zhang, Y.N., Gao, Z.W., Rabczuk, T., Du, X.L.   *
!* Material characteristic length insensitive nonlocal modelling: A computationally *
!* efficient scaled nonlocal integral method. Computers and Geotechnics.            *
!* 2025, 188: 107587." is cited                                                     *     
!************************************************************************************

!************************************************************************************
!*                              Very important!!!                                   * 
!*(1) root%x_min, root%x_max, root%y_min, root%y_max, root%z_min, root%z_max  need  * 
!* to be set according to the size and position of the analysis object when the     *
!* octree algorithm is used (index_oct = 1)                                         *
!*(2) The octree algorithm presented here is only applicable to 3D examples         *
!************************************************************************************

      module GlobalVar
       implicit none
         integer, parameter :: Nel = 2583, Nip = 8  ! NEED TO BE CHANGED FOR EACH SIMULATIONS
         double precision, target :: GPold(Nel,Nip,4),GPnew(Nel,Nip,4),
     & Volint(Nel,Nip,1)
         double precision incem,Iteration_U
      end module        
      
CC#########################Octree module################################         
      module octree_mod
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
      
      
      if (x_max>node%x_min .and. x_min<node%x_max .and. 
     &    y_max>node%y_min .and. y_min<node%y_max .and.
     &    z_max>node%z_min .and. z_min<node%z_max) then

          if (node%leaf) then
              
          do i = 1,size(node%points)
               
          dis = distance (node%points(i)%tx, node%points(i)%ty, 
     &                    node%points(i)%tz)
          
          if ((dis-radius) .LE. 1.0d-13) then
              
             weight = exp(-dis/cLt)/(2.d0*cLt)
     &                * max(1.D-15, node%points(i)%tv)              
              sumwt  = sumwt + weight  
           if(dis .LE. 1.d-13 )then 
             ca2 = weight
           else
             grad = grad +  node%points(i)%tb * weight   
             
           endif     
           
          end if

          end do
          
          else
          do i = 1,size(node%children)
              
          call search_octree(node%children(i), x_min, x_max, y_min, 
     &                            y_max, z_min, z_max)
          end do
          
          
          end if
          end if
     
      end subroutine search_octree
     
      
     
      subroutine allocate_points_in_range(points, x_mid, y_mid, z_mid, 
     &      points_1, points_2, points_3, points_4,points_5,points_6,
     &      points_7,points_8, num_np)
      type(point),dimension(:), pointer :: points, points_1, points_2, 
     &     points_3, points_4, points_5, points_6, points_7, points_8
      real(kind = 8) :: x_mid, y_mid, z_mid
      integer :: num_points, i, j, num_np(8)

      num_np = 0
      
      do i = 1,size(points)           
          
      if (points(i)%tx <= x_mid) then
          
          if (points(i)%ty <= y_mid) then

              if (points(i)%tz <= z_mid) then
                  num_np(1) = num_np(1) + 1
              else
                  num_np(5) = num_np(5) + 1
              end if       
            
          else
              
              if (points(i)%tz <= z_mid) then
                  num_np(2) = num_np(2) + 1
                  
              else
                  num_np(6) = num_np(6) + 1
              end if
            
          end if
          
      else
          
          if (points(i)%ty <= y_mid) then
              
              if (points(i)%tz <= z_mid) then
                  num_np(3) = num_np(3) + 1
              else
                  num_np(7) = num_np(7) + 1
              end if

          else
              if (points(i)%tz <= z_mid) then
                  num_np(4) = num_np(4) + 1
              else
                  num_np(8) = num_np(8) + 1
              end if
          end if
      end if

      end do 
      
      
      allocate (points_1(num_np(1)))
      allocate (points_2(num_np(2)))
      allocate (points_3(num_np(3)))
      allocate (points_4(num_np(4)))
      allocate (points_5(num_np(5)))
      allocate (points_6(num_np(6)))
      allocate (points_7(num_np(7)))
      allocate (points_8(num_np(8)))
      
      num_np = 0
      do i = 1,size(points)
              
      if (points(i)%tx <= x_mid) then
          
          if (points(i)%ty <= y_mid) then
              
              if (points(i)%tz <= z_mid) then
                  num_np(1) = num_np(1) + 1
                  points_1(num_np(1)) = points(i)
              else
                  
                  num_np(5) = num_np(5) + 1
                  points_5(num_np(5)) = points(i)
              end if
              
          else
              
              if (points(i)%tz <= z_mid) then
                  num_np(2) = num_np(2) + 1
                  points_2(num_np(2)) = points(i)
              else
                  num_np(6) = num_np(6) + 1
                  points_6(num_np(6)) = points(i)
              end if
              
          end if
          
      else
          if (points(i)%ty <= y_mid) then
              
              if (points(i)%tz <= z_mid) then
                  num_np(3) = num_np(3) + 1
                  points_3(num_np(3)) = points(i)
              else
                  num_np(7) = num_np(7) + 1
                  points_7(num_np(7)) = points(i)
              end if
          else
              if (points(i)%tz <= z_mid) then
                  num_np(4) = num_np(4) + 1
                  points_4(num_np(4)) = points(i)
              else
                  num_np(8) = num_np(8) + 1
                  points_8(num_np(8)) = points(i)
              end if
            
          end if
      end if

      end do

      
      end 
          
     
      
      function distance(x, y, z) result(dis)
      real(kind = 8), intent(in) :: x, y, z
      real(kind = 8) :: dis
      
      dis = sqrt ((center%tx - x)**2.0d0 + (center%ty - y)**2.0d0+ 
     &           (center%tz - z)**2.0d0)

      end function distance

      end module octree_mod      
   

      
      
      module ModelParam
       implicit none
       double precision EK,miu,fphi,pi,c0,Epsilon,betaS,Psi,As
       double precision Ftol,Newton,rho,zeta,Maxit
       integer index_s,index_nl,index_oct

      
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
      use octree_mod
      use GlobalVar
      use ModelParam

       implicit none
       
       CHARACTER*80 CMNAME

      Integer NTENS, NDI, NSHR, NSTATV, NPROPS, NOEL, NPT,
     1 LAYER, KSPT, KSTEP, KINC, INITTENSION,DTIME
      DOUBLE PRECISION  stress(ntens),statev(nstatv),
     1  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),jstep(4),
     2  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     3  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), DFGRD1(3,3),
     4  SSE,SPD,SCD,RPL,DRPLDT,TEMP,DTEMP,PNEWDT,CELENT

       double precision BK,GK,dphi,gi,normg,ft,alpha,
     $ damage,damage_i,dp,epn,epnw,p,q,theta,f_ep
       
       double precision Ivol(6,6),Isym(6,6),delta(6,1),De(6,6),
     $ CTO(6,6),g(8,1),xnew(8,1),xold(8,1),Jacob(8,8),IJacob(8,8),
     $ s(6,1),dstra(6,1),st(6,1),dx(8,1),A21(1,6),f_s(6,1)
           
      double precision mw,ca1
      integer K1,K2,K3,Iteration
      
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
      Ftol = 1.D-6
! Read model parameters
      EK = props(1)   
      miu = props(2)
      c0 = props(3)
      fphi = props(4)*pi/180.d0
      Psi = props(5)*pi/180.d0
      As = props(6)
      cLt = props(7) 

      index_s = props(8) ! index_s = 0, 1, 2 represent the linear, exponential, and rational fraction softening law, respectively
      index_nl = props(9) ! index_nl = 0, 1, 2 represent the local model, traditional non-local model, and scaled non-local model, respectively
      index_oct = props(10) ! index_oct = 0, 1  represent traditional- and octree- nonlocal computation methods, respectively. 
! Note that the octree algorithm here is only applicable to 3D conditions
      
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
      call yieldfun(xnew,ft,f_ep,f_s)

      xnew = xold
      call cgfun(xnew,xold,De,dstra,9,gi,g)
      normg = norm2(g)  
C##***********************************Nonlocal computation*************************************##C
      if (index_nl.eq.0) then 
          mw = 0.d0
          ca1 = 0.d0
      else 
          if (index_oct.eq.0) then 
              call nlfunTra(mw,ca1,ntens,noel,npt,time,
     & nstatv,statev)
          elseif (index_oct.eq.1) then 
             call nlfunOct(mw,ca1,statev,nstatv,COORDS,time)
          else 
             write(6,*)"x can only be equal to 0 or 1"
              call XIT
          endif      
      endif       
      
      epnw = (1.d0-mw)*statev(7) + ca1
      epnw = max(epnw,1.D-8)
      call dam_fun(damage,dp,epnw)
!      damage = statev(11)

      if (ft.LE.Ftol) then 
        CTO = De*(1.d0-damage)
        normg = 0.d0
        statev(1:6) = st(1:6,1)
        xnew(1:6,1) = st(1:6,1)*(1.d0-damage)
      else    
      do while ((normg.GT.Ftol).and.(Iteration .LT. Newton))
         call Jacobfun(xnew,xold,De,dstra,Jacob)
         call inverse(Jacob,IJacob,8)
         dx = -matmul(IJacob,g)
         call LSM(xnew,xold,dx,De,dstra,alpha,g)
         xnew = xnew + alpha*dx 
         call cgfun(xnew,xold,De,dstra,9,gi,g)
         normg = norm2(g)
         Iteration = Iteration +1
 !        pause
      enddo
      
      call Jacobfun(xnew,xold,De,dstra,Jacob)
      call inverse(Jacob,IJacob,8)
      CTO = 0.0d0
      CTO(1:6,1:6) = IJacob(1:6,1:6)
      CTO = matmul(CTO,De) 
      A21(1,1:6) = IJacob(7,1:6)
      s(1:6,1) = xnew(1:6,1)
      
! Calculate nonlocal damage internal variable and damage variable
        epnw = (1.d0-mw)*xnew(7,1) + ca1
        epnw = max(epnw,1.D-8)
        call dam_fun(damage,dp,epnw)
! Calculate the consistent tangent stiffness matrix and nominal stress, and store the effective stress        
        dp = (1.d0-mw)*dp 
        CTO = (1.d0-damage)*CTO - dp*matmul(s,A21)
        statev(1:6) = xnew(1:6,1)
        xnew(1:6,1) = xnew(1:6,1)*(1.d0-damage)     
      end if 

      call dam_fun(damage_i,dp,xnew(7,1))
      statev(7) = max(xnew(7,1),1.D-8)
      statev(8) = Iteration
      statev(9) = normg         
      statev(10) = damage_i
      statev(11) = damage
! Record the damage/plastic internal variable of the integration point 
      GPnew(Noel,Npt,4) = statev(7)      
      call DimOut(DDSDDE,stress,CTO,xnew,NTENS)   
      
      if ((Iteration .GE. Newton).or.isnan(normg).or.isnan(norm2(CTO))
     $.or.isnan(norm2(xnew)).or.(Iteration_U .GE. 50.d0)) then
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

      
! The subroutine computing the fourth unit order tensors in matrix form     
	subroutine TenMat1(Ivol,Isym,delta)
	implicit none
      double precision Ivol(6,6),Isym(6,6),delta(6,1)
	integer K1,K2

      Ivol = 0.d0
      Isym = 0.d0
      delta = 0.d0
      do  K1 = 1,3
        delta(K1,1) = 1.d0  
        Isym(K1,K1) = 1.d0
        Isym(K1+3,K1+3) = 1.d0/2.d0          
        do  K2 = 1,3
           Ivol(K1,K2) = 1.d0/3.d0
        end do
      end do

      return
      end

! Convert input stress and strain
	subroutine DimIn(stress,dstran,s,dstra,NTENS)
	implicit none
      double precision stress(NTENS),dstran(NTENS),s(6,1),dstra(6,1)
      integer NTENS

      if (NTENS.EQ.6) then
         s(1:6,1) = stress(1:6)
         dstra(1:6,1) = dstran(1:6)
      elseif (NTENS.EQ.4) then
          s(1:4,1) = stress(1:4)
         dstra(1:4,1) = dstran(1:4)
      elseif (NTENS.EQ.3) then   
         s(1:2,1) = stress(1:2)
         s(4,1) = stress(3)
         dstra(1:2,1) = dstran(1:2)  
         dstra(4,1) = dstran(3) 
      end if         

!       write(*,*)"DDSDDE",DDSDDE
      return
      end 

! Convert output stress and strain
	subroutine DimOut(DDSDDE,stress,CTO,xnew,NTENS)
	implicit none
      double precision DDSDDE(NTENS,NTENS),stress(NTENS),CTO(6,6)
     & ,xnew(8,1)
      integer NTENS,K1,K2

      if (NTENS.EQ. 6) then
         stress(1:6) = xnew(1:6,1)
         DDSDDE = 0.0d0
         DDSDDE = CTO
      elseif (NTENS.EQ. 4) then
         stress(1:4) = xnew(1:4,1)
         DDSDDE = 0.0d0  
          do K1 = 1,4
              do K2 = 1,4
              DDSDDE(K1,K2)= CTO(K1,K2)
              end do
          end do
      elseif (NTENS.EQ.3) then    
         stress(1:2) = xnew(1:2,1)
         stress(3) = xnew(4,1)
         DDSDDE = 0.0d0 

         DDSDDE(1:2,1:2)= CTO(1:2,1:2)
         DDSDDE(1:2,3)= CTO(1:2,4)
         DDSDDE(3,1:2)= CTO(4,1:2)
         DDSDDE(3,3)= CTO(4,4)
        
      end if         

!       write(*,*)"DDSDDE",DDSDDE
      return
      end  


! line search method
	subroutine LSM(xnew,xold,dx,De,dstra,alpha,g)
      use ModelParam
	implicit none
      double precision xnew(8,1),xold(8,1),De(6,6),
     & dstra(6,1),dx(8,1),g(8,1),alpha,F0,F1,normg,Iteration,gi
      
      Iteration = 0
      alpha = 1.0d0

      normg = norm2(g)
      F0 = 0.5d0*normg*normg

         do while (Iteration.LT.Maxit) 

           call cgfun(xnew+alpha*dx,xold,De,dstra,9,gi,g)
            normg = norm2(g)            
            F1 = 0.5d0*normg*normg
            alpha = max(zeta*alpha,F0/(F0+F1)) ! 不管怎样，先执行一次
 
            if (F1.LT.((1.0d0-2.0d0*rho*alpha)*F0)) then
                go to 2008
            else
                alpha = max(zeta*alpha,F0/(F0+F1))
            endif
            Iteration = Iteration + 1
         end do    

2008     continue       

      return
      end 


	subroutine pqsd(s,p,q,sd,theta)
!      include 'ABA_PARAM.INC'
	implicit none
      double precision s(6,1),sd(6,1),delta(6,1),p,q,theta,J2,J3,
     & Dum

      p = (s(1,1)+s(2,1)+s(3,1))/3.d0
      q = dsqrt(((s(1,1)-s(2,1))**2.d0+(s(2,1)-s(3,1))**2.d0
     & +(s(3,1)-s(1,1))**2.d0+6.d0*(s(4,1)**2.d0+s(5,1)**2.d0
     & +s(6,1)**2.d0))/2.d0)

      delta = 0.d0
      delta(1:3,1) = 1.d0
      sd = s-p*delta     

      
      J2 = 1.d0/3.d0*q**2.d0
      J3 = sd(1,1)*sd(2,1)*sd(3,1)+2.d0*sd(4,1)*sd(5,1)*sd(6,1) -
     & sd(1,1)*sd(6,1)*sd(6,1)-sd(2,1)*sd(5,1)*sd(5,1)-
     & sd(3,1)*sd(4,1)*sd(4,1)
      
      if (J2.LE.1.d-12) then
          theta = 0.0d0
      else
          Dum = 3.d0*dsqrt(3.d0)/2.d0*J3/(J2**(3.d0/2.00))
          if (Dum.GE.1.d0) then
              Dum = 1.d0 - 2.d-15
          endif
          if (Dum.LE.-1.d0) then
              Dum = -1.d0 + 2.d-15
          endif
          theta = 1.d0/3.d0*dacos(Dum)
      endif      

      return
      end  

! The subroutine computing Jacobian matrix
	subroutine Jacobfun(xnew,xold,De,dstra,Jacob)
!      include 'ABA_PARAM.INC'
      use ModelParam
	implicit none
      double precision xnew(8,1),xold(8,1),De(6,6),
     $ dstra(6,1),g(8,1),Jacob(8,8),f,f_ep,f_s(6,1),
     & g_2s(6,6),g_qs(6,1),r(6,1),rq,rqfun,I1(6,6),
     & snew(6,1),sold(6,1),epnnew,epnold,dphi
      double precision ga,gb,tem11(6,6),tem21(1,6),tem13(6,1),
     & t_tem21(6,1)
	integer K1,K2,K3
      Jacob = 0.d0
      I1 = 0.d0
      
      do K1 = 1,6
          I1(K1,K1) = 1.d0
      enddo
      snew(1:6,1) = xnew(1:6,1)
      epnnew = xnew(7,1)
      dphi = xnew(8,1)
      sold(1:6,1) = xold(1:6,1)
      epnold = xold(7,1)
           
      call analy2nd(xnew,g_2s,g_qs)
      call rfun(xnew,r)
      rq = rqfun(xnew)
      

      tem11 = I1 + dphi*matmul(De,g_2s)
      t_tem21 = -dphi*g_qs
      tem21(1,1:6) = t_tem21(1:6,1)
      
      tem13 =  matmul(De,r)

      call yieldfun(xnew,f,f_ep,f_s)

      Jacob(1:6,1:6) = tem11(1:6,1:6) ! O K !
      Jacob(1:6,7) = 0.d0             ! O K !
      Jacob(1:6,8) = tem13(1:6,1)     ! O K !
      Jacob(7,1:6) = tem21(1,1:6)     ! O K !
      Jacob(7,7) = 1.d0               ! O K !
      Jacob(7,8) = -rq  ! -rqfun_1    ! O K !
      Jacob(8,1:6) = f_s(1:6,1)       ! O K !
      Jacob(8,7) = f_ep               ! O K !
      Jacob(8,8) = 0.d0               ! O K !

      end 

      ! 
      subroutine inverse(A,IA,N)
      implicit none
      integer :: i,j,N
      real*8   :: A(N,N), IA(N,N),B(N,N),HIA(N),HA(N),Co_number
         ! 
      forall(i=1:N,j=1:N,i==j) IA(i,j) = 1.0D0
      forall(i=1:N,j=1:N,i/=j) IA(i,j)=0.0D0
         ! 
         B=A
         ! 
      call Upper(B,IA,N) !
      call Lower(B,IA,N) !
         ! 
      forall(i=1:N) IA(i,:) = IA(i,:)/B(i,i)
      forall(i=1:N) HIA(i) = sum(abs(IA(i,:)))
      forall(i=1:N) HA(i) = sum(abs(A(i,:)))
      Co_number = dlog((maxval(HA))*(maxval(HIA)))
      if (Co_number.GT.12) then
       ! 
      endif    

      return
      end subroutine

      ! 
      subroutine Upper(M,S,N)
      implicit none
      integer :: N
      real*8  :: M(N,N)
      real*8  :: S(N,N)
      integer :: i,j
      real*8 :: E
      do i=1,N-1
      do j=i+1,N
         E = M(j,i)/M(i,i)
         M(j,i:N)=M(j,i:N)-M(i,i:N)*E
         S(j,:)=S(j,:)-S(i,:)*E
      enddo
      enddo
      return
      end subroutine Upper
      
      !
      subroutine Lower(M,S,N)
      implicit none
      integer :: N
      real*8  :: M(N,N)
      real*8  :: S(N,N)
      integer :: i,j
      real*8 :: E
      do i=N,2,-1
      do j=i-1,1,-1
         E = M(j,i)/M(i,i)
         M(j,1:N)=M(j,1:N)-M(i,1:N)*E
         S(j,:)=S(j,:)-S(i,:)*E
      enddo
      enddo
      return
      end subroutine Lower  


! The subroutine nonlinear equations 
	subroutine cgfun(xnew,xold,De,dstra,K1,gi,g)
      use ModelParam
	implicit none
      double precision snew(6,1),sold(6,1),r(6,1),
     & sd(6,1),qs(6,1),ssfun(6,1),qtemp,rqfun
      double precision p,q,dphi,epnnew,epnold,rq,epnfun,f,Ffun
      double precision xnew(8,1),xold(8,1),De(6,6),
     & dstra(6,1),g(8,1),f_ep,f_s(6,1)
      double precision gi
	integer K1,K2
      gi  = 0.0d0  
      g = 0.0d0       
! read variables
      snew(1:6,1) = xnew(1:6,1)
      epnnew = xnew(7,1)
      dphi = xnew(8,1)
      sold(1:6,1) = xold(1:6,1)
      epnold = xold(7,1)
! calculate plastic flow direction
      call rfun(xnew,r)
! calculate nonliear equations 1-6
      ssfun = 0.d0
      ssfun = snew - sold - matmul(De,(dstra-dphi*r))
! calculate nonliear equation 7
      rq = rqfun(xnew)
      epnfun = epnnew -epnold - rq*dphi 
! calculate nonliear equation 8
      call yieldfun(xnew,f,f_ep,f_s)
      Ffun=f
! Output the values of nonlinear equation systems
      g(1:6,1) = ssfun(1:6,1)
      g(7,1) = epnfun
      g(8,1) = Ffun
      if (K1.LE.8) gi = g(K1,1)

      return
      end 

      function rqfun(x)
      use ModelParam
	implicit none
      double precision x(8,1),s(6,1),sd(6,1)
      double precision rqfun,ep,Rmc,Rnu,Rde,Rmw,p,q,theta,e

      s(1:6,1) = x(1:6,1)
      call pqsd(s,p,q,sd,theta)
         
      e = (3.d0-dsin(fphi))/(3.d0+dsin(fphi))
      Rmc = 1.d0/(2.d0*dcos(fphi)) - dtan(fphi)/6.d0
      Rnu =(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0+(2.d0*e-1.d0)**2.d0)*Rmc

      Rde = 2.d0*(1.d0-e*e)*(dcos(theta)) + 
     $ (2.d0*e-1.d0)*dsqrt(4.d0*(1.d0-e*e)*
     $  (dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)
        Rmw = Rnu/Rde

       rqfun=Rmw*Rmw*q/dsqrt((Epsilon*c0*dtan(Psi))**2.d0+(Rmw*q)**2.d0)
      
      Return
      end function      
      
      subroutine rfun(x,r)
      use ModelParam
	implicit none
      double precision x(8,1),s(6,1),sd(6,1),r(6,1),ps(6,1)
      double precision gfun,ep,Rmc,Rnu,Rde,Rmw,p,q,theta,e,Rnud,Rded
      double precision c1s3,qtemp,J2,J3,Rmws(6,1)
      double precision J2s(6,1),J3s(6,1),qs(6,1),delta(6,1),c3s(6,1)

      s(1:6,1) = x(1:6,1)
      delta = 0.d0
      delta(1:3,1) = 1.d0      
      call J2J3(s,p,q,sd,theta,J2,J3)
         
      e = (3.d0-dsin(fphi))/(3.d0+dsin(fphi))
      
      Rmc = 1.d0/(2.d0*dcos(fphi)) - dtan(fphi)/6.d0
      Rnu =(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0+(2.d0*e-1.d0)**2.d0)*Rmc
      Rnud = 8.d0*Rmc*(1.d0-e*e)*dcos(theta)

      Rde = 2.d0*(1.d0-e*e)*(dcos(theta)) + 
     $ (2.d0*e-1.d0)*dsqrt(4.d0*(1.d0-e*e)*
     $  (dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)
      
      Rded =2.d0*(1.d0-e*e) + (2.d0*e-1.d0)*4.d0*(1.d0-e*e)*dcos(theta)/
     $ dsqrt(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)
      
      c1s3 = 1.d0/(4.d0*(dcos(theta))**2.d0-1.d0)/3.d0        
      
      J2s = sd
      J2s(4:6,1) = 2.d0*sd(4:6,1)
      
      J3s(1,1) = sd(2,1)*sd(3,1) - sd(6,1)*sd(6,1)
      J3s(2,1) = sd(1,1)*sd(3,1) - sd(5,1)*sd(5,1)
      J3s(3,1) = sd(1,1)*sd(2,1) - sd(4,1)*sd(4,1)
      J3s(4,1) = 2.d0*(sd(6,1)*sd(5,1) - sd(3,1)*sd(4,1))
      J3s(5,1) = 2.d0*(sd(4,1)*sd(6,1) - sd(2,1)*sd(5,1))
      J3s(6,1) = 2.d0*(sd(5,1)*sd(4,1) - sd(1,1)*sd(6,1))
      J3s = J3s + J2*delta/3.d0
      
      
      ps(1:3,1) = 1.d0/3.d0
      ps(4:6,1) = 0.d0      
      qs = 3.d0*sd/(2.d0*q)
      qs(4:6,1) = 2.d0*qs(4:6,1)
      if (J2.LE.1.d-12)  then
         qtemp = 1.d-12
         qs = 3.d0*sd/(2.d0*qtemp)
      endif      

      if (J2.GE.1.d-12)  then     ! cos(3theta)
        c3s =3.d0*dsqrt(3.d0)/4.d0*(2.d0*J2*J3s-3.d0*J3*J2s)/(J2**2.5d0)
      else
        J2 = 1.d-12  
        c3s = 0.d0
      endif  

      
      Rmws = (Rnud*Rde-Rnu*Rded)*c1s3*c3s/Rde**2.d0
      Rmw = Rnu/Rde
       r=(Rmw*qs+q*Rmws)*Rmw*q/dsqrt((Epsilon*c0*dtan(Psi))**2.d0+
     $ (Rmw*q)**2.d0)+ps*dtan(Psi)

      Return
      end subroutine 
      


! The subroutine yield function  
	subroutine J2J3(s,p,q,sd,theta,J2,J3)
!      include 'ABA_PARAM.INC'
	implicit none
      double precision s(6,1),sd(6,1),delta(6,1),p,q,theta,J2,J3,
     & Dum
      p = (s(1,1)+s(2,1)+s(3,1))/3.d0
      q = dsqrt(((s(1,1)-s(2,1))**2.d0+(s(2,1)-s(3,1))**2.d0
     & +(s(3,1)-s(1,1))**2.d0+6.d0*(s(4,1)**2.d0+s(5,1)**2.d0
     & +s(6,1)**2.d0))/2.d0)

      delta = 0.d0
      delta(1:3,1) = 1.d0
      sd = s-p*delta     

      
      J2 = 1.d0/3.d0*q**2.d0
      J3 = sd(1,1)*sd(2,1)*sd(3,1)+2.d0*sd(4,1)*sd(5,1)*sd(6,1) -
     & sd(1,1)*sd(6,1)*sd(6,1)-sd(2,1)*sd(5,1)*sd(5,1)-
     & sd(3,1)*sd(4,1)*sd(4,1)
      
      if (J2.LE.1.d-12) then
          theta = 0.0d0
      else
          Dum = 3.d0*dsqrt(3.d0)/2.d0*J3/(J2**(3.d0/2.d0))
          if (Dum.GE.(1.d0)) then
              Dum = 1.d0 - 1.d-15
          endif
          if (Dum.LE.-(1.d0)) then
              Dum = -1.d0 + 1.d-15
          endif
          theta = 1.d0/3.d0*dacos(Dum)
      endif      

      return   
         
      end  
 


! The subroutine yield function  
	subroutine yieldfun(x,fun,f_ep,f_s)
!      include 'ABA_PARAM.INC'
      use ModelParam
	implicit none
      double precision, intent(in)::x(8,1)
      double precision s(6,1),sd(6,1),p,q,theta,epn,c,J2,J3
      double precision MS,KS,aw,gamw,gam,Rmc,J3s(6,1)
      double precision p_s(6,1),q_s(6,1),c3s(6,1),J2s(6,1)
      double precision R1,R2,Rmc_R1,R1_R2,R2_c3s,Rmc_s(6,1)
      double precision Hs,Hsp,damage,dp,Hw,Hwp
      double precision, intent(out)::fun,f_ep,f_s(6,1)
	integer K1,K2
      
!      write(*,*)"x",x
      s(1:6,1) = x(1:6,1)
      epn = x(7,1)
      
      call Hs_fun(Hs,Hsp,epn)
      call dam_fun(damage,dp,epn)
      Hw = Hs/(1.d0-damage)                 ! 根据损伤函数与名义函数确定有效函数
      Hwp = ((1.d0-damage)*Hsp + Hs*dp)/(1.d0-damage)**2.d0      
      
      c = c0*Hw
      call Diff_of_Args(x,p,q,sd,theta,p_s,q_s,c3s,J2s,
     & J3s,J2,J3)
      MS = 6.d0*dsin(fphi)/dsqrt(3.d0)/(3.d0-dsin(fphi))
      KS = 6.d0*c*dcos(fphi)/dsqrt(3.d0)/(3.d0-dsin(fphi))
      
      gamw = 6.d0/pi*datan(dsin(fphi)/dsqrt(3.d0))
      gam = 1.d0-gamw
      aw = 1.d0/(dcos((gamw+1.d0)*pi/6.d0))
      theta = theta - pi/6.d0
      Rmc = aw*dcos(1.d0/3.d0*dacos(-betaS*dsin(3.d0*theta))
     $ -pi/6.d0*gam)


      fun = MS*p - KS + Rmc*q*dsqrt(3.d0)/3.d0
      
      f_ep = -6.d0*dcos(fphi) / ((3.d0 - dsin(fphi))*dsqrt(3.d0))*
     & c0*Hwp
      

      R2 = -betaS*dsin(3.d0*theta)
      R1 = 1.d0/3.d0*dacos(R2) - pi/6.d0*gam
      Rmc_R1 = -aw*dsin(R1)
      R1_R2 = -1.d0/(3.d0*dsqrt(1.d0-R2**2.d0))

      R2_c3s = betaS!*dcos(3.d0*theta)/dsin(3.d0*theta)
      if (J2.GE.1.d-12)  then     ! cos(3theta)
        c3s =3.d0*dsqrt(3.d0)/4.d0*(2.d0*J2*J3s-3.d0*J3*J2s)/(J2**2.5d0)
      else
        J2 = 1.d-12  
        c3s = 0.d0
      endif  

      Rmc_s = Rmc_R1*R1_R2*R2_c3s*c3s

      f_s =  MS*p_s + dsqrt(3.d0)/3.d0*(Rmc*q_s + q*Rmc_s)

         
      return
      end  


      subroutine Diff_of_Args(x,p,q,sd,theta,ps,qs,c3s,J2s,
     & J3s,J2,J3)
      use ModelParam
      implicit none
      double precision x(8,1),s(6,1),sd(6,1),ps(6,1)
      double precision p,q,theta,c1s3,qtemp,J2,J3,Rmws(6,1),Dum
      double precision J2s(6,1),J3s(6,1),qs(6,1),delta(6,1),c3s(6,1)
      
      s(1:6,1) = x(1:6,1)
      delta = 0.d0
      delta(1:3,1) = 1.d0      
            

      p = (s(1,1)+s(2,1)+s(3,1))/3.d0
      q = dsqrt(((s(1,1)-s(2,1))**2.d0+(s(2,1)-s(3,1))**2.d0
     & +(s(3,1)-s(1,1))**2.d0+6.d0*(s(4,1)**2.d0+s(5,1)**2.d0
     & +s(6,1)**2.d0))/2.d0)

      delta = 0.d0
      delta(1:3,1) = 1.d0
      sd = s-p*delta     

      J2 = 1.d0/3.d0*q**2.d0
      J3 = sd(1,1)*sd(2,1)*sd(3,1)+2.d0*sd(4,1)*sd(5,1)*sd(6,1) -
     & sd(1,1)*sd(6,1)*sd(6,1)-sd(2,1)*sd(5,1)*sd(5,1)-
     & sd(3,1)*sd(4,1)*sd(4,1)
      
      if (J2.LE.1.d-12) then
          theta = 0.0d0
      else
          Dum = 3.d0*dsqrt(3.d0)/2.d0*J3/(J2**(3.d0/2.d0))
          if (Dum.GE.(1.d0)) then
              Dum = 1.d0 - 1.d-15
          endif
          if (Dum.LE.-(1.d0)) then
              Dum = -1.d0 + 1.d-15
          endif
          theta = 1.d0/3.d0*dacos(Dum)
      endif      
      !write(*,*)"theta",theta
      c1s3 = 1.d0/(4.d0*(dcos(theta))**2.d0-1.d0)/3.d0 
      
      J2s = sd
      J2s(4:6,1) = 2.d0*sd(4:6,1)
      
      J3s(1,1) = sd(2,1)*sd(3,1) - sd(6,1)*sd(6,1)
      J3s(2,1) = sd(1,1)*sd(3,1) - sd(5,1)*sd(5,1)
      J3s(3,1) = sd(1,1)*sd(2,1) - sd(4,1)*sd(4,1)
      J3s(4,1) = 2.d0*(sd(6,1)*sd(5,1) - sd(3,1)*sd(4,1))
      J3s(5,1) = 2.d0*(sd(4,1)*sd(6,1) - sd(2,1)*sd(5,1))
      J3s(6,1) = 2.d0*(sd(5,1)*sd(4,1) - sd(1,1)*sd(6,1))
      J3s = J3s + J2*delta/3.d0
      
      ps(1:3,1) = 1.d0/3.d0
      ps(4:6,1) = 0.d0      
      qs = 3.d0*sd/(2.d0*q)
      qs(4:6,1) = 2.d0*qs(4:6,1)
      if (J2.LE.1.d-12)  then
         qtemp = 1.d-12
         qs = 3.d0*sd/(2.d0*qtemp)
      endif      

      if (J2.GE.1.d-12)  then
        ! pian cos(3theta)/ pian s
        c3s = 3.d0*dsqrt(3.d0)/4.d0*(2.d0*J2*J3s-3.d0*J3*J2s)/
     &    (J2**2.5d0)
      else
        J2 = 1.d-12  
        c3s = 0.d0
!        write(*,*)"c3s=0.d0"
      endif  
      
      end

      subroutine analy2nd(x,g_2s,g_qs)
      use ModelParam
      implicit none
      double precision x(8,1),s(6,1),sd(6,1),g_2s(6,6)
      double precision g_qs(6,1),sba
      double precision p,q,theta,c1s3,qtemp,J2,J3,ps(6,1),qs(6,1)
      double precision J2s(6,1),J3s(6,1),delta(6,1),c3s(6,1),tm_qs(1,6)
      double precision sba_2s(6,6),J3_2s(6,6),J2_2s(6,6),q_2s(6,6)
      double precision tm_sba_2s(6,6),tm_J3_2s(6,6),Rmw_cost_s(6,1),
     & Rmw_2cost,cost_cs3,Rmc,Rnu,Rde,Rmw,e,Rnud,Rded,Upie,Vpie,bigW,
     & Upp,Vpp,Rmw_cost,Rmw_2s(6,6),cs_c3,cs2_c3_s(6,1),c3_2s(6,6),
     & cs2_c3_cs,c32_s2(6,6),A1_s(6,1),Rmw_s(6,1),A1,A2_s(6,1),A2,
     & tem1,tm_A1s(1,6),tm_A2s(1,6),tm_J2s(1,6),Jup(6,1),tm_Jup(1,6),
     & tm_R2(1,6),tm_J3s(1,6),tem_Rc(6,1),tem_cs2(1,6),
     & Rmw_2s_a,Rmw_2s_b,q_2s_a,q_2s_b,Rmw_2s_t,q_2s_t,
     & cos3a,cos3b,cos3,J3_2sa,J3_2sb
      integer K1,K2
      
      s(1:6,1) = x(1:6,1)
      delta = 0.d0
      delta(1:3,1) = 1.d0      
      
      call Diff_of_Args(x,p,q,sd,theta,ps,qs,c3s,J2s,
     & J3s,J2,J3)
      
!*********************** begin J2_2s & J3_2s ***********************
      if (q.LE.1.d-10) then
          q = 1.d-10
      endif
      q_2s(1,1) = -9.d0/(4.d0*q**3.d0)*sd(1,1)*sd(1,1) + 1.d0/q
      
      q_2s(2,1) = -9.d0/(4.d0*q**3.d0)*sd(2,1)*sd(1,1) - 1.d0/2.d0/q
      q_2s(2,2) = -9.d0/(4.d0*q**3.d0)*sd(2,1)*sd(2,1) + 1.d0/q
      q_2s(1,2) = q_2s(2,1)
      
      q_2s(3,1) = -9.d0/(4.d0*q**3.d0)*sd(3,1)*sd(1,1) - 1.d0/2.d0/q
      q_2s(3,2) = -9.d0/(4.d0*q**3.d0)*sd(3,1)*sd(2,1) - 1.d0/2.d0/q
      q_2s(3,3) = -9.d0/(4.d0*q**3.d0)*sd(3,1)*sd(3,1) + 1.d0/q
      q_2s(1,3) = q_2s(3,1)
      q_2s(2,3) = q_2s(3,2)
      
      q_2s(4,1) = -9.d0/(4.d0*q**3.d0)*sd(4,1)*sd(1,1)*2.d0
      q_2s(4,2) = -9.d0/(4.d0*q**3.d0)*sd(4,1)*sd(2,1)*2.d0
      q_2s(4,3) = -9.d0/(4.d0*q**3.d0)*sd(4,1)*sd(3,1)*2.d0
      q_2s(4,4) = -9.d0/(4.d0*q**3.d0)*sd(4,1)*sd(4,1)*4.d0 + 3.d0/q
      q_2s(1,4) = q_2s(4,1)
      q_2s(2,4) = q_2s(4,2)
      q_2s(3,4) = q_2s(4,3)
      
      q_2s(5,1) = -9.d0/(4.d0*q**3.d0)*sd(5,1)*sd(1,1)*2.d0
      q_2s(5,2) = -9.d0/(4.d0*q**3.d0)*sd(5,1)*sd(2,1)*2.d0
      q_2s(5,3) = -9.d0/(4.d0*q**3.d0)*sd(5,1)*sd(3,1)*2.d0
      q_2s(5,4) = -9.d0/(4.d0*q**3.d0)*sd(5,1)*sd(4,1)*4.d0
      q_2s(5,5) = -9.d0/(4.d0*q**3.d0)*sd(5,1)*sd(5,1)*2.d0 + 3.d0/q
      q_2s(1,5) = q_2s(5,1)
      q_2s(2,5) = q_2s(5,2)
      q_2s(3,5) = q_2s(5,3)
      q_2s(4,5) = q_2s(5,4)
      
      q_2s(6,1) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(1,1)*2.d0
      q_2s(6,2) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(2,1)*2.d0
      q_2s(6,3) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(3,1)*2.d0
      q_2s(6,4) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(4,1)*4.d0
      q_2s(6,5) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(5,1)*4.d0
      q_2s(6,6) = -9.d0/(4.d0*q**3.d0)*sd(6,1)*sd(6,1)*4.d0 + 3.d0/q
      q_2s(1,6) = q_2s(6,1)
      q_2s(2,6) = q_2s(6,2)
      q_2s(3,6) = q_2s(6,3)
      q_2s(4,6) = q_2s(6,4)
      q_2s(5,6) = q_2s(6,5)

      
      tm_qs(1,1:6) = qs(1:6,1)
      J2_2s = 2.d0/3.d0*matmul(qs,tm_qs) + 2.d0/3.d0*q*q_2s
      
      tm_J3_2s(1,1) = (4.d0*s(1,1)-2.d0*s(2,1)-2.d0*s(3,1))/9.d0
      tm_J3_2s(2,1) = (4.d0*s(3,1)-2.d0*s(2,1)-2.d0*s(1,1))/9.d0
      tm_J3_2s(2,2) = (4.d0*s(2,1)-2.d0*s(3,1)-2.d0*s(1,1))/9.d0
      tm_J3_2s(1,2) = tm_J3_2s(2,1)
      
      tm_J3_2s(3,1) = (4.d0*s(2,1)-2.d0*s(3,1)-2.d0*s(1,1))/9.d0
      tm_J3_2s(3,2) = (4.d0*s(1,1)-2.d0*s(2,1)-2.d0*s(3,1))/9.d0
      tm_J3_2s(3,3) = (-2.d0*s(1,1)-2.d0*s(2,1)+4.d0*s(3,1))/9.d0
      tm_J3_2s(1,3) = tm_J3_2s(3,1)
      tm_J3_2s(2,3) = tm_J3_2s(3,2)
      
      tm_J3_2s(4,1) = 2.d0*s(4,1)/3.d0
      tm_J3_2s(4,2) = 2.d0*s(4,1)/3.d0
      tm_J3_2s(4,3) = -4.d0*s(4,1)/3.d0
      tm_J3_2s(4,4) = (2.d0*s(1,1)+2.d0*s(2,1)-4.d0*s(3,1))/3.d0
      tm_J3_2s(1,4) = tm_J3_2s(4,1)
      tm_J3_2s(2,4) = tm_J3_2s(4,2)
      tm_J3_2s(3,4) = tm_J3_2s(4,3)
      
      tm_J3_2s(5,1) = 2.d0*s(5,1)/3.d0
      tm_J3_2s(5,2) = -4.d0*s(5,1)/3.d0
      tm_J3_2s(5,3) = 2.d0*s(5,1)/3.d0
      tm_J3_2s(5,4) = 2.d0*s(6,1)
      tm_J3_2s(5,5) = (2.d0*s(1,1)+2.d0*s(3,1)-4.d0*s(2,1))/3.d0
      tm_J3_2s(1,5) = tm_J3_2s(5,1)
      tm_J3_2s(2,5) = tm_J3_2s(5,2)
      tm_J3_2s(3,5) = tm_J3_2s(5,3)
      tm_J3_2s(4,5) = tm_J3_2s(5,4)
      
      tm_J3_2s(6,1) = -4.d0*s(6,1)/3.d0  !-4.d0*sd(6,1)
      tm_J3_2s(6,2) = 2.d0*s(6,1)/3.d0
      tm_J3_2s(6,3) = 2.d0*s(6,1)/3.d0
      tm_J3_2s(6,4) = 2.d0*s(5,1)
      tm_J3_2s(6,5) = 2.d0*s(4,1)
      tm_J3_2s(6,6) = (2.d0*s(2,1)+2.d0*s(3,1)-4.d0*s(1,1))/3.d0
      tm_J3_2s(1,6) = tm_J3_2s(6,1)
      tm_J3_2s(2,6) = tm_J3_2s(6,2)
      tm_J3_2s(3,6) = tm_J3_2s(6,3)
      tm_J3_2s(4,6) = tm_J3_2s(6,4)
      tm_J3_2s(5,6) = tm_J3_2s(6,5)
      
!      J3_2s = 1.d0/3.d0*tm_J3_2s
      J3_2s = tm_J3_2s
!***************** end   J2_2s & J3_2s ************************
      
!*********************  begin Rmw_2s  *********************************
!********cs_c3 & c3s & Rmw_cost_s & Rmw_cost & cs2_c3_s & c32_s2*******

      e = (3.d0-dsin(fphi))/(3.d0+dsin(fphi))
      Rmc = 1.d0/(2.d0*dcos(fphi)) - dtan(fphi)/6.d0
      Rnu =(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0+(2.d0*e-1.d0)**2.d0)*Rmc
      Rnud = 8.d0*Rmc*(1.d0-e*e)*dcos(theta)
      Rde = 2.d0*(1.d0-e*e)*(dcos(theta)) + 
     $ (2.d0*e-1.d0)*dsqrt(4.d0*(1.d0-e*e)*
     $  (dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)  ! V
      
      Rded =2.d0*(1.d0-e*e) + (2.d0*e-1.d0)*4.d0*(1.d0-e*e)*dcos(theta)/
     $ dsqrt(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)
      
      Rmw = Rnu/Rde ! U / V
      Upie = 8.d0*Rmc*(1.d0-e*e)*dcos(theta)  !(√)
      Vpie = 2.d0*(1.d0-e*e) + (2.d0*e-1.d0)*4.d0*(1.d0-e*e)*dcos(theta)
     & /(dsqrt(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e))!(√)
      
      bigW =dsqrt(4.d0*(1.d0-e*e)*(dcos(theta))**2.d0 + 5.d0*e*e-4.d0*e)
      Upp = 8.d0*Rmc*(1.d0-e*e)
      Vpp = (2.d0*e-1.d0)*(4.d0*(1.d0-e*e)*(bigW**2.d0) -
     & (4.d0*(1.d0-e*e)*dcos(theta))**2.d0) / (bigW**3.d0)!(√)
      

!    Rmw_2s --> Rmw_cost_s*cs_c3*c3s
!             + Rmw_cost*c3s*cs2_c3_s
!             + Rmw_cost*cs_c3*c32_s2
!    Rmw_cost_s ----> Rmw/costheta/s ---> Rmw_2cost & cs_c3 & c3s
      Rmw_cost = (Upie*Rde - Vpie*Rnu)/(Rde**2.d0) !(√)
      
      Rmw_2cost = ((Upp*Rde-Rnu*Vpp)*Rde**2.d0 - 2.d0*Rde*Vpie*
     & (Upie*Rde-Rnu*Vpie)) / (Rde**4.d0)               !(√)
      cs_c3 = 1.d0/3.d0/(4.d0*(dcos(theta)**2.d0)-1.d0) !(√)
              !1.d0/3.d0/(4.d0*(dcos(theta))**2.d0-1.d0)
      Rmw_cost_s = Rmw_2cost*cs_c3*c3s                  !(√)
      
!    cs2_c3_s ---> cs2_c3_cs*cs_c3*c3s
      cs2_c3_cs = -8.d0*dcos(theta)/ 3.d0 / ( (4.d0*dcos(theta)**2.d0
     &  - 1.d0) **2.d0)
      cs2_c3_s = cs2_c3_cs*cs_c3*c3s
      
      tm_J2s(1,1:6) = J2s(1:6,1)
      tm_J3s(1,1:6) = J3s(1:6,1)
      Jup = 2.d0*J2*J3s - 3.d0*J3*J2s
      tm_Jup(1,1:6) = Jup(1:6,1)
      
                              !(√)
        c32_s2 = 3.d0*dsqrt(3.d0)/4.d0*
     & ( (2.d0*matmul(J2s,tm_J3s) + 2.d0*J2*J3_2s
     & - 3.d0*matmul(J3s,tm_J2s) - 3.d0*J3*J2_2s)*
     & (J2**2.5d0) - 
     & 2.5d0*(J2**1.5d0)*matmul(J2s,tm_Jup) ) / 
     & (J2**5.d0)                                  !(√)
      

      
      tem_Rc = Rmw_cost_s*cs_c3
      tem_cs2(1,1:6) = cs2_c3_s(1:6,1)
      Rmw_2s = cs_c3*matmul(c3s,transpose(Rmw_cost_s)) + 
     & Rmw_cost*matmul(c3s,transpose(cs2_c3_s)) +Rmw_cost*cs_c3*c32_s2
      
 
      

!***********************   end Rmw_2s  ******************************
      
!*************************  begin g_2s  *****************************
!**** A1_s & Rmw_s & A1 & Rmw_2s(√) & A2_s & qs(√) & A2 & q_2s(√) ***
      Rmw_s = Rmw_cost*cs_c3*c3s  !(√)

      
      A1 = Rmw*q**2.d0 / dsqrt((Epsilon*c0*dtan(Psi))**2.d0+
     & (Rmw*q)**2.d0)
      A2 = Rmw**2.d0*q / dsqrt((Epsilon*c0*dtan(Psi))**2.d0+
     & (Rmw*q)**2.d0)
      
      tem1 = (Epsilon*c0*dtan(Psi))**2.d0+(Rmw*q)**2.d0
      
      A1_s = ( (Rmw_s*q**2.d0 + 2.d0*Rmw*q*qs)*dsqrt(tem1) -
     & Rmw*q**2.d0*Rmw*q*(q*Rmw_s+Rmw*qs)/dsqrt(tem1) ) / tem1
      
      A2_s = ( (2.d0*Rmw_s*Rmw*q + Rmw**2.d0*qs)*dsqrt(tem1) -
     & Rmw**2.d0*q*Rmw*q*(q*Rmw_s+Rmw*qs)/dsqrt(tem1) ) / tem1

      
      tm_A1s(1,1:6) = A1_s(1:6,1)
      tm_A2s(1,1:6) = A2_s(1:6,1)
      
      g_2s = matmul(Rmw_s,tm_A1s) + A1*Rmw_2s + matmul(qs,tm_A2s) +
     &  A2*q_2s

!************************************  end g_2s  ******************
      g_qs = A2_s
      
      return
      end

      
      subroutine dam_fun(damage,dp,epn)
      use ModelParam
	implicit none
      double precision damage,dp,epn,Hs,Hsp,Hw,Hwp,Hsd,Ad,
     & Hspd

      
      call Hs_fun(Hs,Hsp,epn)
      call Hw_fun(Hw,Hwp,epn)
     

      damage = 1.d0 - Hs/Hw             
      dp = -(Hsp*Hw-Hs*Hwp)/Hw**2.d0        
      
      if (damage.GT.0.999999d0) then
          damage = 0.999999d0
          dp = 0.d0
      elseif (damage.LT.0.0d0) then
          damage = 0.d0
      endif      

      Return
      end subroutine      
    

      
      subroutine Hs_fun(Hs,Hsp,epn)
      use octree_mod
      use ModelParam
	implicit none
      double precision Hs,Hsp,epn,epnd,AL

      
      if (index_nl.EQ.2) then
          AL = cLt/0.12d0
      else
          AL = 1.d0
      endif
      epnd = epn*AL    
      
      if (index_s.EQ.0) then
      Hs=((1.d0-As*epnd)+dsqrt((1.d0-As*epnd)**2.d0+4.d0*1.d-3))/2.d0 
      Hsp=0.5d0*((As*epnd-1.d0)/
     $ dsqrt((1.d0-As*epnd)**2.d0+4.d0*1.d-3)-As)   
        
      elseif (index_s.EQ.1) then
      Hs = dexp(-As*epnd)             
      Hsp =  -As*dexp(-As*epnd)  
         
      elseif (index_s.EQ.2) then  
      Hs = (epnd+1.d0)/(As*epnd*epnd+epnd+1.d0)      
      Hsp = ((As*epnd*epnd+epnd+1.d0)-(2.d0*As*epnd+1.d0)*(epnd+1.d0))/
     $ (As*epnd*epnd+epnd+1.d0)**2.d0        
          
      endif   
      
       Hsp = Hsp*AL 
      
      Return
      end subroutine    
      
      
      subroutine Hw_fun(Hw,Hwp,epn)
      use ModelParam
	implicit none
      double precision Hw,Hwp,epn

 
      Hw = 1.d0         
      Hwp = 0.d0
      
      Return
      end subroutine        
      

      subroutine nlfunTra(mw,ca1,ntens,noel,npt,time,
     & nstatv,statev)
      use octree_mod
      use ModelParam
      use GlobalVar
      
	implicit none
      double precision ca1,cutoff,mw,dist,weight,time(2),statev(nstatv)
      integer K1,K2,noel,npt,ntens,nstatv
      
      grad = 0.0D0    
      sumwt = 0.0D0   !sumwt = sum(w_j)
      ca1 = 0.d0      !ca1 = mw0*sum(w_j*kappa_j)/sum(w_j)
      ca2 = 0.0d0     !ca2 = mw0*kappa_i//sum(w_j)
      
      
      cutoff = 2.d0   
      mw0  = 1.1d0 


      do K1 = 1, Nel
          do K2 = 1, Nip

             dist = sqrt((GPold(noel,npt,1)-GPold(K1,K2,1))**2.0d0+
     &       (GPold(noel,npt,2)-GPold(K1,K2,2))**2.0d0         
     &       +(GPold(noel,npt,3)-GPold(K1,K2,3))**2.0d0)
 
            if(dist .LE. (cutoff*cLt) )then   ! Cut off at 2-3*cLt to save computation time, this causes little error. See the reference papers. If we cut off at 3, the results will be very strange.
! Note that the coefficient 2.d0*cLt does not affect the calculation result, as it appears both in the numerator and denominator of the non-local formula

              weight=exp(-dist/cLt)/(2.d0*cLt)  
     $            * max(1.D-15, Volint(K1,K2,1))   ! VOLUME !              
              
              sumwt = sumwt + weight    ! Using uniform mesh in the averaging area saves a lot of compuation time !              
              if(dist .LE. 1.d-13 )then 
                  ca2 = weight
              else
                  grad = grad + GPold(K1,K2,4) * weight     ! The volume of each element is not considered because we use uniform element size in the area which needs averaging.
              endif     
              
            endif   

          end do
          end do 

      
      ca2 = mw0*ca2/sumwt
      if (time(2).LT.1.D-13) then
          ca1 = 0.d0
          ca2 = 0.0d0
      else    
          ca1 = mw0*grad/sumwt
      end if  
      mw = mw0 - ca2  

 
      
      Return
      end subroutine        
  
      subroutine nlfunOct(mw,ca1,statev,nstatv,COORDS,time)
      use GlobalVar
      use octree_mod
      use ModelParam
      
	implicit none
      double precision time(2),statev(nstatv),ca1,mw
      integer nstatv,di,K1,K2,i
      double precision, target :: coords(3)
       
      
      grad = 0.0D0   
      sumwt = 0.0D0   !sumwt = sum(w_j)
      ca1 = 0.d0      !ca1 = mw0*sum(w_j*kappa_j)/sum(w_j)
      ca2 = 0.0d0     !ca2 = mw0*kappa_i//sum(w_j)      

      
      if ((statev(7).GE.1.d-4)) then
      
      center%tx => COORDS(1)
      center%ty => COORDS(2)
      center%tz => COORDS(3)
      
      
      radius = 2.d0*cLt
      mw0 = 1.1d0  
            
            
      if (time(2).LT.1.D-7) then
          
          ca1 = 0.d0
          ca2 = 0.0d0
          
      elseif (associated(root)) then
          
          call search_octree(root, center%tx-radius, center%tx+radius, 
     &            center%ty-radius, center%ty+radius,center%tz-radius,
     &            center%tz+radius)
      
      else
          
          di = NEL * NIP                                    !!!!!!

      if (associated(root)) then 
          deallocate(root)
      end if

      allocate(root, root%points(di))
          
      do K1 = 1, NEL                                         !!!!!!
         do K2 = 1, NIP
          i = (K1 -1)*8+K2
          
          root%points(i)%tx => GPold(K1,K2,1)
          root%points(i)%ty => GPold(K1,K2,2)
          root%points(i)%tz => GPold(K1,K2,3)
          root%points(i)%tv => VOLINT(K1,K2,1)
          root%points(i)%tb => GPold(K1,K2,4)

         end do
      end do
! we need to set the vertex coordinates of a box containing the analysis object      
      root%x_min = 0.0d0                             
      root%x_max = 12.0d0
      root%y_min = 0.0d0
      root%y_max = 3.50d0
      root%z_min = 0.0d0
      root%z_max = 2.0d0
      
      call build_octree(root)
      
      call search_octree(root, center%tx-radius, center%tx+radius, 
     &            center%ty-radius, center%ty+radius,center%tz-radius,
     &            center%tz+radius)

      end if

      
      ca2 = mw0*ca2/sumwt    

      
      if (time(2).LT.1.D-7) then
          ca1 = 0.d0
          ca2 = 0.0d0
      else    
          ca1 = mw0*grad/sumwt
      end if  
      mw = mw0 - ca2  

      
      else
         mw =0.d0 
         mw0 = 0.d0
         ca1 = 0.d0
      end if 
      
      Return
      end subroutine               
      