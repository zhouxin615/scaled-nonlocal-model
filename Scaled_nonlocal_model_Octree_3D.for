      module GlobalVariables
       implicit none
         integer, parameter :: Nel = 17459, Nip = 8  ! NEED TO BE CHANGED FOR EACH SIMULATIONS
         double precision, target :: Encd(Nel,Nip,4),Encdnew(Nel,Nip,4),
     & Volint(Nel,Nip,1)
         double precision incem,Iteration_U
      end module        
      
CC#########################�˲���ģ��################################            
CC#########################�˲���ģ��################################       
      module octree_mod
      implicit none
!���������
      type point
      real(kind = 8), pointer :: tx, ty, tz, tv, tb
      end type point
      

!����˲����ڵ�����
      type octree_node
      type(point), dimension(:), pointer :: points
      type(octree_node), dimension(:), pointer :: children
      real(kind = 8) :: x_min, x_max, y_min, y_max, z_min, z_max
      logical leaf
      end type octree_node
      
      
      type(point) :: center
      real(kind = 8) ::  grad, sumwt, radius, cLt, mw0, ca2
      type(octree_node), pointer :: root
  
!�ݹ鹹���˲���
      Contains
      

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
      
!      if (sum(num_np) . NE. size(points)) 
!     &   write(*,*) "There is unallocated data present"
      
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
CC#########################�˲���ģ��################################            
CC#########################�˲���ģ��################################     
      
      
      
      
      
      module ModelParam
       implicit none
! ģ�Ͳ���      
       double precision EK,miu,fphi,pi,c0,Epsilon,betaS,Psi,As
! �㷨����
       double precision Ftol,Newton,rho,zeta,Maxit,cLt
       integer index_s,index_nl
! index_s = 0,1,2 �ֱ��ʾ����,ָ��,�Լ������ʽ������
! index_nl = 0,1,2 �ֱ��ʾ�ֲ�ģ��,�ɷǾֲ�ģ��,�·Ǿֲ�ģ��     
       
       contains 
          subroutine Initialize(para,npara)
        !===================================
            double precision para(npara)
            integer npara
            !********************************************

! �㷨������ģ�Ͳ���
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
      Use GlobalVariables
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
      use GlobalVariables
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
       
! �ͻ��ַǾֲ������йصĲ���       
      double precision mw,ca1,index_s,index_nl,Amw,Aca1
      integer K1,K2,K3,Iteration,Newton,npara
      
      Encdnew(noel,npt,1:3) = Coords(1:3)
      if (time(2).LT.1.D-10) then
          statev = 0.d0
          statev(1:6) = stress(1:6)
          statev(7) = 1.d-8
          Encdnew(noel,npt,4) = statev(7)
          Encd=Encdnew 
! ע:��n��n+1���������У�Encd�������£�Encdnew���£��ڼ����У��õ���Encd
! ����������ʱ����������µ�����������Encdnew�е������ǿ��ŵģ���Encd=Encdnew
! ������������ʱ��Encdnew�е������Ѿ�����Ⱦ����Encdnew=Encd��������ǰ�������µ�ֵ
      endif

! �����ж��Ƿ�ʼ���µ��������������Encd = Encdnew��Encdnew���ڴ洢δ����Ⱦ����Ϣ
      if ((Noel.Eq.1).and.(Npt.eq.1)) then
          if  (incem.LT.kinc) then  
               Encd = Encdnew
               incem = kinc
               Iteration_U = 0.d0
          else
               Iteration_U = Iteration_U + 1.d0
          endif     
      endif
      
      
! �����㷨����(�⻬Ħ�����ס���������ţ�ٵ���)
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
! ��ȡ���ϲ���
      EK = props(1)   
      miu = props(2)
      c0 = props(3)
      fphi = props(4)*pi/180.d0
      Psi = props(5)*pi/180.d0
      As = props(6)
      cLt = props(7)     !��������
      index_s = props(8)
      index_nl = props(9)
! �����ϲ�������para��
      para(1:6) = props(1:6)
      para(4) = fphi       ! �Ƕ��򻡶�ת��
      para(5) = Psi       ! �Ƕ��򻡶�ת��
      para(7) = Epsilon
      para(8) = betaS
      para(9) = pi
      para(10) = Maxit
      para(11) = rho
      para(12) = zeta
      para(13) = cLt 
      para(14) = index_s
      para(15) = index_nl  !�Ǿֲ�ϵ��
c �ſɱȾ���DDSDDE��ʼ��Ĭ��Ϊ0,���嵯�ԸնȾ���
      De = 0.d0
      BK = EK/3.d0/(1.d0-2.d0*miu)
      GK = EK/2.d0/(1.d0+miu)
      call TenMat1(Ivol,Isym,delta)
      De = (3.0d0*BK - 2.0d0*GK)*Ivol+2.0d0*GK*Isym
! ��statev(1:6)�ж�ȡ��һ������ЧӦ��
      stress(1:6) = statev(1:6) 
      s = 0.0d0
      dstra = 0.0d0
      call DimIn(stress,dstran,s,dstra,NTENS)    
      epn = statev(7)
! ���㵯����̽Ӧ�������м�ж���ж�
      st = s + matmul(De,dstra)
      xold(1:6,1)=s(1:6,1)
      xold(7,1) = epn
      xold(8,1) = 0.d0
      xnew = xold
      xnew(1:6,1) = st(1:6,1)
      ft = 0.d0
! ������̽��������
      call yieldfun(xnew,para,npara,ft,f_ep,f_s)
! ����һ����Ӧ������Ϊ�����ĳ�ʼ��
      xnew = xold
      call cgfun(xnew,xold,para,npara,De,dstra,9,gi,g)
      normg = norm2(g)  
C##***********************************�Ǿֲ�����*************************************##C
!      call nlfunTra(mw,ca1,ntens,noel,npt,para,npara,time,
!     & nstatv,statev)

      call nlfunOct(mw,ca1,statev,nstatv,COORDS,time)

      
      !if ((statev(7).GE.1.d-5)) then
      !    write(*,*)"Amw-mw,Aca1-ca1",(Amw-mw),(Aca1-ca1)
      !    write(*,*)"mw,ca1",mw,ca1
      !    write(*,*)"Amw,Aca1",Amw,Aca1
      !endif
      !write(*,*)"UMAT,Volint(noel,npt,1)",Volint(noel,npt,1) 
! ca1 ��ʾ�������ֵ�ԷǾֲ�����Ĺ���
      epnw = (1.d0-mw)*statev(7) + ca1
      epnw = max(epnw,1.D-8)
      call dam_fun(damage,dp,epnw,para,npara)
!      damage = statev(11)
C##***********************************�Ǿֲ�����*************************************##C      

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
      
! ����Ǿֲ������ڱ���epnw���Ǿֲ����˱���damage
        epnw = (1.d0-mw)*xnew(7,1) + ca1
        epnw = max(epnw,1.D-8)
        call dam_fun(damage,dp,epnw,para,npara)
! ����һ�������߸նȾ���������Ӧ�����洢��ЧӦ��        
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
! ��¼��ǰ������_��ǰ�������¸����ֵ�������ڱ�����Ҳ�������ڱ����� 
      Encdnew(Noel,Npt,4) = statev(7)      
      call DimOut(DDSDDE,stress,CTO,xnew,NTENS)   
      
      if ((Iteration .GE. Newton).or.isnan(normg).or.isnan(norm2(CTO))
     2.or.isnan(norm2(xnew)).or.(Iteration_U .GE. 50.d0)) then
         PNEWDT = 0.25d0
! ��һ���������Ӧ���͸ն�
         CTO = -1.d60
         stress = 1.d60
         Encdnew  = Encd
         Iteration_U = 0.d0
!         call XIT
      end if        

      Return
      end 

 
      
      