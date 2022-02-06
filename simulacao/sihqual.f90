
program sihqualv2
 
! ************************************************
! 
! Algortimo de solução das eq. de Saint-Venant e Advecção-Dispersão-Reação: SihQual
! Autora: Danieli M. Ferreira | contato: danielimaraferreira@gmail.com
! Versão: Tutorial simplificado (10/03/2021) | Projeto Enquadramento Rio Paranapanema - ANA/UFPR
!
! ************************************************

implicit none

! Declaration

             
integer:: il                    ! auxiliares para leitura de arquivos
integer:: i, k                  ! índices espaço (i) e tempo (k)
integer, parameter:: J=108, NN=1576801
real, parameter:: dx=500., dt=20., DD=200.
real, dimension (J):: cq       
real, dimension(NN)::up
real, dimension(J)::bet4, Porg, sg2, alf2, mi 
real, dimension(J):: c2, Wtot, c1
real, dimension(NN):: concig5
integer, dimension(1)::pos_imp_wq=59

real, dimension(J):: A1, v1, y1
              
! Open files
open(unit=31, file="initial_c.txt")
open(unit=35, file="upboundary.txt")
open(unit=37, file="carga.txt")

! Create file to print the result
open(unit=32, file="res_qualidade.txt")

write(*,*) "Executando..."

! ........... READING DATA ...........

!  Initial conditions
do il=1,J
  read(31,*)c1(il)
end do

!  upstream boundary (up)
do il=1,NN
  read(35,*)up(il)
end do

! cargas
do il=1,J
  read(37,*)Wtot(il)
end do

! ........... READING DATA END ...........

! Calling hydrodynamic subroutine
call hidrodin(A1,v1,y1)

! Kinetic rates
bet4(1:J)=0.01/86400.
Porg(1:J)=0.5
sg2(1:J)=0.
alf2(1:J)=0.01
mi(1:J)=(1./86400.)*(4.*0.01*0.001) 

do k=1,(NN-1)

    cq(1:J)=Wtot(1:J)/(86400.*A1(1:J)*dx) 
    
   do i=2,(J-1)
           
     c2(i)=c1(i)-(dt/(2.*dx))*v1(i)*(c1(i+1)-c1(i-1))+(DD*dt/A1(i))*(A1(i+1)-A1(i-1))/(2.*dx)*(c1(i+1)-c1(i-1))/(2.*dx)   &
     +(DD*dt/(dx**2.))*(c1(i+1)-2.*c1(i)+c1(i-1))+(bet4(i)*Porg(i)*c1(i)+sg2(i)/y1(i)-alf2(i)*mi(i)*0.)*dt*(c1(i))+cq(i)*dt
    
   end do
     
   if (v1(i)*dt/dx>1.0)  then 
          write(*,*)'erro'
   end  if

    ! Boundary conditions
    c2(1)=up(k)/1000. 
    c2(J)=c2(J-1)  

   do i=1,J   
     c1(i)=c2(i)    
   end do
   
! saving results for output file
   concig5(k)=c1(pos_imp_wq(1))*1000.
   concig5(k+1)=c2(pos_imp_wq(1))*1000.
   
end do

do k=1,NN
   write(32,*)concig5(k)
end do

contains


! ********************************
!
!    Subprogram - hydrodynamic
!
! ********************************

subroutine hidrodin(A1,v1,y1)

implicit none

! Declaration

integer:: il                           ! auxiliares para leitura de arquivos
integer:: i, k                         ! índices espaço (i) e tempo (k)
integer, dimension(1):: pos_imp=35
integer, parameter:: J=108, NN=1576801
real, parameter:: dx=500., dt=20., g=9.81, alfa=0.1
real, dimension(J):: y1, v1, b1, So, n, mm
real, dimension(NN):: Qaam, yU, Qjc1
real, dimension(J):: A1, Bt1, Rh1, Sf1 
real, dimension(J):: y2, v2, b2, A2, Rh2, Sf2, Bt2!, yy, vv, SSf, AA, BBt
real, dimension(NN,J):: q

! Open files
open(unit=19, file="vazao-m.txt")  
open(unit=3, file="vel.txt")     
open(unit=5, file="mm.txt") 
open(unit=7, file="cota-m.txt")
open(unit=9, file="bb.txt") 
open(unit=11, file="so.txt") 
open(unit=13, file="manning.txt") 
open(unit=15, file="prof.txt")  
open(unit=17, file="hm_lateralcontrib.txt")

! Create file to print the result
open(unit=4, file="res_hidrodinamico.txt")

! ........... READING DATA ...........

read(17,*)q(:,:)

! upstream hydrograph
do il=1,NN
  read(19,*)Qaam(il)
end do

! upstream water depths
do il=1,NN
  read(7,*)yU(il)
end do

! Initial water depth and velocity
do il=1,J
  read(3,*)v1(il)
end do
do il=1,J
  read(15,*)y1(il)
end do

! Other input data
do il=1,J
  read(5,*)mm(il)
end do
do il=1,J
  read(9,*)b1(il)
end do
do il=1,J
  read(11,*)So(il)
end do
do il=1,J
  read(13,*)n(il)
end do

! ........... READING DATA END ...........

! INITIAL CONDITIONS

do i=1,J
	Bt1(i)=b1(i)+2.d0*mm(i)*y1(i)
	A1(i)=b1(i)*y1(i)+mm(i)*(y1(i)**2.d0)
	Rh1(i)=A1(i)/(b1(i)+2.d0*y1(i)*((1.d0+(mm(i)**2.d0))**0.5d0))
	Sf1(i)=((n(i)*v1(i))**2.d0)/(Rh1(i)**(4.d0/3.d0)) 
end do


do k=1,(NN-1)                        
   
    do i=2,(J-1)  
         
    y2(i)=alfa*y1(i)+(1.d0-alfa)*(0.5d0*(y1(i-1)+y1(i+1)))-0.5d0*(dt/dx)*(0.5d0*(v1(i-1)+v1(i+1)))*(y1(i+1)-y1(i-1))   &
   			-0.5d0*(0.5d0*(v1(i-1)+v1(i+1)))*(dt/dx)*(A1(i+1)-A1(i-1))/(0.5d0*(Bt1(i-1)+Bt1(i+1)))       &  
            -0.5d0*(dt/dx)*((0.5d0*(A1(i-1)+A1(i+1)))/(0.5d0*(Bt1(i-1)+Bt1(i+1))))*(v1(i+1)-v1(i-1))     &
            +(q(k,i)*dt)/(0.5d0*(Bt1(i-1)+Bt1(i+1)))
     
    v2(i)=alfa*v1(i)+(1.d0-alfa)*(0.5d0*(v1(i-1)+v1(i+1)))-0.5d0*(0.5d0*(v1(i-1)+v1(i+1)))*(dt/dx)*(v1(i+1)-v1(i-1))   &
              -0.5d0*g*(dt/dx)*(y1(i+1)-y1(i-1))+g*dt*(So(i)-0.5d0*(Sf1(i-1)+Sf1(i+1)))    &
              +q(k,i)*(v1(i)-0.5*(v1(i-1)+v1(i+1)))*dt/(0.5*(A1(i-1)+A1(i+1)))
      
    b2(i)=b1(i)            
    Bt2(i)=b2(i)+2.d0*mm(i)*y2(i)
    A2(i)=b2(i)*y2(i)+mm(i)*(y2(i)**2.d0)
    Rh2(i)=A2(i)/(b2(i)+2.d0*y2(i)*((1.d0+(mm(i)**2.d0))**0.5d0))
    Sf2(i)=((n(i)*v2(i))**2.d0)/(Rh2(i)**(4.d0/3.d0))
    
    end do
 
    ! Upstream boundary condition
    y2(1)=yU(k)   
        
    b2(1)=b1(1)
    Bt2(1)=b2(1)+2.d0*mm(1)*y2(1)
    A2(1)=b2(1)*y2(1)+mm(1)*(y2(1)**2.d0)
    v2(1)=Qaam(k+1)/A2(1)
    Rh2(1)=A2(1)/(b2(1)+2.d0*y2(1)*((1.d0+(mm(1)**2.d0))**0.5d0))
    Sf2(1)=((n(1)*v2(1))**2.d0)/(Rh2(1)**(4.d0/3.d0))

    ! Downstream boundary condition
    y2(J)=y2(J-1)    

    b2(J)=b1(J)
    Bt2(J)=b2(J)+2.d0*mm(J)*y2(J)
    A2(J)=b2(J)*y2(J)+mm(i)*(y2(J)**2.d0)
    v2(J)=v2(J-1)
    Rh2(J)=A2(J)/(b2(J)+2.d0*y2(J)*((1.d0+(mm(J)**2.d0))**0.5d0))
    Sf2(J)=((n(J)*v2(J))**2.d0)/(Rh2(J)**(4.d0/3.d0))

! Variables redefinition
    do i=1,J                     
    	y1(i)=y2(i)     
    	v1(i)=v2(i)
    	Bt1(i)=Bt2(i)
    	A1(i)=A2(i)
    	Rh1(i)=Rh2(i)
    	Sf1(i)=Sf2(i)
    end do
  
    Qjc1(1)=Qaam(1)     
    
    Qjc1(k+1)=v2(pos_imp(1))*A2(pos_imp(1))

end do
 
do k=1,NN
  write(4,*)Qjc1(k)
end do

return
end subroutine


end program sihqualv2
