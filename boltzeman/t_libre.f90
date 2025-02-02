!-----------------------------------------------------------------------------------------!
!  Code résolvant l'équation du transport libre 1D/1D pour le problème de sod
!  - partie 1 : modules (mod_precision)
!  - partie 2 : programme principal (nom transport)
!-----------------------------------------------------------------------------------------!


!---------------------------------------- modules ----------------------------------------!

module mod_precision

  implicit none

  integer, parameter :: pr=selected_real_kind(15,307)

end module mod_precision


!----------------------------------- programme principal ---------------------------------!
program transport

  use mod_precision

  implicit none

  integer :: imax,i
  real(pr), parameter :: pi = acos(-1._pr)
  real(pr) :: a, b, dx
  real(pr) :: tg,ug,td,ud,x,rhog,rhod,pg,pd
  real(pr) :: xg, xd
  real(pr) :: tmax
  real(pr) :: rho,u,p,t,e,rhou

    

  !--- bornes du domaine
  a = -1._pr ; b = 1._pr

  !--- valeurs aux bords
  rhog = 1._pr ; rhod = 0.125_pr
  ug = 0._pr ; ud = 0._pr
  pg = 1._pr ; pd = 0.1_pr
  tg=pg/rhog ; td=pd/rhod

  !--- maillage en espace pour la sortie fichier
  imax = 300
  dx = (b-a)/imax

  !--- temps maximum
  tmax =  0.2_pr

 
  !--- tracé de la solution finale
  open(unit=12,file='rho_transport.dat')
  open(unit=22,file='u_transport.dat')
  open(unit=32,file='t_transport.dat')
  open(unit=42,file='p_transport.dat')

  do i=1,imax

     x = a+dx/2 + (i-1)*dx

     xg = x/(tmax*sqrt(2._pr*tg))
     xd = x/(tmax*sqrt(2._pr*td))

     rho = 0.5_pr * (rhog*erfc(xg) + rhod*erfc(-xd)) 
         
     rhou = 0.5_pr * ( rhog*sqrt(2._pr*tg)/sqrt(pi) *exp(-xg**2)            &
          &            - rhod*sqrt(2._pr*td)/sqrt(pi) *exp(-xd**2) )
     u = rhou/rho

     e = 0.5_pr * ( rhog*tg*( xg/sqrt(pi)*exp(-xg**2) + 0.5_pr*erfc(xg)  )    &
          &          + rhod*td*( -xd/sqrt(pi)*exp(-xd**2) + 0.5_pr*erfc(-xd)  )    )
     t = ( e - 0.5_pr*rho*u**2) / (0.5_pr*rho)

     p = rho * t
     
     write(12,*) x,rho
     write(22,*) x,u
     write(32,*) x,t
     write(42,*) x,p

  end do

  close(12)
  close(22)
  close(32)
  close(42)

end program transport


