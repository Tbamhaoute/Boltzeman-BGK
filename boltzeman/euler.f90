!---------------------------------------------------------------------------------------------!
!  Code résolvant les équations d'Euler avec gamma=3 (problème de Sod, schéma de Roe, ordre 1)
!  - partie 1 : modules ()
!  - partie 2 : programme principal (nom euler)
!---------------------------------------------------------------------------------------------!


!---------------------------------------- modules ----------------------------------------!

module mod_precision

  implicit none

  integer, parameter :: qr=selected_real_kind(15,70)

end module mod_precision

module mod_gamma

  use mod_precision

  implicit none

  real(qr), parameter :: gamma = 3._qr

end module mod_gamma

module mod_decode

  use mod_precision
  use mod_gamma

contains

  subroutine decode(mts,ro,u,p,t,c)

    implicit none

    !--- entrees
    real(qr), dimension(:,:), intent(in) :: mts

    !--- sorties
    real(qr), dimension(1:size(mts,2)), intent(out) :: ro,u,p,t,c

    !--- locales
    integer :: i
    real(qr) :: e


    do i=1,size(mts,2)
       ro(i)=mts(1,i)

       if (ro(i)<=0._qr) then
          print*,'mod_decode : densite <=0 : densite=', ro(i)
          stop
       end if

       e=mts(3,i)
       u(i)=mts(2,i)/ro(i)
       p(i)=(e-ro(i)*(u(i)**2)/2._qr)*(gamma-1._qr)
       t(i)=p(i)/ro(i)
       c(i)=sqrt(gamma*p(i)/ro(i));

    end do

  end subroutine decode

end module mod_decode

module mod_flux_roe

contains

  function flux_roe(mtsl,mtsr) result(flux)

    use mod_precision
    use mod_gamma

    implicit none

    !--- entrees
    real(qr), dimension(1:3), intent(in) :: mtsl,mtsr

    !--- sorties
    real(qr), dimension(1:3) :: flux

    !--- locales
    real(qr) :: rol,ul,pl,El,Hl
    real(qr) :: ror,ur,pr,Er,Hr
    real(qr) :: ut,Ht,at
    real(qr), dimension(1:3) :: l,K1,K2,K3,du,alpha
    real(qr), dimension(1:3) :: fluxl,fluxr


    !--- quantites ro,u,H gauche et droite
    call macros(mtsl,rol,ul,pl,El,Hl)
    call macros(mtsr,ror,ur,pr,Er,Hr)

    !--- test sur la validite de l'entree
    if ((rol<=0._qr).or.(ror<=0._qr)) then
       print*,'flux_roe : rol ou ror incorrects : rol,ror=',rol,ror
    end if

    !--- quantites moyennes (toro, 11.60)
    ut=(sqrt(rol)*ul+sqrt(ror)*ur)/(sqrt(rol)+sqrt(ror));
    Ht=(sqrt(rol)*Hl+sqrt(ror)*Hr)/(sqrt(rol)+sqrt(ror));
    at=sqrt((gamma-1)*(Ht-0.5*ut**2));
    if (at<=0._qr) then
       print*,'flux_roe : at incorrecte : at=',at
    end if

    !--- valeurs propres (toro 11.58)
    l(1)=ut-at;
    l(2)=ut;
    l(3)=ut+at;

    !--- vecteurs Ki (toro 11.59)
    K1=(/1._qr , ut-at , Ht-ut*at/)
    K2=(/1._qr , ut , 0.5*ut**2/)
    K3=(/1._qr , ut+at , Ht+ut*at/)

    !--- alphai (toro 11.68-11.70)
    du=mtsr-mtsl

    alpha(2)=(gamma-1)/(at**2)*(du(1)*(Ht-ut**2)+ut*du(2)-du(3))
    alpha(1)=(du(1)*(ut+at)-du(2)-at*alpha(2))/(2*at)
    alpha(3)=du(1)-(alpha(1)+alpha(2))

    !--- flux
    fluxl=(/rol*ul , rol*ul**2+pl , ul*(El+pl)/)
    fluxr=(/ror*ur , ror*ur**2+pr , ur*(Er+pr)/)

    flux=0.5_qr*(fluxl+fluxr)                       &
         &          -0.5*(alpha(1)*abs(l(1))*K1+alpha(2)*abs(l(2))*K2+alpha(3)*abs(l(3))*K3);

  end function flux_roe

  subroutine macros(mts,ro,u,p,E,H)

    use mod_precision
    use mod_gamma

    implicit none

    !--- entrees
    real(qr), dimension(1:3), intent(in) :: mts

    !--- sorties
    real(qr), intent(out) :: ro,u,p,e,h

    ro=mts(1)
    e=mts(3)
    if (ro<=0._qr) then  
       print*,'macros : densite <=0 : densite=',ro
    end if
    u=mts(2)/ro
    p=(e-ro*(u**2)/2._qr)*(gamma-1._qr)
    h=(e+p)/ro

  end subroutine macros

end module mod_flux_roe

!---------------------------------------- programme principal ----------------------------------------!
program euler

  use mod_precision
  use mod_gamma
  use mod_decode
  use mod_flux_roe

  implicit none

  integer :: N,i
  real(qr) :: rol,ul,pl,ror,ur,pr
  real(qr) :: tmax,cfl,temps,dt
  real(qr) :: dx,xi,a,b
  real(qr), dimension(1:3) :: fluxp,fluxm
  real(qr), dimension(:), allocatable :: x
  real(qr), dimension(:), allocatable :: ro,u,p,t,c
  real(qr), dimension(:,:), allocatable :: mts,mtsn


  !--- tube a choc de Sod, euler, gamma=3, euler, roe ordre 1

  !--- parametres physiques (état gauche, état droit, tmax)
  rol = 1.0_qr ; ul = 0._qr ; pl = 1.0_qr 
  ror = 0.125_qr ; ur = 0._qr ; pr = 0.1_qr 
  tmax = 0.2_qr

  !--- parametres numeriques
  N = 300
  a=-1._qr ; b=1._qr
  dx = (b-a)/N
  allocate(x(1:N))
  x = (/(a+dx/2+(i-1)*dx,i=1,N)/)
  cfl = 0.9_qr


  !--- donnee initiale
  temps = 0._qr;
  allocate(mts(1:3,1:N),mtsn(1:3,1:N))
  mts(1,1:N/2) = rol ; mts(2,1:N/2) = rol*ul ; mts(3,1:N/2) = rol*ul**2/2+1/(gamma-1)*pl 
  mts(1,N/2+1:N) = ror ; mts(2,N/2+1:N) = ror*ur ; mts(3,N/2+1:N) = ror*ur**2/2+1/(gamma-1)*pr 
  allocate(ro(1:N),u(1:N),p(1:N),t(1:N),c(1:N))
  call decode(mts,ro,u,p,t,c)

  !--- boucle en temps
  do while (temps<tmax)

     !--- pas de temps
     dt=cfl*dx/maxval(abs(u)+c)

     !--- boucle espace
     do i=1,N

        if (i==1) then
           fluxp=flux_roe(mts(1:3,i),mts(1:3,i+1))
           fluxm=flux_roe(mts(1:3,i),mts(1:3,i))
        elseif (i==N) then
           fluxp=flux_roe(mts(1:3,i),mts(1:3,i))
           fluxm=flux_roe(mts(1:3,i-1),mts(1:3,i))
        else
           fluxp=flux_roe(mts(1:3,i),mts(1:3,i+1))
           fluxm=flux_roe(mts(1:3,i-1),mts(1:3,i))
        end if

        mtsn(1:3,i) = mts(1:3,i) -dt/dx*(fluxp-fluxm)

     end do

     mts=mtsn

     call decode(mts,ro,u,p,t,c)

     temps=temps+dt

  end do


  !--- tracé de la solution finale
  open(unit=12,file='rho_euler.dat')
  open(unit=22,file='u_euler.dat')
  open(unit=32,file='t_euler.dat')
  open(unit=42,file='p_euler.dat')

  do i=1,n

     xi = a+dx/2 + (i-1)*dx

     write(12,*) xi,ro(i)
     write(22,*) xi,u(i)
     write(32,*) xi,t(i)
     write(42,*) xi,p(i)

  end do

  close(12)
  close(22)
  close(32)
  close(42)


end program euler