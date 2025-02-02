program bgk
  use omp_lib  ! Inclure la bibliothèque OpenMP
  implicit none
  integer, parameter :: pr=selected_real_kind(15,70)
  real(pr), parameter :: pi = acos(-1._pr)
  integer,parameter :: nx=500 , nv=100
  real(pr) :: a, b, dx, dv, vmin, vmax, tcount, epsilon
  real(pr) :: tg, ug, td, ud, x, rhog, rhod, pg, pd
  real, dimension(6) :: epsilon_values = [1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6]
  real(pr) :: tmax, cfl
  real(pr), allocatable :: fn(:,:), fnpun(:,:),rho(:), rhou(:), T(:), E(:), v(:)
  integer :: i,j
  real(pr) :: start_time, end_time, duration

  ! --- Définir le nombre de threads OpenMP ---
  call omp_set_num_threads(4)  ! Définir le nombre de proc à 4

  ! --- Conditions aux bords ---
  rhog = 1._pr
  rhod = 0.125_pr
  ug = 0._pr
  ud = 0._pr
  pg = 1._pr
  pd = 0.1_pr
  tg = pg / rhog
  td = pd / rhod
  vmin = -14 * sqrt(tg)
  vmax = 14 * sqrt(tg)
  dv = (vmax - vmin) / nv
  cfl = 0.4_pr
  !--- temps maximum
  tmax =  0.2_pr  

  allocate(fn(0:nx,0:nv), fnpun(0:nx,0:nv),rho(0:nx), rhou(0:nx), T(0:nx), E(0:nx), v(0:nv))

  print*,"nx =",nx
  print*,"nv =",nv
    ! Simulation sans collision
   call init_velocity(vmin, vmax, nv, v)
   call init_fn(nx, nv, v, rhog, ug, tg, rhod, ud, td, fn)
   call sans_collision(nx, nv, tmax, tg, td, rhog, rhod,ug ,ud, vmax, cfl,fn)

     ! Simulation avec collision
      ! Boucle sur les différentes valeurs d'epsilon
    do i = 1, 6
        epsilon = epsilon_values(i)
        call cpu_time(start_time)
        call avec_collision(nx, nv, tmax, epsilon, tg, td, rhog, rhod, ug, ud, vmax, cfl, fn)
        call cpu_time(end_time)
        duration = end_time - start_time
        print *, "Temps d'exécution pour epsilon = ", epsilon, ": ", duration, " secondes"
    end do
   deallocate(fn, fnpun,rho, rhou, T, E, v)


 contains

  subroutine init_velocity(vmin, vmax, nv, v)
    real(pr), intent(in) :: vmin, vmax
    integer, intent(in) :: nv
    real(pr), intent(out) :: v(0:nv)
    integer :: j

    do j = 0, nv 
      v(j) = vmin + j * (vmax - vmin) / nv
    end do
  end subroutine init_velocity

  subroutine init_fn(nx, nv, v, rhog, ug, tg, rhod, ud, td, fn)
    integer, intent(in) :: nx, nv
    real(pr), intent(in) :: rhog, ug, tg, rhod, ud, td
    real(pr), dimension(0:nv), intent(in)::v
    real(pr), dimension(0:nx, 0:nv), intent(inout) :: fn
    integer :: i, j

    !$OMP PARALLEL DO PRIVATE(i, j) SHARED(fn, v, rhog, ug, tg, rhod, ud, td)
    do i = 0, nx
      do j = 0, nv
        if (i < nx / 2) then
          fn(i, j) = Maxwellienne(rhog, 0._pr, v(j), tg)
        else
          fn(i, j) = Maxwellienne(rhod, 0._pr, v(j), td)
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine init_fn

  function Maxwellienne(rho, u, v, T) result(f)
    real(pr), intent(in) :: rho, u, v, T
    real(pr) :: f
    f = (rho / sqrt(2._pr * pi * T)) * exp(-((v - u)**2) / (2 * T))
  end function Maxwellienne

  subroutine macro(nx, nv,dv, v, fn, rho, rhou, T, E)
    integer, intent(in) :: nx, nv
    real(pr), dimension(0:nv), intent(in)::v
    real(pr), dimension(0:nx, 0:nv), intent(inout) :: fn
    real(pr), dimension(0:nx), intent(out) :: rho, rhou, T, E
    real(pr), intent(in) :: dv
    integer :: i

    !$OMP PARALLEL DO PRIVATE(i) SHARED(fn, v, rho, rhou, T, E, dv)
    
    do i = 0, nx
      rho(i) = dv * (sum(fn(i, 0:nv - 1)) + 0.5 * (fn(i, 0) + fn(i, nv)))
      rhou(i) = dv * (sum(fn(i, 1:nv - 1) * v(1:nv - 1)) + &
                      0.5 * (fn(i, 0) * v(0) + fn(i, nv) * v(nv)))
      E(i) = dv * (sum(fn(i, 1:nv - 1) * v(1:nv - 1)**2) + &
                   0.5 * (fn(i, 0) * v(0)**2 + fn(i, nv) * v(nv)**2))
      if (rho(i) > 0._pr) then
        T(i) = E(i) / rho(i) - (rhou(i) / rho(i))**2
      else
        T(i) = 0._pr
      end if
    end do
    !$OMP END PARALLEL DO
  end subroutine macro

  subroutine sans_collision(nx, nv, tmax, tg, td, rhog, rhod,ug ,ud, vmax, cfl,fn)
    integer, intent(in) :: nx, nv
    real(pr), intent(in) :: tmax,tg,td,rhog,rhod,vmax,cfl,ug,ud
    real(pr), dimension(0:nv)::v
    real(pr), dimension(0:nx, 0:nv), intent(inout) :: fn
    real(pr), dimension(0:nx, 0:nv) :: fnpun
    real(pr), dimension(0:nx) :: rho, rhou, T, E
    real(pr) :: tcount, dx, a, b, dt
    integer :: i, j, iter

    tcount = 0.0_pr
    a = -1._pr
    b = 1._pr
    dx = (b-a)/nx
    dt = cfl * dx / vmax   !ici dt est fixe ne varie pas donc on peux le sortir de la boucle 
    call init_velocity(vmin, vmax, nv, v)
    call init_fn(nx, nv, v, rhog, ug, tg, rhod, ud, td, fn)
    iter = 0
    do while (tcount < tmax)
      ! Transport sans collision
      !$OMP PARALLEL DO PRIVATE(i, j) SHARED(fn, fnpun, v, tg, td, rhog, rhod, ug, ud, dx, dt)
     do j = 0, nv
          !if (v(j)>0) then
          fn(0, j) = Maxwellienne(rhog, ug, v(j), tg)   ! Bord gauche
          fnpun(0, j) = Maxwellienne(rhog, ug, v(j), tg)   ! Bord gauche

          !end if 
          !if (v(j)<0) then
            fn(nx, j) = Maxwellienne(rhod, ud, v(j), td)  ! Bord droit 
            fnpun(nx, j) = Maxwellienne(rhod, ud, v(j), td)  ! Bord droit 

          !end if  
      do i = 1, nx - 1
          fnpun(i, j) = fn(i, j) - &
                        (dt / dx) * max(v(j), 0._pr) * (fn(i, j) - fn(i - 1, j)) - &
                        (dt / dx) * min(v(j), 0._pr) * (fn(i + 1, j) - fn(i, j))
        end do
      end do
     ! !$OMP END PARALLEL DO

      ! Mise à jour des variables pour le prochain pas de temps
      fn = fnpun
      call macro(nx, nv, dv, v, fn, rho, rhou, T, E)
      tcount = tcount + dt
    iter = iter + 1 
      
    end do
    print*,"nombre d'itération sans collision :",iter

    open(unit=12, file='rho_transport_bgk.dat')
    open(unit=22, file='u_transport_bgk.dat')
    open(unit=32, file='t_transport_bgk.dat')

    do i = 0, nx
      write(12, *) a + i * dx, rho(i)
      write(22, *) a + i * dx, rhou(i) / rho(i)
      write(32, *) a + i * dx, T(i)
    end do
  end subroutine sans_collision

subroutine avec_collision(nx, nv, tmax, epsilon, tg, td, rhog, rhod,ug ,ud, vmax, cfl,fn)
    integer, intent(in) :: nx, nv
    real(pr), intent(in) :: tmax, epsilon, tg, td, rhog, rhod, vmax, cfl, ug, ud
    real(pr), dimension(0:nx, 0:nv), intent(inout) :: fn
    real(pr), dimension(0:nx, 0:nv) :: fnpun,wn
    real(pr), dimension(0:nx) :: rho, rhou, T, E
    real(pr), dimension(0:nv) :: v 
    real(pr) :: tcount, dt
    integer :: i, j,iter
    tcount = 0.0_pr
    a = -1._pr
    b = 1._pr
    dx = (b-a)/nx
    call init_velocity(vmin, vmax, nv, v)
    call init_fn(nx, nv, v, rhog, ug, tg, rhod, ud, td, fn)
    iter=0
    call macro(nx, nv, dv, v, fn, rho, rhou, T, E)
    do while (tcount < tmax)  
    dt = cfl/((vmax/dx)+(1/(epsilon*maxval(rho))))
    !dt=cfl/((vmax/dx))   !pour lancer le schéma implicite
      ! Transport et collision
      !$OMP PARALLEL DO PRIVATE(i, j) SHARED(fn, fnpun, v, tg, td, rhog, rhod, ug, ud, dt)
    do j = 0, nv
        !if (v(j)>0) then
          fn(0, j) = Maxwellienne(rhog, ug, v(j), tg)   ! Bord gauche
          fnpun(0, j) = Maxwellienne(rhog, ug, v(j), tg)   ! Bord gauche
        !end if 
        !if (v(j)<0) then
          fn(nx, j) = Maxwellienne(rhod, ud, v(j), td)  ! Bord droit 
          fnpun(nx, j) = Maxwellienne(rhod, ud, v(j), td)  ! Bord droit 

        !end if  
        do i = 1, nx - 1
        wn(i,j) = fn(i, j) - &
                        (dt / dx) * max(v(j), 0._pr) * (fn(i, j) - fn(i - 1, j)) - &
                        (dt / dx) * min(v(j), 0._pr) * (fn(i + 1, j) - fn(i, j)) 
        !call macro(nx, nv, dv, v, wn, rho, rhou, T, E)   !pour lancer le schéma implicite
        fnpun(i, j) =( wn(i,j) + &
                        dt * (rho(i) / epsilon) * &
                        (Maxwellienne(rho(i), rhou(i)/rho(i), v(j), T(i)))) / (1+dt * (rho(i) / epsilon))

        end do
    end do

      !$OMP END PARALLEL DO

      ! Mise à jour des variables
      fn = fnpun
      tcount = tcount + dt
      iter = iter + 1 
     call macro(nx, nv, dv, v, fn, rho, rhou, T, E)
    end do
    
    call macro(nx, nv, dv, v, fn, rho, rhou, T, E)
    print*,"nombre d'itération avec collision epsilon=",epsilon,":",iter
    call write_file_with_epsilon('rho_transport_col', epsilon, 1, rho, nx, dx)
    call write_file_with_epsilon('u_transport_col', epsilon, 2, rhou / rho, nx, dx)
    call write_file_with_epsilon('t_transport_col', epsilon, 3, T, nx, dx)
  end subroutine avec_collision
  subroutine write_file_with_epsilon(base_name, epsilon_val, unit_num, data, nx, dx)

    character(len=*) :: base_name
    real(pr),intent(in) :: epsilon_val,dx
    integer,intent(in) :: unit_num,nx
    real(pr), dimension(0:nx) :: data
    integer :: i
    real(pr) :: x
    character(len=100) :: filename, epsilon_str

    ! Formatage de epsilon pour l'inclure dans le nom du fichier
    write(epsilon_str, '(1pE14.5)') epsilon_val
    filename = trim(adjustl(base_name)) // "_" // trim(adjustl(epsilon_str)) // '.dat'

    ! Ouverture du fichier
    open(unit=unit_num, file=filename)

    ! Écriture des données dans le fichier
    do i = 0, nx
      x = -1_pr + i * dx
      write(unit_num, *) x, data(i)
    end do

    ! Fermeture du fichier
    close(unit_num)
  end subroutine write_file_with_epsilon

end program bgk