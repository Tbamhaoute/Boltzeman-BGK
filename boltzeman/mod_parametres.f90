module mod_parametres

  implicit none

  !--- double precision
  integer, parameter :: pr = selected_real_kind(15,307)

  !--- nombre pi
  real(pr), parameter :: pi = acos(-1._pr)
  
end module mod_parametres

