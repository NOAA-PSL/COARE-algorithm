module machine
  implicit none
  private

  integer, parameter :: kind_phys = KIND(1.0)
  integer, parameter :: kind_rad = KIND(1.0)

  public :: kind_phys, kind_rad

end module machine
