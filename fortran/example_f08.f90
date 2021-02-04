module kspclib_example
  implicit none

  contains

  subroutine snf3x3(A)
    use kspclib_f08, only: ksp_get_snf3x3

    integer(8), intent(in) :: A(3, 3)
    integer(4) :: succeeded
    integer(8) :: D_diag(3), P(3, 3), Q(3, 3)


    succeeded = ksp_get_snf3x3(D_diag, P, Q, A)

    write(*,*) "D_diag", D_diag
    write(*,*) ""
    write(*,*) "P", P(1, :)
    write(*,*) "P", P(2, :)
    write(*,*) "P", P(3, :)
    write(*,*) ""
    write(*,*) "Q", Q(1, :)
    write(*,*) "Q", Q(2, :)
    write(*,*) "Q", Q(3, :)
    write(*,*) ""
  end subroutine snf3x3

  subroutine kspclib_version()
    use kspclib_f08, only: ksp_get_major_version, ksp_get_minor_version, &
         ksp_get_micro_version

    integer :: major, minor, micro

    major = ksp_get_major_version()
    minor = ksp_get_minor_version()
    micro = ksp_get_micro_version()
    print '("Kspclib version ", i0, ".", i0, ".", i0)', major, minor, micro
  end subroutine kspclib_version

end module kspclib_example

program kspclib_example_f08

  use kspclib_example, only: snf3x3, kspclib_version

  integer(8) :: A(3, 3)

  A(1, 1) = -1
  A(2, 1) = 1
  A(3, 1) = 1
  A(1, 2) = 1
  A(2, 2) = -1
  A(3, 2) = 1
  A(1, 3) = 1
  A(2, 3) = 1
  A(3, 3) = -1

  call kspclib_version()

  call snf3x3(A)

end program kspclib_example_f08
