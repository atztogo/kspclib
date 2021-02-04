module kspclib_example
  implicit none

  contains

  subroutine kspclib_version()
    use kspclib_f08, only: ksp_get_major_version, ksp_get_minor_version, &
         ksp_get_micro_version

    integer :: major, minor, micro

    major = ksp_get_major_version()
    minor = ksp_get_minor_version()
    micro = ksp_get_micro_version()
    print '("ksp_get_*_version example")'
    print '("Kspclib version ", i0, ".", i0, ".", i0)', major, minor, micro
    print *
  end subroutine kspclib_version


  subroutine get_all_grid_addresses(mesh)
    use kspclib_f08, only: ksp_get_all_grid_addresses

    integer(4), intent(in) :: mesh(3)
    integer(4), allocatable :: grid_address(:, :)
    integer :: i

    allocate(grid_address(3, product(mesh)))
    call ksp_get_all_grid_addresses(grid_address, mesh)

    print '("ksp_get_all_grid_addresses example")'
    print '("All grid addresses for sampling mesh ", i0, "x", i0, "x", i0)', mesh(:)
    do i = 1, product(mesh)
       print *, grid_address(:, i)
    end do
  end subroutine get_all_grid_addresses


  subroutine get_double_grid_address(address, mesh, is_shift)
    use kspclib_f08, only: ksp_get_double_grid_address

    integer(4), intent(in) :: address(3), mesh(3), is_shift(3)
    integer(4) :: address_double(3)

    call ksp_get_double_grid_address(address_double, address, mesh, is_shift)

    print '("ksp_get_double_grid_address")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Half-shift", 3i3)', is_shift(:)
    print '("Double-grid address (address * 2 + is_shift)", 3i3)', address_double(:)
  end subroutine get_double_grid_address


  subroutine snf3x3(A)
    use kspclib_f08, only: ksp_get_snf3x3

    integer(8), intent(in) :: A(3, 3)
    integer(4) :: succeeded
    integer(8) :: D_diag(3), P(3, 3), Q(3, 3)


    succeeded = ksp_get_snf3x3(D_diag, P, Q, A)

    print '("ksp_get_snf3x3 example")'
    print *, "D_diag", D_diag
    print *, ""
    print *, "P", P(1, :)
    print *, "P", P(2, :)
    print *, "P", P(3, :)
    print *, ""
    print *, "Q", Q(1, :)
    print *, "Q", Q(2, :)
    print *, "Q", Q(3, :)
    print *, ""
  end subroutine snf3x3
end module kspclib_example

program kspclib_example_f08

  use kspclib_example, only: snf3x3, kspclib_version, get_all_grid_addresses, &
       get_double_grid_address

  integer(8) :: A(3, 3)
  integer(4) :: mesh(3), address(3), is_shift(3)

  A(1, 1) = -1
  A(2, 1) = 1
  A(3, 1) = 1
  A(1, 2) = 1
  A(2, 2) = -1
  A(3, 2) = 1
  A(1, 3) = 1
  A(2, 3) = 1
  A(3, 3) = -1

  mesh(:) = 4
  address(:) = 1
  is_shift(:) = 1

  call kspclib_version()
  call snf3x3(A)
  call get_all_grid_addresses(mesh)
  call get_double_grid_address(address, mesh, is_shift)

end program kspclib_example_f08
