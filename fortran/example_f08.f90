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
       print '(i0, 3i4)', i - 1, grid_address(:, i)
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
    print *, ""
  end subroutine get_double_grid_address


  subroutine get_double_grid_index(address_double, mesh)
    use kspclib_f08, only: ksp_get_double_grid_index

    integer(4), intent(in) :: address_double(3), mesh(3)
    integer(8) :: grid_point

    grid_point = ksp_get_double_grid_index(address_double, mesh)

    print '("ksp_get_double_grid_index")'
    print '("(Note that the index starts with 0.)")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Double-grid address", 3i3)', address_double(:)
    print '("Grid point ", i0)', grid_point
    print *, ""
  end subroutine get_double_grid_index


  subroutine get_thm_relative_grid_addresses(rec_lattice)
    use kspclib_f08, only: ksp_get_thm_relative_grid_addresses

    real(8), intent(in) :: rec_lattice(3, 3)
    integer(4) :: relative_grid_addresses(3, 4, 24)
    integer :: i, j

    call ksp_get_thm_relative_grid_addresses(relative_grid_addresses, &
         rec_lattice)

    print '("ksp_get_thm_relative_grid_addresses")'
    print '("rec_lattice:")'
    print '("  a: ", 3f11.8)', rec_lattice(1, :)
    print '("  b: ", 3f11.8)', rec_lattice(2, :)
    print '("  c: ", 3f11.8)', rec_lattice(3, :)
    print *, ""
    do i = 1, 24
       do j = 1, 4
          print '(3i3)', relative_grid_addresses(:, j, i)
       end do
       print *, ""
    end do
    print *, ""
  end subroutine get_thm_relative_grid_addresses


  subroutine get_thm_integration_weight
    use kspclib_f08, only: ksp_get_thm_integration_weight, &
         ksp_get_thm_relative_grid_addresses, ksp_get_all_grid_addresses, &
         ksp_get_double_grid_address, ksp_get_double_grid_index

    integer :: i, j, k, bi, fi
    integer(4) :: ga_d(3), ga(3)
    integer(8) :: gp
    real(8) :: buf(6000)
    real(8) :: rec_lattice(3, 3)
    integer(4) :: relative_grid_addresses(3, 4, 24)
    integer(4) :: mesh(3), shift(3)
    integer(4), allocatable :: grid_address(:, :)
    integer(4) :: tetrahedra_ga(3, 4, 24)
    integer(8) :: tetrahedra_gps(4, 24)
    real(8) :: tetrahedra_freqs(4, 24)
    real(8) :: f
    real(8) :: dos(8), acc(8)

    dos(:) = 0
    acc(:) = 0

    mesh(:) = 10
    shift(:) = 0
    allocate(grid_address(3, product(mesh)))
    call ksp_get_all_grid_addresses(grid_address, mesh)

    rec_lattice(:, :) = 0.1757376132331179
    rec_lattice(1, 1) = -0.1757376132331179
    rec_lattice(2, 2) = -0.1757376132331179
    rec_lattice(3, 3) = -0.1757376132331179
    call ksp_get_thm_relative_grid_addresses(relative_grid_addresses, &
         rec_lattice)

    open(17, file='../python/test/NaCl-freqs-101010.dat', status='old')
    do i = 1, 6000
       read (17, *) buf(i)
    end do
    close(17)

    ! do i = 1, 1000
    !    print '(6f12.8)', buf((i - 1)*6 + 1:  i * 6)
    ! end do

    do i = 1, product(mesh)
       do fi = 1, 8
          do bi = 1, 6
             do j = 1, 24
                do k = 1, 4
                   ga = relative_grid_addresses(:, k, j) + grid_address(:, i)
                   call ksp_get_double_grid_address(ga_d, ga, mesh, shift)
                   gp = ksp_get_double_grid_index(ga_d, mesh)
                   tetrahedra_freqs(k, j) = buf(gp * 6 + bi)
                end do
             end do
             f = fi - 1
             dos(fi) = dos(fi) + &
                  ksp_get_thm_integration_weight(f, tetrahedra_freqs, 'J')
          end do
       end do
    end do

    do fi = 1, 8
       print *, dos(fi) / 1000
    end do

  end subroutine get_thm_integration_weight


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
       get_double_grid_address, get_double_grid_index, &
       get_thm_relative_grid_addresses, get_thm_integration_weight

  integer(8) :: A(3, 3)
  integer(4) :: mesh(3), address(3), is_shift(3), address_double(3)
  real(8) :: rec_lattice(3, 3)

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
  address_double(:) = 3

  rec_lattice(:, :) = 0
  rec_lattice(1, 1) = 3.23752762
  rec_lattice(2, 1) = -1.61876380
  rec_lattice(2, 2) = 2.80378116
  rec_lattice(3, 3) = 5.22655942

  call kspclib_version()
  call snf3x3(A)
  call get_all_grid_addresses(mesh)
  call get_double_grid_address(address, mesh, is_shift)
  call get_double_grid_index(address_double, mesh)
  call get_thm_relative_grid_addresses(rec_lattice)
  call get_thm_integration_weight()

end program kspclib_example_f08
