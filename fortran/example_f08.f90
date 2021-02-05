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


  subroutine get_all_grid_addresses()
    use kspclib_f08, only: ksp_get_all_grid_addresses

    integer(4) :: mesh(3)
    integer(4), allocatable :: grid_address(:, :)
    integer :: i

    mesh(:) = 4

    allocate(grid_address(3, product(mesh)))
    call ksp_get_all_grid_addresses(grid_address, mesh)

    print '("ksp_get_all_grid_addresses example")'
    print '("All grid addresses for sampling mesh ", i0, "x", i0, "x", i0)', mesh(:)
    do i = 1, product(mesh)
       print '(i0, 3i4)', i - 1, grid_address(:, i)
    end do
    deallocate(grid_address)
  end subroutine get_all_grid_addresses


  subroutine get_double_grid_address()
    use kspclib_f08, only: ksp_get_double_grid_address

    integer(4) :: address(3), mesh(3), is_shift(3)
    integer(4) :: address_double(3)

    mesh(:) = 4
    address(:) = 1
    is_shift(:) = 1

    call ksp_get_double_grid_address(address_double, address, mesh, is_shift)

    print '("ksp_get_double_grid_address")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Half-shift", 3i3)', is_shift(:)
    print '("Double-grid address (address * 2 + is_shift)", 3i3)', address_double(:)
    print *, ""
  end subroutine get_double_grid_address


  subroutine get_double_grid_index()
    use kspclib_f08, only: ksp_get_double_grid_index

    integer(4) :: address_double(3), mesh(3)
    integer(8) :: grid_point

    mesh(:) = 4
    address_double(:) = 3
    grid_point = ksp_get_double_grid_index(address_double, mesh)

    print '("ksp_get_double_grid_index")'
    print '("(Note that the index starts with 0.)")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Double-grid address", 3i3)', address_double(:)
    print '("Grid point ", i0)', grid_point
    print *, ""
  end subroutine get_double_grid_index


  subroutine get_thm_relative_grid_addresses()
    use kspclib_f08, only: ksp_get_thm_relative_grid_addresses

    real(8) :: rec_lattice(3, 3)
    integer(4) :: relative_grid_addresses(3, 4, 24)
    integer :: i, j

    rec_lattice(:, :) = 0
    rec_lattice(1, 1) = 3.23752762
    rec_lattice(2, 1) = -1.61876380
    rec_lattice(2, 2) = 2.80378116
    rec_lattice(3, 3) = 5.22655942

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


  !> @brief
  !> Fortran implementation of test_get_thm_integration_weight
  !> in test_tetrahdron_method.py.
  subroutine get_thm_integration_weight()
    use kspclib_f08, only: ksp_get_thm_integration_weight, &
         ksp_get_thm_relative_grid_addresses, ksp_get_all_grid_addresses, &
         ksp_get_double_grid_address, ksp_get_double_grid_index

    integer(8) :: i, j, k, gi, bi, fi, num_freq_points, num_bands
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
    real(8), allocatable :: freq_points(:), dos(:), acc(:)

    print '("ksp_get_thm_integration_weight")'

    mesh(:) = 10
    shift(:) = 0

    num_bands = 6
    num_freq_points = 8

    allocate(grid_address(3, product(mesh)))
    allocate(freq_points(num_freq_points))
    allocate(dos(num_freq_points))
    allocate(acc(num_freq_points))

    dos(:) = 0
    acc(:) = 0
    do fi = 1, num_freq_points
       freq_points(fi) = fi - 1
    end do

    open(17, file='../python/test/NaCl-freqs-101010.dat', status='old')
    do i = 1, product(mesh) * num_bands
       read (17, *) buf(i)
    end do
    close(17)

    rec_lattice(:, :) = 0.1757376132331179
    rec_lattice(1, 1) = -0.1757376132331179
    rec_lattice(2, 2) = -0.1757376132331179
    rec_lattice(3, 3) = -0.1757376132331179

    call ksp_get_all_grid_addresses(grid_address, mesh)
    call ksp_get_thm_relative_grid_addresses(relative_grid_addresses, &
         rec_lattice)

    do gi = 1, product(mesh)
       do fi = 1, 8
          do bi = 1, num_bands
             do j = 1, 24
                do k = 1, 4
                   ga = relative_grid_addresses(:, k, j) + grid_address(:, gi)
                   call ksp_get_double_grid_address(ga_d, ga, mesh, shift)
                   gp = ksp_get_double_grid_index(ga_d, mesh)
                   tetrahedra_freqs(k, j) = buf(gp * 6 + bi)
                end do
             end do
             dos(fi) = dos(fi) + ksp_get_thm_integration_weight(freq_points(fi), &
                  tetrahedra_freqs, "I")
             acc(fi) = acc(fi) + ksp_get_thm_integration_weight(freq_points(fi), &
                  tetrahedra_freqs, "J")
          end do
       end do
    end do

    do fi = 1, num_freq_points
       print *, freq_points(fi), dos(fi) / product(mesh), acc(fi) / product(mesh)
    end do

    print *, ""

    deallocate(acc)
    deallocate(dos)
    deallocate(grid_address)
  end subroutine get_thm_integration_weight


  subroutine snf3x3()
    use kspclib_f08, only: ksp_get_snf3x3

    integer(8) :: A(3, 3)
    integer(4) :: succeeded
    integer(8) :: D_diag(3), D_diag_ref(3), P_ref(3, 3), Q_ref(3, 3), P(3, 3), Q(3, 3)

    A(:, :) = reshape([0, 5, 5, 5, 0, 5, 2, 2, 0], [3, 3])
    D_diag_ref(:) = [1, 5, 20]
    P_ref(:, :) = reshape([0, 1, -2, 1, 0, 0, -2 ,2, -5], [3, 3])
    Q_ref(:, :) = reshape([1, -5, -9, 0, 0, -1, 0, 1, 1], [3, 3])

    succeeded = ksp_get_snf3x3(D_diag, P, Q, A)

    print '("ksp_get_snf3x3")'
    print '("D_diag ref", 3i5)', D_diag_ref
    print '("D_diag ret", 3i5)', D_diag
    print *, ""
    print '("P_ref", 9i5)', P_ref
    print '("P_ret", 9i5)', P
    print *, ""
    print '("Q_ref", 9i5)', Q_ref
    print '("Q_ret", 9i5)', Q
    print *, ""
  end subroutine snf3x3

  !> @brief
  !> Fortran implementation of test_snf_transform_rotations
  !> in test_grgrid.py (transformed_rots_2).
  subroutine snf_transform_rotations()
    use kspclib_f08, only: ksp_snf_transform_rotations

    integer :: i, j, k
    integer(8) :: transformed_rots(3, 3, 16), transformed_rots_ret(3, 3, 16)
    integer(4) :: rotations(3, 3, 16)
    integer(4) :: num_rot, succeeded
    integer(8) :: D_diag(3)
    integer(8) :: Q(3, 3)


    num_rot = 16
    rotations(:, :, :) = reshape([ &
         1, 0, 0, 0, 1, 0, 0, 0, 1, &
         0, 0, -1, 1, 1, 1, 0, -1, 0, &
         0, 1, 0, 1, 0, 0, -1, -1, -1, &
         1, 1, 1, 0, 0, -1, -1, 0, 0, &
         -1, -1, -1, 0, 0, 1, 0, 1, 0, &
         0, -1, 0, -1, 0, 0, 0, 0, -1, &
         0, 0, 1, -1, -1, -1, 1, 0, 0, &
         -1, 0, 0, 0, -1, 0, 1, 1, 1, &
         -1, 0, 0, 0, -1, 0, 0, 0, -1, &
         0, 0, 1, -1, -1, -1, 0, 1, 0, &
         0, -1, 0, -1, 0, 0, 1, 1, 1, &
         -1, -1, -1, 0, 0, 1, 1, 0, 0, &
         1, 1, 1, 0, 0, -1, 0, -1, 0, &
         0, 1, 0, 1, 0, 0, 0, 0, 1, &
         0, 0, -1, 1, 1, 1, -1, 0, 0, &
         1, 0, 0, 0, 1, 0, -1, -1, -1], [3, 3, 16])

    transformed_rots(:, :, :) = reshape([ &
         1, 0, 0, 0, 1, 0, 0, 0, 1, &
         -4, 3, 2, 5, -4, -2, -20, 16, 9, &
         -9, 8, 4, 0, -1, 0, -20, 20, 9, &
         -4, 5, 2, -5, 4, 2, 0, 4, 1, &
         -1, 0, 0, 0, 1, 0, 0, -4, -1, &
         4, -5, -2, -5, 4, 2, 20, -20, -9, &
         9, -8, -4, 0, -1, 0, 20, -16, -9, &
         4, -3, -2, 5, -4, -2, 0, 0, -1, &
         -1, 0, 0, 0, -1, 0, 0, 0, -1, &
         4, -3, -2, -5, 4, 2, 20, -16, -9, &
         9, -8, -4, 0, 1, 0, 20, -20, -9, &
         4, -5, -2, 5, -4, -2, 0, -4, -1, &
         1, 0, 0, 0, -1, 0, 0, 4, 1, &
         -4, 5, 2, 5, -4, -2, -20, 20, 9, &
         -9, 8, 4, 0, 1, 0, -20, 16, 9, &
         -4, 3, 2, -5, 4, 2, 0, 0, 1], [3, 3, 16])

    D_diag(:) = [1, 5, 20]
    Q(:, :) = reshape([1, -5, -9, 0, 0, -1, 0, 1, 1], [3, 3])

    succeeded = ksp_snf_transform_rotations(transformed_rots_ret, rotations, &
         num_rot, D_diag, Q)

    print '("ksp_snf_transform_rotations")'
    do i = 1, 16
       print '("ref", 9i5)', transformed_rots(:, :, i)
       print '("ret", 9i5)', transformed_rots_ret(:, :, i)
       print *, ""
    end do
    print *, ""

  end subroutine snf_transform_rotations


end module kspclib_example


program kspclib_example_f08

  use kspclib_example, only: kspclib_version, get_all_grid_addresses, &
       get_double_grid_address, get_double_grid_index, &
       get_thm_relative_grid_addresses, get_thm_integration_weight, &
       snf3x3, snf_transform_rotations

  call kspclib_version()
  call get_all_grid_addresses()
  call get_double_grid_address()
  call get_double_grid_index()
  call get_thm_relative_grid_addresses()
  call get_thm_integration_weight()
  call snf3x3()
  call snf_transform_rotations()

end program kspclib_example_f08
