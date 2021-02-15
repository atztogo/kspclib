module kspclib_example
  use iso_c_binding, only: c_long, c_double
  implicit none
  private
  public kspclib_version, get_all_grid_addresses, get_double_grid_address, &
       get_double_grid_index, get_thm_relative_grid_addresses, &
       get_thm_integration_weight, snf3x3, snf_transform_rotations, &
       get_all_grgrid_addresses, get_double_grgrid_address, &
       get_grgrid_index, get_double_grgrid_index, &
       get_grgrid_address_from_index, rotate_grgrid_index, &
       get_ir_grgrid_map, get_reciprocal_point_group
  contains

  subroutine kspclib_version()
    use kspclib_f08, only: ksp_get_major_version, ksp_get_minor_version, &
         ksp_get_micro_version

    integer(c_long) :: major, minor, micro

    major = ksp_get_major_version()
    minor = ksp_get_minor_version()
    micro = ksp_get_micro_version()
    print '("[ksp_get_*_version]")'
    print '("Kspclib version ", i0, ".", i0, ".", i0)', major, minor, micro
    print *
  end subroutine kspclib_version


  subroutine get_all_grid_addresses()
    use kspclib_f08, only: ksp_get_all_grid_addresses
    integer(c_long) :: mesh(3)
    integer(c_long), allocatable :: grid_address(:, :)
    integer :: i

    mesh(:) = 4

    allocate(grid_address(3, product(mesh)))
    call ksp_get_all_grid_addresses(grid_address, mesh)

    print '("[ksp_get_all_grid_addresses]")'
    print '("All grid addresses for sampling mesh ", i0, "x", i0, "x", i0)', mesh(:)
    do i = 1, product(mesh)
       print '(i0, 3i4)', i - 1, grid_address(:, i)
    end do
    print *, ""
    deallocate(grid_address)
  end subroutine get_all_grid_addresses


  subroutine get_double_grid_address()
    use kspclib_f08, only: ksp_get_double_grid_address

    integer(c_long) :: address(3), mesh(3), is_shift(3)
    integer(c_long) :: address_double(3)

    mesh(:) = 4
    address(:) = 1
    is_shift(:) = 1

    call ksp_get_double_grid_address(address_double, address, mesh, is_shift)

    print '("[ksp_get_double_grid_address]")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Half-shift", 3i3)', is_shift(:)
    print '("Double-grid address (address * 2 + is_shift)", 3i3)', address_double(:)
    print *, ""
  end subroutine get_double_grid_address


  subroutine get_double_grid_index()
    use kspclib_f08, only: ksp_get_double_grid_index

    integer(c_long) :: address_double(3), mesh(3)
    integer(c_long) :: grid_index

    mesh(:) = 4
    address_double(:) = 3
    grid_index = ksp_get_double_grid_index(address_double, mesh)

    print '("[ksp_get_double_grid_index]")'
    print '("(Note that the index starts with 0.)")'
    print '("Mesh", 3i3)', mesh(:)
    print '("Double-grid address", 3i3)', address_double(:)
    print '("Grid index ", i0)', grid_index
    print *, ""
  end subroutine get_double_grid_index


  subroutine get_thm_relative_grid_addresses()
    use kspclib_f08, only: ksp_get_thm_relative_grid_addresses

    real(c_double) :: rec_lattice(3, 3)
    integer(c_long) :: relative_grid_addresses(3, 4, 24)
    integer :: i, j

    rec_lattice(:, :) = 0
    rec_lattice(1, 1) = 3.23752762
    rec_lattice(2, 1) = -1.61876380
    rec_lattice(2, 2) = 2.80378116
    rec_lattice(3, 3) = 5.22655942

    call ksp_get_thm_relative_grid_addresses(relative_grid_addresses, &
         rec_lattice)

    print '("[ksp_get_thm_relative_grid_addresses]")'
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

    integer(c_long) :: i, j, k, gi, bi, fi, num_freq_points, num_bands
    integer(c_long) :: ga_d(3), ga(3)
    integer(c_long) :: gp
    real(c_double) :: buf(6000)
    real(c_double) :: rec_lattice(3, 3)
    integer(c_long) :: relative_grid_addresses(3, 4, 24)
    integer(c_long) :: mesh(3), shift(3)
    integer(c_long), allocatable :: grid_address(:, :)
    integer(c_long) :: tetrahedra_ga(3, 4, 24)
    integer(c_long) :: tetrahedra_gps(4, 24)
    real(c_double) :: tetrahedra_freqs(4, 24)
    real(c_double), allocatable :: freq_points(:), dos(:), acc(:)

    print '("[ksp_get_thm_integration_weight]")'

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

    integer(c_long) :: A(3, 3)
    integer(c_long) :: succeeded
    integer(c_long) :: D_diag(3), D_diag_ref(3), P_ref(3, 3), Q_ref(3, 3), P(3, 3), Q(3, 3)

    A(:, :) = reshape([0, 5, 5, 5, 0, 5, 2, 2, 0], [3, 3])
    D_diag_ref(:) = [1, 5, 20]
    P_ref(:, :) = reshape([0, 1, -2, 1, 0, 0, -2 ,2, -5], [3, 3])
    Q_ref(:, :) = reshape([1, -5, -9, 0, 0, -1, 0, 1, 1], [3, 3])

    succeeded = ksp_get_snf3x3(D_diag, P, Q, A)

    print '("[ksp_get_snf3x3]")'
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
    integer(c_long) :: transformed_rots(3, 3, 16), transformed_rots_ret(3, 3, 16)
    integer(c_long) :: rotations(3, 3, 16)
    integer(c_long) :: num_rot, succeeded
    integer(c_long) :: D_diag(3)
    integer(c_long) :: Q(3, 3)


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

    print '("[ksp_snf_transform_rotations]")'
    do i = 1, 16
       print '("ref", 9i5)', transformed_rots(:, :, i)
       print '("ret", 9i5)', transformed_rots_ret(:, :, i)
       print *, ""
    end do
    print *, ""

  end subroutine snf_transform_rotations


  subroutine get_all_grgrid_addresses()
    use kspclib_f08, only: ksp_get_all_grgrid_addresses

    integer(c_long) :: D_diag(3)
    integer(c_long), allocatable :: grid_address(:, :)
    integer :: i

    D_diag(:) = 4

    allocate(grid_address(3, product(D_diag)))
    call ksp_get_all_grgrid_addresses(grid_address, D_diag)

    print '("[ksp_get_all_grgrid_addresses]")'
    print '("All generalized grid addresses for D_diag ", i0, "x", i0, "x", i0)', &
         D_diag(:)
    do i = 1, product(D_diag)
       print '(i0, 3i4)', i - 1, grid_address(:, i)
    end do
    print *, ""
    deallocate(grid_address)
  end subroutine get_all_grgrid_addresses


  subroutine get_double_grgrid_address()
    use kspclib_f08, only: ksp_get_double_grgrid_address

    integer(c_long) :: address(3), D_diag(3), PS(3)
    integer(c_long) :: address_double(3)

    D_diag(:) = 4
    address(:) = 1
    PS(:) = 1

    call ksp_get_double_grgrid_address(address_double, address, D_diag, PS)

    print '("[ksp_get_double_grid_address]")'
    print '("D_diag", 3i3)', D_diag(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Half-shift", 3i3)', PS(:)
    print '("Double-grid address (address * 2 + PS)", 3i3)', address_double(:)
    print *, ""
  end subroutine get_double_grgrid_address


  subroutine get_grgrid_index()
    use kspclib_f08, only: ksp_get_grgrid_index

    integer(c_long) :: address(3), D_diag(3)
    integer(c_long) :: grid_index

    D_diag(:) = 4
    address(:) = 3
    grid_index = ksp_get_grgrid_index(address, D_diag)

    print '("[ksp_get_grgrid_index]")'
    print '("(Note that the index starts with 0.)")'
    print '("D_diag", 3i3)', D_diag(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Grid index ", i0)', grid_index
    print *, ""
  end subroutine get_grgrid_index


  subroutine get_double_grgrid_index()
    use kspclib_f08, only: ksp_get_double_grgrid_index, &
         ksp_get_double_grgrid_address

    integer(c_long) :: address(3), address_double(3), D_diag(3), PS(3)
    integer(c_long) :: grid_index

    D_diag(:) = 4
    address(:) = 1
    PS(:) = 1

    call ksp_get_double_grgrid_address(address_double, address, D_diag, PS)
    grid_index = ksp_get_double_grgrid_index(address_double, D_diag, PS)

    print '("[ksp_get_double_grgrid_index]")'
    print '("(Note that the index starts with 0.)")'
    print '("D_diag", 3i3)', D_diag(:)
    print '("Single-grid address", 3i3)', address(:)
    print '("Double-grid address", 3i3)', address_double(:)
    print '("PS", 3i3)', PS(:)
    print '("Grid index ", i0)', grid_index
    print *, ""
  end subroutine get_double_grgrid_index


  subroutine get_grgrid_address_from_index
    use kspclib_f08, only: ksp_get_grgrid_address_from_index

    integer(c_long) :: address(3), D_diag(3)
    integer(c_long) :: grid_index

    grid_index = 21
    D_diag(:) = 4
    call ksp_get_grgrid_address_from_index(address, grid_index, D_diag)

    print '("[ksp_get_grgrid_address_from_index]")'
    print '("(Note that the index starts with 0.)")'
    print '("D_diag", 3i3)', D_diag(:)
    print '("Grid index ", i0)', grid_index
    print '("Single-grid address", 3i3)', address(:)
    print *, ""
  end subroutine get_grgrid_address_from_index


  subroutine rotate_grgrid_index
    use kspclib_f08, only: ksp_rotate_grgrid_index, &
         ksp_get_grgrid_address_from_index, ksp_get_double_grgrid_address

    integer(c_long) :: D_diag(3), PS(3)
    integer(c_long) :: rotation(3, 3)
    integer(c_long) :: grid_index, gp_rot
    integer(c_long) :: P(3, 3)
    integer(c_long) :: shift(3)
    integer(c_long) :: address(3), address_rot(3)
    integer(c_long) :: address_double(3), address_double_rot(3), address_double_matmul(3)
    integer :: i

    shift(:) = 1
    D_diag(:) = [1, 5, 20]
    P(:, :) = reshape([0, 1, -2, 1, 0, 0, -2, 2, -5], [3, 3])
    rotation(:, :) = reshape([4, -5, -2, -5, 4, 2, 20, -20, -9], [3, 3])
    PS(:) = matmul(transpose(P), shift)
    grid_index = 21
    gp_rot = ksp_rotate_grgrid_index(grid_index, rotation, D_diag, PS)

    call ksp_get_grgrid_address_from_index(address, grid_index, D_diag)
    call ksp_get_double_grgrid_address(address_double, address, D_diag, PS)
    call ksp_get_grgrid_address_from_index(address_rot, gp_rot, D_diag)
    call ksp_get_double_grgrid_address(address_double_rot, address_rot, &
         D_diag, PS)

    print '("[ksp_rotate_grgrid_index]")'
    print '("(Note that the index starts with 0.)")'
    print '("D_diag", 3i3)', D_diag
    print '("P", 9i5)', P
    print '("PS", 3i3)', PS
    print '("rotation", 9i5)', rotation
    print '("Grid index ", i0)', grid_index
    print '("Grid index after rotation ", i0)', gp_rot
    print '("Double-grid address", 3i3)', address_double
    print '("Double-grid address after rotation", 3i3)', address_double_rot
    address_double_matmul = matmul(transpose(rotation), address_double)
    do i = 1, 3
       address_double_matmul(i) = modulo(address_double_matmul(i), D_diag(i) * 2)
    end do
    print '("Double-grid address after rotation (matmul)", 3i3)', address_double_matmul
    print *, ""
  end subroutine rotate_grgrid_index


  subroutine get_ir_grgrid_map
    use kspclib_f08, only: ksp_get_ir_grgrid_map

    integer(c_long), allocatable :: ir_grid_indices(:), ir_grid_indices_ref(:)
    integer(c_long) :: rotations(3, 3, 16)
    integer(c_long) :: num_rot
    integer(c_long) :: D_diag(3), PS(3), shift(3)
    integer(c_long) :: P(3, 3)
    integer :: i

    rotations(:, :, :) = reshape([ &
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
    shift(:) = 1
    P(:, :) = reshape([0, 1, -2, 1, 0, 0, -2, 2, -5], [3, 3])
    PS(:) = matmul(transpose(P), shift)
    num_rot = 16

    allocate(ir_grid_indices(product(D_diag)))
    call ksp_get_ir_grgrid_map(ir_grid_indices, rotations, num_rot, D_diag, PS)

    print '("[ksp_get_ir_grgrid_map]")'
    print '("(Note that the index starts with 0.)")'
    print '("D_diag", 3i3)', D_diag
    print '("P", 9i5)', P
    print '("PS", 3i3)', PS

    print '("Grid index mapping to ir-grid indices")'
    do i = 1, 5
       print '(20i4)', ir_grid_indices((i - 1) * 20 + 1 : i * 20)
    end do

    allocate(ir_grid_indices_ref(product(D_diag)))
    ir_grid_indices_ref(:) = [ &
         0, 1, 2, 1, 0, 5, 5, 7, 8, 7, 5, 11, 11, 5, 0, 0, 5, 11, 11, 5, &
         7, 8, 7, 5, 5, 0, 1, 2, 1, 0, 30, 30, 7, 1, 7, 5, 11, 11, 5, 0, &
         40, 30, 11, 11, 30, 7, 8, 7, 5, 5, 40, 8, 2, 8, 40, 30, 30, 7, 1, 7, &
         30, 11, 11, 30, 40, 40, 30, 11, 11, 30, 7, 1, 7, 30, 30, 40, 8, 2, 8, 40, &
         5, 5, 7, 8, 7, 30, 11, 11, 30, 40, 0, 5, 11, 11, 5, 7, 1, 7, 30, 30]

    print '("Grid index mapping to ir-grid indices (ref)")'
    do i = 1, 5
       print '(20i4)', ir_grid_indices_ref((i - 1) * 20 + 1 : i * 20)
    end do

    print *, ""

    deallocate(ir_grid_indices)
    deallocate(ir_grid_indices_ref)
  end subroutine get_ir_grgrid_map


  subroutine get_reciprocal_point_group()
    use kspclib_f08, only: ksp_get_reciprocal_point_group

    integer(c_long) :: rotations_tipn3(3, 3, 4)
    integer(c_long) :: rotations_tio2(3, 3, 16)
    integer(c_long) :: rec_rotations(3, 3, 48)
    integer(c_long) :: num_rot
    integer(c_long) :: is_time_reversal = 1
    integer :: num_rot_ret, i

    rotations_tipn3(:, :, :) = reshape( &
         [1, 0, 0, 0, 1, 0, 0, 0, 1, &
         -1, 0, 0, 0, 0, 1, 0, 1, 0, &
         -1, 0, 0, 0, 1, 0, 0, 0, 1, &
         1, 0, 0, 0, 0, 1, 0, 1, 0], [3, 3, 4])
    num_rot = 4
    num_rot_ret = ksp_get_reciprocal_point_group(rec_rotations, &
         rotations_tipn3, num_rot, is_time_reversal)

    print '("[ksp_get_reciprocal_point_group]")'
    print *, ""
    print '("TiPN3 rotations in direct space")'
    do i = 1, num_rot
       print '(9i3)', rotations_tipn3(:, :, i)
    end do
    print '("TiPN3 rotations in reciprocal space (is_time_reversal=1)")'
    do i = 1, num_rot_ret
       print '(9i3)', rec_rotations(:, :, i)
    end do
    print *, ""

    rotations_tio2(:, :, :) = reshape( &
         [1, 0, 0, 0, 1, 0, 0, 0, 1, &
         0, 1, 0, 0, 1, -1, -1, 1, 0, &
         0, 1, -1, 1, 0, -1, 0, 0, -1, &
         1, 0, -1, 1, 0, 0, 1, -1, 0, &
         -1, 0, 0, -1, 0, 1, -1, 1, 0, &
         0, -1, 0, -1, 0, 0, 0, 0, -1, &
         0, -1, 1, 0, -1, 0, 1, -1, 0, &
         -1, 0, 1, 0, -1, 1, 0, 0, 1, &
         -1, 0, 0, 0, -1, 0, 0, 0, -1, &
         0, -1, 0, 0, -1, 1, 1, -1, 0, &
         0, -1, 1, -1, 0, 1, 0, 0, 1, &
         -1, 0, 1, -1, 0, 0, -1, 1, 0, &
         1, 0, 0, 1, 0, -1, 1, -1, 0, &
         0, 1, 0, 1, 0, 0, 0, 0, 1, &
         0, 1, -1, 0, 1, 0, -1, 1, 0, &
         1, 0, -1, 0, 1, -1, 0, 0, -1], [3, 3, 16])
    num_rot = 16
    num_rot_ret = ksp_get_reciprocal_point_group(rec_rotations, &
         rotations_tio2, num_rot, is_time_reversal)

    print '("TiO2 (anatase) rotations in direct space")'
    do i = 1, num_rot
       print '(9i3)', rotations_tio2(:, :, i)
    end do
    print '("TiO2 (anatase) rotations in reciprocal space (is_time_reversal=1)")'
    do i = 1, num_rot_ret
       print '(9i3)', rec_rotations(:, :, i)
    end do
    print *, ""
  end subroutine get_reciprocal_point_group

end module kspclib_example


program kspclib_example_f08

  use kspclib_example, only: kspclib_version, get_all_grid_addresses, &
       get_double_grid_address, get_double_grid_index, &
       get_thm_relative_grid_addresses, get_thm_integration_weight, &
       snf3x3, snf_transform_rotations, get_all_grgrid_addresses, &
       get_double_grgrid_address, get_grgrid_index, get_double_grgrid_index, &
       get_grgrid_address_from_index, rotate_grgrid_index, get_ir_grgrid_map, &
       get_reciprocal_point_group

  call kspclib_version()
  call get_all_grid_addresses()
  call get_double_grid_address()
  call get_double_grid_index()
  call get_thm_relative_grid_addresses()
  call get_thm_integration_weight()
  call snf3x3()
  call snf_transform_rotations()
  call get_all_grgrid_addresses()
  call get_double_grgrid_address()
  call get_grgrid_index()
  call get_double_grgrid_index()
  call get_grgrid_address_from_index()
  call rotate_grgrid_index()
  call get_ir_grgrid_map()
  call get_reciprocal_point_group()

end program kspclib_example_f08
