program noCaps_to_withCaps
!==============================================================================
!   Smooth Voronoi CAP integral generator + APO vector
!
!   Reference: Sommerfeld & Ehara, J. Chem. Theory Comput. 11, 4627 (2015),
!   "Complex Absorbing Potentials with Voronoi Isosurfaces Wrapping Perfectly
!   around Molecules".  Smooth Voronoi-CAP (eqs. 6-8):
!       r_WA(r) = sqrt( sum_a w_a r_a^2 / sum_a w_a )
!       w_a     = 1 / (r_a^2 - r_nearest^2 + 1 a0^2)^2
!       W(r)    = 0                  if r_WA <= r_cut
!       W(r)    = (r_WA - r_cut)^2   if r_WA  > r_cut
!
!   Grid strategy (Sommerfeld & Ehara, sec. 2.1, last paragraph):
!     For every primitive-pair integral, place a single large grid centered
!     at the Gaussian product center P = (a*A + b*B)/(a+b), using:
!       - Lebedev angular grid (paper: 3470 points; default here: 1454)
!       - Equidistant radial grid (paper: 500 points; default here: 500)
!     The radial quadrature is the leading (trapezoidal) term of the
!     Euler-Maclaurin formula; the Bernoulli-number corrections involving
!     derivatives at the endpoints are dropped, as is standard for QC grids.
!
!   The CAP magnitude prefactor -i*eta is applied downstream (e.g. by
!   SNOWMAN); this generator produces eta-independent integrals.
!
!   Authors: L. M. Thompson, A. M. Kinyua, & P. D. Simpson, Louisville KY, 2023-2026
!==============================================================================

    use mqc_gaussian
    use iso_fortran_env, only: int64, real64
    use omp_lib

    implicit none

    !---------------------------------------------------------------------------
    ! MQC Types (required for Gaussian matrix file interface)
    !---------------------------------------------------------------------------
    type(mqc_gaussian_unformatted_matrix_file) :: temp_file
    type(mqc_molecule_data) :: moleculeInfo
    type(mqc_wavefunction) :: common_wave
    type(mqc_scf_integral) :: overlap, density, core_ham
    type(mqc_matrix) :: sh2AtMp, shlTyp, nPrmSh, prmExp, conCoef, conCoTwo, shCoor
    type(mqc_matrix) :: numerical_overlap, numerical_overlap_cap
    type(mqc_vector) :: coords, nPrmSh_list, prmExp_list, conCoef_list
    type(mqc_vector) :: shell_list, shell2AtomMap_list, conCoefP_list
    type(mqc_vector) :: apo_vector  ! APO vector W_p for each orbital
    type(mqc_scalar) :: numerical_trace, numerical_cap_trace, analytical_trace

    !---------------------------------------------------------------------------
    ! Command line and file handling
    !---------------------------------------------------------------------------
    character(len=:), allocatable :: command, fileName, sphereLebedev
    character(len=256), dimension(:), allocatable :: fileList
    integer(kind=int64) :: iOut = 6, numFile = 0, stat_num, iUnit, io_stat_number

    !---------------------------------------------------------------------------
    ! System parameters.  Defaults are chosen to match the Sommerfeld-Ehara
    ! JCTC 2015 reference values (3470 angular, 500 radial points) where
    ! applicable, with a generous radial extent so the CAP tail is captured
    ! out to ~20 a0 regardless of basis-set diffuseness.
    !---------------------------------------------------------------------------
    integer(kind=int64) :: nBasis, nAtoms, nShells
    integer(kind=int64) :: n_lebedev = 1454_int64
    integer(kind=int64) :: n_radial  = 500_int64
    real(kind=real64) :: r_min = 0.0_real64
    real(kind=real64) :: r_max = 20.0_real64
    real(kind=real64) :: r_cut = 3.0_real64
    real(kind=real64) :: cutoff_threshold = 1.0e-12_real64
    logical :: compute_validation = .false.

    !---------------------------------------------------------------------------
    ! Loop indices and counters
    !---------------------------------------------------------------------------
    integer(kind=int64) :: i, j, k, v, l, u, q, m, n, b, p
    integer(kind=int64) :: nShellEndI, nShellEndJ, shell_pair_idx, n_shell_pairs
    integer(kind=int64) :: total_grid_pts, tid, max_threads
    integer(kind=int64) :: max_n_prim

    !---------------------------------------------------------------------------
    ! Lebedev grid arrays
    !---------------------------------------------------------------------------
    real(kind=real64), dimension(:), allocatable :: xx, yy, zz, wweight  ! Raw input
    real(kind=real64), dimension(:), allocatable :: leb_x, leb_y, leb_z, leb_w  ! Preprocessed
    real(kind=real64) :: theta_deg, phi_deg, phi_rad, theta_rad

    !---------------------------------------------------------------------------
    ! Coordinate caches (for fast access in parallel loops)
    !---------------------------------------------------------------------------
    real(kind=real64), dimension(:,:), allocatable :: atom_coords_cache   ! (3, nAtoms)
    real(kind=real64), dimension(:,:), allocatable :: shell_coords_cache  ! (3, nShells)
    real(kind=real64), dimension(:), allocatable :: shell_type_cache, nPrmSh_cache
    integer(kind=int64), dimension(:), allocatable :: shellOffsets, basis_offsets
    integer(kind=int64), dimension(:,:), allocatable :: shell_pairs

    !---------------------------------------------------------------------------
    ! Thread-local result accumulators
    !---------------------------------------------------------------------------
    real(kind=real64), dimension(:,:,:), allocatable :: overlap_local, cap_local
    real(kind=real64), dimension(:,:), allocatable :: apo_local  ! Thread-local APO accumulator

    !---------------------------------------------------------------------------
    ! Calculation variables
    !---------------------------------------------------------------------------
    integer(kind=int64), parameter :: max_shell_size = 28  ! i cartesian
    real(kind=real64) :: x_r, y_r, z_r, x_i, y_i, z_i, x_j, y_j, z_j
    real(kind=real64) :: px, py, pz, R_shell, r_1, r_2
    real(kind=real64) :: mPrExp_val, nPrExp_val, normPA
    real(kind=real64) :: exp_arg, w_full, W_r, X_1
    real(kind=real64) :: dx, dy, dz, d_sq, min_d_sq, r_WA_distance
    real(kind=real64) :: apo_contrib, prod, gto_sum
    integer(kind=int64) :: l_val_i, l_val_j
    character(len=7) :: primGaussType_i, primGaussType_j

    ! Radial grid variables (equidistant, trapezoidal weights)
    real(kind=real64) :: radial_step, r_val, r_weight, trap_weight
    integer(kind=int64) :: radial_idx, angular_idx

    ! Per-grid-point smooth Voronoi temporaries (thread-private, squared distances)
    real(kind=real64), dimension(:), allocatable :: dsq_atoms_tp, weights_atoms_tp

    ! Per-shell-pair scratch arrays (thread-private, sized to max_shell_size)
    character(len=7), dimension(:), allocatable :: prim_types_i_tp, prim_types_j_tp
    integer(kind=int64), dimension(:), allocatable :: l_vals_i_tp, l_vals_j_tp
    real(kind=real64), dimension(:), allocatable :: mconCoef_arr_tp, nconCoef_arr_tp
    real(kind=real64), dimension(:), allocatable :: normPA_arr_tp, normPB_arr_tp
    real(kind=real64), dimension(:), allocatable :: X1_arr_tp, X2_arr_tp
    real(kind=real64), dimension(:), allocatable :: apo_block_tp
    real(kind=real64), dimension(:,:), allocatable :: gv1_pair_tp, cap_pair_tp
    real(kind=real64), dimension(:,:), allocatable :: gv1_uq_tp, cap_uq_tp
    ! APO precomputed (per shell): contraction coefficients and norms indexed by (u, m)
    real(kind=real64), dimension(:), allocatable :: mPrExp_arr_tp
    real(kind=real64), dimension(:,:), allocatable :: mconCoef_um_tp, normPA_um_tp

    !---------------------------------------------------------------------------
    ! Timing
    !---------------------------------------------------------------------------
    real(kind=real64) :: t_start, t_end, t_apo_start, t_apo_end

!==============================================================================
! PROGRAM START
!==============================================================================

    write(iOut,'(*(A))') NEW_LINE('a'), ' ', repeat('*',73), NEW_LINE('a'), &
        '  Smooth Voronoi CAP integral generator', NEW_LINE('a'), &
        '  trapezoidal radial + Lebedev angular (per-integral grid)', NEW_LINE('a'), &
        '  with APO vector W_p output', NEW_LINE('a'), &
        ' ', repeat('*',73), NEW_LINE('a'), NEW_LINE('a'), &
        repeat(' ',30), 'Version 26.05-SmoothVoronoi', NEW_LINE('a'), NEW_LINE('a'), &
        'Thompson, Kinyua, & Simpson, Louisville KY, 2023-2026.', NEW_LINE('a'), &
        'Reference: Sommerfeld & Ehara, JCTC 11, 4627 (2015).', NEW_LINE('a')

    !---------------------------------------------------------------------------
    ! Parse command line arguments
    !---------------------------------------------------------------------------
    j = 1
    do i = 1, command_argument_count()
        if (i /= j) cycle
        call mqc_get_command_argument(i, command)
        
        select case(command)
        case('-f')
            call mqc_get_command_argument(i+1, fileName)
            j = i + 2
        case('--lebedev_grid')
            call mqc_get_command_argument(i+1, sphereLebedev)
            j = i + 2
        case('--r_min')
            call mqc_get_command_argument(i+1, command)
            read(command, '(F12.6)') r_min
            j = i + 2
        case('--r_max')
            call mqc_get_command_argument(i+1, command)
            read(command, '(F12.6)') r_max
            j = i + 2
        case('--n_lebedev')
            call mqc_get_command_argument(i+1, command)
            read(command, '(I20)') n_lebedev
            j = i + 2
        case('--n_radial')
            call mqc_get_command_argument(i+1, command)
            read(command, '(I20)') n_radial
            j = i + 2
        case('--r_cut')
            call mqc_get_command_argument(i+1, command)
            read(command, '(F12.6)') r_cut
            j = i + 2
        case('--validate')
            compute_validation = .true.
            j = i + 1
        case default
            call mqc_error_a('Unknown command-line argument', 6, &
                'argument', trim(command))
        end select
    end do

    write(6,'(A)') 'Using Smooth Voronoi CAP (Sommerfeld & Ehara, JCTC 2015) + APO Vector'

    !---------------------------------------------------------------------------
    ! Validate command-line parameters.  Each value has a default; we only
    ! sanity-check that they are physically reasonable.  Note: the strength
    ! parameter eta is intentionally not handled here -- the integrals
    ! generated by this program are eta-independent (W(r) = (r_WA - r_cut)^2);
    ! the caller (SNOWMAN) multiplies by -i*eta when assembling H - i*eta*W.
    !---------------------------------------------------------------------------
    if (r_min < 0.0_real64) &
        call mqc_error('r_min must be >= 0', 6)
    if (r_max <= r_min) &
        call mqc_error('r_max must be > r_min', 6)
    if (r_cut < 0.0_real64) &
        call mqc_error('r_cut must be >= 0', 6)
    if (n_radial < 2_int64) &
        call mqc_error('n_radial must be >= 2', 6)
    if (n_lebedev < 1_int64) &
        call mqc_error('n_lebedev must be >= 1', 6)
    write(6,'(A)') 'Grid parameters (after defaults + validation):'
    write(6,'(A,F10.4)') '  r_min      (a0) = ', r_min
    write(6,'(A,F10.4)') '  r_max      (a0) = ', r_max
    write(6,'(A,F10.4)') '  r_cut      (a0) = ', r_cut
    write(6,'(A,I10)')   '  n_radial        = ', n_radial
    write(6,'(A,I10)')   '  n_lebedev       = ', n_lebedev
    if (compute_validation) &
        write(6,'(A)') 'Validation mode: numerical-vs-analytical overlap + symmetry checks will be reported'

    !---------------------------------------------------------------------------
    ! Read input file list
    !---------------------------------------------------------------------------
    if (.not. allocated(fileName)) call mqc_error('No input file provided', iOut)
    
    open(newunit=iUnit, file=fileName, status='old', iostat=stat_num)
    if (stat_num /= 0) call mqc_error_a('Error opening file', iOut, 'fileName', fileName)
    
    read(unit=iUnit, fmt='(i20)', iostat=stat_num) numFile
    if (stat_num /= 0) call mqc_error('Error reading file number', iOut)
    
    allocate(fileList(numFile))
    do i = 1, numFile
        read(unit=iUnit, fmt='(A)', iostat=stat_num) fileList(i)
        if (stat_num < 0 .and. i <= numFile) call mqc_error('File EOF reached early', iOut)
    end do
    close(unit=iUnit)

    !---------------------------------------------------------------------------
    ! Read Gaussian matrix file data
    !---------------------------------------------------------------------------
    call temp_file%getArray('SHELL TO ATOM MAP', sh2AtMp, fileName=fileList(1))
    call temp_file%getArray('SHELL TYPES', shlTyp, fileName=fileList(1))
    call temp_file%getArray('NUMBER OF PRIMITIVES PER SHELL', nPrmSh, fileName=fileList(1))
    call temp_file%getArray('PRIMITIVE EXPONENTS', prmExp, fileName=fileList(1))
    call temp_file%getArray('CONTRACTION COEFFICIENTS', conCoef, fileName=fileList(1))
    call temp_file%getArray('P(S=P) CONTRACTION COEFFICIENTS', conCoTwo, fileName=fileList(1))
    call temp_file%getArray('COORDINATES OF EACH SHELL', shCoor, fileName=fileList(1))

    call temp_file%getESTObj('overlap', est_integral=overlap, fileName=fileList(1))
    call temp_file%getESTObj('wavefunction', common_wave, fileName=fileList(1))
    call temp_file%getESTObj('density', est_integral=density, fileName=fileList(1))
    call temp_file%getESTObj('core hamiltonian', est_integral=core_ham)
    call temp_file%getMolData(moleculeInfo)

    !---------------------------------------------------------------------------
    ! Extract system parameters
    !---------------------------------------------------------------------------
    nBasis = common_wave%nbasis%ival()
    nAtoms = MQC_Scalar_Get_Intrinsic_Integer(moleculeInfo%getNumAtoms())
    
    coords = shCoor%diagonal()
    nPrmSh_list = nPrmSh%diagonal()
    prmExp_list = prmExp%diagonal()
    conCoef_list = conCoef%diagonal()
    shell_list = shlTyp%diagonal()
    shell2AtomMap_list = sh2AtMp%diagonal()
    conCoefP_list = conCoTwo%diagonal()
    
    nShells = size(nPrmSh_list)

    call mqc_print(nBasis, 6, 'nBasis')
    call sh2AtMp%print(iOut, 'Shell to Atom Map', Blank_At_Bottom=.true.)
    call shlTyp%print(iOut, 'Shell Type', Blank_At_Bottom=.true.)
    call moleculeInfo%print(iOut)

    !---------------------------------------------------------------------------
    ! Build shell offsets for primitive indexing
    !---------------------------------------------------------------------------
    allocate(shellOffsets(nShells))
    shellOffsets(1) = 0
    do i = 2, nShells
        shellOffsets(i) = shellOffsets(i-1) + int(nPrmSh_list%at(i-1))
    end do

    !---------------------------------------------------------------------------
    ! Build basis function offsets
    !---------------------------------------------------------------------------
    allocate(basis_offsets(nShells + 1))
    basis_offsets(1) = 0
    do v = 1, nShells
        basis_offsets(v+1) = basis_offsets(v) + get_shell_size(int(int(shell_list%at(v)), int64))
    end do

    !---------------------------------------------------------------------------
    ! Read Lebedev grid
    !---------------------------------------------------------------------------
    if (.not. allocated(sphereLebedev)) call mqc_error('No Lebedev file provided', 6)
    
    open(newunit=iUnit, file=sphereLebedev, status='old', iostat=io_stat_number)
    if (io_stat_number /= 0) call mqc_error('Error opening Lebedev file', 6)

    allocate(xx(n_lebedev), yy(n_lebedev), zz(n_lebedev), wweight(n_lebedev))
    
    do k = 1, n_lebedev
        read(unit=iUnit, fmt=*, iostat=io_stat_number) theta_deg, phi_deg, wweight(k)
        if (io_stat_number /= 0) then
            write(*,*) 'Error reading Lebedev data at line ', k
            stop
        end if
        
        theta_rad = theta_deg * (pi / 180.0_real64)
        phi_rad = phi_deg * (pi / 180.0_real64)
        
        xx(k) = cos(theta_rad) * sin(phi_rad)
        yy(k) = sin(theta_rad) * sin(phi_rad)
        zz(k) = cos(phi_rad)
    end do
    close(iUnit)

    !---------------------------------------------------------------------------
    ! Precompute Lebedev grid with 4*pi factor
    !---------------------------------------------------------------------------
    allocate(leb_x(n_lebedev), leb_y(n_lebedev), leb_z(n_lebedev), leb_w(n_lebedev))
    do k = 1, n_lebedev
        leb_x(k) = xx(k)
        leb_y(k) = yy(k)
        leb_z(k) = zz(k)
        leb_w(k) = wweight(k) * 4.0_real64 * pi
    end do
    deallocate(xx, yy, zz, wweight)  ! No longer needed

    !---------------------------------------------------------------------------
    ! Cache coordinate data for fast access
    !---------------------------------------------------------------------------
    write(6,'(A)') 'Caching coordinate data...'

    allocate(atom_coords_cache(3, nAtoms))
    do i = 1, nAtoms
        atom_coords_cache(1, i) = moleculeInfo%Cartesian_Coordinates%at(1, i)
        atom_coords_cache(2, i) = moleculeInfo%Cartesian_Coordinates%at(2, i)
        atom_coords_cache(3, i) = moleculeInfo%Cartesian_Coordinates%at(3, i)
    end do

    allocate(shell_coords_cache(3, nShells))
    do v = 1, nShells
        shell_coords_cache(1, v) = coords%at((v-1)*3 + 1)
        shell_coords_cache(2, v) = coords%at((v-1)*3 + 2)
        shell_coords_cache(3, v) = coords%at((v-1)*3 + 3)
    end do

    allocate(shell_type_cache(nShells), nPrmSh_cache(nShells))
    do v = 1, nShells
        shell_type_cache(v) = shell_list%at(v)
        nPrmSh_cache(v) = nPrmSh_list%at(v)
    end do

    !---------------------------------------------------------------------------
    ! Print grid parameters
    !---------------------------------------------------------------------------
    write(6,'(A)') ''
    write(6,'(A)') 'Grid summary:'
    total_grid_pts = n_radial * n_lebedev
    write(6,'(A,I0)')    '  Grid points per integral:   ', total_grid_pts

    !---------------------------------------------------------------------------
    ! Equidistant radial step.  Trapezoidal weights below (leading EM term).
    !---------------------------------------------------------------------------
    radial_step = (r_max - r_min) / dble(n_radial - 1)

    !---------------------------------------------------------------------------
    ! Build shell pair list (upper triangle only for symmetry)
    !---------------------------------------------------------------------------
    n_shell_pairs = (nShells * (nShells + 1)) / 2
    allocate(shell_pairs(2, n_shell_pairs))
    
    shell_pair_idx = 0
    do v = 1, nShells
        do l = v, nShells
            shell_pair_idx = shell_pair_idx + 1
            shell_pairs(1, shell_pair_idx) = v
            shell_pairs(2, shell_pair_idx) = l
        end do
    end do
    write(6,'(A,I0,A)') 'Processing ', n_shell_pairs, ' unique shell pairs'

    !---------------------------------------------------------------------------
    ! Initialize output matrices and APO vector
    !---------------------------------------------------------------------------
    call numerical_overlap%init(int(nBasis), int(nBasis))
    call numerical_overlap_cap%init(int(nBasis), int(nBasis))
    call apo_vector%init(int(nBasis))

    !---------------------------------------------------------------------------
    ! Get thread count and allocate thread-local storage
    !---------------------------------------------------------------------------
    max_threads = 1
    !$OMP PARALLEL
        !$OMP MASTER
            max_threads = omp_get_num_threads()
            write(6,'(A,I0,A)') 'Running with ', max_threads, ' OpenMP threads'
        !$OMP END MASTER
    !$OMP END PARALLEL

    allocate(overlap_local(nBasis, nBasis, max_threads))
    allocate(cap_local(nBasis, nBasis, max_threads))
    allocate(apo_local(nBasis, max_threads))
    overlap_local = 0.0_real64
    cap_local = 0.0_real64
    apo_local = 0.0_real64

!==============================================================================
! MAIN PARALLEL COMPUTATION - CAP MATRIX
!
! Grid strategy (Sommerfeld & Ehara, JCTC 2015, sec. 2.1):
!   For each primitive pair (m in shell v, n in shell l):
!     1. Compute Gaussian product center P = (a*A + b*B)/(a+b)
!     2. Place a single large grid at P:
!        - Lebedev angular grid (n_lebedev points)
!        - Equidistant radial grid (n_radial points, trapezoidal weights)
!     3. At each grid point: evaluate X_1[u] for all u in shell v, X_2[q] for
!        all q in shell l, and W(r) once. Accumulate gv1[u,q] and cap[u,q].
!     4. After all primitives, fold in normalization + contraction.
!
!   Loop order is (shell-pair) -> (m,n primitives) -> grid -> (u,q basis).
!   This computes W(r), exp(-alpha r^2), and the Voronoi geometry once per
!   grid point per primitive pair instead of once per (u,q) basis pair.
!==============================================================================

    t_start = omp_get_wtime()

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(shell_pairs, n_shell_pairs, shell_type_cache, nPrmSh_cache, &
    !$OMP        shellOffsets, shell_coords_cache, atom_coords_cache, &
    !$OMP        prmExp_list, conCoef_list, conCoefP_list, &
    !$OMP        leb_x, leb_y, leb_z, leb_w, &
    !$OMP        n_radial, n_lebedev, r_min, r_max, r_cut, radial_step, &
    !$OMP        nAtoms, nBasis, basis_offsets, &
    !$OMP        overlap_local, cap_local, cutoff_threshold) &
    !$OMP PRIVATE(shell_pair_idx, v, l, u, q, m, n, i, j, k, b, tid, &
    !$OMP         nShellEndI, nShellEndJ, primGaussType_i, primGaussType_j, &
    !$OMP         R_shell, prod, w_full, &
    !$OMP         x_r, y_r, z_r, x_i, y_i, z_i, x_j, y_j, z_j, &
    !$OMP         px, py, pz, r_1, r_2, r_val, r_weight, trap_weight, &
    !$OMP         radial_idx, angular_idx, &
    !$OMP         dx, dy, dz, d_sq, min_d_sq, r_WA_distance, W_r, &
    !$OMP         dsq_atoms_tp, weights_atoms_tp, &
    !$OMP         prim_types_i_tp, prim_types_j_tp, l_vals_i_tp, l_vals_j_tp, &
    !$OMP         mconCoef_arr_tp, nconCoef_arr_tp, normPA_arr_tp, normPB_arr_tp, &
    !$OMP         X1_arr_tp, X2_arr_tp, gv1_pair_tp, cap_pair_tp, gv1_uq_tp, cap_uq_tp, &
    !$OMP         exp_arg, mPrExp_val, nPrExp_val, l_val_i, l_val_j)

    tid = omp_get_thread_num() + 1

    ! Allocate thread-private scratch
    allocate(dsq_atoms_tp(nAtoms), weights_atoms_tp(nAtoms))
    allocate(prim_types_i_tp(max_shell_size), prim_types_j_tp(max_shell_size))
    allocate(l_vals_i_tp(max_shell_size), l_vals_j_tp(max_shell_size))
    allocate(mconCoef_arr_tp(max_shell_size), nconCoef_arr_tp(max_shell_size))
    allocate(normPA_arr_tp(max_shell_size), normPB_arr_tp(max_shell_size))
    allocate(X1_arr_tp(max_shell_size), X2_arr_tp(max_shell_size))
    allocate(gv1_pair_tp(max_shell_size, max_shell_size))
    allocate(cap_pair_tp(max_shell_size, max_shell_size))
    allocate(gv1_uq_tp(max_shell_size, max_shell_size))
    allocate(cap_uq_tp(max_shell_size, max_shell_size))

    !$OMP DO SCHEDULE(DYNAMIC, 1)
    do shell_pair_idx = 1, n_shell_pairs
        v = shell_pairs(1, shell_pair_idx)
        l = shell_pairs(2, shell_pair_idx)

        nShellEndI = get_shell_size(nint(shell_type_cache(v), int64))
        nShellEndJ = get_shell_size(nint(shell_type_cache(l), int64))

        ! Shell coordinates and inter-shell distance
        x_i = shell_coords_cache(1, v)
        y_i = shell_coords_cache(2, v)
        z_i = shell_coords_cache(3, v)
        x_j = shell_coords_cache(1, l)
        y_j = shell_coords_cache(2, l)
        z_j = shell_coords_cache(3, l)
        R_shell = sqrt((x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2)

        ! Per-shell-pair: precompute primitive labels and l-values per basis function
        do u = 1, nShellEndI
            prim_types_i_tp(u) = get_primitive_type(nint(shell_type_cache(v), int64), int(u, int64))
            l_vals_i_tp(u)     = get_l_value(prim_types_i_tp(u))
        end do
        do q = 1, nShellEndJ
            prim_types_j_tp(q) = get_primitive_type(nint(shell_type_cache(l), int64), int(q, int64))
            l_vals_j_tp(q)     = get_l_value(prim_types_j_tp(q))
        end do

        ! Reset shell-pair accumulators (already-normalized contributions)
        gv1_uq_tp(1:nShellEndI, 1:nShellEndJ) = 0.0_real64
        cap_uq_tp(1:nShellEndI, 1:nShellEndJ) = 0.0_real64

        ! Loop over primitives in shell v
        do m = 1, int(nPrmSh_cache(v))
            mPrExp_val = prmExp_list%at(m + shellOffsets(v))
            do u = 1, nShellEndI
                if (index(prim_types_i_tp(u), 'sp') > 0) then
                    mconCoef_arr_tp(u) = conCoefP_list%at(m + shellOffsets(v))
                else
                    mconCoef_arr_tp(u) = conCoef_list%at(m + shellOffsets(v))
                end if
                normPA_arr_tp(u) = normalization_factor(mPrExp_val, prim_types_i_tp(u), l_vals_i_tp(u))
            end do

            ! Loop over primitives in shell l
            do n = 1, int(nPrmSh_cache(l))
                nPrExp_val = prmExp_list%at(n + shellOffsets(l))
                do q = 1, nShellEndJ
                    if (index(prim_types_j_tp(q), 'sp') > 0) then
                        nconCoef_arr_tp(q) = conCoefP_list%at(n + shellOffsets(l))
                    else
                        nconCoef_arr_tp(q) = conCoef_list%at(n + shellOffsets(l))
                    end if
                    normPB_arr_tp(q) = normalization_factor(nPrExp_val, prim_types_j_tp(q), l_vals_j_tp(q))
                end do

                ! Schwarz-like screening on primitive overlap exponential
                exp_arg = -(mPrExp_val * nPrExp_val) / (mPrExp_val + nPrExp_val) * R_shell * R_shell
                if (exp(exp_arg) < cutoff_threshold) cycle

                ! Gaussian product center: P = (a A + b B) / (a + b)
                call gaussian_product_center(mPrExp_val, nPrExp_val, &
                     x_i, y_i, z_i, x_j, y_j, z_j, px, py, pz)

                ! Per-(m,n) accumulators (no normalization yet)
                gv1_pair_tp(1:nShellEndI, 1:nShellEndJ) = 0.0_real64
                cap_pair_tp(1:nShellEndI, 1:nShellEndJ) = 0.0_real64

                do radial_idx = 1, n_radial
                    r_val = r_min + dble(radial_idx - 1) * radial_step
                    ! Trapezoidal weight: h/2 at endpoints, h interior (leading EM term)
                    if (radial_idx == 1 .or. radial_idx == n_radial) then
                        trap_weight = 0.5_real64 * radial_step
                    else
                        trap_weight = radial_step
                    end if
                    r_weight = trap_weight * r_val * r_val
                    if (r_weight < 1.0e-30_real64) cycle

                    do angular_idx = 1, n_lebedev
                        x_r = px + r_val * leb_x(angular_idx)
                        y_r = py + r_val * leb_y(angular_idx)
                        z_r = pz + r_val * leb_z(angular_idx)

                        ! Evaluate every basis function in shell v at this grid point
                        r_1 = (x_r - x_i)**2 + (y_r - y_i)**2 + (z_r - z_i)**2
                        do u = 1, nShellEndI
                            X1_arr_tp(u) = eval_spherical_GTO(x_r, y_r, z_r, x_i, y_i, z_i, r_1, &
                                                              prim_types_i_tp(u), mPrExp_val, l_vals_i_tp(u))
                        end do
                        ! ... and shell l
                        r_2 = (x_r - x_j)**2 + (y_r - y_j)**2 + (z_r - z_j)**2
                        do q = 1, nShellEndJ
                            X2_arr_tp(q) = eval_spherical_GTO(x_r, y_r, z_r, x_j, y_j, z_j, r_2, &
                                                              prim_types_j_tp(q), nPrExp_val, l_vals_j_tp(q))
                        end do

                        ! Smooth Voronoi W(r) at this point (one computation, shared by all u,q)
                        ! Use squared distances throughout to avoid nAtoms sqrts.
                        min_d_sq = huge(1.0_real64)
                        do b = 1, nAtoms
                            dx = x_r - atom_coords_cache(1, b)
                            dy = y_r - atom_coords_cache(2, b)
                            dz = z_r - atom_coords_cache(3, b)
                            d_sq = dx*dx + dy*dy + dz*dz
                            dsq_atoms_tp(b) = d_sq
                            if (d_sq < min_d_sq) min_d_sq = d_sq
                        end do
                        do b = 1, nAtoms
                            weights_atoms_tp(b) = 1.0_real64 / &
                                ((dsq_atoms_tp(b) - min_d_sq + 1.0_real64)**2)
                        end do
                        r_WA_distance = weighted_average_distance_sq( &
                            dsq_atoms_tp, weights_atoms_tp, nAtoms)
                        W_r = cap_weight(r_WA_distance, r_cut)

                        ! Combined quadrature weight at this point
                        w_full = leb_w(angular_idx) * r_weight

                        ! Outer-product accumulation
                        do q = 1, nShellEndJ
                            do u = 1, nShellEndI
                                prod = X1_arr_tp(u) * X2_arr_tp(q) * w_full
                                gv1_pair_tp(u, q) = gv1_pair_tp(u, q) + prod
                                cap_pair_tp(u, q) = cap_pair_tp(u, q) + prod * W_r
                            end do
                        end do
                    end do  ! angular_idx
                end do  ! radial_idx

                ! Fold in normalization * contraction for this primitive pair
                do q = 1, nShellEndJ
                    do u = 1, nShellEndI
                        prod = normPA_arr_tp(u) * normPB_arr_tp(q) * &
                               mconCoef_arr_tp(u) * nconCoef_arr_tp(q)
                        gv1_uq_tp(u, q) = gv1_uq_tp(u, q) + prod * gv1_pair_tp(u, q)
                        cap_uq_tp(u, q) = cap_uq_tp(u, q) + prod * cap_pair_tp(u, q)
                    end do
                end do

            end do  ! n (primitives in l)
        end do  ! m (primitives in v)

        ! Stash into thread-local matrix accumulators (with symmetry)
        do u = 1, nShellEndI
            i = basis_offsets(v) + u
            do q = 1, nShellEndJ
                j = basis_offsets(l) + q
                if (v == l .and. j < i) cycle

                overlap_local(i, j, tid) = overlap_local(i, j, tid) + gv1_uq_tp(u, q)
                cap_local(i, j, tid)     = cap_local(i, j, tid)     + cap_uq_tp(u, q)
                if (i /= j) then
                    overlap_local(j, i, tid) = overlap_local(j, i, tid) + gv1_uq_tp(u, q)
                    cap_local(j, i, tid)     = cap_local(j, i, tid)     + cap_uq_tp(u, q)
                end if
            end do
        end do
    end do  ! shell_pair_idx
    !$OMP END DO

    deallocate(dsq_atoms_tp, weights_atoms_tp)
    deallocate(prim_types_i_tp, prim_types_j_tp, l_vals_i_tp, l_vals_j_tp)
    deallocate(mconCoef_arr_tp, nconCoef_arr_tp, normPA_arr_tp, normPB_arr_tp)
    deallocate(X1_arr_tp, X2_arr_tp)
    deallocate(gv1_pair_tp, cap_pair_tp, gv1_uq_tp, cap_uq_tp)

    !$OMP END PARALLEL

    t_end = omp_get_wtime()
    write(6,'(A,F10.2,A)') 'CAP matrix loop completed in ', t_end - t_start, ' seconds'

!==============================================================================
! APO VECTOR COMPUTATION
! W_p = integral phi_p(r) * W(r) dV for each basis function p.
!
! For each shell, one grid is centered on the shell center; every basis
! function in the shell is evaluated and accumulated against W(r) in the
! same pass.  Contraction coefficients and normalizations are hoisted out
! of the (radial x angular) grid loop because they are grid-independent.
!==============================================================================

    write(6,'(A)') ''
    write(6,'(A)') repeat('=', 60)
    write(6,'(A)') 'Computing APO Vector (W_p for each orbital)'
    write(6,'(A)') repeat('=', 60)

    t_apo_start = omp_get_wtime()

    ! Bound for primitive count across shells (used to size thread-private
    ! (u,m) tables of contraction coefficients + GTO normalizations).
    max_n_prim = 1_int64
    do v = 1, nShells
        if (int(nPrmSh_cache(v), int64) > max_n_prim) &
            max_n_prim = int(nPrmSh_cache(v), int64)
    end do

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(nShells, shell_type_cache, nPrmSh_cache, shellOffsets, &
    !$OMP        shell_coords_cache, atom_coords_cache, &
    !$OMP        prmExp_list, conCoef_list, conCoefP_list, &
    !$OMP        leb_x, leb_y, leb_z, leb_w, &
    !$OMP        n_radial, n_lebedev, r_min, r_max, r_cut, radial_step, &
    !$OMP        nAtoms, nBasis, basis_offsets, max_n_prim, &
    !$OMP        apo_local) &
    !$OMP PRIVATE(v, u, m, p, b, tid, nShellEndI, &
    !$OMP         x_r, y_r, z_r, x_i, y_i, z_i, r_1, X_1, w_full, prod, &
    !$OMP         mPrExp_val, normPA, gto_sum, &
    !$OMP         radial_idx, angular_idx, r_val, r_weight, trap_weight, &
    !$OMP         dx, dy, dz, d_sq, min_d_sq, r_WA_distance, W_r, &
    !$OMP         dsq_atoms_tp, weights_atoms_tp, &
    !$OMP         prim_types_i_tp, l_vals_i_tp, &
    !$OMP         mconCoef_arr_tp, X1_arr_tp, apo_block_tp, &
    !$OMP         mPrExp_arr_tp, mconCoef_um_tp, normPA_um_tp)

    tid = omp_get_thread_num() + 1

    allocate(dsq_atoms_tp(nAtoms), weights_atoms_tp(nAtoms))
    allocate(prim_types_i_tp(max_shell_size), l_vals_i_tp(max_shell_size))
    allocate(mconCoef_arr_tp(max_shell_size))
    allocate(X1_arr_tp(max_shell_size), apo_block_tp(max_shell_size))
    allocate(mPrExp_arr_tp(max_n_prim))
    allocate(mconCoef_um_tp(max_shell_size, max_n_prim))
    allocate(normPA_um_tp(max_shell_size, max_n_prim))

    !$OMP DO SCHEDULE(DYNAMIC, 1)
    do v = 1, nShells
        nShellEndI = get_shell_size(nint(shell_type_cache(v), int64))

        x_i = shell_coords_cache(1, v)
        y_i = shell_coords_cache(2, v)
        z_i = shell_coords_cache(3, v)

        ! Per-shell: precompute primitive labels and l-values per basis function
        do u = 1, nShellEndI
            prim_types_i_tp(u) = get_primitive_type(nint(shell_type_cache(v), int64), int(u, int64))
            l_vals_i_tp(u)     = get_l_value(prim_types_i_tp(u))
        end do

        ! Per-shell: hoist grid-independent primitive data out of the (r,Omega)
        ! loop.  mPrExp_arr_tp(m), mconCoef_um_tp(u,m), normPA_um_tp(u,m) all
        ! depend only on (shell v, primitive m, basis fn u within the shell).
        do m = 1, int(nPrmSh_cache(v))
            mPrExp_arr_tp(m) = prmExp_list%at(m + shellOffsets(v))
            do u = 1, nShellEndI
                if (index(prim_types_i_tp(u), 'sp') > 0) then
                    mconCoef_um_tp(u, m) = conCoefP_list%at(m + shellOffsets(v))
                else
                    mconCoef_um_tp(u, m) = conCoef_list%at(m + shellOffsets(v))
                end if
                normPA_um_tp(u, m) = normalization_factor( &
                    mPrExp_arr_tp(m), prim_types_i_tp(u), l_vals_i_tp(u))
            end do
        end do

        apo_block_tp(1:nShellEndI) = 0.0_real64

        do radial_idx = 1, n_radial
            r_val = r_min + dble(radial_idx - 1) * radial_step
            if (radial_idx == 1 .or. radial_idx == n_radial) then
                trap_weight = 0.5_real64 * radial_step
            else
                trap_weight = radial_step
            end if
            r_weight = trap_weight * r_val * r_val
            if (r_weight < 1.0e-30_real64) cycle

            ! Grid is centered on the shell center and Lebedev points are unit
            ! vectors, so |r - r_shell|^2 = r_val^2 exactly.
            r_1 = r_val * r_val

            do angular_idx = 1, n_lebedev
                x_r = x_i + r_val * leb_x(angular_idx)
                y_r = y_i + r_val * leb_y(angular_idx)
                z_r = z_i + r_val * leb_z(angular_idx)

                ! Contracted phi_p[u] at this grid point, for all u in shell v
                X1_arr_tp(1:nShellEndI) = 0.0_real64
                do m = 1, int(nPrmSh_cache(v))
                    do u = 1, nShellEndI
                        X_1 = eval_spherical_GTO(x_r, y_r, z_r, x_i, y_i, z_i, r_1, &
                                                 prim_types_i_tp(u), mPrExp_arr_tp(m), &
                                                 l_vals_i_tp(u))
                        X1_arr_tp(u) = X1_arr_tp(u) + &
                            mconCoef_um_tp(u, m) * normPA_um_tp(u, m) * X_1
                    end do
                end do

                ! Smooth Voronoi W(r) once for all u
                min_d_sq = huge(1.0_real64)
                do b = 1, nAtoms
                    dx = x_r - atom_coords_cache(1, b)
                    dy = y_r - atom_coords_cache(2, b)
                    dz = z_r - atom_coords_cache(3, b)
                    d_sq = dx*dx + dy*dy + dz*dz
                    dsq_atoms_tp(b) = d_sq
                    if (d_sq < min_d_sq) min_d_sq = d_sq
                end do
                do b = 1, nAtoms
                    weights_atoms_tp(b) = 1.0_real64 / &
                        ((dsq_atoms_tp(b) - min_d_sq + 1.0_real64)**2)
                end do
                r_WA_distance = weighted_average_distance_sq( &
                    dsq_atoms_tp, weights_atoms_tp, nAtoms)
                W_r = cap_weight(r_WA_distance, r_cut)

                w_full = W_r * leb_w(angular_idx) * r_weight
                do u = 1, nShellEndI
                    apo_block_tp(u) = apo_block_tp(u) + X1_arr_tp(u) * w_full
                end do
            end do  ! angular_idx
        end do  ! radial_idx

        ! Store per-u APO contributions into thread-local accumulator
        do u = 1, nShellEndI
            p = basis_offsets(v) + u
            apo_local(p, tid) = apo_local(p, tid) + apo_block_tp(u)
        end do
    end do  ! v (shells)
    !$OMP END DO

    deallocate(dsq_atoms_tp, weights_atoms_tp)
    deallocate(prim_types_i_tp, l_vals_i_tp)
    deallocate(mconCoef_arr_tp, X1_arr_tp, apo_block_tp)
    deallocate(mPrExp_arr_tp, mconCoef_um_tp, normPA_um_tp)

    !$OMP END PARALLEL

    t_apo_end = omp_get_wtime()
    write(6,'(A,F10.2,A)') 'APO vector loop completed in ', t_apo_end - t_apo_start, ' seconds'

!==============================================================================
! REDUCE THREAD-LOCAL RESULTS
!==============================================================================

    write(6,'(A)') 'Reducing thread-local results...'
    ! Sum over threads into plain real(:,:) buffers first; the mqc_matrix
    ! accessors are heavy enough that doing nBasis*nBasis %at+%put calls
    ! per thread dominates the reduction wall time.  One %put per element,
    ! after the sum, is dramatically faster.
    block
      real(kind=real64), allocatable :: ovlp_sum(:,:), cap_sum(:,:), apo_sum(:)
      allocate(ovlp_sum(nBasis, nBasis), cap_sum(nBasis, nBasis), apo_sum(nBasis))
      ovlp_sum = 0.0_real64
      cap_sum  = 0.0_real64
      apo_sum  = 0.0_real64
      do tid = 1, max_threads
          ovlp_sum = ovlp_sum + overlap_local(:, :, tid)
          cap_sum  = cap_sum  + cap_local(:, :, tid)
          apo_sum  = apo_sum  + apo_local(:, tid)
      end do
      do j = 1, nBasis
          do i = 1, nBasis
              call numerical_overlap%put(ovlp_sum(i, j), i, j)
              call numerical_overlap_cap%put(cap_sum(i, j), i, j)
          end do
      end do
      do i = 1, nBasis
          call apo_vector%put(apo_sum(i), i)
      end do
      deallocate(ovlp_sum, cap_sum, apo_sum)
    end block

!==============================================================================
! OUTPUT RESULTS
!==============================================================================

    call overlap%print(iOut, 'Analytical Overlap', Blank_At_Bottom=.true.)
    call numerical_overlap%print(iOut, 'Numerical Overlap')
    call numerical_overlap_cap%print(iOut, 'CAP matrix')
    
    write(6,'(A)') ''
    write(6,'(A)') repeat('=', 60)
    write(6,'(A)') 'APO Vector (W_p for each basis function)'
    write(6,'(A)') repeat('=', 60)
    call apo_vector%print(iOut, 'APO Vector W_p')

    numerical_trace = numerical_overlap%trace()
    numerical_cap_trace = numerical_overlap_cap%trace()
    analytical_trace = overlap%trace()

    write(*,'(A,F20.6)') 'Numerical Trace     = ', numerical_trace%rval()
    write(*,'(A,F20.6)') 'Numerical CAP Trace = ', numerical_cap_trace%rval()
    write(*,'(A,F20.6)') 'Analytical Trace    = ', analytical_trace%rval()
    write(*,'(A,F20.6)') 'APO Vector Sum      = ', sum_vector(apo_vector)

    if (compute_validation) call validate_results(overlap, numerical_overlap, &
                                                  numerical_overlap_cap, apo_vector, &
                                                  nBasis)

    call write_matrix_file(temp_file, 'CAP_Mat', overlap, numerical_overlap, numerical_overlap_cap, apo_vector, nBasis)

!==============================================================================
! CLEANUP
!==============================================================================

    deallocate(overlap_local, cap_local, apo_local)
    deallocate(atom_coords_cache, shell_coords_cache, shell_type_cache, nPrmSh_cache)
    deallocate(leb_x, leb_y, leb_z, leb_w)
    deallocate(shell_pairs, shellOffsets, basis_offsets)

contains

!==============================================================================
! HELPER FUNCTIONS
!==============================================================================

!------------------------------------------------------------------------------
! Grid validation: numerical-vs-analytical overlap and CAP/overlap symmetry.
! Run when --validate is passed.  No physical check on the APO vector beyond
! "no NaN/Inf" since the smooth Voronoi integral has no closed form.
!------------------------------------------------------------------------------
subroutine validate_results(analytical, numerical, cap_mat, apo_vec, nb)
    type(mqc_scf_integral), intent(in) :: analytical
    type(mqc_matrix), intent(in) :: numerical, cap_mat
    type(mqc_vector), intent(in) :: apo_vec
    integer(kind=int64), intent(in) :: nb
    real(kind=real64) :: diff, max_diff_ovlp, max_asym_ovlp, max_asym_cap
    real(kind=real64) :: max_apo, val
    real(kind=real64) :: num_ij, num_ji, ana_ij, cap_ij, cap_ji, apo_i
    integer(kind=int64) :: ii, jj
    type(mqc_matrix) :: analytical_mat

    write(6,'(/,A)') repeat('=', 60)
    write(6,'(A)')   'Validation checks'
    write(6,'(A)')   repeat('=', 60)

    ! Analytical overlap, materialized as mqc_matrix for elementwise compare.
    analytical_mat = analytical%getBlock('full')

    max_diff_ovlp = 0.0_real64
    max_asym_ovlp = 0.0_real64
    max_asym_cap  = 0.0_real64
    do jj = 1, nb
        do ii = 1, nb
            num_ij = numerical%at(ii, jj)
            ana_ij = analytical_mat%at(ii, jj)
            diff = abs(num_ij - ana_ij)
            if (diff > max_diff_ovlp) max_diff_ovlp = diff
            if (ii < jj) then
                num_ji = numerical%at(jj, ii)
                diff = abs(num_ij - num_ji)
                if (diff > max_asym_ovlp) max_asym_ovlp = diff
                cap_ij = cap_mat%at(ii, jj)
                cap_ji = cap_mat%at(jj, ii)
                diff = abs(cap_ij - cap_ji)
                if (diff > max_asym_cap) max_asym_cap = diff
            end if
        end do
    end do

    max_apo = 0.0_real64
    do ii = 1, nb
        apo_i = apo_vec%at(ii)
        val = abs(apo_i)
        if (val > max_apo) max_apo = val
    end do

    write(6,'(A,ES12.4)') '  max |S_num - S_analytical|     = ', max_diff_ovlp
    write(6,'(A,ES12.4)') '  max |S_ij - S_ji| (numerical)  = ', max_asym_ovlp
    write(6,'(A,ES12.4)') '  max |CAP_ij - CAP_ji|          = ', max_asym_cap
    write(6,'(A,ES12.4)') '  max |W_p| (APO vector)         = ', max_apo

    if (max_diff_ovlp > 1.0e-5_real64) &
        write(6,'(A)') '  WARNING: numerical overlap deviates from analytical by > 1e-5 -- increase n_radial/n_lebedev?'
    if (max_asym_ovlp > 1.0e-10_real64 .or. max_asym_cap > 1.0e-10_real64) &
        write(6,'(A)') '  WARNING: symmetric-copy logic may be broken (off-diagonal asymmetry > 1e-10)'
    write(6,'(A)') ''
end subroutine validate_results

!------------------------------------------------------------------------------
! Sum vector elements (helper for output)
!------------------------------------------------------------------------------
function sum_vector(vec) result(total)
    type(mqc_vector), intent(in) :: vec
    real(kind=real64) :: total
    integer(kind=int64) :: i, n
    
    n = size(vec)
    total = 0.0_real64
    do i = 1, n
        total = total + vec%at(i)
    end do
end function sum_vector

!------------------------------------------------------------------------------
! CAP weight function W(r) - Smooth Voronoi (Eq. 9 of Ehara et al.)
!   W(r_WA) = 0                      if r_WA <= r_cut
!   W(r_WA) = (r_WA - r_cut)^2       if r_WA >  r_cut
!------------------------------------------------------------------------------
pure function cap_weight(r_WA, r_cut) result(W_r)
    real(kind=real64), intent(in) :: r_WA, r_cut
    real(kind=real64) :: W_r
    
    if (r_WA <= r_cut) then
        W_r = 0.0_real64
    else
        W_r = (r_WA - r_cut)**2
    end if
end function cap_weight

!------------------------------------------------------------------------------
! Gaussian product center P = (a*A + b*B)/(a+b)
!------------------------------------------------------------------------------
pure subroutine gaussian_product_center(a, b, xi, yi, zi, xj, yj, zj, px, py, pz)
    real(kind=real64), intent(in) :: a, b, xi, yi, zi, xj, yj, zj
    real(kind=real64), intent(out) :: px, py, pz
    real(kind=real64) :: sab
    
    sab = a + b
    px = (a * xi + b * xj) / sab
    py = (a * yi + b * yj) / sab
    pz = (a * zi + b * zj) / sab
end subroutine gaussian_product_center

!------------------------------------------------------------------------------
! Weighted average distance for smooth Voronoi (Eq. 8a of Ehara et al.)
!   r_WA = sqrt( sum_a(omega_a * r_a^2) / sum_a(omega_a) )
! where omega_a = 1 / (r_a^2 - r_nearest^2 + 1)^2   (Eq. 8b)
!
! _sq variant takes squared distances directly to avoid nAtoms sqrt's per pt.
!------------------------------------------------------------------------------
pure function weighted_average_distance_sq(d_sq, weights, n_atoms) result(r_WA)
    integer(kind=int64), intent(in) :: n_atoms
    real(kind=real64), dimension(n_atoms), intent(in) :: d_sq, weights
    real(kind=real64) :: r_WA
    real(kind=real64) :: sum_weighted_sq, sum_weights
    integer(kind=int64) :: i

    sum_weighted_sq = 0.0_real64
    sum_weights = 0.0_real64
    do i = 1, n_atoms
        sum_weighted_sq = sum_weighted_sq + weights(i) * d_sq(i)
        sum_weights = sum_weights + weights(i)
    end do

    if (sum_weights > 1.0e-12_real64) then
        r_WA = sqrt(sum_weighted_sq / sum_weights)
    else
        r_WA = sqrt(d_sq(1))
    end if
end function weighted_average_distance_sq

!------------------------------------------------------------------------------
! GTO normalization factor
!------------------------------------------------------------------------------
function normalization_factor(alpha, prim_type, l_val) result(norm)
    real(kind=real64), intent(in) :: alpha
    character(len=7), intent(in) :: prim_type
    integer(kind=int64), intent(in) :: l_val
    real(kind=real64) :: norm
    integer(kind=int64), dimension(3) :: lmn
    real(kind=real64) :: df_l, df_m, df_n
    integer(kind=int64) :: sum_lmn
    
    if (prim_type(1:1) == 'S') then
        ! Spherical GTO
        norm = sqrt((2.0_real64**(l_val+2)) * ((2.0_real64*alpha)**(l_val+1.5_real64)) / &
                   (double_factorial(2*l_val + 1) * sqrt(pi)))
    else
        ! Cartesian GTO
        lmn = get_lmn(prim_type)
        sum_lmn = sum(lmn)
        
        norm = (2.0_real64 * alpha / pi)**(0.75_real64)
        
        if (sum_lmn > 0) then
            df_l = double_factorial(2*lmn(1) - 1)
            df_m = double_factorial(2*lmn(2) - 1)
            df_n = double_factorial(2*lmn(3) - 1)
            norm = norm * sqrt((4.0_real64 * alpha)**sum_lmn / (df_l * df_m * df_n))
        end if
    end if
end function normalization_factor

!------------------------------------------------------------------------------
! Double factorial n!!
!------------------------------------------------------------------------------
pure recursive function double_factorial(n) result(df)
    integer(kind=int64), intent(in) :: n
    integer(kind=int64) :: df
    
    if (n <= 1) then
        df = 1
    else
        df = n * double_factorial(n - 2)
    end if
end function double_factorial

!------------------------------------------------------------------------------
! Get angular momentum quantum number l from orbital type
!------------------------------------------------------------------------------
pure function get_l_value(orbital_type) result(l)
    character(len=*), intent(in) :: orbital_type
    integer(kind=int64) :: l
    
    if (orbital_type(1:1) == 'S') then
        select case(orbital_type(2:2))
        case('s'); l = 0
        case('p'); l = 1
        case('d'); l = 2
        case('f'); l = 3
        case('g'); l = 4
        case('h'); l = 5
        case('i'); l = 6
        case default; l = 0
        end select
    else
        l = sum(get_lmn(orbital_type))
    end if
end function get_l_value

!------------------------------------------------------------------------------
! Get number of basis functions in a shell
!------------------------------------------------------------------------------
pure function get_shell_size(shell_type) result(n)
    integer(kind=int64), intent(in) :: shell_type
    integer(kind=int64) :: n
    
    select case(shell_type)
    case(0);   n = 1    ! s
    case(1);   n = 3    ! p (cartesian)
    case(-1);  n = 4    ! sp
    case(2);   n = 6    ! d (cartesian)
    case(-2);  n = 5    ! d (spherical)
    case(3);   n = 10   ! f (cartesian)
    case(-3);  n = 7    ! f (spherical)
    case(4);   n = 15   ! g (cartesian)
    case(-4);  n = 9    ! g (spherical)
    case(5);   n = 21   ! h (cartesian)
    case(-5);  n = 11   ! h (spherical)
    case(6);   n = 28   ! i (cartesian)
    case(-6);  n = 13   ! i (spherical)
    case default; n = 1
    end select
end function get_shell_size

!------------------------------------------------------------------------------
! Get primitive orbital type label
!------------------------------------------------------------------------------
pure function get_primitive_type(shell_type, idx) result(prim_type)
    integer(kind=int64), intent(in) :: shell_type, idx
    character(len=7) :: prim_type
    
    select case(shell_type)
    case(0)
        prim_type = 's'
        
    case(1)  ! p cartesian
        select case(idx)
        case(1); prim_type = 'Px'
        case(2); prim_type = 'Py'
        case(3); prim_type = 'Pz'
        case default; prim_type = 'Px'
        end select
        
    case(-1)  ! sp
        select case(idx)
        case(1); prim_type = 's'
        case(2); prim_type = 'spPx'
        case(3); prim_type = 'spPy'
        case(4); prim_type = 'spPz'
        case default; prim_type = 's'
        end select
        
    case(2)  ! d cartesian
        select case(idx)
        case(1); prim_type = 'Dxx'
        case(2); prim_type = 'Dyy'
        case(3); prim_type = 'Dzz'
        case(4); prim_type = 'Dxy'
        case(5); prim_type = 'Dxz'
        case(6); prim_type = 'Dyz'
        case default; prim_type = 'Dxx'
        end select
        
    case(-2)  ! d spherical
        select case(idx)
        case(1); prim_type = 'Sd0'
        case(2); prim_type = 'Sd1c'
        case(3); prim_type = 'Sd1s'
        case(4); prim_type = 'Sd2c'
        case(5); prim_type = 'Sd2s'
        case default; prim_type = 'Sd0'
        end select
        
    case(3)  ! f cartesian
        select case(idx)
        case(1);  prim_type = 'Fxxx'
        case(2);  prim_type = 'Fyyy'
        case(3);  prim_type = 'Fzzz'
        case(4);  prim_type = 'Fxyy'
        case(5);  prim_type = 'Fxxy'
        case(6);  prim_type = 'Fxxz'
        case(7);  prim_type = 'Fxzz'
        case(8);  prim_type = 'Fyzz'
        case(9);  prim_type = 'Fyyz'
        case(10); prim_type = 'Fxyz'
        case default; prim_type = 'Fxxx'
        end select
        
    case(-3)  ! f spherical
        select case(idx)
        case(1); prim_type = 'Sf0'
        case(2); prim_type = 'Sf1c'
        case(3); prim_type = 'Sf1s'
        case(4); prim_type = 'Sf2c'
        case(5); prim_type = 'Sf2s'
        case(6); prim_type = 'Sf3c'
        case(7); prim_type = 'Sf3s'
        case default; prim_type = 'Sf0'
        end select
        
    case(4)  ! g cartesian (15 components, Gaussian z-major ordering)
        select case(idx)
        case(1);  prim_type = 'Gzzzz'
        case(2);  prim_type = 'Gyzzz'
        case(3);  prim_type = 'Gxzzz'
        case(4);  prim_type = 'Gyyzz'
        case(5);  prim_type = 'Gxyzz'
        case(6);  prim_type = 'Gxxzz'
        case(7);  prim_type = 'Gyyyz'
        case(8);  prim_type = 'Gxyyz'
        case(9);  prim_type = 'Gxxyz'
        case(10); prim_type = 'Gxxxz'
        case(11); prim_type = 'Gyyyy'
        case(12); prim_type = 'Gxyyy'
        case(13); prim_type = 'Gxxyy'
        case(14); prim_type = 'Gxxxy'
        case(15); prim_type = 'Gxxxx'
        case default; prim_type = 'Gzzzz'
        end select
        
    case(-4)  ! g spherical
        select case(idx)
        case(1); prim_type = 'Sg0'
        case(2); prim_type = 'Sg1c'
        case(3); prim_type = 'Sg1s'
        case(4); prim_type = 'Sg2c'
        case(5); prim_type = 'Sg2s'
        case(6); prim_type = 'Sg3c'
        case(7); prim_type = 'Sg3s'
        case(8); prim_type = 'Sg4c'
        case(9); prim_type = 'Sg4s'
        case default; prim_type = 'Sg0'
        end select
        
    case(5)  ! h cartesian (21 components, Gaussian z-major ordering)
        select case(idx)
        case(1);  prim_type = 'Hzzzzz'
        case(2);  prim_type = 'Hyzzzz'
        case(3);  prim_type = 'Hxzzzz'
        case(4);  prim_type = 'Hyyzzz'
        case(5);  prim_type = 'Hxyzzz'
        case(6);  prim_type = 'Hxxzzz'
        case(7);  prim_type = 'Hyyyzz'
        case(8);  prim_type = 'Hxyyzz'
        case(9);  prim_type = 'Hxxyzz'
        case(10); prim_type = 'Hxxxzz'
        case(11); prim_type = 'Hyyyyz'
        case(12); prim_type = 'Hxyyyz'
        case(13); prim_type = 'Hxxyyz'
        case(14); prim_type = 'Hxxxyz'
        case(15); prim_type = 'Hxxxxz'
        case(16); prim_type = 'Hyyyyy'
        case(17); prim_type = 'Hxyyyy'
        case(18); prim_type = 'Hxxyyy'
        case(19); prim_type = 'Hxxxyy'
        case(20); prim_type = 'Hxxxxy'
        case(21); prim_type = 'Hxxxxx'
        case default; prim_type = 'Hzzzzz'
        end select
        
    case(-5)  ! h spherical
        select case(idx)
        case(1);  prim_type = 'Sh0'
        case(2);  prim_type = 'Sh1c'
        case(3);  prim_type = 'Sh1s'
        case(4);  prim_type = 'Sh2c'
        case(5);  prim_type = 'Sh2s'
        case(6);  prim_type = 'Sh3c'
        case(7);  prim_type = 'Sh3s'
        case(8);  prim_type = 'Sh4c'
        case(9);  prim_type = 'Sh4s'
        case(10); prim_type = 'Sh5c'
        case(11); prim_type = 'Sh5s'
        case default; prim_type = 'Sh0'
        end select
        
    case(6)  ! i cartesian (28 components, Gaussian z-major ordering)
        select case(idx)
        case(1);  prim_type = 'Izzzzzz'
        case(2);  prim_type = 'Iyzzzzz'
        case(3);  prim_type = 'Ixzzzzz'
        case(4);  prim_type = 'Iyyzzzz'
        case(5);  prim_type = 'Ixyzzzz'
        case(6);  prim_type = 'Ixxzzzz'
        case(7);  prim_type = 'Iyyyzzz'
        case(8);  prim_type = 'Ixyyzzz'
        case(9);  prim_type = 'Ixxyzzz'
        case(10); prim_type = 'Ixxxzzz'
        case(11); prim_type = 'Iyyyyzz'
        case(12); prim_type = 'Ixyyyzz'
        case(13); prim_type = 'Ixxyyzz'
        case(14); prim_type = 'Ixxxyzz'
        case(15); prim_type = 'Ixxxxzz'
        case(16); prim_type = 'Iyyyyyz'
        case(17); prim_type = 'Ixyyyyz'
        case(18); prim_type = 'Ixxyyyz'
        case(19); prim_type = 'Ixxxyyz'
        case(20); prim_type = 'Ixxxxyz'
        case(21); prim_type = 'Ixxxxxz'
        case(22); prim_type = 'Iyyyyyy'
        case(23); prim_type = 'Ixyyyyy'
        case(24); prim_type = 'Ixxyyyy'
        case(25); prim_type = 'Ixxxyyy'
        case(26); prim_type = 'Ixxxxyy'
        case(27); prim_type = 'Ixxxxxy'
        case(28); prim_type = 'Ixxxxxx'
        case default; prim_type = 'Izzzzzz'
        end select
        
    case(-6)  ! i spherical
        select case(idx)
        case(1);  prim_type = 'Si0'
        case(2);  prim_type = 'Si1c'
        case(3);  prim_type = 'Si1s'
        case(4);  prim_type = 'Si2c'
        case(5);  prim_type = 'Si2s'
        case(6);  prim_type = 'Si3c'
        case(7);  prim_type = 'Si3s'
        case(8);  prim_type = 'Si4c'
        case(9);  prim_type = 'Si4s'
        case(10); prim_type = 'Si5c'
        case(11); prim_type = 'Si5s'
        case(12); prim_type = 'Si6c'
        case(13); prim_type = 'Si6s'
        case default; prim_type = 'Si0'
        end select
        
    case default
        prim_type = 's'
    end select
end function get_primitive_type

!------------------------------------------------------------------------------
! Get Cartesian exponents (l,m,n) from orbital label
!------------------------------------------------------------------------------
pure function get_lmn(label) result(lmn)
    character(len=*), intent(in) :: label
    integer(kind=int64), dimension(3) :: lmn
    
    select case(trim(label))
    case('s');                lmn = [0, 0, 0]
    case('Px', 'spPx');       lmn = [1, 0, 0]
    case('Py', 'spPy');       lmn = [0, 1, 0]
    case('Pz', 'spPz');       lmn = [0, 0, 1]
    case('Dxx');              lmn = [2, 0, 0]
    case('Dyy');              lmn = [0, 2, 0]
    case('Dzz');              lmn = [0, 0, 2]
    case('Dxy');              lmn = [1, 1, 0]
    case('Dxz');              lmn = [1, 0, 1]
    case('Dyz');              lmn = [0, 1, 1]
    case('Fxxx');             lmn = [3, 0, 0]
    case('Fyyy');             lmn = [0, 3, 0]
    case('Fzzz');             lmn = [0, 0, 3]
    case('Fxyy');             lmn = [1, 2, 0]
    case('Fxxy');             lmn = [2, 1, 0]
    case('Fxxz');             lmn = [2, 0, 1]
    case('Fxzz');             lmn = [1, 0, 2]
    case('Fyzz');             lmn = [0, 1, 2]
    case('Fyyz');             lmn = [0, 2, 1]
    case('Fxyz');             lmn = [1, 1, 1]
    ! g cartesian
    case('Gzzzz'); lmn = [0, 0, 4]
    case('Gyzzz'); lmn = [0, 1, 3]
    case('Gxzzz'); lmn = [1, 0, 3]
    case('Gyyzz'); lmn = [0, 2, 2]
    case('Gxyzz'); lmn = [1, 1, 2]
    case('Gxxzz'); lmn = [2, 0, 2]
    case('Gyyyz'); lmn = [0, 3, 1]
    case('Gxyyz'); lmn = [1, 2, 1]
    case('Gxxyz'); lmn = [2, 1, 1]
    case('Gxxxz'); lmn = [3, 0, 1]
    case('Gyyyy'); lmn = [0, 4, 0]
    case('Gxyyy'); lmn = [1, 3, 0]
    case('Gxxyy'); lmn = [2, 2, 0]
    case('Gxxxy'); lmn = [3, 1, 0]
    case('Gxxxx'); lmn = [4, 0, 0]
    ! h cartesian
    case('Hzzzzz'); lmn = [0, 0, 5]
    case('Hyzzzz'); lmn = [0, 1, 4]
    case('Hxzzzz'); lmn = [1, 0, 4]
    case('Hyyzzz'); lmn = [0, 2, 3]
    case('Hxyzzz'); lmn = [1, 1, 3]
    case('Hxxzzz'); lmn = [2, 0, 3]
    case('Hyyyzz'); lmn = [0, 3, 2]
    case('Hxyyzz'); lmn = [1, 2, 2]
    case('Hxxyzz'); lmn = [2, 1, 2]
    case('Hxxxzz'); lmn = [3, 0, 2]
    case('Hyyyyz'); lmn = [0, 4, 1]
    case('Hxyyyz'); lmn = [1, 3, 1]
    case('Hxxyyz'); lmn = [2, 2, 1]
    case('Hxxxyz'); lmn = [3, 1, 1]
    case('Hxxxxz'); lmn = [4, 0, 1]
    case('Hyyyyy'); lmn = [0, 5, 0]
    case('Hxyyyy'); lmn = [1, 4, 0]
    case('Hxxyyy'); lmn = [2, 3, 0]
    case('Hxxxyy'); lmn = [3, 2, 0]
    case('Hxxxxy'); lmn = [4, 1, 0]
    case('Hxxxxx'); lmn = [5, 0, 0]
    ! i cartesian
    case('Izzzzzz'); lmn = [0, 0, 6]
    case('Iyzzzzz'); lmn = [0, 1, 5]
    case('Ixzzzzz'); lmn = [1, 0, 5]
    case('Iyyzzzz'); lmn = [0, 2, 4]
    case('Ixyzzzz'); lmn = [1, 1, 4]
    case('Ixxzzzz'); lmn = [2, 0, 4]
    case('Iyyyzzz'); lmn = [0, 3, 3]
    case('Ixyyzzz'); lmn = [1, 2, 3]
    case('Ixxyzzz'); lmn = [2, 1, 3]
    case('Ixxxzzz'); lmn = [3, 0, 3]
    case('Iyyyyzz'); lmn = [0, 4, 2]
    case('Ixyyyzz'); lmn = [1, 3, 2]
    case('Ixxyyzz'); lmn = [2, 2, 2]
    case('Ixxxyzz'); lmn = [3, 1, 2]
    case('Ixxxxzz'); lmn = [4, 0, 2]
    case('Iyyyyyz'); lmn = [0, 5, 1]
    case('Ixyyyyz'); lmn = [1, 4, 1]
    case('Ixxyyyz'); lmn = [2, 3, 1]
    case('Ixxxyyz'); lmn = [3, 2, 1]
    case('Ixxxxyz'); lmn = [4, 1, 1]
    case('Ixxxxxz'); lmn = [5, 0, 1]
    case('Iyyyyyy'); lmn = [0, 6, 0]
    case('Ixyyyyy'); lmn = [1, 5, 0]
    case('Ixxyyyy'); lmn = [2, 4, 0]
    case('Ixxxyyy'); lmn = [3, 3, 0]
    case('Ixxxxyy'); lmn = [4, 2, 0]
    case('Ixxxxxy'); lmn = [5, 1, 0]
    case('Ixxxxxx'); lmn = [6, 0, 0]
    ! Spherical harmonics use l value for normalization
    case('Sp0', 'Sp1c', 'Sp1s');                                          lmn = [1, 0, 0]
    case('Sd0', 'Sd1c', 'Sd1s', 'Sd2c', 'Sd2s');                          lmn = [2, 0, 0]
    case('Sf0', 'Sf1c', 'Sf1s', 'Sf2c', 'Sf2s', 'Sf3c', 'Sf3s');          lmn = [3, 0, 0]
    case('Sg0', 'Sg1c', 'Sg1s', 'Sg2c', 'Sg2s', 'Sg3c', 'Sg3s', 'Sg4c', 'Sg4s')
        lmn = [4, 0, 0]
    case('Sh0', 'Sh1c', 'Sh1s', 'Sh2c', 'Sh2s', 'Sh3c', 'Sh3s', 'Sh4c', 'Sh4s', 'Sh5c', 'Sh5s')
        lmn = [5, 0, 0]
    case('Si0', 'Si1c', 'Si1s', 'Si2c', 'Si2s', 'Si3c', 'Si3s', 'Si4c', 'Si4s', 'Si5c', 'Si5s', 'Si6c', 'Si6s')
        lmn = [6, 0, 0]
    case default
        lmn = [0, 0, 0]
    end select
end function get_lmn

!------------------------------------------------------------------------------
! Evaluate spherical or Cartesian GTO at a point
!------------------------------------------------------------------------------
function eval_spherical_GTO(x_r, y_r, z_r, x_c, y_c, z_c, r_sq, orbital_type, alpha, l_val) result(value)
    real(kind=real64), intent(in) :: x_r, y_r, z_r, x_c, y_c, z_c, r_sq, alpha
    character(len=*), intent(in) :: orbital_type
    integer(kind=int64), intent(in) :: l_val
    real(kind=real64) :: value
    real(kind=real64) :: dx, dy, dz, r, theta, phi, exp_part, angular_part
    integer(kind=int64), dimension(3) :: lmn
    
    dx = x_r - x_c
    dy = y_r - y_c
    dz = z_r - z_c
    r = sqrt(r_sq)
    
    exp_part = exp(-alpha * r_sq)
    
    if (orbital_type(1:1) == 'S') then
        ! Spherical harmonic GTO
        if (r > 1.0e-10_real64) then
            theta = acos(dz / r)
            phi = atan2(dy, dx)
        else
            theta = 0.0_real64
            phi = 0.0_real64
        end if
        
        angular_part = real_spherical_harmonic(orbital_type, theta, phi)
        value = angular_part * (r**l_val) * exp_part
    else
        ! Cartesian GTO
        lmn = get_lmn(orbital_type)
        value = (dx**lmn(1)) * (dy**lmn(2)) * (dz**lmn(3)) * exp_part
    end if
end function eval_spherical_GTO

!------------------------------------------------------------------------------
! Real spherical harmonics with proper normalization
!------------------------------------------------------------------------------
function real_spherical_harmonic(orbital_type, theta, phi) result(value)
    character(len=*), intent(in) :: orbital_type
    real(kind=real64), intent(in) :: theta, phi
    real(kind=real64) :: value
    real(kind=real64) :: ct, st, cp, sp      ! cos/sin of theta/phi
    real(kind=real64) :: ct2, st2, c2p, s2p  ! squares and double angles
    real(kind=real64) :: ct3, st3, c3p, s3p
    real(kind=real64) :: ct4, st4, c4p, s4p
    real(kind=real64) :: ct5, st5, c5p, s5p
    real(kind=real64) :: c6p, s6p
    real(kind=real64) :: N  ! 1/sqrt(4*pi)
    
    N = 1.0_real64 / sqrt(4.0_real64 * pi)
    
    ct = cos(theta);  st = sin(theta)
    cp = cos(phi);    sp = sin(phi)
    ct2 = ct * ct;    st2 = st * st
    c2p = cos(2.0_real64 * phi);  s2p = sin(2.0_real64 * phi)
    ct3 = ct * ct2;   st3 = st * st2
    c3p = cos(3.0_real64 * phi);  s3p = sin(3.0_real64 * phi)
    ct4 = ct2 * ct2;  st4 = st2 * st2
    c4p = cos(4.0_real64 * phi);  s4p = sin(4.0_real64 * phi)
    ct5 = ct * ct4;   st5 = st * st4
    c5p = cos(5.0_real64 * phi);  s5p = sin(5.0_real64 * phi)
    c6p = cos(6.0_real64 * phi);  s6p = sin(6.0_real64 * phi)
    
    select case(trim(orbital_type))
    
    ! s (l=0)
    case('s')
        value = N
    
    ! p (l=1)
    case('Sp0');  value = N * sqrt(3.0_real64) * ct
    case('Sp1c'); value = N * sqrt(3.0_real64) * st * cp
    case('Sp1s'); value = N * sqrt(3.0_real64) * st * sp
    
    ! d (l=2)
    case('Sd0');  value = N * sqrt(5.0_real64) * 0.5_real64 * (3.0_real64*ct2 - 1.0_real64)
    case('Sd1c'); value = N * sqrt(15.0_real64) * st * ct * cp
    case('Sd1s'); value = N * sqrt(15.0_real64) * st * ct * sp
    case('Sd2c'); value = N * sqrt(15.0_real64) * 0.5_real64 * st2 * c2p
    case('Sd2s'); value = N * sqrt(15.0_real64) * 0.5_real64 * st2 * s2p
    
    ! f (l=3)
    case('Sf0');  value = N * sqrt(7.0_real64) * 0.5_real64 * (5.0_real64*ct3 - 3.0_real64*ct)
    case('Sf1c'); value = N * sqrt(21.0_real64/2.0_real64) * 0.5_real64 * st * (5.0_real64*ct2 - 1.0_real64) * cp
    case('Sf1s'); value = N * sqrt(21.0_real64/2.0_real64) * 0.5_real64 * st * (5.0_real64*ct2 - 1.0_real64) * sp
    case('Sf2c'); value = N * sqrt(105.0_real64) * 0.5_real64 * st2 * ct * c2p
    case('Sf2s'); value = N * sqrt(105.0_real64) * 0.5_real64 * st2 * ct * s2p
    case('Sf3c'); value = N * sqrt(35.0_real64/2.0_real64) * 0.5_real64 * st3 * c3p
    case('Sf3s'); value = N * sqrt(35.0_real64/2.0_real64) * 0.5_real64 * st3 * s3p
    
    ! g (l=4)
    case('Sg0');  value = N * 3.0_real64 * 0.375_real64 * (35.0_real64*ct4 - 30.0_real64*ct2 + 3.0_real64)
    case('Sg1c'); value = N * sqrt(10.0_real64) * 0.75_real64 * st * ct * (7.0_real64*ct2 - 3.0_real64) * cp
    case('Sg1s'); value = N * sqrt(10.0_real64) * 0.75_real64 * st * ct * (7.0_real64*ct2 - 3.0_real64) * sp
    case('Sg2c'); value = N * sqrt(5.0_real64) * 0.75_real64 * st2 * (7.0_real64*ct2 - 1.0_real64) * c2p
    case('Sg2s'); value = N * sqrt(5.0_real64) * 0.75_real64 * st2 * (7.0_real64*ct2 - 1.0_real64) * s2p
    case('Sg3c'); value = N * sqrt(70.0_real64) * 0.75_real64 * st3 * ct * c3p
    case('Sg3s'); value = N * sqrt(70.0_real64) * 0.75_real64 * st3 * ct * s3p
    case('Sg4c'); value = N * sqrt(35.0_real64) * 0.375_real64 * st4 * c4p
    case('Sg4s'); value = N * sqrt(35.0_real64) * 0.375_real64 * st4 * s4p
    
    ! h (l=5)
    case('Sh0');  value = N * sqrt(11.0_real64) * 0.25_real64 * (63.0_real64*ct5 - 70.0_real64*ct3 + 15.0_real64*ct)
    case('Sh1c'); value = N * sqrt(165.0_real64) * 0.25_real64 * st * (21.0_real64*ct4 - 14.0_real64*ct2 + 1.0_real64) * cp
    case('Sh1s'); value = N * sqrt(165.0_real64) * 0.25_real64 * st * (21.0_real64*ct4 - 14.0_real64*ct2 + 1.0_real64) * sp
    case('Sh2c'); value = N * sqrt(1155.0_real64) * 0.5_real64 * st2 * ct * (3.0_real64*ct2 - 1.0_real64) * c2p
    case('Sh2s'); value = N * sqrt(1155.0_real64) * 0.5_real64 * st2 * ct * (3.0_real64*ct2 - 1.0_real64) * s2p
    case('Sh3c'); value = N * sqrt(385.0_real64) * 0.25_real64 * st3 * (9.0_real64*ct2 - 1.0_real64) * c3p
    case('Sh3s'); value = N * sqrt(385.0_real64) * 0.25_real64 * st3 * (9.0_real64*ct2 - 1.0_real64) * s3p
    case('Sh4c'); value = N * sqrt(385.0_real64) * 0.75_real64 * st4 * ct * c4p
    case('Sh4s'); value = N * sqrt(385.0_real64) * 0.75_real64 * st4 * ct * s4p
    case('Sh5c'); value = N * sqrt(154.0_real64) * 0.375_real64 * st5 * c5p
    case('Sh5s'); value = N * sqrt(154.0_real64) * 0.375_real64 * st5 * s5p
    
    ! i (l=6)
    case('Si0')
        value = N * sqrt(13.0_real64) * 0.0625_real64 * &
                (231.0_real64*ct**6 - 315.0_real64*ct4 + 105.0_real64*ct2 - 5.0_real64)
    case('Si1c')
        value = N * sqrt(273.0_real64) * 0.125_real64 * st * &
                (33.0_real64*ct5 - 30.0_real64*ct3 + 5.0_real64*ct) * cp
    case('Si1s')
        value = N * sqrt(273.0_real64) * 0.125_real64 * st * &
                (33.0_real64*ct5 - 30.0_real64*ct3 + 5.0_real64*ct) * sp
    case('Si2c')
        value = N * sqrt(1365.0_real64) * 0.25_real64 * st2 * &
                (11.0_real64*ct4 - 3.0_real64*ct2) * c2p
    case('Si2s')
        value = N * sqrt(1365.0_real64) * 0.25_real64 * st2 * &
                (11.0_real64*ct4 - 3.0_real64*ct2) * s2p
    case('Si3c')
        value = N * sqrt(1365.0_real64) * 0.125_real64 * st3 * ct * &
                (11.0_real64*ct2 - 3.0_real64) * c3p
    case('Si3s')
        value = N * sqrt(1365.0_real64) * 0.125_real64 * st3 * ct * &
                (11.0_real64*ct2 - 3.0_real64) * s3p
    case('Si4c')
        value = N * sqrt(91.0_real64) * 0.375_real64 * st4 * &
                (11.0_real64*ct2 - 1.0_real64) * c4p
    case('Si4s')
        value = N * sqrt(91.0_real64) * 0.375_real64 * st4 * &
                (11.0_real64*ct2 - 1.0_real64) * s4p
    case('Si5c')
        value = N * sqrt(1001.0_real64) * 0.375_real64 * st5 * ct * c5p
    case('Si5s')
        value = N * sqrt(1001.0_real64) * 0.375_real64 * st5 * ct * s5p
    case('Si6c')
        value = N * sqrt(3003.0_real64) * 0.1875_real64 * st**6 * c6p
    case('Si6s')
        value = N * sqrt(3003.0_real64) * 0.1875_real64 * st**6 * s6p
    
    case default
        value = 0.0_real64
    end select
end function real_spherical_harmonic

!------------------------------------------------------------------------------
! Write output matrix file (updated to include APO vector)
!------------------------------------------------------------------------------
subroutine write_matrix_file(file_handle, file_name, analytical, numerical, cap, apo_vec, nb)
    type(mqc_gaussian_unformatted_matrix_file), intent(inout) :: file_handle
    character(len=*), intent(in) :: file_name
    type(mqc_scf_integral), intent(in) :: analytical
    type(mqc_matrix), intent(in) :: numerical, cap
    type(mqc_vector), intent(in) :: apo_vec
    integer(kind=int64), intent(in) :: nb
    type(mqc_matrix) :: apo_matrix
    integer(kind=int64) :: i
    real(kind=real64), dimension(:,:), allocatable :: apo_array
    
    ! Convert vector to column matrix for file writing
    ! Using explicit 2D array to ensure proper dimensions
    allocate(apo_array(nb, 1))
    do i = 1, nb
        apo_array(i, 1) = apo_vec%at(i)
    end do
    apo_matrix = apo_array
    deallocate(apo_array)
    
    ! Debug: print dimensions to verify
    write(6,'(A,I0,A,I0)') 'APO matrix dimensions: ', size(apo_matrix, 1), ' x ', size(apo_matrix, 2)
    
    if (file_handle%isOpen()) call Close_MatF(file_handle%UnitNumber)
    
    call file_handle%create(trim(file_name) // '.mat')
    call file_handle%writeESTObj('analytical overlap', est_integral=analytical)
    call file_handle%writeArray('numerical overlap', numerical)
    call file_handle%writeArray('Complex Absorbing Potential', cap)
    call file_handle%writeArray('Alpha Complex Absorbing Potential', cap)
    call file_handle%writeArray('Beta Complex Absorbing Potential', cap)
    call file_handle%writeArray('APO Vector', apo_matrix)
    call Close_MatF(file_handle%UnitNumber)
    
    write(6,'(A)') ''
    write(6,'(A)') 'Written to matrix file:'
    write(6,'(A)') '  - Analytical overlap'
    write(6,'(A)') '  - Numerical overlap'
    write(6,'(A)') '  - Complex Absorbing Potential (matrix)'
    write(6,'(A)') '  - Alpha/Beta Complex Absorbing Potential'
    write(6,'(A)') '  - APO Vector (W_p as column matrix, nBasis x 1)'
end subroutine write_matrix_file

end program noCaps_to_withCaps