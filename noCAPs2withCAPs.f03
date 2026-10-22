! Positive blocked quadrature with MQCPack basis evaluation (2026-10-02).
program noCaps_to_withCaps
!==============================================================================
!   Smooth Voronoi CAP integral generator + APO vector
!
!   Reference: Sommerfeld & Ehara, J. Chem. Theory Comput. 11, 4627 (2015),
!   "Complex Absorbing Potentials with Voronoi Isosurfaces Wrapping Perfectly
!   around Molecules".  Smooth Voronoi-CAP (eqs. 6-8):
!       r_WA(r) = sqrt( sum_a w_a r_a^2 / sum_a w_a )
!       w_a     = 1 / (r_a^2 - r_nearest^2 + a^2)^2, a = 1 a0 in paper
!       W(r)    = 0                  if r_WA <= r_cut
!       W(r)    = (r_WA - r_cut)^2   if r_WA  > r_cut
!
!   Positivity-preserving grid strategy:
!     All contracted AOs use the same positive quadrature, the default
!     Becke atom-centered grid or optional --centroid_grid:
!       - Lebedev angular grid (paper: 3470 points; default here: 1454)
!       - Gauss-Legendre radial grid (default: 500 points), selected after the
!         original equidistant leading-EM rule failed tight-core overlap tests.
!     --equidistant_radial retains a composite-Simpson equidistant comparison.
!
!   The CAP magnitude prefactor -i*eta is applied downstream (e.g. by
!   SNOWMAN); this generator produces eta-independent integrals.
!
!   Authors: L. M. Thompson, A. M. Kinyua, & P. D. Simpson, Louisville KY, 2023-2026
!==============================================================================

    use mqc_gaussian
    use mqc_integrals, only: mqc_gtoBasisSet, basisSetValuesList1, setup_quadrature_gauss_legendre1d, &
        molecular_quadrature_grid, gridBeckeMolecular, becke1988_radii_from_atomic_numbers
    use iso_fortran_env, only: int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib

    implicit none

    type(mqc_gaussian_unformatted_matrix_file) :: temp_file
    type(mqc_molecule_data) :: moleculeInfo
    type(mqc_scf_integral) :: overlap
    type(mqc_gtoBasisSet) :: basis
    type(molecular_quadrature_grid) :: becke_grid
    logical :: use_becke=.true., auto_rmax=.false., explicit_rmax=.false.
    integer(int64), allocatable :: atomic_numbers(:)
    real(real64), allocatable :: atomic_radii(:), angular_points(:,:), block_points(:,:)
    type(mqc_matrix) :: numerical_overlap, numerical_overlap_cap, density_matrix, analytical_matrix, density_beta
    logical :: automation_data=.false.
    integer(int64) :: data_unit, profile_unit, n_electrons
    real(real64), allocatable :: density_array(:,:), analytical_array(:,:), profile(:,:)
    real(real64) :: density_trace, rwa_point
    type(mqc_vector) :: apo_vector
    type(mqc_scalar) :: numerical_trace, numerical_cap_trace, analytical_trace
    character(len=:), allocatable :: command, fileName, sphereLebedev
    character(len=4096) :: matrix_path
    character(len=400) :: cutoff_text
    character(len=:), allocatable :: output_name, input_stem
    integer(int64) :: iOut=6, numFile, stat_num, iUnit, io_stat_number
    integer(int64) :: nBasis, nAtoms, i, j, k, shell_index, angular_momentum
    integer(int64), allocatable :: cap_atoms(:)
    real(real64) :: tail_radius, center_offset, alpha_min
    integer(int64) :: n_lebedev=1454, n_radial=500
    real(real64) :: r_min=0, r_max=20, r_cut=3, smoothing_length=1
    logical :: compute_validation=.false., unit_potential=.false.
    logical :: use_gauss_legendre_radial=.true.
    real(real64), allocatable :: xx(:), yy(:), zz(:), wweight(:)
    real(real64), allocatable :: leb_x(:), leb_y(:), leb_z(:), leb_w(:)
    real(real64), allocatable :: atom_coords_cache(:,:), radial_nodes(:), radial_weights(:)
    real(real64) :: theta_deg, phi_deg, theta_rad, phi_rad, w_full, radial_step
    ! Bounded grid blocks; no full-grid AO cache or dense per-thread matrices.
    integer(int64), parameter :: block_size=256
    integer(int64) :: first, count_points, point, global_point, ir, ia, total_grid_pts
    real(real64) :: centroid(3), xyz(3), radius, potential, t_start
    real(real64), allocatable :: values(:,:), weighted(:,:), values_point(:)
    real(real64), allocatable :: measure(:), cap_measure(:), s_sum(:,:), w_sum(:,:), v_sum(:)

!==============================================================================
! PROGRAM START
!==============================================================================

    write(iOut,'(*(A))') NEW_LINE('a'), ' ', repeat('*',73), NEW_LINE('a'), &
        '  Smooth Voronoi CAP integral generator', NEW_LINE('a'), &
        '  Gauss-Legendre radial + Lebedev angular (default: atom-centered Becke grid)', NEW_LINE('a'), &
        ' ', repeat('*',73), NEW_LINE('a'), NEW_LINE('a'), &
        repeat(' ',30), 'Version 26.10-SmoothVoronoi', NEW_LINE('a'), NEW_LINE('a'), &
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
        case('-h', '--help')
            call print_usage(iOut)
            stop
        case('-f')
            call mqc_get_command_argument(i+1, fileName)
            j = i + 2
        case('--output')
            call mqc_get_command_argument(i+1, output_name)
            k=len_trim(output_name)
            if (k>4) then
                if (output_name(k-3:k)=='.mat') output_name=output_name(:k-4)
            end if
            if (len_trim(output_name)==0) call mqc_error('Empty output filename',iOut)
            j = i + 2
        case('--lebedev_grid')
            call mqc_get_command_argument(i+1, sphereLebedev)
            j = i + 2
        case('--r_min')
            call mqc_get_command_argument(i+1, command)
            read(command, *) r_min
            j = i + 2
        case('--auto_rmax')
            auto_rmax=.true.
            j=i+1
        case('--r_max')
            explicit_rmax=.true.
            call mqc_get_command_argument(i+1, command)
            read(command, *) r_max
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
            read(command, *) r_cut
            j = i + 2
        case('--smoothing_length')
            call mqc_get_command_argument(i+1, command)
            read(command, *) smoothing_length
            j = i + 2
        case('--automation_data')
            automation_data=.true.
            j=i+1
        case('--centroid_grid')
            use_becke=.false.
            j=i+1
        case('--becke_grid')
            use_becke=.true.
            j=i+1
        case('--unit_potential')
            unit_potential = .true.
            j = i + 1
        case('--equidistant_radial')
            use_gauss_legendre_radial = .false.
            j = i + 1
        case('--validate')
            compute_validation = .true.
            j = i + 1
        case default
            call mqc_error_a('Unknown command-line argument', 6, &
                'argument', trim(command))
        end select
    end do

    write(6,'(A)') 'Using Smooth Voronoi CAP (Sommerfeld & Ehara, JCTC 2015)'

    !---------------------------------------------------------------------------
    ! Validate command-line parameters.  Each value has a default; we only
    ! sanity-check that they are physically reasonable.  Note: the strength
    ! parameter eta is intentionally not handled here -- the integrals
    ! generated by this program are eta-independent (W(r) = (r_WA - r_cut)^2);
    ! the caller (SNOWMAN) multiplies by -i*eta when assembling H - i*eta*W.
    !---------------------------------------------------------------------------
    if (auto_rmax.and.explicit_rmax) call mqc_error('Choose --auto_rmax or --r_max, not both',iOut)
    if (.not.all(ieee_is_finite([r_min,r_max,r_cut,smoothing_length]))) &
        call mqc_error('All length parameters must be finite', iOut)
    if (r_min < 0.0_real64) &
        call mqc_error('r_min must be >= 0', 6)
    if (.not.auto_rmax.and.r_max <= r_min) &
        call mqc_error('r_max must be > r_min', 6)
    if (r_cut < 0.0_real64) &
        call mqc_error('r_cut must be >= 0', 6)
    if (smoothing_length <= 0.0_real64) &
        call mqc_error('smoothing_length a must be > 0 bohr', 6)
    if (n_radial < 2_int64) &
        call mqc_error('n_radial must be >= 2', 6)
    if (n_lebedev < 1_int64) &
        call mqc_error('n_lebedev must be >= 1', 6)
    !---------------------------------------------------------------------------
    ! Read input file list
    !---------------------------------------------------------------------------
    if (.not. allocated(fileName)) call mqc_error('No input file provided', iOut)
    
    open(newunit=iUnit, file=fileName, status='old', iostat=stat_num)
    if (stat_num /= 0) call mqc_error_a('Error opening file', iOut, 'fileName', fileName)
    
    read(unit=iUnit, fmt='(i20)', iostat=stat_num) numFile
    if (stat_num /= 0) call mqc_error('Error reading file number', iOut)
    
    if (numFile /= 1) call mqc_error('Input list must contain exactly one matrix file', iOut)
    read(unit=iUnit, fmt='(A)', iostat=stat_num) matrix_path
    if (stat_num /= 0 .or. len_trim(matrix_path) == 0) call mqc_error('Missing matrix file path', iOut)
    close(unit=iUnit)

    !---------------------------------------------------------------------------
    ! Read Gaussian matrix file data
    !---------------------------------------------------------------------------
    call temp_file%getESTObj('overlap', est_integral=overlap, fileName=trim(matrix_path))
    analytical_matrix=overlap%getBlock('alpha')
    call temp_file%getMolData(moleculeInfo)
    call temp_file%loadGTOBasisSet(basis)
    nBasis = basis%nBasis
    nAtoms = MQC_Scalar_Get_Intrinsic_Integer(moleculeInfo%getNumAtoms())
    if (nAtoms < 1) call mqc_error('At least one atom is required', iOut)
    call mqc_print(nBasis, iOut, 'nBasis')
    call moleculeInfo%print(iOut)

    if (automation_data) then
        if (temp_file%isGeneral()) call mqc_error('Automation currently requires real RHF/ROHF/UHF density',iOut)
        call temp_file%getArray('ALPHA DENSITY MATRIX',density_matrix)
        if (MQC_Matrix_HaveComplex(density_matrix)) call mqc_error('Automation requires a real density',iOut)
        if (temp_file%isUnrestricted()) then
            call temp_file%getArray('BETA DENSITY MATRIX',density_beta)
            if (MQC_Matrix_HaveComplex(density_beta)) call mqc_error('Automation requires a real beta density',iOut)
            density_matrix=density_matrix+density_beta
        end if
        allocate(density_array(nBasis,nBasis),analytical_array(nBasis,nBasis),profile(2,block_size))
        do j=1,nBasis
            do i=1,nBasis
                density_array(i,j)=density_matrix%at(i,j)
                analytical_array(i,j)=analytical_matrix%at(i,j)
            end do
        end do
        n_electrons=temp_file%getVal('nElectrons')
        density_trace=sum(density_array*transpose(analytical_array))
        ! Gaussian restricted files may expose a one-spin or spin-summed record.
        ! Accept only these explicit conventions, with an electron-count check.
        if (temp_file%isRestricted().and.abs(2*density_trace-real(n_electrons,real64)) < 1e-7_real64) &
            density_array=2*density_array
        density_trace=sum(density_array*transpose(analytical_array))
        if (abs(density_trace-real(n_electrons,real64)) > 1e-6_real64) &
            call mqc_error('Density convention/electron count mismatch',iOut)
        if (.not.all(ieee_is_finite(density_array))) call mqc_error('Nonfinite density',iOut)
        write(iOut,'(A,ES24.16)') 'Automation density electron count = ',density_trace
        open(newunit=profile_unit,file='density_samples.bin',access='stream',form='unformatted',status='replace')
    end if

    !---------------------------------------------------------------------------
    ! Read Lebedev grid
    !---------------------------------------------------------------------------
    if (.not. allocated(sphereLebedev)) call mqc_error('No Lebedev file provided', 6)
    if (.not. is_supported_lebedev_size(n_lebedev)) &
        call mqc_error('Unsupported Lebedev point count', 6)
    
    open(newunit=iUnit, file=sphereLebedev, status='old', iostat=io_stat_number)
    if (io_stat_number /= 0) call mqc_error('Error opening Lebedev file', 6)

    allocate(xx(n_lebedev), yy(n_lebedev), zz(n_lebedev), wweight(n_lebedev))
    
    do k = 1, n_lebedev
        read(unit=iUnit, fmt=*, iostat=io_stat_number) theta_deg, phi_deg, wweight(k)
        if (io_stat_number /= 0) then
            write(*,*) 'Error reading Lebedev data at line ', k
            error stop 1
        end if
        
        theta_rad = theta_deg * (pi / 180.0_real64)
        phi_rad = phi_deg * (pi / 180.0_real64)
        
        xx(k) = cos(theta_rad) * sin(phi_rad)
        yy(k) = sin(theta_rad) * sin(phi_rad)
        zz(k) = cos(phi_rad)
    end do
    read(unit=iUnit, fmt=*, iostat=io_stat_number) theta_deg, phi_deg, w_full
    if (io_stat_number == 0) call mqc_error('Lebedev file has more points than --n_lebedev', 6)
    if (io_stat_number > 0) call mqc_error('Malformed trailing Lebedev data', 6)
    close(iUnit)

    if (.not.all(ieee_is_finite(xx)).or..not.all(ieee_is_finite(yy)).or..not.all(ieee_is_finite(zz))) &
        call mqc_error('Non-finite Lebedev coordinate', iOut)
    if (.not. all(ieee_is_finite(wweight))) call mqc_error('Non-finite Lebedev weight', 6)
    if (minval(wweight) <= 0.0_real64) call mqc_error('Common-grid CAP requires positive weights', 6)
    if (abs(sum(wweight)-1.0_real64) > 1.0e-10_real64) &
        call mqc_error('Lebedev weights must sum to one in supplied-file convention', 6)

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
    if (compute_validation) write(6,'(A,ES24.16)') 'Lebedev solid-angle weight sum = ', sum(leb_w)
    deallocate(xx, yy, zz, wweight)  ! No longer needed

    !---------------------------------------------------------------------------
    ! Cache coordinate data for fast access
    !---------------------------------------------------------------------------
    write(6,'(A)') 'Caching coordinate data...'

    ! Ghost centers (Bq, Z = 0) only carry basis functions: they are not CAP atoms.
    allocate(atom_coords_cache(3, nAtoms), cap_atoms(nAtoms))
    j = 0
    do i = 1, nAtoms
        if (int(moleculeInfo%Atomic_Numbers%at(i)) == 0) cycle
        j = j + 1
        cap_atoms(j) = i
        atom_coords_cache(1, j) = moleculeInfo%Cartesian_Coordinates%at(1, i)
        atom_coords_cache(2, j) = moleculeInfo%Cartesian_Coordinates%at(2, i)
        atom_coords_cache(3, j) = moleculeInfo%Cartesian_Coordinates%at(3, i)
    end do
    if (j < 1) call mqc_error('At least one non-ghost atom is required', iOut)
    if (j < nAtoms) write(iOut,'(A,I0,A)') 'Ghost centers skipped for the CAP: ', nAtoms-j, ' (Z = 0)'
    nAtoms = j
    atom_coords_cache = atom_coords_cache(:, :nAtoms)

    centroid = sum(atom_coords_cache, dim=2)/real(nAtoms,real64)
    if (auto_rmax) then
        ! Basis-tail heuristic, NOT a quadrature error bound or a physical CAP width.
        ! Enclose every shell's estimated tail about every integration center.
        r_max=0.0_real64
        do shell_index=1,basis%nShells
            alpha_min=minval(basis%shells(shell_index)%alpha)
            if (.not.ieee_is_finite(alpha_min).or.alpha_min <= 0) &
                call mqc_error('Automatic extent requires finite positive exponents',iOut)
            angular_momentum=maxval(sum(basis%shells(shell_index)%lArrays,dim=1))
            tail_radius=sqrt((50.0_real64+2.0_real64*real(angular_momentum,real64))/(2.0_real64*alpha_min))
            center_offset=0.0_real64
            if (use_becke) then
                do i=1,nAtoms
                    center_offset=max(center_offset,sqrt(sum((basis%shells(shell_index)%center &
                        -atom_coords_cache(:,i))**2)))
                end do
            else
                center_offset=sqrt(sum((basis%shells(shell_index)%center-centroid)**2))
            end if
            r_max=max(r_max,tail_radius+center_offset)
        end do
        if (.not.ieee_is_finite(r_max).or.r_max <= r_min) &
            call mqc_error('Automatic extent is nonfinite or does not exceed r_min',iOut)
        write(iOut,'(A,ES20.12)') 'Automatic initial r_max (bohr) = ',r_max
        write(iOut,'(A)') 'Heuristic only: converge radial/angular resolution and extent; r_cut is unchanged.'
    end if
    write(6,'(A)') 'Grid parameters (after defaults + validation):'
    write(6,'(A,F10.4)') '  r_min      (a0) = ', r_min
    write(6,'(A,F10.4)') '  r_max      (a0) = ', r_max
    write(6,'(A,F10.4)') '  r_cut      (a0) = ', r_cut
    write(6,'(A,F10.4)') '  smoothing a(a0) = ', smoothing_length
    write(6,'(A,I10)')   '  n_radial        = ', n_radial
    write(6,'(A,I10)')   '  n_lebedev       = ', n_lebedev
    if (unit_potential) write(6,'(A)') '  TEST MODE: scalar potential W(r)=1'
    if (use_gauss_legendre_radial) then
        write(6,'(A)') '  radial rule      = Gauss-Legendre'
    else
        write(6,'(A)') '  radial rule      = equidistant composite Simpson'
    end if
    if (compute_validation) &
        write(6,'(A)') 'Validation mode: numerical-vs-analytical overlap + symmetry checks will be reported'

    allocate(radial_nodes(n_radial), radial_weights(n_radial))
    if (use_gauss_legendre_radial) then
        call setup_quadrature_gauss_legendre1d(n_radial, r_min, r_max, radial_nodes, radial_weights)
    else
        radial_step = (r_max-r_min)/real(n_radial-1,real64)
        do i=1,n_radial
            radial_nodes(i) = r_min+real(i-1,real64)*radial_step
            radial_weights(i) = radial_quadrature_weight(i,n_radial,radial_step)
        end do
    end if
    if (n_radial > huge(total_grid_pts)/n_lebedev) call mqc_error('Grid count overflow',iOut)
    total_grid_pts = n_radial*n_lebedev
    if (use_becke) then
        allocate(atomic_numbers(nAtoms),angular_points(3,n_lebedev),block_points(3,block_size))
        do i=1,nAtoms
            atomic_numbers(i)=int(moleculeInfo%Atomic_Numbers%at(cap_atoms(i)))
        end do
        call becke1988_radii_from_atomic_numbers(atomic_numbers,atomic_radii)
        angular_points(1,:)=leb_x; angular_points(2,:)=leb_y; angular_points(3,:)=leb_z
        call gridBeckeMolecular(becke_grid,atom_coords_cache,atomic_radii, &
            radial_nodes,radial_weights,angular_points,leb_w)
        total_grid_pts=becke_grid%pointCount()
        write(iOut,'(A)') 'Grid: MQCPack atom-centered Becke partition (same smooth Voronoi potential)'
    else
        write(iOut,'(A)') 'Grid: molecular centroid'
    end if
    write(iOut,'(A,I0)') 'Grid points: ', total_grid_pts
    write(iOut,'(A,I0)') 'Maximum OpenMP threads: ', omp_get_max_threads()
    write(iOut,'(A)') 'Basis evaluation: MQCPack loadGTOBasisSet / basisSetValuesList1'
    allocate(values(nBasis,block_size), weighted(nBasis,block_size))
    allocate(measure(block_size), cap_measure(block_size))
    allocate(s_sum(nBasis,nBasis), w_sum(nBasis,nBasis), v_sum(nBasis))
    s_sum=0; w_sum=0; v_sum=0
    t_start=omp_get_wtime()
    ! Point evaluation is parallel; block reductions always occur in grid order.
    ! This keeps results independent of OpenMP scheduling and thread count.
    do first=1,total_grid_pts,block_size
        count_points=min(block_size,total_grid_pts-first+1)
        if (use_becke) then
            call becke_grid%copyBlock(first,block_points(:,:count_points),measure(:count_points))
            if (minval(measure(:count_points)) < 0) call mqc_error('Negative Becke quadrature weight',iOut)
        end if
        !$omp parallel default(shared) private(point,global_point,ir,ia,xyz,radius,potential,values_point,rwa_point)
        allocate(values_point(nBasis))
        !$omp do schedule(static)
        do point=1,count_points
            if (use_becke) then
                xyz=block_points(:,point)
            else
                global_point=first+point-2
                ir=global_point/n_lebedev+1
                ia=mod(global_point,n_lebedev)+1
                radius=radial_nodes(ir)
                xyz=centroid+radius*[leb_x(ia),leb_y(ia),leb_z(ia)]
                measure(point)=radial_weights(ir)*radius**2*leb_w(ia)
            end if
            call basisSetValuesList1(basis,xyz,values_point)
            values(:,point)=values_point
            potential=1.0_real64
            if (automation_data) then
                potential=smooth_voronoi(xyz,rwa_point)
                profile(1,point)=rwa_point
                if (unit_potential) potential=1.0_real64
            else if (.not.unit_potential) then
                potential=smooth_voronoi(xyz)
            end if
            cap_measure(point)=measure(point)*potential
        end do
        !$omp end do
        deallocate(values_point)
        !$omp end parallel
        if (.not.all(ieee_is_finite(values(:,:count_points)))) call mqc_error('Nonfinite AO values',iOut)
        if (.not.all(ieee_is_finite(cap_measure(:count_points)))) call mqc_error('Nonfinite CAP weights',iOut)
        if (automation_data) then
            weighted(:,:count_points)=matmul(density_array,values(:,:count_points))
            profile(2,:count_points)=sum(values(:,:count_points)*weighted(:,:count_points),dim=1)
            if (.not.all(ieee_is_finite(profile(:,:count_points)))) call mqc_error('Nonfinite density sample',iOut)
            if (minval(profile(2,:count_points)) < -1e-10_real64) call mqc_error('Negative sampled density',iOut)
            profile(2,:count_points)=max(0.0_real64,profile(2,:count_points))*measure(:count_points)
            write(profile_unit) profile(:,:count_points)
        end if
        v_sum=v_sum+matmul(values(:,:count_points),cap_measure(:count_points))
        do point=1,count_points
            weighted(:,point)=values(:,point)*sqrt(measure(point))
        end do
        s_sum=s_sum+matmul(weighted(:,:count_points),transpose(weighted(:,:count_points)))
        do point=1,count_points
            weighted(:,point)=values(:,point)*sqrt(cap_measure(point))
        end do
        w_sum=w_sum+matmul(weighted(:,:count_points),transpose(weighted(:,:count_points)))
    end do
    if (automation_data) then
        close(profile_unit)
        open(newunit=data_unit,file='cap_data.bin',access='stream',form='unformatted',status='replace')
        write(data_unit) nBasis,n_electrons
        write(data_unit) analytical_array,density_array,s_sum,w_sum
        close(data_unit)
    end if
    numerical_overlap=s_sum
    numerical_overlap_cap=w_sum
    apo_vector=v_sum
    write(iOut,'(A,F12.4,A)') 'CAP integration completed in ',omp_get_wtime()-t_start,' seconds'

!==============================================================================
! OUTPUT RESULTS
!==============================================================================

    call overlap%print(iOut, 'Analytical Overlap', Blank_At_Bottom=.true.)
    call numerical_overlap%print(iOut, 'Numerical Overlap')
    call numerical_overlap_cap%print(iOut, 'CAP matrix')
    

    numerical_trace = numerical_overlap%trace()
    numerical_cap_trace = numerical_overlap_cap%trace()
    analytical_trace = analytical_matrix%trace()

    write(*,'(A,F20.6)') 'Numerical Trace     = ', numerical_trace%rval()
    write(*,'(A,F20.6)') 'Numerical CAP Trace = ', numerical_cap_trace%rval()
    write(*,'(A,F20.6)') 'Analytical Trace    = ', analytical_trace%rval()

    if (compute_validation) call validate_results(overlap, numerical_overlap, &
                                                  numerical_overlap_cap, apo_vector, &
                                                  nBasis)

    if (.not.allocated(output_name)) then
      output_name='CAP_Mat'
      if (.not.automation_data) then
        input_stem=trim(matrix_path(scan(trim(matrix_path),'/',back=.true.)+1:))
        i=scan(input_stem,'.',back=.true.)
        if (i>1) input_stem=input_stem(:i-1)
        write(cutoff_text,'(F0.6)') r_cut
        cutoff_text=adjustl(cutoff_text)
        i=len_trim(cutoff_text)
        do while (i>1)
            if (cutoff_text(i:i)/='0') exit
            i=i-1
        end do
        if (cutoff_text(i:i)=='.') i=i-1
        if (i==0) then
            output_name='CAP_'//input_stem//'_rcut_0'
        else if (cutoff_text(1:1)=='.') then
            output_name='CAP_'//input_stem//'_rcut_0'//cutoff_text(:i)
        else
            output_name='CAP_'//input_stem//'_rcut_'//cutoff_text(:i)
        end if
      end if
    end if
    call write_matrix_file(temp_file, output_name, overlap, numerical_overlap, numerical_overlap_cap, apo_vector, nBasis)

contains

!==============================================================================
! HELPER FUNCTIONS
!==============================================================================

subroutine print_usage(unit_number)
    integer(kind=int64), intent(in) :: unit_number
    write(unit_number,'(A)') 'Usage: noCAPs2withCAPs -f FILELIST --lebedev_grid GRID [options]'
    write(unit_number,'(A)') 'Required:'
    write(unit_number,'(A)') '  -f FILELIST              text file: 1, then one Gaussian .mat path'
    write(unit_number,'(A)') '  --lebedev_grid FILE      theta(deg), polar(deg), normalized-weight file'
    write(unit_number,'(A)') 'Options (lengths in bohr):'
    write(unit_number,'(A)') '  --output FILE           override the generated matrix filename'
    write(unit_number,'(A)') '  --r_min X                radial lower bound (default 0; recommended 0)'
    write(unit_number,'(A)') '  --r_max X                radial upper bound (default 20)'
    write(unit_number,'(A)') '  --auto_rmax              estimate initial extent from Gaussian basis (still converge)'
    write(unit_number,'(A)') '  --r_cut X                CAP onset distance (default 3)'
    write(unit_number,'(A)') '  --smoothing_length X     Voronoi smoothing a (default 1)'
    write(unit_number,'(A)') '  --n_radial N             radial points (default 500)'
    write(unit_number,'(A)') '  --n_lebedev N            must match a supported grid and file (default 1454)'
    write(unit_number,'(A)') '  --equidistant_radial     use composite Simpson instead of Gauss-Legendre'
    write(unit_number,'(A)') '  --validate               overlap, Hermiticity, PSD, and closure diagnostics'
    write(unit_number,'(A)') '  --becke_grid             MQCPack atom-centered Becke integration grid (default)'
    write(unit_number,'(A)') '  --automation_data        export checked density samples and binary arrays for cap_auto.py'
    write(unit_number,'(A)') '  --centroid_grid          use the older common molecular-centroid grid'
    write(unit_number,'(A)') '  --unit_potential         test mode W(r)=1; CAP matrix must reproduce overlap'
    write(unit_number,'(A)') '  -h, --help               print this help'
    write(unit_number,'(A)') 'Output: CAP_<input-stem>_rcut_<cutoff>.mat in the current directory (unscaled W; no eta).'
end subroutine print_usage

pure function is_supported_lebedev_size(n_points) result(supported)
    integer(kind=int64), intent(in) :: n_points
    logical :: supported
    integer(kind=int64), parameter :: sizes(32) = [ &
        6_int64, 14_int64, 26_int64, 38_int64, 50_int64, 74_int64, 86_int64, 110_int64, &
        146_int64, 170_int64, 194_int64, 230_int64, 266_int64, 302_int64, 350_int64, &
        434_int64, 590_int64, 770_int64, 974_int64, 1202_int64, 1454_int64, 1730_int64, &
        2030_int64, 2354_int64, 2702_int64, 3074_int64, 3470_int64, 3890_int64, &
        4334_int64, 4802_int64, 5294_int64, 5810_int64 ]
    supported = any(sizes == n_points)
end function is_supported_lebedev_size

pure function radial_quadrature_weight(index, n_points, step) result(weight)
    integer(kind=int64), intent(in) :: index, n_points
    real(kind=real64), intent(in) :: step
    real(kind=real64) :: weight
    integer(kind=int64) :: tail_start

    weight = 0.0_real64
    if (n_points == 2) then
        weight = 0.5_real64*step
    else if (mod(n_points,2_int64) == 1) then
        if (index == 1 .or. index == n_points) then
            weight = step/3.0_real64
        else if (mod(index,2_int64) == 0) then
            weight = 4.0_real64*step/3.0_real64
        else
            weight = 2.0_real64*step/3.0_real64
        end if
    else if (n_points == 4) then
        if (index == 1 .or. index == 4) then
            weight = 3.0_real64*step/8.0_real64
        else
            weight = 9.0_real64*step/8.0_real64
        end if
    else
        ! Simpson 1/3 through point n-3, plus Simpson 3/8 on the last
        ! three intervals. The shared endpoint receives both endpoint weights.
        tail_start = n_points-3
        if (index <= tail_start) then
            if (index == 1 .or. index == tail_start) then
                weight = weight + step/3.0_real64
            else if (mod(index,2_int64) == 0) then
                weight = weight + 4.0_real64*step/3.0_real64
            else
                weight = weight + 2.0_real64*step/3.0_real64
            end if
        end if
        if (index >= tail_start) then
            if (index == tail_start .or. index == n_points) then
                weight = weight + 3.0_real64*step/8.0_real64
            else
                weight = weight + 9.0_real64*step/8.0_real64
            end if
        end if
    end if
end function radial_quadrature_weight

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
    type(mqc_matrix) :: cap_eigenvectors
    type(mqc_vector) :: cap_eigenvalues
    real(kind=real64), allocatable :: reconstructed(:,:)
    real(kind=real64) :: min_eval, max_eval, most_negative, scale, diff_norm, ref_norm
    real(kind=real64), parameter :: eta_diagnostic = 0.01_real64
    integer(kind=int64) :: kk, n_negative

    write(6,'(/,A)') repeat('=', 60)
    write(6,'(A)')   'Validation checks'
    write(6,'(A)')   repeat('=', 60)

    ! Analytical overlap, materialized as mqc_matrix for elementwise compare.
    analytical_mat = analytical%getBlock('alpha')

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

    call cap_mat%diag(EVals=cap_eigenvalues, EVecs=cap_eigenvectors)
    min_eval = huge(1.0_real64)
    max_eval = -huge(1.0_real64)
    most_negative = 0.0_real64
    n_negative = 0
    do kk = 1, nb
        val = cap_eigenvalues%at(kk)
        min_eval = min(min_eval, val)
        max_eval = max(max_eval, val)
        if (val < 0.0_real64) then
            n_negative = n_negative + 1
            most_negative = min(most_negative, val)
        end if
    end do
    write(6,'(A,ES12.4)') '  CAP minimum eigenvalue         = ', min_eval
    write(6,'(A,ES12.4)') '  CAP maximum eigenvalue         = ', max_eval
    write(6,'(A,I0,A,ES12.4)') '  CAP negative eigenvalues       = ', n_negative, &
        '; most negative = ', most_negative

    ! Diagnostic construction B=sqrt(2 eta Lambda_+) V^T. Reconstructing
    ! B^T B tests both the factor of two and placement of eta. Negative values
    ! are not silently altered in the CAP matrix; projection is confined to
    ! this explicitly reported diagnostic factorization.
    allocate(reconstructed(nb,nb))
    reconstructed = 0.0_real64
    do jj = 1, nb
        do ii = 1, nb
            do kk = 1, nb
                scale = cap_eigenvalues%at(kk)
                scale = max(0.0_real64, scale)
                reconstructed(ii,jj) = reconstructed(ii,jj) + 2.0_real64*eta_diagnostic*scale * &
                    cap_eigenvectors%at(ii,kk)*cap_eigenvectors%at(jj,kk)
            end do
        end do
    end do
    diff_norm = 0.0_real64
    ref_norm = 0.0_real64
    do jj = 1, nb
        do ii = 1, nb
            diff_norm = diff_norm + (reconstructed(ii,jj) - &
                2.0_real64*eta_diagnostic*cap_mat%at(ii,jj))**2
            ref_norm = ref_norm + (2.0_real64*eta_diagnostic*cap_mat%at(ii,jj))**2
        end do
    end do
    if (ref_norm > 0.0_real64) then
        write(6,'(A,ES12.4)') '  ||B^T B-2 eta W||/||2 eta W|| = ', sqrt(diff_norm/ref_norm)
    else
        write(6,'(A)') '  factorization closure          = zero CAP matrix'
    end if
    deallocate(reconstructed)

    if (max_diff_ovlp > 1.0e-5_real64) &
        write(6,'(A)') '  WARNING: numerical overlap deviates from analytical by > 1e-5 -- increase n_radial/n_lebedev?'
    if (max_asym_ovlp > 1.0e-10_real64 .or. max_asym_cap > 1.0e-10_real64) &
        write(6,'(A)') '  WARNING: symmetric-copy logic may be broken (off-diagonal asymmetry > 1e-10)'
    if (min_eval < -1.0e-10_real64*max(1.0_real64,abs(max_eval))) &
        call mqc_error('CAP has a significant negative eigenvalue',6)
    write(6,'(A)') ''
end subroutine validate_results

! Stable rescaling of Sommerfeld--Ehara eqs. 6--8. All distances in bohr.
function smooth_voronoi(point_xyz,rwa_out) result(w)
    real(real64), intent(in) :: point_xyz(3)
    real(real64), optional, intent(out) :: rwa_out
    real(real64) :: w, d2(nAtoms), scaled(nAtoms), delta(nAtoms), nearest, rwa
    integer(int64) :: atom
    do atom=1,nAtoms
        d2(atom)=sum((point_xyz-atom_coords_cache(:,atom))**2)
    end do
    nearest=minval(d2)
    delta=d2-nearest
    ! Divide differences by a twice, avoiding overflow/underflow of a**2.
    scaled=(1.0_real64/(1.0_real64+(delta/smoothing_length)/smoothing_length))**2
    rwa=sqrt(dot_product(scaled,d2)/sum(scaled))
    if (present(rwa_out)) rwa_out=rwa
    w=max(0.0_real64,rwa-r_cut)**2
end function smooth_voronoi

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
    write(6,'(A)') 'Written matrix file: '//trim(file_name)//'.mat'
end subroutine write_matrix_file

end program noCaps_to_withCaps
