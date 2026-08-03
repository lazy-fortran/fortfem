program sheet_current_surface_gallery
    !! Physical-first slab/cylinder/sphere/torus sheet-current gallery.
    !!
    !! Every geometry supplies only physical points, oriented unit normals,
    !! positive surface measures, and plus/minus traces to the canonical
    !! surface-parity contract.  The direct sums below are an independent
    !! ledger oracle; the fixture deliberately carries no plasma closure.
    use fortfem_interop, only: compare_sheet_current_surface_representations
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_parametric_surface, add_scatter, &
        colorbar, figure, savefig, title, view_init, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: thickness = 0.02_dp
    real(dp), parameter :: half_width = 5.0_dp*thickness
    real(dp), parameter :: axis_a(3) = [0.35_dp, -0.45_dp, 0.8_dp]
    real(dp), parameter :: axis_b(3) = [-0.55_dp, 0.2_dp, 0.65_dp]
    real(dp), parameter :: major_radius = 2.4_dp, minor_radius = 0.65_dp
    integer, parameter :: profile_count = 401
    integer, parameter :: slab_u = 8, slab_v = 8
    integer, parameter :: cylinder_phi = 16, cylinder_z = 8
    integer, parameter :: sphere_theta = 10, sphere_phi = 16
    integer, parameter :: torus_theta = 10, torus_phi = 16
    character(*), parameter :: output_directory = &
        "output/example/sheet_current_surface_gallery"

    real(dp) :: signed_distance(profile_count), normal_weights(profile_count)
    real(dp) :: spacing, start_time, end_time, elapsed_seconds
    real(dp) :: ampere_errors(4), maximum_ampere_error
    real(dp), allocatable :: slab_points(:, :), slab_normals(:, :), slab_weights(:)
    real(dp), allocatable :: slab_current(:, :), slab_plus(:, :), slab_minus(:, :)
    real(dp), allocatable :: cylinder_points(:, :), cylinder_normals(:, :)
    real(dp), allocatable :: cylinder_weights(:), cylinder_current(:, :)
    real(dp), allocatable :: cylinder_plus(:, :), cylinder_minus(:, :)
    real(dp), allocatable :: sphere_points(:, :), sphere_normals(:, :)
    real(dp), allocatable :: sphere_weights(:), sphere_current(:, :)
    real(dp), allocatable :: sphere_plus(:, :), sphere_minus(:, :)
    real(dp), allocatable :: torus_points(:, :), torus_normals(:, :)
    real(dp), allocatable :: torus_weights(:), torus_current(:, :)
    real(dp), allocatable :: torus_plus(:, :), torus_minus(:, :)
    integer :: sequence_unit, surface_unit, ledger_unit, metadata_unit, sample
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
        status="replace", action="write")
    open (newunit=surface_unit, file=output_directory//"/surface_current.csv", &
        status="replace", action="write")
    write (surface_unit, "(a)") &
        "geometry,point,x,y,z,nx,ny,nz,weight,kx,ky,kz"
    open (newunit=ledger_unit, file=output_directory//"/geometry_ledger.csv", &
        status="replace", action="write")
    write (ledger_unit, "(a)") &
        "geometry,points,area,fitted_x,fitted_y,fitted_z,resolved_x,"// &
        "resolved_y,resolved_z,relative_error"

    spacing = 2.0_dp*half_width/real(profile_count - 1, dp)
    do sample = 1, profile_count
        signed_distance(sample) = -half_width + spacing*real(sample - 1, dp)
        normal_weights(sample) = spacing
    end do
    normal_weights([1, profile_count]) = 0.5_dp*spacing

    call cpu_time(start_time)
    call build_slab(slab_points, slab_normals, slab_weights, slab_current, &
        slab_plus, slab_minus)
    call build_cylinder(cylinder_points, cylinder_normals, cylinder_weights, &
        cylinder_current, cylinder_plus, cylinder_minus)
    call build_sphere(sphere_points, sphere_normals, sphere_weights, sphere_current, &
        sphere_plus, sphere_minus)
    call build_torus(torus_points, torus_normals, torus_weights, torus_current, &
        torus_plus, torus_minus)

    call render_case("slab", slab_points, slab_normals, slab_weights, slab_current, &
        slab_plus, slab_minus, surface_unit, ledger_unit, ampere_errors(1), status)
    call write_stage(sequence_unit, "physical_solution_slab")
    call render_case("cylinder", cylinder_points, cylinder_normals, cylinder_weights, &
        cylinder_current, cylinder_plus, cylinder_minus, surface_unit, ledger_unit, &
        ampere_errors(2), status)
    call write_stage(sequence_unit, "physical_solution_cylinder")
    call render_case("sphere", sphere_points, sphere_normals, sphere_weights, sphere_current, &
        sphere_plus, sphere_minus, surface_unit, ledger_unit, ampere_errors(3), status)
    call write_stage(sequence_unit, "physical_solution_sphere")
    call render_case("torus", torus_points, torus_normals, torus_weights, torus_current, &
        torus_plus, torus_minus, surface_unit, ledger_unit, ampere_errors(4), status)
    call write_stage(sequence_unit, "physical_solution_torus")
    call write_stage(sequence_unit, "diagnostics")
    close (surface_unit)
    close (ledger_unit)

    call cpu_time(end_time)
    elapsed_seconds = end_time - start_time
    maximum_ampere_error = maxval(ampere_errors)
    if (.not. (elapsed_seconds > 0.0_dp)) &
        error stop "surface-current timing must be positive"
    if (.not. (maximum_ampere_error >= 0.0_dp)) &
        error stop "surface-current Ampere error must be finite"
    if (maximum_ampere_error >= 4.0e-12_dp) &
        error stop "surface-current Ampere error exceeds tolerance"
    close (sequence_unit)
    open (newunit=metadata_unit, file=output_directory//"/benchmark.json", &
        status="replace", action="write")
    write (metadata_unit, "(a)") "{"
    write (metadata_unit, "(a)") &
        '  "schema": "fortfem-sheet-current-surface-gallery-v2",'
    write (metadata_unit, "(a)") '  "physical_solution_first": true,'
    write (metadata_unit, "(a)") &
        '  "geometries": ["slab", "cylinder", "sphere", "torus"],'
    write (metadata_unit, "(a,i0,a)") &
        '  "total_surface_points": ', &
        slab_u*slab_v + cylinder_phi*cylinder_z + &
        sphere_theta*sphere_phi + torus_theta*torus_phi, ","
    write (metadata_unit, "(a,es24.16,a)") &
        '  "elapsed_seconds": ', elapsed_seconds, ","
    write (metadata_unit, "(a,es24.16)") &
        '  "max_integrated_ampere_relative_error": ', maximum_ampere_error
    write (metadata_unit, "(a)") "}"
    close (metadata_unit)

    deallocate(slab_points, slab_normals, slab_weights, slab_current, slab_plus, slab_minus)
    deallocate(cylinder_points, cylinder_normals, cylinder_weights, cylinder_current, &
        cylinder_plus, cylinder_minus)
    deallocate(sphere_points, sphere_normals, sphere_weights, sphere_current, sphere_plus, &
        sphere_minus)
    deallocate(torus_points, torus_normals, torus_weights, torus_current, torus_plus, &
        torus_minus)

contains

    subroutine build_slab(points, normals, weights, current, plus, minus)
        real(dp), allocatable, intent(out) :: points(:, :), normals(:, :), weights(:)
        real(dp), allocatable, intent(out) :: current(:, :), plus(:, :), minus(:, :)
        integer :: i, j, point
        real(dp) :: y, z, du, dv, normal(3), sheet(3)

        allocate(points(slab_u*slab_v, 3), normals(slab_u*slab_v, 3), &
            weights(slab_u*slab_v), current(slab_u*slab_v, 3), &
            plus(slab_u*slab_v, 3), minus(slab_u*slab_v, 3))
        du = 2.0_dp/real(slab_u, dp)
        dv = 2.0_dp/real(slab_v, dp)
        normal = [1.0_dp, 0.0_dp, 0.0_dp]
        point = 0
        do j = 1, slab_v
            z = -1.0_dp + (real(j, dp) - 0.5_dp)*dv
            do i = 1, slab_u
                y = -1.0_dp + (real(i, dp) - 0.5_dp)*du
                point = point + 1
                points(point, :) = [0.0_dp, y, z]
                normals(point, :) = normal
                weights(point) = du*dv
                call manufactured_current(normal, y, z, sheet)
                current(point, :) = sheet
                plus(point, :) = cross(sheet, normal)
                minus(point, :) = 0.0_dp
            end do
        end do
    end subroutine build_slab

    subroutine build_cylinder(points, normals, weights, current, plus, minus)
        real(dp), allocatable, intent(out) :: points(:, :), normals(:, :), weights(:)
        real(dp), allocatable, intent(out) :: current(:, :), plus(:, :), minus(:, :)
        integer :: i, j, point
        real(dp) :: phi, z, dphi, dz, normal(3), sheet(3)

        allocate(points(cylinder_phi*cylinder_z, 3), normals(cylinder_phi*cylinder_z, 3), &
            weights(cylinder_phi*cylinder_z), current(cylinder_phi*cylinder_z, 3), &
            plus(cylinder_phi*cylinder_z, 3), minus(cylinder_phi*cylinder_z, 3))
        dphi = 2.0_dp*pi/real(cylinder_phi, dp)
        dz = 2.0_dp/real(cylinder_z, dp)
        point = 0
        do j = 1, cylinder_z
            z = -1.0_dp + (real(j, dp) - 0.5_dp)*dz
            do i = 1, cylinder_phi
                phi = (real(i, dp) - 0.5_dp)*dphi
                point = point + 1
                normal = [cos(phi), sin(phi), 0.0_dp]
                points(point, :) = [cos(phi), sin(phi), z]
                normals(point, :) = normal
                weights(point) = dphi*dz
                call manufactured_current(normal, z, phi, sheet)
                current(point, :) = sheet
                plus(point, :) = cross(sheet, normal)
                minus(point, :) = 0.0_dp
            end do
        end do
    end subroutine build_cylinder

    subroutine build_sphere(points, normals, weights, current, plus, minus)
        real(dp), allocatable, intent(out) :: points(:, :), normals(:, :), weights(:)
        real(dp), allocatable, intent(out) :: current(:, :), plus(:, :), minus(:, :)
        integer :: i, j, point
        real(dp) :: theta, phi, dtheta, dphi, normal(3), sheet(3)

        allocate(points(sphere_theta*sphere_phi, 3), normals(sphere_theta*sphere_phi, 3), &
            weights(sphere_theta*sphere_phi), current(sphere_theta*sphere_phi, 3), &
            plus(sphere_theta*sphere_phi, 3), minus(sphere_theta*sphere_phi, 3))
        dtheta = pi/real(sphere_theta, dp)
        dphi = 2.0_dp*pi/real(sphere_phi, dp)
        point = 0
        do j = 1, sphere_theta
            theta = (real(j, dp) - 0.5_dp)*dtheta
            do i = 1, sphere_phi
                phi = (real(i, dp) - 0.5_dp)*dphi
                point = point + 1
                normal = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
                points(point, :) = normal
                normals(point, :) = normal
                weights(point) = sin(theta)*dtheta*dphi
                call manufactured_current(normal, theta, phi, sheet)
                current(point, :) = sheet
                plus(point, :) = cross(sheet, normal)
                minus(point, :) = 0.0_dp
            end do
        end do
    end subroutine build_sphere

    subroutine build_torus(points, normals, weights, current, plus, minus)
        real(dp), allocatable, intent(out) :: points(:, :), normals(:, :), weights(:)
        real(dp), allocatable, intent(out) :: current(:, :), plus(:, :), minus(:, :)
        integer :: i, j, point
        real(dp) :: theta, phi, dtheta, dphi, radius, normal(3), sheet(3)

        allocate(points(torus_theta*torus_phi, 3), normals(torus_theta*torus_phi, 3), &
            weights(torus_theta*torus_phi), current(torus_theta*torus_phi, 3), &
            plus(torus_theta*torus_phi, 3), minus(torus_theta*torus_phi, 3))
        dtheta = 2.0_dp*pi/real(torus_theta, dp)
        dphi = 2.0_dp*pi/real(torus_phi, dp)
        point = 0
        do j = 1, torus_theta
            theta = (real(j, dp) - 0.5_dp)*dtheta
            do i = 1, torus_phi
                phi = (real(i, dp) - 0.5_dp)*dphi
                point = point + 1
                radius = major_radius + minor_radius*cos(theta)
                normal = [cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]
                points(point, :) = [radius*cos(phi), radius*sin(phi), minor_radius*sin(theta)]
                normals(point, :) = normal
                weights(point) = minor_radius*radius*dtheta*dphi
                call manufactured_current(normal, theta, phi, sheet)
                current(point, :) = sheet
                plus(point, :) = cross(sheet, normal)
                minus(point, :) = 0.0_dp
            end do
        end do
    end subroutine build_torus

    subroutine manufactured_current(normal, first, second, current)
        real(dp), intent(in) :: normal(3), first, second
        real(dp), intent(out) :: current(3)
        real(dp) :: envelope

        envelope = 0.14_dp + 0.035_dp*cos(first - 2.0_dp*second) + &
            0.02_dp*sin(2.0_dp*first + second)
        current = envelope*cross(normal, axis_a)*dot_product(normal, axis_b)
    end subroutine manufactured_current

    subroutine render_case(name, points, normals, weights, current, plus, minus, &
            surface_unit, ledger_unit, relative_error, local_status)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: points(:, :), normals(:, :), weights(:), current(:, :)
        real(dp), intent(in) :: plus(:, :), minus(:, :)
        integer, intent(in) :: surface_unit, ledger_unit
        real(dp), intent(out) :: relative_error
        type(fortsparse_status_t), intent(out) :: local_status
        real(dp) :: fitted(3), resolved(3)
        real(dp) :: direct(3), jump(3), area, arrow_scale, maximum_current
        real(dp) :: point_extent(3), reference_extent
        real(dp) :: x2(2), y2(2), z2(2), direction(3)
        integer :: point, stride

        call compare_sheet_current_surface_representations(plus, minus, normals, weights, &
            signed_distance, normal_weights, thickness, fitted, resolved, relative_error, &
            local_status)
        if (local_status%code /= 0) error stop "surface parity call failed"
        direct = 0.0_dp
        area = sum(weights)
        do point = 1, size(weights)
            jump = plus(point, :) - minus(point, :)
            direct = direct + weights(point)*cross(normals(point, :), jump)
            if (abs(dot_product(normals(point, :), current(point, :))) > 2.0e-12_dp) &
                error stop "manufactured current is not tangential"
            write (surface_unit, '(a,",",i0,10(",",es24.16))') trim(name), point, &
                points(point, :), normals(point, :), weights(point), current(point, :)
        end do
        if (maxval(abs(direct - fitted)) > 2.0e-12_dp) &
            error stop "independent surface-current ledger disagrees"
        if (relative_error > 4.0e-12_dp) error stop "surface-current layer parity failed"
        write (ledger_unit, '(a,",",i0,8(",",es24.16))') trim(name), size(weights), &
            area, fitted, resolved, relative_error

        call figure(figsize=[8.8_dp, 7.0_dp])
        call add_geometry_surface(name)
        call add_scatter(points(:, 1), points(:, 2), points(:, 3), c=sqrt(sum(current**2, dim=2)), &
            cmap="viridis", marker="o", markersize=4.0_dp, label="surface samples")
        maximum_current = maxval(sqrt(sum(current**2, dim=2)))
        point_extent = maxval(points, dim=1) - minval(points, dim=1)
        reference_extent = minval(pack(point_extent, point_extent > 1.0e-12_dp))
        arrow_scale = 0.22_dp*reference_extent/max(1.0e-14_dp, maximum_current)
        stride = max(1, size(weights)/32)
        do point = 1, size(weights), stride
            direction = arrow_scale*current(point, :)
            x2 = [points(point, 1), points(point, 1) + direction(1)]
            y2 = [points(point, 2), points(point, 2) + direction(2)]
            z2 = [points(point, 3), points(point, 3) + direction(3)]
            call add_3d_plot(x2, y2, z2, color="black", linewidth=1.0_dp)
        end do
        call colorbar(label="|K| [A/m]")
        call xlabel("x [m]")
        call ylabel("y [m]")
        call title("Physical surface current on "//trim(name)//" (arrows show K)")
        call view_init(elev=27.0_dp, azim=-52.0_dp)
        call savefig(output_directory//"/surface_current_solutions_"//trim(name)//"_3d.png")
    end subroutine render_case

    subroutine add_geometry_surface(name)
        character(len=*), intent(in) :: name
        integer :: i, j
        real(dp), allocatable :: x(:, :), y(:, :), z(:, :)
        real(dp) :: theta, phi, radius, du, dv

        select case (trim(name))
        case ("slab")
            allocate(x(slab_u + 1, slab_v + 1), y(slab_u + 1, slab_v + 1), &
                z(slab_u + 1, slab_v + 1))
            do j = 1, slab_v + 1
                do i = 1, slab_u + 1
                    du = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(slab_u, dp)
                    dv = -1.0_dp + 2.0_dp*real(j - 1, dp)/real(slab_v, dp)
                    x(i, j) = 0.0_dp
                    y(i, j) = du
                    z(i, j) = dv
                end do
            end do
        case ("cylinder")
            allocate(x(cylinder_phi + 1, cylinder_z + 1), y(cylinder_phi + 1, cylinder_z + 1), &
                z(cylinder_phi + 1, cylinder_z + 1))
            do j = 1, cylinder_z + 1
                do i = 1, cylinder_phi + 1
                    phi = 2.0_dp*pi*real(i - 1, dp)/real(cylinder_phi, dp)
                    x(i, j) = cos(phi)
                    y(i, j) = sin(phi)
                    z(i, j) = -1.0_dp + 2.0_dp*real(j - 1, dp)/real(cylinder_z, dp)
                end do
            end do
        case ("sphere")
            allocate(x(sphere_phi + 1, sphere_theta + 1), y(sphere_phi + 1, sphere_theta + 1), &
                z(sphere_phi + 1, sphere_theta + 1))
            do j = 1, sphere_theta + 1
                theta = pi*real(j - 1, dp)/real(sphere_theta, dp)
                do i = 1, sphere_phi + 1
                    phi = 2.0_dp*pi*real(i - 1, dp)/real(sphere_phi, dp)
                    x(i, j) = sin(theta)*cos(phi)
                    y(i, j) = sin(theta)*sin(phi)
                    z(i, j) = cos(theta)
                end do
            end do
        case ("torus")
            allocate(x(torus_phi + 1, torus_theta + 1), y(torus_phi + 1, torus_theta + 1), &
                z(torus_phi + 1, torus_theta + 1))
            do j = 1, torus_theta + 1
                theta = 2.0_dp*pi*real(j - 1, dp)/real(torus_theta, dp)
                do i = 1, torus_phi + 1
                    phi = 2.0_dp*pi*real(i - 1, dp)/real(torus_phi, dp)
                    radius = major_radius + minor_radius*cos(theta)
                    x(i, j) = radius*cos(phi)
                    y(i, j) = radius*sin(phi)
                    z(i, j) = minor_radius*sin(theta)
                end do
            end do
        end select
        call add_parametric_surface(x, y, z, color="lightsteelblue", alpha=0.24_dp, &
            linewidth=0.0_dp, filled=.true., label="oriented surface")
        deallocate(x, y, z)
    end subroutine add_geometry_surface

    subroutine write_stage(unit, stage)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: stage
        write (unit, "(a)") trim(stage)
    end subroutine write_stage

    pure function cross(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)
        product = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross

end program sheet_current_surface_gallery
