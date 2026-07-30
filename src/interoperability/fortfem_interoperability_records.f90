module fortfem_interoperability_records
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: interoperability_record_t
        character(len=64) :: mesh_id = ""
        integer :: cells = 0
        integer :: dofs = 0
        real(dp) :: h = 0.0_dp
        real(dp) :: l2_error = 0.0_dp
        real(dp) :: graph_error = 0.0_dp
        real(dp) :: assembly_seconds = 0.0_dp
        real(dp) :: solve_seconds = 0.0_dp
        real(dp) :: total_seconds = 0.0_dp
    end type interoperability_record_t

    public :: read_interoperability_records

    character(len=*), parameter :: poisson_header = &
        "mesh_id,cells,dofs,h,l2_error,h1_error,"// &
        "assembly_seconds,solve_seconds,total_seconds"
    character(len=*), parameter :: ampere_header = &
        "mesh_id,cells,dofs,h,l2_error,hcurl_error,"// &
        "assembly_seconds,solve_seconds,total_seconds"

contains

    subroutine read_interoperability_records(path, records, status)
        character(len=*), intent(in) :: path
        type(interoperability_record_t), allocatable, intent(out) :: records(:)
        integer, intent(out) :: status

        character(len=512) :: line
        integer :: count, index, io_status, unit

        status = 0
        open (newunit=unit, file=path, status="old", action="read", &
            iostat=io_status)
        if (io_status /= 0) then
            allocate (records(0))
            status = 1
            return
        end if

        read (unit, "(a)", iostat=io_status) line
        if (io_status /= 0) then
            close (unit)
            allocate (records(0))
            status = 2
            return
        end if
        if (.not. valid_header(trim(line))) then
            close (unit)
            allocate (records(0))
            status = 2
            return
        end if

        count = 0
        do
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0) exit
            if (len_trim(line) > 0) count = count + 1
        end do
        rewind (unit)
        read (unit, "(a)") line
        allocate (records(count))

        index = 0
        do while (index < count)
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0) then
                status = 3
                exit
            end if
            if (len_trim(line) == 0) cycle
            index = index + 1
            read (line, *, iostat=io_status) records(index)%mesh_id, &
                records(index)%cells, records(index)%dofs, records(index)%h, &
                records(index)%l2_error, records(index)%graph_error, &
                records(index)%assembly_seconds, records(index)%solve_seconds, &
                records(index)%total_seconds
            if (io_status /= 0) then
                status = 3
                exit
            end if
            if (.not. valid_record(records(index))) then
                status = 4
                exit
            end if
        end do
        close (unit)
    end subroutine read_interoperability_records

    pure logical function valid_header(header)
        character(len=*), intent(in) :: header

        valid_header = header == poisson_header .or. header == ampere_header
    end function valid_header

    pure logical function valid_record(record)
        type(interoperability_record_t), intent(in) :: record
        real(dp) :: values(6)

        values = [record%h, record%l2_error, record%graph_error, &
            record%assembly_seconds, record%solve_seconds, &
            record%total_seconds]
        valid_record = len_trim(record%mesh_id) > 0
        if (.not. valid_record) return
        valid_record = record%cells > 0
        if (.not. valid_record) return
        valid_record = record%dofs > 0
        if (.not. valid_record) return
        valid_record = all(ieee_is_finite(values))
        if (.not. valid_record) return
        valid_record = all(values >= 0.0_dp)
        if (.not. valid_record) return
        valid_record = record%h > 0.0_dp
    end function valid_record

end module fortfem_interoperability_records
