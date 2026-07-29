module fortfem_codegen_provenance
    implicit none

    private

    public :: fortsym_revision
    public :: generated_path

contains

    function fortsym_revision() result(revision)
        character(:), allocatable :: revision

        character(128) :: line
        integer :: ios, unit

        open( &
            newunit=unit, file="fortsym.lock", status="old", action="read", &
            iostat=ios)
        if (ios /= 0) error stop "cannot read fortsym.lock"
        read(unit, "(a)", iostat=ios) line
        close(unit)
        if (ios /= 0) error stop "cannot parse fortsym.lock"
        if (len_trim(line) /= 40) then
            error stop "fortsym.lock must contain a full SHA"
        end if
        revision = "fortsym@"//trim(line)
    end function fortsym_revision

    function generated_path(filename) result(path)
        character(*), intent(in) :: filename
        character(:), allocatable :: path

        character(4096) :: output_directory
        integer :: length, status

        call get_environment_variable( &
            "FORTFEM_CODEGEN_OUTPUT_DIR", output_directory, length=length, &
            status=status)
        if (status /= 0 .or. length == 0) then
            path = "../../src/generated/"//filename
            return
        end if
        if (output_directory(length:length) == "/") then
            path = output_directory(:length)//filename
        else
            path = output_directory(:length)//"/"//filename
        end if
    end function generated_path

end module fortfem_codegen_provenance
