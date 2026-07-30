program verify_polar_complex
    use fortsym_arena, only: arena_t
    use fortsym_check, only: &
        check_identity, suite_begin, suite_end, suite_t
    use fortsym_engine_symengine, only: &
        make_symengine_engine, symengine_engine_t
    use fortsym_expr, only: expr_t
    use fortsym_parse, only: parse_expr
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: identity
    type(suite_t) :: suite
    character(:), allocatable :: message
    logical :: parsed

    call arena%init()
    engine = make_symengine_engine(arena)
    call suite_begin(suite, "FortFEM polar FEEC identities")

    identity = parse_expr(arena, &
        "(vn-(1-n2-n3)*p1-n2*p2-n3*p3)"// &
        "-(v-(1-c2-c3)*p1-c2*p2-c3*p3)"// &
        "-(vn-v)+(n2-c2)*(p2-p1)+(n3-c3)*(p3-p1)", parsed, message)
    if (.not. parsed) error stop message
    call check_identity(suite, engine, &
        "polar-cap exterior derivative squares to zero", identity)

    identity = parse_expr(arena, &
        "(b-a)+(d-b)-(d-c)-(c-a)", parsed, message)
    if (.not. parsed) error stop message
    call check_identity(suite, engine, &
        "regular-cell exterior derivative squares to zero", identity)

    call suite_end(suite, "results.json")
end program verify_polar_complex
