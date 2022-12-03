program main

    use slbm_2d_kinds, only : wp, ip
    use slbm_2d,       only : sim_test
    use slbm_2d,       only : line_seg_test

    implicit none

    !call sim_test
    call line_seg_test

end program main
