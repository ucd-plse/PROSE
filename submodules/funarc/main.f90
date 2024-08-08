program main
    use target_module, only: funarc
    implicit none

    real        :: result, start_t, end_t
    integer         :: dummy=0

    call CPU_TIME(start_t)
    call funarc(result)
    call CPU_TIME(end_t)

    print *, "out: ", result
    print *, "time: ", (end_t-start_t)
    
end program main