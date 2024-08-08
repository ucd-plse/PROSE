module target_module
  implicit none

    contains
    
    subroutine funarc(result)

        real :: result
        integer :: i, j, k, l, n
        real :: h, t1, t2, dppi
        real :: s1

        n = 1000000

        do l=1,100
            t1 = -1.0
            dppi = acos(t1)
            s1 = 0.0;
            t1 = 0.0
            h = dppi / n

            do i=1,n
                t2 = fun(i*h) 
                s1 = s1 + sqrt (h*h + (t2 - t1)*(t2 - t1))
                t1 = t2
            enddo
        enddo

        result = s1

    end subroutine funarc


    function fun(x) result(t1)
        real, intent(in) :: x
        integer       :: k, n
        real  :: t1, d1

        n = 5
        d1 = 1.0

        t1 = x
        do k=1,n
            d1 = 2.0 * d1
            t1 = t1 + sin (d1 * x) / d1
        end do

    end function fun

end module target_module