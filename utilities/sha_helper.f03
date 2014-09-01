module sha_helper
    implicit none
    private
    public :: xfact, mod_plgndr, ylmr
    double precision, public, parameter :: pi = 3.1415926536
    !This module contains functions that we need for spherical harmonic analysis
    !This includes an efficient method to calculate associated Legendre polynomials
    !and the double precision spherical harmonics
    !These functions are adapted from a proposal thread by Jozef Vesely on the SciPy development
    !pages which are in turn an improvement on the functions in Numerical Recipes in Fortran 77
    !https://github.com/scipy/scipy/issues/1280

contains
    !computes (2m-1)!!/sqrt((2m)!) 
    function xfact(m) result(xfact_result)
        integer, intent(in) :: m
        integer             :: i
        double precision                :: xfact_result
        xfact_result = 1.0D+0
        do i = 1, 2*m
            if (mod(i,2).eq.1) then
                xfact_result = xfact_result*i !(2m-1)!!
            end if
            xfact_result = xfact_result/sqrt(dble(i)) !sqrt((2m)!)
        end do
    end function xfact
    
    function mod_plgndr(l,m,x) result(mod_plgndr_result)
        integer, intent(in) :: l, m
        integer             :: ll
        double precision, intent(in)    :: x
        !double precision, parameter :: pi = 3.1415926536
        double precision :: norm, pmm, pmmp1, pll, mod_plgndr_result
        pll = 0.0D+0 !get rid of compiler issues
        if (m.lt.0.or.m.gt.l.or.abs(x).gt.1) then 
            write (*, *) "bad arguments to plgndr, aborting", m, x
            mod_plgndr_result=-10e6 !return a ridiculous value
        else
            norm = sqrt(2.0D+0*l+1.0D+0)/sqrt(4*pi)
            if (m.eq.0) then
                pmm = norm
            else
                pmm = (-1)**m*norm*xfact(m)*(1.0D+0-x*x)**(m/2.0D+0)
            end if
            if (l.eq.m) then
                mod_plgndr_result = pmm
            else 
                pmmp1 = x*pmm*sqrt(2.0D+0*m+1.0D+0)
                if (l.eq.m+1) then
                    mod_plgndr_result = pmmp1
                else
                    do ll = m+2, l
                        pll = (x*(2.0D+0*ll-1)*pmmp1 - sqrt((ll-1.0D+0)**2 - m*m)*pmm)/sqrt(ll**2.0D+0-m**2.0D+0)
                        pmm = pmmp1
                        pmmp1 = pll
                    end do
                    mod_plgndr_result = pll
                end if
            end if
        end if
        end function mod_plgndr
                
        function ylmr(l,m,phi,theta) result(ylmr_result)
            integer, intent(in) :: l, m
            double precision, intent(in)    :: phi, theta 
            double precision                :: ylmr_result
            if (m.gt.0) then
                ylmr_result = mod_plgndr(l,m,cos(theta))*cos(m*phi)*sqrt(2.0D+0)
            else if (m.lt.0) then
                ylmr_result = (-1)**m*mod_plgndr(l,-m,cos(theta))*sin(-m*phi)*sqrt(2.0D+0)
            else
                ylmr_result = mod_plgndr(l,m,cos(theta))
            end if
        end function ylmr

end module sha_helper     