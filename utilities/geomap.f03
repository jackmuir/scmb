program geomap
    !computes average and standard deviation maps from SHA
    use sha_helper
    implicit none
    integer, parameter :: no_paramsets = 40000 !can burn just by dropping the no here as scala results printed in reverse
    integer, parameter :: harm_order = 8
    integer, parameter :: harm_array_length = (harm_order+1)**2
    integer, parameter :: deg_spacing = 1
    integer, parameter :: lat_no = 180/deg_spacing+1, lon_no = 360/deg_spacing+1
    integer, parameter :: alength = lat_no*lon_no
    integer            :: counter, i, j, k
    integer            :: l_index_list(harm_array_length), m_index_list(harm_array_length)
    double precision               :: lat(alength), lon(alength), phi(alength), theta(alength), &
                                      mean(alength), std(alength), temp(alength)
    double precision               :: parameters(harm_array_length)
    double precision               :: shmatrix(alength,harm_array_length)
    counter = 1
    mean = 0
    std = 0
    temp = 0
    do i = 0, harm_order
        do j = -i, i
            l_index_list(counter) = i
            m_index_list(counter) = j
            counter = counter+1
        end do 
    end do
    counter = 1
    do i = 0, lon_no-1
        do j = 0, lat_no-1
            lat(counter) = deg_spacing*j-90.0
            lon(counter) = deg_spacing*i-180.0
            counter = counter+1
        end do
    end do
    theta = -(lat-90)*(pi/180.0)
    phi   = lon*pi/180.0
    print *, 'Setting up matrix'
    do k = 1, harm_array_length
        do j = 1, alength
            shmatrix(j,k) = ylmr(l_index_list(k),m_index_list(k),phi(j),theta(j))
        end do
    end do
    print *, 'Matrix set up'
    print *, 'Calculating mean'
    open(10, file = 'sha_parameters.dat', status = 'old')
    do i = 1, no_paramsets
        read(10, *) parameters
        do k = 1, harm_array_length
            mean = mean+parameters(k)*shmatrix(:,k)
        end do
    end do
    mean = mean/no_paramsets
    close(10)
    print *, 'Calculating std. dev.'
    open(10, file = 'sha_parameters.dat', status = 'old')
    do i = 1, no_paramsets
        read(10, *) parameters
        temp = 0
        do k = 1, harm_array_length
            temp = temp+parameters(k)*shmatrix(:,k)
        end do
        std = std + (mean-temp)**2
    end do
    std = sqrt(std/(no_paramsets-1)) !-1 to correct for bias slightly
    close(10)
    print *, 'Writing to file'
    open(11, file = 'geomaps.dat', status = 'replace')
    do j = 1, alength
        write(11, *) lat(j), lon(j), mean(j), std(j) 
    end do
end program geomap
