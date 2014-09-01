BLASPATH = /opt/OpenBlas/lib
F03 = gfortran

geomap : geomap.f03
	$(F03) -L$(BLASPATH) -lblas sha_helper.f03 geomap.f03 -O3 -mno-sse4.2 -Wall -Wextra -o geomap.a