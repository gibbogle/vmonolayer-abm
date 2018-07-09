module envelope

use global
implicit none

integer, parameter :: MAXN = 5000
contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function getCentre() result (C)
integer :: kcell, n
real(REAL_KIND) :: csum(3), C(3)
type(cell_type), pointer :: cp

n = 0
csum = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
	cp => cell_list(kcell)
	csum = csum + cp%site
enddo
C = csum/n
end function

!-----------------------------------------------------------------------------------------
! Estimate blob radius (units sites) from max z-slice area.
!-----------------------------------------------------------------------------------------
function getRadius() result(R)
real(REAL_KIND) :: R
real(REAL_KIND) :: vol, area

!call getBlobCentreRange(cntr,rng, R)
call getVolume(vol, area)
R = sqrt(area/PI)
end function

!---------------------------------------------------------------------------------------
! Returns total blob volume, max z-slice area in units of sites
!---------------------------------------------------------------------------------------
subroutine getVolume(volume,maxarea)
real(REAL_KIND) :: volume, maxarea
real(REAL_KIND) :: area(NZ), cntr(2)
integer :: iz, iz0, izmin, izmax

iz0 = blob_centre(3)
do iz = iz0,NZ
	call getArea(iz,area(iz),cntr)
	if (area(iz) == 0) exit
	izmax = iz
enddo
do iz = iz0-1,1,-1
	call getArea(iz,area(iz),cntr)
	if (area(iz) == 0) exit
	izmin = iz
enddo
volume = 0
maxarea = 0
do iz = izmin,izmax
	volume = volume + area(iz)
	if (area(iz) > maxarea) maxarea = area(iz)
enddo
!write(*,'(a,i6,f10.1,f8.3)') 'getVolume: Ncells,volume: ',Ncells,volume,volume/Ncells
end subroutine

!---------------------------------------------------------------------------------------
! nsmax = dimension of rad(:,2)
!---------------------------------------------------------------------------------------
subroutine getSlices(nslices,dzslice,nsmax,rad)
integer :: nslices, nsmax
real(REAL_KIND) :: dzslice, rad(:,:)
real(REAL_KIND) :: area, cntr(2), r_inner(NZ), r_outer(NZ)
integer :: iz, iz1, iz2, i

!iz0 = blob_centre(3)
!do iz = iz0,NZ
!	call getArea(iz,area(iz),cntr)
!	if (area(iz) == 0) exit
!	izmax = iz
!enddo
!do iz = iz0-1,1,-1
!	call getArea(iz,area(iz),cntr)
!	if (area(iz) == 0) exit
!	izmin = iz
!enddo

do iz = 1,NZ
	call getArea(iz,area,cntr)
    r_outer(iz) = sqrt(area/PI)	! cm
    if (r_outer(iz) > 0) then
        call getInnerRadius(iz,cntr,r_inner(iz))
    else
        r_inner(iz) = 0
    endif
enddo
iz1 = 0
do iz = 1,NZ
    if (r_outer(iz) > 0) then
        if (iz1 == 0) iz1 = iz
        iz2 = iz
    endif
enddo
nslices = iz2 - iz1 + 1
nslices = min(nslices,nsmax)
do i = 1,nslices
    iz = iz1 + i - 1
    rad(i,1) = DELTA_X*r_inner(iz)
    rad(i,2) = DELTA_X*r_outer(iz)
enddo
dzslice = DELTA_X

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine getInnerRadius(iz,cntr,r_in)
integer :: iz
real(REAL_KIND) :: r_in, cntr(2)
integer :: ix, iy
real(REAL_KIND) :: r2, r2min, dx, dy

r2min = 1.0e10
do ix = 1,NX
	do iy = 1,NY
		if (occupancy(ix,iy,iz)%indx(1) > 0) then
	        dx = ix - cntr(1)
	        dy = iy - cntr(2)
	        r2 = dx**2 + dy**2
	        r2min = min(r2,r2min)
	    endif
	enddo
enddo
if (r2min > 0) then
    r_in = sqrt(r2min)
else
    r_in = 0
endif
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine getArea(iz, area, cntr)
integer :: iz
real(REAL_KIND) :: area, cntr(2)
real(REAL_KIND), allocatable :: x(:), y(:), xv(:), yv(:)
real(REAL_KIND) :: ave(2), xc, yc, d2, d2max, dx(4), dy(4), dx_um
real(REAL_KIND) :: delta = 0.20
integer, allocatable :: iwk(:), vertex(:)
integer :: ix, iy, n, nvert, i, i1, i2, k, kmax

allocate(x(MAXN))
allocate(y(MAXN))
allocate(iwk(MAXN))
allocate(vertex(1000))
allocate(xv(1000))
allocate(yv(1000))

! First test by looking at the z-slice approximately through the centre.
!iz = Centre(3)
n = 0
do ix = 1,NX
	do iy = 1,NY
		if (occupancy(ix,iy,iz)%indx(1) > 0) then
			n = n+1
			x(n) = ix
			y(n) = iy
		endif
	enddo
enddo
if (n <= 3) then
	area = n
	return
endif
call enveloper(x,y,n,vertex,nvert,iwk)
!write(*,*) 'get_diameter: nvert: ',nvert
ave = 0
do i = 1,nvert
	xv(i) = x(vertex(i))
	yv(i) = y(vertex(i))
!	write(*,'(i6,2f6.1)') i,xv(i),yv(i)
	ave(1) = ave(1) + xv(i)
	ave(2) = ave(2) + yv(i)
enddo
cntr = ave/nvert	! this is the approximate centre
! Now adjust (xv,yv) to the site corner most remote from the centre
dx = [-delta, delta, delta, -delta]
dy = [-delta, -delta, delta, delta]
do i = 1,nvert
	d2max = 0
	do k = 1,4
		xc = xv(i) + dx(k)
		yc = yv(i) + dy(k)
		d2 = (xc - cntr(1))**2 + (yc-cntr(2))**2
		if (d2 > d2max) then
			d2max = d2
			kmax = k
		endif
	enddo
	xv(i) = xv(i) + dx(kmax)
	yv(i) = yv(i) + dy(kmax)
!	write(*,'(i6,2f6.1)') i,xv(i),yv(i)
enddo
area = 0
do i1 = 1,nvert
	i2 = i1+1
	if (i2 > nvert) i2 = 1
	area = area + xv(i1)*yv(i2) - xv(i2)*yv(i1)
enddo
area = -area/2

deallocate(xv)
deallocate(yv)
deallocate(x)
deallocate(y)
deallocate(iwk)
deallocate(vertex)

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
SUBROUTINE enveloper(x, y, n, vertex, nvert, iwk)

!  Find the vertices (in clockwise order) of a polygon enclosing
!  the points (x(i), y(i), i=1, ..., n.

!  On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
!  iwk() is an integer work array which must have dimension at least n
!  in the calling program.

!  There is a limit of 1000 vertices imposed by the dimension of array next.

!  Programmer: Alan Miller
!  Latest revision - 12 September 1987
!  Fortran 90 version - 8 August 1996

INTEGER :: n, vertex(n), nvert, iwk(n)
real(REAL_KIND)    :: x(n), y(n)

!       Local variables

INTEGER :: next(1000), i, i1, i2, j, jp1, jp2, i2save, i3, i2next
real(REAL_KIND)    :: xmax, xmin, ymax, ymin, dist, dmax, dmin, x1, y1, dx, dy, x2, y2, &
           dx1, dx2, dmax1, dmax2, dy1, dy2, temp, zero = 0.0

IF (n < 2) RETURN

!  Choose the points with smallest & largest x- values as the
!  first two vertices of the polygon.

IF (x(1) > x(n)) THEN
  vertex(1) = n
  vertex(2) = 1
  xmin = x(n)
  xmax = x(1)
ELSE
  vertex(1) = 1
  vertex(2) = n
  xmin = x(1)
  xmax = x(n)
END IF

DO i = 2, n-1
  temp = x(i)
  IF (temp < xmin) THEN
    vertex(1) = i
    xmin = temp
  ELSE IF (temp > xmax) THEN
    vertex(2) = i
    xmax = temp
  END IF
END DO

!       Special case, xmax = xmin.

IF (xmax == xmin) THEN
  IF (y(1) > y(n)) THEN
    vertex(1) = n
    vertex(2) = 1
    ymin = y(n)
    ymax = y(1)
  ELSE
    vertex(1) = 1
    vertex(2) = n
    ymin = y(1)
    ymax = y(n)
  END IF

  DO i = 2, n-1
    temp = y(i)
    IF (temp < ymin) THEN
      vertex(1) = i
      ymin = temp
    ELSE IF (temp > ymax) THEN
      vertex(2) = i
      ymax = temp
    END IF
  END DO

  nvert = 2
  IF (ymax == ymin) nvert = 1
  RETURN
END IF

!  Set up two initial lists of points; those points above & those below the
!  line joining the first two vertices.    next(i) will hold the pointer to the
!  point furthest from the line joining vertex(i) to vertex(i+1) on the left
!  hand side.

i1 = vertex(1)
i2 = vertex(2)
iwk(i1) = -1
iwk(i2) = -1
dx = xmax - xmin
y1 = y(i1)
dy = y(i2) - y1
dmax = zero
dmin = zero
next(1) = -1
next(2) = -1

DO i = 1, n
  IF (i == vertex(1) .OR. i == vertex(2)) CYCLE
  dist = (y(i) - y1)*dx - (x(i) - xmin)*dy
  IF (dist > zero) THEN
    iwk(i1) = i
    i1 = i
    IF (dist > dmax) THEN
      next(1) = i
      dmax = dist
    END IF
  ELSE IF (dist < zero) THEN
    iwk(i2) = i
    i2 = i
    IF (dist < dmin) THEN
      next(2) = i
      dmin = dist
    END IF
  END IF
END DO

!  Ends of lists are indicated by pointers to -ve positions.

iwk(i1) = -1
iwk(i2) = -1
nvert = 2

j = 1

!  Start of main process.

!  Introduce new vertex between vertices j & j+1, if one has been found.
!  Otherwise increase j.   Exit if no more vertices.

40 IF (next(j) < 0) THEN
IF (j == nvert) RETURN
j = j + 1
GO TO 40
END IF

jp1 = j + 1
DO i = nvert, jp1, -1
  vertex(i+1) = vertex(i)
  next(i+1) = next(i)
END DO
jp2 = jp1 + 1
nvert = nvert + 1
IF (jp2 > nvert) jp2 = 1
i1 = vertex(j)
i2 = next(j)
i3 = vertex(jp2)
vertex(jp1) = i2

!  Process the list of points associated with vertex j.   New list at vertex j
!  consists of those points to the left of the line joining it to the new
!  vertex (j+1).   Similarly for the list at the new vertex.
!  Points on or to the right of these lines are dropped.

x1 = x(i1)
x2 = x(i2)
y1 = y(i1)
y2 = y(i2)
dx1 = x2 - x1
dx2 = x(i3) - x2
dy1 = y2 - y1
dy2 = y(i3) - y2
DMAX1 = zero
dmax2 = zero
next(j) = -1
next(jp1) = -1
i2save = i2
i2next = iwk(i2)
i = iwk(i1)
iwk(i1) = -1
iwk(i2) = -1

60 IF (i /= i2save) THEN
  dist = (y(i) - y1)*dx1 - (x(i) - x1)*dy1
  IF (dist > zero) THEN
    iwk(i1) = i
    i1 = i
    IF (dist > DMAX1) THEN
      next(j) = i
      DMAX1 = dist
    END IF
  ELSE
    dist = (y(i) - y2)*dx2 - (x(i) - x2)*dy2
    IF (dist > zero) THEN
      iwk(i2) = i
      i2 = i
      IF (dist > dmax2) THEN
        next(jp1) = i
        dmax2 = dist
      END IF
    END IF
  END IF
  i = iwk(i)
ELSE
  i = i2next
END IF

!  Get next point from old list at vertex j.

IF (i > 0) GO TO 60

!  End lists with -ve values.

iwk(i1) = -1
iwk(i2) = -1

GO TO 40
END SUBROUTINE enveloper

end module
