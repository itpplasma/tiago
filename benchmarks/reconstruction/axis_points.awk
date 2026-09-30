# Points on the magnetic axis guess R(phi) = sum_n RAXIS(n) cos(n nfp phi),
# Z(phi) = -sum_n ZAXIS(n) sin(n nfp phi) (VMEC's RAXIS_CC, ZAXIS_CS), as
# free-boundary consistency points (x y z rows) inside the plasma.
# Usage: awk -v nfp=3 -v npts=12 -v raxis="1.5 0.1 ..." -v zaxis="0 -0.05 ..." -f axis_points.awk
BEGIN {
    pi = atan2(0, -1)
    nr = split(raxis, rc, " "); nz = split(zaxis, zs, " ")
    print "# x y z [m]: points on the magnetic axis guess"
    for (k = 0; k < npts; k++) {
        phi = 2 * pi * k / npts; r = 0; z = 0
        for (n = 0; n < nr; n++) r += rc[n + 1] * cos(n * nfp * phi)
        for (n = 0; n < nz; n++) z -= zs[n + 1] * sin(n * nfp * phi)
        printf "% .12E % .12E % .12E\n", r * cos(phi), r * sin(phi), z
    }
}
