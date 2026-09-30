# Synthetic magnetic diagnostics on a torus (major radius R0, minor radius rs)
# around the plasma, in DIAGNO format: diamagnetic loops (idia = 1), toroidal
# loops, saddle loops and segmented Rogowskis on nphi x nth patches, one
# Ampere Rogowski around the plasma, and nprobe B-probes with pseudo-random
# orientations (fixed seed).
# Usage: awk -v R0=1.42 -v rs=0.85 -v nphi=6 -v nth=4 -v nprobe=20 -v dir=DIR -f sensors.awk
function torus(r0, r, phi, th) {
    X = (r0 + r * cos(th)) * cos(phi); Y = (r0 + r * cos(th)) * sin(phi); Z = r * sin(th)
}
function point(file, area) {
    if (area == "") printf " % .12E % .12E % .12E\n", X, Y, Z > file
    else printf " % .12E % .12E % .12E % .12E\n", X, Y, Z, area > file
}
function header(file, n, idia, label) { printf "%6d%6d%6d %-48s\n", n, 0, idia, label > file }
BEGIN {
    pi = atan2(0, -1); side = 6; seg_area = 3.4e-4
    flux = dir "/flux.diagno"; seg = dir "/seg.diagno"; probes = dir "/probes.diagno"
    printf "%6d\n", nphi + 5 + nphi * nth > flux
    for (k = 0; k < nphi; k++) {
        header(flux, 49, 1, sprintf("DIA_%02d", k))
        for (i = 0; i < 49; i++) { torus(R0, rs, 2 * pi * k / nphi, 2 * pi * i / 48); point(flux) }
    }
    for (j = 0; j < 5; j++) {
        th = -pi / 2 + pi * j / 4
        header(flux, 97, 0, sprintf("TOR_%02d", j))
        for (i = 0; i < 97; i++) { torus(R0, rs, 2 * pi * i / 96, th); point(flux) }
    }
    dph = 0.8 * 2 * pi / nphi; dth = 0.8 * 2 * pi / nth
    printf "%6d\n", nphi * nth + 1 > seg
    for (k = 0; k < nphi; k++) for (j = 0; j < nth; j++) {
        p0 = 2 * pi * k / nphi; t0 = 2 * pi * j / nth
        header(flux, 4 * (side - 1) + 1, 0, sprintf("SAD_%02d_%02d", k, j))
        for (i = 0; i < side - 1; i++) { u = i / (side - 1); torus(R0, rs, p0 + dph * u, t0); point(flux) }
        for (i = 0; i < side - 1; i++) { u = i / (side - 1); torus(R0, rs, p0 + dph, t0 + dth * u); point(flux) }
        for (i = 0; i < side - 1; i++) { u = i / (side - 1); torus(R0, rs, p0 + dph * (1 - u), t0 + dth); point(flux) }
        for (i = 0; i < side - 1; i++) { u = i / (side - 1); torus(R0, rs, p0, t0 + dth * (1 - u)); point(flux) }
        torus(R0, rs, p0, t0); point(flux)
        header(seg, 8, 0, sprintf("SEG_%02d_%02d", k, j))
        for (i = 0; i < 8; i++) {
            torus(R0, rs, 2 * pi * (k + 0.5) / nphi, t0 + dth * i / 7); point(seg, seg_area / 7)
        }
    }
    header(seg, 101, 0, "AMPERE")
    for (i = 0; i < 101; i++) {
        t = 2 * pi * i / 100; X = R0 + 0.9 * cos(t); Y = 0; Z = 0.9 * sin(t); point(seg, seg_area / 100)
    }
    srand(1)
    printf "%d\n", nprobe > probes
    for (k = 0; k < nprobe; k++) {
        torus(R0, rs, 2 * pi * k / nprobe, 2 * pi * rand())
        printf " % .12E % .12E % .12E % .12E % .12E % .12E\n", X, Y, Z, 360 * rand(), 180 * rand(), 1.0e-3 > probes
    }
}
