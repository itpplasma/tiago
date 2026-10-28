"""Synthetic SIMPLE-shaped files using only the standard library and libnetcdf."""
import ctypes as c
import ctypes.util
import pathlib
import sys

nc = c.CDLL(ctypes.util.find_library('netcdf'))
nc.nc_create.argtypes = [c.c_char_p, c.c_int, c.POINTER(c.c_int)]
nc.nc_def_dim.argtypes = [c.c_int, c.c_char_p, c.c_size_t, c.POINTER(c.c_int)]
nc.nc_def_var.argtypes = [c.c_int, c.c_char_p, c.c_int, c.c_int, c.POINTER(c.c_int), c.POINTER(c.c_int)]
nc.nc_put_att_text.argtypes = [c.c_int, c.c_int, c.c_char_p, c.c_size_t, c.c_char_p]
nc.nc_put_att_int.argtypes = [c.c_int, c.c_int, c.c_char_p, c.c_int, c.c_size_t, c.POINTER(c.c_int)]
nc.nc_put_att_double.argtypes = [c.c_int, c.c_int, c.c_char_p, c.c_int, c.c_size_t, c.POINTER(c.c_double)]
nc.nc_put_var_double.argtypes = [c.c_int, c.c_int, c.POINTER(c.c_double)]
nc.nc_put_var_int.argtypes = [c.c_int, c.c_int, c.POINTER(c.c_int)]

def ok(result):
    if result != 0:
        raise RuntimeError(f'netCDF status {result}')

def fixture(path, legacy=False, unknown=False, bad=False):
    ident = c.c_int()
    ok(nc.nc_create(str(path).encode(), 0x1000, c.byref(ident)))
    ident = ident.value
    dims = []
    for name, size in [('particle', 4), ('phase', 5)]:
        dim = c.c_int()
        ok(nc.nc_def_dim(ident, name.encode(), size, c.byref(dim)))
        dims.append(dim.value)
    def text(var, name, value):
        ok(nc.nc_put_att_text(ident, var, name.encode(), len(value), value.encode()))
    def integer(var, name, value):
        ok(nc.nc_put_att_int(ident, var, name.encode(), 4, 1, c.byref(c.c_int(value))))
    variables = []
    def var(name, values, integer_type=False, units=None, shape=None):
        v = c.c_int()
        shape = dims[:1] if shape is None else shape
        ok(nc.nc_def_var(ident, name.encode(), 4 if integer_type else 6,
                         len(shape), (c.c_int * len(shape))(*shape), c.byref(v)))
        if units:
            text(v.value, 'units', units)
        variables.append((v.value, values, integer_type))
        return v.value
    var('times_lost', [0.25, 1.5, 2.0, 0.0], units='s')
    fv = var('wall_hit_triangle', [0, 1, -1, 99 if bad else -1], True)
    integer(fv, 'start_index', 0)
    var('wall_hit', [1, 1, 0, 1], True)
    var('orbit_exit_code', [2, 2, 1, 2 if bad else 3], True)
    if legacy:
        ok(nc.nc_put_att_double(ident, -1, b'particle_energy_eV', 6, 1, c.byref(c.c_double(100))))
        var('zend', [x for p in [0.5, 0.8, 1.0, 1.0] for x in [0, 0, 0, p, 0]], shape=dims)
    else:
        text(-1, 'results_schema_version', '1.1')
        text(-1, 'cartesian_frame', 'synthetic-Cartesian')
        var('kinetic_energy_eV', [25, 64, 100, 100], units='eV')
    wv = var('marker_weight', [1, 3, 2, 4], units='1')
    text(wv, 'source', 'equal_marker_default; source sampling weights unknown; external weights may override'
         if unknown else 'supplied_relative_weights')
    integer(-1, 'stl_facet_count', 2)
    text(-1, 'stl_file', 'synthetic-two-triangles.stl')
    ok(nc.nc_enddef(ident))
    for v, values, integral in variables:
        typ = c.c_int if integral else c.c_double
        put = nc.nc_put_var_int if integral else nc.nc_put_var_double
        ok(put(ident, v, (typ * len(values))(*values)))
    ok(nc.nc_close(ident))

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
fixture(out / 'events.nc')
fixture(out / 'legacy.nc', legacy=True)
fixture(out / 'unknown.nc', unknown=True)

fixture(out / 'bad-late.nc', bad=True)
