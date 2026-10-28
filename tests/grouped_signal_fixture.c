#include <netcdf.h>
#include <string.h>

#define NC(call) do { int status = (call); if (status) return status; } while (0)

int grouped_create(const char *path, int *id)
{
    NC(nc_create(path, NC_NETCDF4 | NC_CLOBBER, id));
    return nc_enddef(*id);
}

int grouped_close(int id)
{
    return nc_close(id);
}

static int text(int id, const char *name, const char *value)
{
    return nc_put_att_text(id, NC_GLOBAL, name, strlen(value), value);
}

/* Use system netCDF independently of the Fortran compatibility writer. */
int grouped_write(int parent, const char *name, const char *channel,
                  const char *version, const char *record, const char *metadata,
                  const char *checksum, const double *times, const double *samples,
                  const double *variance, const double *covariance, int correlated,
                  double first, double last)
{
    int id, nested, dim, vars[4];
    int dims[2];
    const char *names[] = {"time", "samples", "variance", "covariance"};
    const double *values[] = {times, samples, variance, covariance};
    size_t chunks[2] = {4, 4};
    NC(nc_redef(parent));
    NC(nc_def_grp(parent, name, &id));
    NC(nc_def_dim(id, "sample", 4, &dim));
    /* The nested channel deliberately inherits its dimension. */
    if (strcmp(name, "channel-b") == 0) {
        NC(nc_def_grp(id, "nested", &nested));
        id = nested;
    }
    dims[0] = dims[1] = dim;
    for (int i = 0; i < 3 + correlated; ++i) {
        NC(nc_def_var(id, names[i], NC_DOUBLE, i == 3 ? 2 : 1, dims, &vars[i]));
        NC(nc_def_var_chunking(id, vars[i], NC_CHUNKED, chunks));
        NC(nc_def_var_deflate(id, vars[i], 1, 1, 1));
    }
    NC(text(id, "format", "tiago-signal-1"));
    NC(text(id, "channel", channel));
    NC(text(id, "source_id", "synthetic-source"));
    NC(text(id, "version", version));
    NC(text(id, "record_id", record));
    NC(text(id, "metadata", metadata));
    NC(text(id, "checksum", checksum));
    NC(text(id, "units", "V"));
    NC(nc_put_att_double(id, NC_GLOBAL, "valid_from", NC_DOUBLE, 1, &first));
    NC(nc_put_att_double(id, NC_GLOBAL, "valid_to", NC_DOUBLE, 1, &last));
    NC(nc_enddef(parent));
    for (int i = 0; i < 3 + correlated; ++i)
        NC(nc_put_var_double(id, vars[i], values[i]));
    return 0;
}
