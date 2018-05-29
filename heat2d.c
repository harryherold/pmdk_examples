#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include <libpmem.h>

#include "layout_heat.h"

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define IDX(y, x, w) ((y) * (w) + (x))

#define SYNC(is_pmem, args...) ((is_pmem) ? pmem_persist(args) : pmem_msync(args) )

static void print_help(const char** argv)
{
    fprintf(stderr, "Usage: %s <N> <M> <MAX STEPS> <OUT-FILE>\n", argv[0]);
    fprintf(stderr, "<N> x <M> - Matrix size\n");
    fprintf(stderr, "<MAX STEPS> - Maximum iteration steps\n");
    fprintf(stderr, "<OUT-FILE> - filename for the output file\n");
}

void printstate(int step, double time, double h, double* vec, int n, int m);

void
compute(double * u, double * v, double * residual, int n, int m)
{
    double ln = 2.0;
    double h = ln / n;
    double hi = 1.0 / h;
    double hi2 = hi * hi;
    double alpha = 22.0;
    double eps = 1.0e-3;
    double u_min = 10.0;
    double u_max = 100.0;
    double dt = h * h / 4 / alpha;

    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            double du = (u[IDX(i, (j - 1), (m + 1))] + u[IDX(i, (j + 1), (m + 1))]
                    + u[IDX((i - 1), j, (m + 1))] + u[IDX((i + 1), j, (m + 1))]
                    - 4.0 * u[IDX(i, j, (m + 1))]) * dt * hi2 * alpha;

            v[IDX(i, j, (m + 1))] = u[IDX(i, j, (m + 1))] + du;

            du = MAX(du, -du);
            *residual = MAX(*residual, du);
        }
    }
}

int main(int argc, const char** argv)
{

    int n = 100;
    int m = 100;
    double ln = 2.0;
    double h = ln / n;
    double hi = 1.0 / h;
    double hi2 = hi * hi;
    double alpha = 22.0;
    double eps = 1.0e-3;
    double u_min = 10.0;
    double u_max = 100.0;
    double dt = h * h / 4 / alpha;
    double* u = NULL;
    double* v = NULL;
    double* t = NULL;

    if (argc < 5)
    {
        print_help(argv);
        return -1;
    }

    n = atoi(argv[1]);
    m = atoi(argv[2]);
    int    steps_max = atoi(argv[3]);
    int    is_pmem;
    size_t mapped_len;
    const size_t matrix_size = (n + 1) * (m + 1) * sizeof(double);
    const size_t file_size = sizeof(heat_info_t) + matrix_size * 2;

    void* pmemaddr = pmem_map_file(argv[4], file_size, PMEM_FILE_CREATE | PMEM_FILE_EXCL, 0666, &mapped_len, &is_pmem);

	if (pmemaddr == NULL)
    {
		perror("pmem_map_file");
		return 1;
	}

    heat_info_t * heat_info = (heat_info_t *) pmemaddr;
    heat_info->m = m;
    heat_info->n = n;
    heat_info->nsteps = 0;
    heat_info->offset_u = sizeof(*heat_info);
    heat_info->offset_v = sizeof(*heat_info) + matrix_size;

    u = (double *) ((char *) pmemaddr + heat_info->offset_u);
    v = (double *) ((char *) pmemaddr + heat_info->offset_v);

    for (int i = 0; i <= n; i++)
    {
        char is_hot = (i < (m / 3));
        for (int j = 0; j <= m; j++)
        {
            if (is_hot)
            {
                u[IDX(i, j, m + 1)] = u_max;
                v[IDX(i, j, m + 1)] = u_max;
            }
            else
            {
                u[IDX(i, j, m + 1)] = u_min;
                v[IDX(i, j, m + 1)] = u_min;
            }
        }
    }

    SYNC(is_pmem, pmemaddr, file_size);

    double residual = 0.0;
    int step = 0;

    do
    {

        residual = 0.0;
        compute(u, v, &residual, n, m);
        t = u;
        u = v;
        v = t;

        step++;

    } while ( steps_max > step);

    printf("step %7u done, residual %10.8f <= %10.8f\n", step, residual, eps);

    heat_info->nsteps = step;
    SYNC(is_pmem, pmemaddr, file_size);

    pmem_unmap(pmemaddr, mapped_len);

    return 0;
}

void printstate(int step, double time, double h, double* vec, int n, int m)
{
    int i, j;
    printf("### time step(%d)  %lf \n", step, time);
    for (i = 0; i <= 6; i++)
    {
        printf("[%d] :", i);
        for (j = 0; j <= m; j++)
        {

            printf("%lf ", vec[i * (m + 1) + j]);
        }
        printf("\n");
    }
    return;
}
