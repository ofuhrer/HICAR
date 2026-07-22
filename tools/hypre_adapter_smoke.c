/* Four-rank smoke test for HICAR's decomposition-aware HYPRE adapter. */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "HYPRE.h"
#include "HYPRE_utilities.h"

int hicar_hypre_build(MPI_Fint, int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,
                      const float*,const float*,const float*,const float*,const float*,
                      const float*,const float*,const float*,const float*,const float*,
                      const float*,const float*,const float*,const float*,const float*);
int hicar_hypre_solve(const double*, double*, int, double, int*, double*);
void hicar_hypre_destroy(void);

int main(int argc, char **argv)
{
  const int xm=8, ym=8, zm=6, mx=16, my=16, mz=6;
  int rank, nprocs, xs, ys, q, i, j, k, its, rc;
  size_t n=(size_t)xm*ym*zm;
  float *a=calloc(n,sizeof(*a)), *b=calloc(n,sizeof(*b)), *c=calloc(n,sizeof(*c));
  float *d=calloc(n,sizeof(*d)), *e=calloc(n,sizeof(*e)), *f=calloc(n,sizeof(*f));
  float *g=calloc(n,sizeof(*g)), *z=calloc(n,sizeof(*z));
  double *rhs=calloc(n,sizeof(*rhs)), *x=calloc(n,sizeof(*x)), rnorm=0.0, rmax=0.0;
  MPI_Init(&argc,&argv); MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  if(nprocs != 4) { if(rank==0) fprintf(stderr,"requires four ranks\n"); MPI_Abort(MPI_COMM_WORLD,2); }
  xs=(rank%2)*xm; ys=(rank/2)*ym;
  for(j=0;j<ym;j++) for(k=0;k<zm;k++) for(i=0;i<xm;i++) {
    int gi=xs+i, gj=ys+j; q=(j*zm+k)*xm+i;
    if(gi==0 || gi==mx-1 || gj==0 || gj==my-1 || k==0 || k==mz-1) { a[q]=1.0f; }
    else { a[q]=-6.0f; b[q]=c[q]=d[q]=e[q]=f[q]=g[q]=1.0f; rhs[q]=1.0; }
  }
  rc=hicar_hypre_build(MPI_Comm_c2f(MPI_COMM_WORLD),xs,ys,0,xm,ym,zm,mx,my,mz,
                       xs,xs+xm-1,0,zm-1,ys,ys+ym-1,
                       a,b,c,d,e,f,g,z,z,z,z,z,z,z,z);
  if(rank==0) printf("adapter build_rc=%d\n",rc);
  if(!rc) rc=hicar_hypre_solve(rhs,x,500,1.e-8,&its,&rnorm);
  for(q=0;q<(int)n;q++) if(!isfinite(x[q])) rc=99;
  MPI_Allreduce(&rc,&q,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&rnorm,&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0) printf("adapter rc=%d iterations=%d relative_residual=%.8e\n",q,its,rmax);
  hicar_hypre_destroy(); HYPRE_Finalize();
  free(a);free(b);free(c);free(d);free(e);free(f);free(g);free(z);free(rhs);free(x);
  MPI_Finalize(); return q ? 1 : 0;
}
