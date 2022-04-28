/* ThinEPU2Pass.c
   Accelerator Toolbox 
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "elempass.h"
#include "r8lib.c"
#include "pwl_interp_2d.c"
#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"

struct elem 
{
    double NumX;
    double NumY;
    double *XGrid;
    double *YGrid;
    double *PxGrid;
    double *PyGrid;
    /* Optional fields */
    double *R1;
    double *R2;
    double *T1;
    double *T2;  
};


void fill(double *pout, double *pin, int m, int n)
{ int i, j;
  for (i=0; i<m; i++) {
     for (j=0; j<n; j++) {
        pout[i+j*m] = pin[i+j*m];
     }
  }
}

void ThinEPU2Pass(double *r, int nx, int ny, double *XEPU, double *YEPU, double *PXEPU, double *PYEPU,
                 double *R1, double *R2, double *T1, double *T2, int num_particles)

{ int c, ix, iy;
  double *r6, *xi, *yi, *xd, *yd, *PYEPUnew, *PXEPUnew;
  double *xkick, *ykick;

  /*
  printf("nP: %d \n",num_particles);
  printf("nx: %d \n", nx);
  printf("ny: %d \n", ny);
  printf("XEPU: %.3e \n",XEPU[0]);
  printf("YEPU: %.3e \n",XEPU[0]);
  */

  xi       = (double*)malloc(1*sizeof(double));
  yi       = (double*)malloc(1*sizeof(double));
  xd       = (double*)malloc(nx*sizeof(double));
  yd       = (double*)malloc(ny*sizeof(double));
  PXEPUnew = (double*)malloc(nx*ny*sizeof(double));
  PYEPUnew = (double*)malloc(nx*ny*sizeof(double));
  
  for (ix = 0;ix<nx;ix++){/* values for iy=0 */
    xd[ix]=XEPU[ny*ix];
  }
  for (iy = 0;iy<ny;iy++){/* values for ix=0 */
    yd[iy]=YEPU[ny-1-iy];
  }
  for (iy = 0;iy<ny;iy++){
    for (ix = 0;ix<nx;ix++){
      PXEPUnew[iy+ix*ny] = PXEPU[ny-1-iy+ix*ny];
      PYEPUnew[iy+ix*ny] = PYEPU[ny-1-iy+ix*ny];
    }
  }
  /*
  printf("xd: %.3e %.3e %.3e ...  %.3e \n",xd[0],xd[1],xd[2],xd[nx-1]);
  printf("yd: %.3e %.3e %.3e ...  %.3e \n",yd[0],yd[1],yd[2],yd[ny-1]);
  
  printf("PXEPUnew: %.3e \n",PXEPU[0]);
  printf("PYEPUnew: %.3e \n",PYEPU[0]);

  printf("T1: %d \n",T1[0]);
  printf("R1: %d \n",R1[0]);
  */



  for (c = 0;c<num_particles;c++){/* Loop over particles */
     r6 = r+c*6;
     /* linearized misalignment transformation at the entrance */
     ATaddvv(r6,T1);
     ATmultmv(r6,R1);
     xi[0] = r6[0]; 
     yi[0] = r6[2]; 


     xkick = pwl_interp_2d(ny,nx,yd,xd,PXEPUnew,1,yi,xi);
     ykick = pwl_interp_2d(ny,nx,yd,xd,PYEPUnew,1,yi,xi);
     /*xkick = pwl_interp_2d(nx,ny,xd,yd,PXEPU,1,xi,yi);
     ykick = pwl_interp_2d(nx,ny,xd,yd,PYEPU,1,xi,yi);*/

     r6[1] += xkick[0];
     r6[3] += ykick[0];
    
     /*
     printf("xi: %.3e \n", xi[0]);
     printf("yi: %.3e \n", yi[0]);
     printf("xd: %.3e \n", xd[0]);
     printf("yd: %.3e \n", yd[0]);
     printf("xkick: %.3e \n", xkick[0]);
     printf("ykick: %.3e \n", ykick[0]);
     */

     /* linearized misalignment transformation at the exit */
     ATmultmv(r6,R2);
     ATaddvv(r6,T2);
  }
  free(yi);
  free(xi);
  free(yd);
  free(xd);
  free(PXEPUnew);
  free(PYEPUnew);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
            double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        int NumX, NumY;
        double *XGrid, *YGrid, *PxGrid, *PyGrid, *R1, *R2, *T1, *T2;
        NumX=atGetLong(ElemData,"NumX");
        NumY=atGetLong(ElemData,"NumY");
        XGrid=atGetDoubleArray(ElemData,"XGrid"); check_error();
        YGrid=atGetDoubleArray(ElemData,"YGrid"); check_error();
        PxGrid=atGetDoubleArray(ElemData,"PxGrid"); check_error();
        PyGrid=atGetDoubleArray(ElemData,"PyGrid"); check_error();
          /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
 
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->NumX=NumX;
        Elem->NumY=NumY;
        Elem->XGrid=XGrid;
        Elem->YGrid=YGrid;
        Elem->PxGrid=PxGrid;
        Elem->PyGrid=PyGrid;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    ThinEPU2Pass(r_in,Elem->NumX,Elem->NumY,Elem->XGrid,Elem->YGrid,
            Elem->PxGrid,Elem->PyGrid,Elem->R1,Elem->R2,
            Elem->T1,Elem->T2,num_particles);
    return Elem;
}

MODULE_DEF(ThinEPU2Pass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/





#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        int NumX, NumY;
        double *XGrid, *YGrid, *PxGrid, *PyGrid, *R1, *R2, *T1, *T2;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        NumX=atGetLong(ElemData,"NumX");
        NumY=atGetLong(ElemData,"NumY");
        XGrid=atGetDoubleArray(ElemData,"XGrid"); check_error();
        YGrid=atGetDoubleArray(ElemData,"YGrid"); check_error();
        PxGrid=atGetDoubleArray(ElemData,"PxGrid"); check_error();
        PyGrid=atGetDoubleArray(ElemData,"PyGrid"); check_error();
          /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        ThinEPU2Pass(r_in, NumX, NumY, XGrid, YGrid, PxGrid, PyGrid, R1, R2, T1, T2, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("NumX"));
        mxSetCell(plhs[0],1,mxCreateString("NumY"));
        mxSetCell(plhs[0],2,mxCreateString("XGrid"));
        mxSetCell(plhs[0],3,mxCreateString("YGrid"));
        mxSetCell(plhs[0],4,mxCreateString("PxGrid"));
        mxSetCell(plhs[0],5,mxCreateString("PyGrid"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1);
            mxSetCell(plhs[1],0,mxCreateString("R1"));
            mxSetCell(plhs[1],1,mxCreateString("R2"));
            mxSetCell(plhs[1],2,mxCreateString("T1"));
            mxSetCell(plhs[1],3,mxCreateString("T2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */