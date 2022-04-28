# include <stdbool.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "pwl_interp_2d.h"
# include "r8lib.h"

/******************************************************************************/

double *pwl_interp_2d ( int nxd, int nyd, double *xd, double *yd, double *zd, 
  int ni, double *xi, double *yi )

/******************************************************************************/
/*
  Purpose:

    PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.

  Discussion:

    Thanks to Adam Hirst for pointing out an error in the formula that
    chooses the interpolation triangle, 04 February 2018.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 February 2018

  Author:

    John Burkardt

  Parameters:

    Input, int NXD, NYD, the number of X and Y data values.

    Input, double XD[NXD], YD[NYD], the sorted X and Y data.

    Input, double ZD[NXD*NYD}, the Z data.

    Input, int NI, the number of interpolation points.

    Input, double XI[NI], YI[NI], the coordinates of the
    interpolation points.

    Output, double PWL_INTERP_2D[NI], the value of the interpolant.
*/
{
  double alpha;
  double beta;
  double det;
  double dxa;
  double dxb;
  double dxi;
  double dya;
  double dyb;
  double dyi;
  double gamma;
  int i;
  int j;
  int k;
  double *zi;

  zi = ( double * ) malloc ( ni * sizeof ( double ) );

  for ( k = 0; k < ni; k++ )
  {
/*
  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:

    xd(i) <= xi(k) <= xd(i+1),
    yd(j) <= yi(k) <= yd(j+1).

  But if the interpolation point is not within a data interval, 
  assign the dummy interpolant value zi(k) = infinity.
*/
    i = r8vec_bracket5 ( nxd, xd, xi[k] );
    if ( i == -1 )
    {
      zi[k] = r8_huge ( );
      continue;
    }

    j = r8vec_bracket5 ( nyd, yd, yi[k] );
    if ( j == -1 )
    {
      zi[k] = r8_huge ( );
      continue;
    }
/*
  The rectangular cell is arbitrarily split into two triangles.
  The linear interpolation formula depends on which triangle 
  contains the data point.

    (I,J+1)--(I+1,J+1)
      |\       |
      | \      |
      |  \     |
      |   \    |
      |    \   |
      |     \  |
    (I,J)---(I+1,J)
*/
    if ( yi[k] < yd[j+1] + ( yd[j] - yd[j+1] ) * ( xi[k] - xd[i] ) / ( xd[i+1] - xd[i] ) )
    {
      dxa = xd[i+1] - xd[i];
      dya = yd[j]   - yd[j];

      dxb = xd[i]   - xd[i];
      dyb = yd[j+1] - yd[j];

      dxi = xi[k]   - xd[i];
      dyi = yi[k]   - yd[j];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

      zi[k] = alpha * zd[i+1+j*nxd] + beta * zd[i+(j+1)*nxd] + gamma * zd[i+j*nxd];
    }
    else
    {
      dxa = xd[i]   - xd[i+1];
      dya = yd[j+1] - yd[j+1];

      dxb = xd[i+1] - xd[i+1];
      dyb = yd[j]   - yd[j+1];

      dxi = xi[k]   - xd[i+1];
      dyi = yi[k]   - yd[j+1];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

      zi[k] = alpha * zd[i+(j+1)*nxd] + beta * zd[i+1+j*nxd] + gamma * zd[i+1+(j+1)*nxd];
    }
  }

  return zi;
}
