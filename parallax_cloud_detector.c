/*----------------------------------------------------------------------------

  Cloud detector by parallax between the channels of push-broom satellite images

  Copyright (c) 2019-2020 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "iio.h" /* CHH */

/*----------------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif /* !M_PI */

/*----------------------------------------------------------------------------*/
#define UNDEF -1000

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard-error output and exit
 */
static void error(char * msg)
{
  fprintf(stderr,"error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, print an error and exit if fail
 */
static void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc: zero size");
  p = malloc(size);
  if( p == NULL ) error("xmalloc: out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
/* normalized angle difference between 'a' and the symmetric of 'b'
   relative to a vertical axis; the result is a value in [0,1]
 */
static double norm_angle_error(double a, double b)
{
  if( a == UNDEF || b == UNDEF ) return 1.0; /* return maximum error if undef */

  a -= b;
  while( a <= -M_PI ) a += 2.0*M_PI;
  while( a >   M_PI ) a -= 2.0*M_PI;

  return fabs(a) / M_PI;
}

/*----------------------------------------------------------------------------*/
/* compute image gradient direction (vectors are normalized to equal norm);
   the output components gx and gy must be already allocated
 */
static void image_gradient_direction( double * image, double * gx, double * gy,
                                      int X, int Y,
                                      unsigned char * input_mask) /* CHH */
{
  int x,y;

  /* initialize borders */
  for(x=0; x<X; x++) gx[x+0*X] = gx[x+(Y-1)*X] = gy[x+0*X] = gy[x+(Y-1)*X] = 0.;
  for(y=0; y<Y; y++) gx[0+y*X] = gx[X-1+y*X]   = gy[0+y*X] = gy[X-1+y*X]   = 0.;

  /* compute image gradient direction at each position */
  for(x=1; x<X-1; x++)
  for(y=1; y<Y-1; y++)
    {
      if(input_mask[x+y*X] > 0) /* CHH: use mask, set to 0. gradient for */
        {                       /*      invalid pixels (else clause)     */
          double norm;

          /* compute image gradient */
          gx[x+y*X] = image[x+1 +  y   *X] - image[x-1 +  y   *X];
          gy[x+y*X] = image[x   + (y+1)*X] - image[x   + (y-1)*X];

          norm = sqrt( gx[x+y*X] * gx[x+y*X] + gy[x+y*X] * gy[x+y*X] );

          /* normalize norm to 1, except for null gradient */
          if( norm > 0.0 )
            {
              gx[x+y*X] /= norm;
              gy[x+y*X] /= norm;
            }
        }
      else
        {
          gx[x+y*X] = 0.;
          gy[x+y*X] = 0.;
        }
    }
}

/*----------------------------------------------------------------------------*/
/* compute a sub-sampled 'optical flow' for each image pair and return
   only the offset angle for the non-null offsets (otherwise UNDEF).

   images : pointer to the images data
            there are 2N images (thus N image pairs)
            the value of pixel x,y in image i is given by:
               image[x + y*X + 2*i*X*Y]
   X,Y    : image domain size
   N      : number of image pairs
   W      : size of the square correlation window
   D      : maximal offset length to be evaluated
 */
static double * offset_angle_per_pair( double * images, int X, int Y, int N,
                                       int W, int D, double T,
                                       unsigned char * input_mask ) /* CHH */
{
  int U = (X - 2 - 2*D) / W;
  int V = (Y - 2 - 2*D) / W;
  int DD = 2 * D + 1; /* size of correlation matrix */
  double * correlation = (double *) xmalloc( DD * DD * sizeof(double) );
  double * angle;
  double * oflow_x; /* CHH */
  double * oflow_y; /* CHH */
  double * gx;
  double * gy;
  int u,v,i,iii;  /* CHH iii */
  char filenameA[128]; /* CHH */
  char filenameB[128]; /* CHH */

  /* get memory */
  gx    = (double *) xmalloc( X * Y * 2 * N * sizeof(double) );
  gy    = (double *) xmalloc( X * Y * 2 * N * sizeof(double) );
  angle = (double *) xmalloc( U * V * N * sizeof(double) );
  oflow_x = (double *) xmalloc( U * V * N * sizeof(double) ); /* CHH */
  oflow_y = (double *) xmalloc( U * V * N * sizeof(double) ); /* CHH */

  /* compute image gradient direction (the norm is normalized to 1) */
  for(i=0; i<2*N; i++)
    {
      image_gradient_direction(images+i*X*Y, gx+i*X*Y, gy+i*X*Y, X, Y, input_mask);
      /* << CHH */
      /* sprintf(filenameA, "/tmp/pcd_gradient_x_image_%d.tif", i); */
      /* sprintf(filenameB, "/tmp/pcd_gradient_y_image_%d.tif", i); */
      /* fprintf(stderr, "Writing gradient x for pair %d\n", i); */
      /* iio_write_image_double(filenameA, gx+i*X*Y, X, Y); */
      /* fprintf(stderr, "Writing gradient y for pair %d\n", i); */
      /* iio_write_image_double(filenameB, gy+i*X*Y, X, Y); */
      /* CHH >> */
    }

  /* compute local offset direction for each pairs of images */
  for(u=0; u<U; u++)  /* loops on spatial coordinates */
  for(v=0; v<V; v++)
    {
      int x = 1 + D + u * W;
      int y = 1 + D + v * W;

      if(input_mask[x+y*X] > 0)
        {
          for(i=0; i<N; i++)  /* loop on image pairs */
            {
              double max = -DBL_MAX;  /* max of gradient correlation */
              int mx,my,dx,dy,xx,yy;

              mx = my = -D;  /*  initialize offset of best correlation */
              for(dx=-D; dx<=D; dx++)  /* loops on relative offset to be evaluated */
              for(dy=-D; dy<=D; dy++)
                {
                  double c = 0.0;  /* initialize gradient correlation */

                  /* compute the correlation of the image gradient direction */
                  for(xx=x; xx<x+W; xx++)  /* loop on correlation window */
                  for(yy=y; yy<y+W; yy++)
                    c += gx[xx + yy*X + 2*i*X*Y] * gx[xx+dx + (yy+dy)*X + (2*i+1)*X*Y]
                       + gy[xx + yy*X + 2*i*X*Y] * gy[xx+dx + (yy+dy)*X + (2*i+1)*X*Y];

                  /* store correlation (to be used for sub-pixel interpolation */
                  correlation[ (dx+D) + (dy+D) * DD ] = c;

                  if( c > max )  /* update max of gradient correlation */
                    {
                      max = c;
                      mx = dx;  /* mx,my is offset of best correlation */
                      my = dy;
                    }
                }

              /* if the maximal correlation is not at the border of the correlation
                 domain perform a second order interpolation of the offset.
                 otherwise, or if the best offset is null, set as undefined
               */
              if( mx > -D && mx < D && my > -D && my < D ) // && (mx != 0 || my != 0) )
                {
                  double ax = correlation[ (mx-1+D) + (my+D) * DD ];
                  double bx = correlation[ (mx  +D) + (my+D) * DD ];
                  double cx = correlation[ (mx+1+D) + (my+D) * DD ];
                  double DX = (double) mx + 0.5 * (ax - cx) / (ax - bx - bx + cx);

                  double ay = correlation[ (mx+D) + (my-1+D) * DD ];
                  double by = correlation[ (mx+D) + (my  +D) * DD ];
                  double cy = correlation[ (mx+D) + (my+1+D) * DD ];
                  double DY = (double) my + 0.5 * (ay - cy) / (ay - by - by + cy);

                  if(DX*DX + DY*DY < T*T)
                    {
                      oflow_x[u + v*U + i*U*V] = NAN; /* CHH */
                      oflow_y[u + v*U + i*U*V] = NAN; /* CHH */
                      angle[u + v*U + i*U*V] = UNDEF; /* UNDEF */
                    }
                  else
                    {
                      oflow_x[u + v*U + i*U*V] = DX; /* CHH */
                      oflow_y[u + v*U + i*U*V] = DY; /* CHH */
                      angle[u + v*U + i*U*V] = atan2(DY,DX);
                    }
                }
              else
                {
                  oflow_x[u + v*U + i*U*V] = NAN; /* CHH */
                  oflow_y[u + v*U + i*U*V] = NAN; /* CHH */
                  angle[u + v*U + i*U*V] = UNDEF; /* UNDEF */
                }
            }
        }
      else /* CHH invalid pixels according to input_mask */
      {
        for(i=0; i<N; i++)
          angle[u + v*U + i*U*V] = UNDEF;
      }
    }

  for(iii=0; iii<N; iii++)  /* loop on image pairs */
    {
      /* << CHH write angle images */
      /* fprintf(stderr, "Writing angle\n"); */
      /* sprintf(filenameA, "/tmp/pcd_angle_image_%d.tif", iii); */
      /* iio_write_image_double(filenameA, angle+iii*U*V, U, V); */
      /* fprintf(stderr, "Writing optical flow\n"); */
      /* sprintf(filenameA, "/tmp/pcd_oflow_x_image_%d.tif", iii); */
      /* iio_write_image_double(filenameA, oflow_x+iii*U*V, U, V); */
      /* sprintf(filenameA, "/tmp/pcd_oflow_y_image_%d.tif", iii); */
      /* iio_write_image_double(filenameA, oflow_y+iii*U*V, U, V); */
      /* CHH >> */
    }

  /* free memory */
  free( (void *) gx );
  free( (void *) gy );
  free( (void *) correlation );
  free( (void *) oflow_x );
  free( (void *) oflow_y );

  return angle;
}

/*----------------------------------------------------------------------------*/
/* compute the number of pairs where the angles at position 'u,v' is defined
   and corresponds to the angle 'ref' up to a tolerance 'p'
 */
static int number_angle_agreements( double * angle, int U, int V, int N,
                                    int u, int v, double ref, double p )
{
  int k = 0;
  int l;

  for(l=0; l<N; l++)  /* loop on number of pairs */
    if( angle[u + v*U + l*U*V] != UNDEF &&
        norm_angle_error(ref, angle[u + v*U + l*U*V]) < p )
      ++k;

  return k;
}

/*----------------------------------------------------------------------------*/
/* compute a cloud mask by grouping regions with the same parallax angle
   between pairs of corresponding satellite images

   images : pointer to the images data
            there are 2N images (thus N image pairs)
            the value of pixel x,y in image i is given by:
               image[x + y*X + 2*i*X*Y]
   X,Y    : image domain size
   N      : number of image pairs
   W      : size of the square correlation window
   D      : maximal offset length to be evaluated
   input_mask: mask for valid pixels (do not apply the detector where
               pixels are not valid). 0 for invalid and 255 for valid.
 */
void parallax_cloud_detector( double * cloud_mask, double * images,
                              int X, int Y, int N, int W, int D, double T,
                              unsigned char * input_mask) /* CHH */
{
  int P = 6;  /* p of size P: list of precision values used in region growing */
  double p[] = { 0.025, 0.05, 0.1, 0.2, 0.3, 0.4 };
  int neigh_u[4] = { -1,  1,  0,  0 };    /* 4-connection neighbors */
  int neigh_v[4] = {  0,  0, -1,  1 };
  int U = (X - 2 - 2*D) / W;
  int V = (Y - 2 - 2*D) / W;
  int * mask  = (int *) xmalloc( U * V * sizeof(int) );
  int * reg_x = (int *) xmalloc( U * V * sizeof(int) );
  int * reg_y = (int *) xmalloc( U * V * sizeof(int) );
  int reg_n;
  double * angle;
  int i,j,k,l,n,u,v,uu,vv;
  double log_nfa;
  int region_label = 1;  /* CHH */

  /* compute 'optical flow' for each of the N pairs of images
     and keep the angle at each position for each pair */
  angle = offset_angle_per_pair(images,X,Y,N,W,D,T,input_mask); /* CHH */

  /* initialize cloud masks */
  for(i=0; i<X*Y; i++) cloud_mask[i] = 0.0;

  /* search for clouds */
  for(l=0; l<P; l++)  /* loop on precision value used in the region growing */
    {
      for(i=0; i<U*V; i++) mask[i] = 0;  /* initialize current precision mask */

      for(u=0; u<U; u++)  /* loops on spatial coordinates */
      for(v=0; v<V; v++)
      for(i=0; i<N; i++)  /* loop on the image pairs */
      if( mask[u+v*U] == 0 && angle[u + v*U + i*U*V] != UNDEF )
        {
          double ref = angle[u + v*U + i*U*V];

          /* all image pairs must agree on the angle, otherwise continue */
          if( number_angle_agreements(angle,U,V,N,u,v,ref,p[l]) < N ) continue;

          /* initialize region */
          mask[u+v*U] = 1;
          reg_x[0] = u;
          reg_y[0] = v;
          reg_n = 1;

          /* number of gradient angle agreements */
          k = N - 1; /* the reference point must not be counted in k, thus -1 */

          /* grow region on neighbors */
          for(n=0; n<reg_n; n++)  /* loop on pixels on the region */
          for(j=0; j<4; j++)      /* loop on the four neighbors */
            {
              uu = reg_x[n] + neigh_u[j];  /* coordinates of the neighbor */
              vv = reg_y[n] + neigh_v[j];

              /* add the new point if it is inside the image domain,
                 it was not already used, and all the angles agree with ref */
              if( uu >= 0 && uu < U && vv >= 0 && vv < V && mask[uu+vv*U] == 0 )
                if( number_angle_agreements(angle,U,V,N,uu,vv,ref,p[l]) == N )
                  {
                    /* add new point to region */
                    mask[uu+vv*U] = 1;
                    reg_x[reg_n] = uu;
                    reg_y[reg_n] = vv;
                    ++reg_n;

                    /* update the number of gradient angle agreements */
                    k += N;  /* the condition is satisfied for all pairs */
                  }
            }

          /* NFA = Nt  *  P( observing k agreements in H0 up to precision p )
                 = (UV)^2 * N * P * b_(reg_n)  *  p^k

             U*V : possible seed pixels
             N   : the angle of each of the N pairs is a possible reference
             P   : number of angular tolerances used (each is one test)
             U*V : number of sizes for the polyomino regions
             b_n : number of polyominoes of size n

             b_n ~ B * tau^n / n    where B ~ 0.316915 and tau ~ 4.062570
                                    (see Jensen & Guttmann 2000)
          */
          log_nfa = 2.0 * log10( (double) U ) + 2.0 * log10( (double) V )
                  + log10( (double) N ) + log10( (double) P )
                  + log10(0.316915) + (double) reg_n * log10(4.062570)
                  - log10( (double) reg_n ) + (double) k * log10( p[l] );

          /* if not meaningful, remove the points from mask */
          if( log_nfa >= 0.0 )
            for(n=0; n<reg_n; n++)
              mask[ reg_x[n] + reg_y[n] * U ] = 0;

          /* if meaningful, add the region to the output cloud mask */
          if( log_nfa < 0.0 )
            {
              for(n=0; n<reg_n; n++)
                {
                  int x = 1 + D + reg_x[n] * W;
                  int y = 1 + D + reg_y[n] * W;
                  int xx,yy;
                  for(xx=x; xx < x+W; xx++)  /* loops on pixels corresponding */
                  for(yy=y; yy < y+W; yy++)  /* to the correlation window     */
                    cloud_mask[xx+yy*X] += 16;
                    /* cloud_mask[xx+yy*X] = (double) region_label;  /1* CHH 255.0; *1/ */
                }
              region_label = region_label + 1; /* CHH */
            }
        }
      /* CHH >>
      sprintf(filenameA, "/tmp/cloud_mask_precision_%d.tif", l);
      fprintf(stderr, "%s\n", filenameA);
      iio_write_image_double(cloud_mask, X, Y);
      CHH << */
    }

  for(int y=0; y<Y; y++)
  for(int x=0; x<X; x++)
  if(cloud_mask[x+y*X] > 0.0)
    cloud_mask[x+y*X] += 159.0; /* 159 * 6*16 = 255 */

  /* free memory */
  free( (void *) angle );
  free( (void *) reg_x );
  free( (void *) reg_y );
  free( (void *) mask  );
}
/*----------------------------------------------------------------------------*/
