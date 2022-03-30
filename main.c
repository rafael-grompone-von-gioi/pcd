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
#include "parallax_cloud_detector.h"
#include "iio.h"

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
/* print usage and exit
 */
static void usage(void)
{
  fprintf(stderr,"error: invalid input\n\n");
  fprintf(stderr,"Parallax Cloud Detector %s\n",
                 PARALLAX_CLOUD_DETECTOR_VERSION);
  fprintf(stderr,"Copyright (c) 2019-2020 rafael grompone von gioi");
  fprintf(stderr," <grompone@gmail.com>\n");
  fprintf(stderr,"Cloud detector by parallax between the channels of");
  fprintf(stderr," push-broom satellite images.\n\n");
  fprintf(stderr,"usage: parallax_cloud_detector W D T <input_mask> <img1A> ");
  fprintf(stderr,"<img1B> [...imgNA imgNB] <cloud_mask>\n");
  fprintf(stderr,"       W              : size of the correlation window\n");
  fprintf(stderr,"       D              : maximal cloud displacement in");
  fprintf(stderr," pixels\n");
  fprintf(stderr,"       T              : threshold on mimimal gradient norm ");
  fprintf(stderr,"in pixels\n");
  fprintf(stderr,"       img1A to imgNB : N pairs of single-channel images");
  fprintf(stderr," (all of size X x Y)\n");
  fprintf(stderr,"       cloud_mask     : output cloud mask (X x Y image,");
  fprintf(stderr," 255=cloud, 0=non-cloud)\n");
  fprintf(stderr,"       input_mask     : valid pixels mask\n\n");
  fprintf(stderr,"All the image pairs must be ordered so that the acquisition");
  fprintf(stderr," order is coherent\namong all pairs.  In other words, the");
  fprintf(stderr," observed parallax must be in the same\ndirection for a");
  fprintf(stderr," given object in all pairs.\n\n");
  fprintf(stderr,"The pairs should not be redundant.  In RGB channels, for");
  fprintf(stderr," example, the pairs\nred/green and green/blue are");
  fprintf(stderr," independent; however the pairs red/green,\ngreen/blue and");
  fprintf(stderr," red/blue include redundant information, which undermines\n");
  fprintf(stderr,"the false detection control mechanism.\n\n");
  fprintf(stderr,"examples: ./parallax_cloud_detector 5 20 1 mask.png ");
  fprintf(stderr,"imageA.tif imageB.tif out.png\n");
  fprintf(stderr,"          ./parallax_cloud_detector 20 20 0.1 mask.png ");
  fprintf(stderr,"R.tif G.tif G.tif B.tif out.png\n");

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  double * images;      /* array to contain the whole image stack */
  double * cloud_mask;  /* output cloud mask */
  int X,Y;              /* image domain size */
  int W;                /* the correlation window size is 2W+1 */
  int D;                /* maximal cloud parallax 'movement' in pixels */
  int N;                /* number of image pairs */
  double T;             /* threshold on gradient norm */
  int n,i;

  /* usage */
  /* if( argc < 6 || (argc%2) != 0 ) usage(); */
  if( argc < 8 || (argc%2) != 0 ) usage();

  /* CHH debug */
  /* for(n=0; n<argc; n++) fprintf(stderr,"argv[%d]=%s\n",n,argv[n]); */

  /* read parameters W and D */
  W = atoi(argv[1]);
  if( W < 2 ) error("W must be at least 2");
  D = atoi(argv[2]);
  if( D < 2 ) error("D must be at least 2");
  T = atof(argv[3]);

  /* read images */
  N = (argc - 6) / 2;  /* number of input image pars from number of arguments */
  for(n=0; n<2*N; n++) /* loop on the number of images to read */
    {
      /* read new image */
      int XX,YY,CC;
      double * new_img = iio_read_image_double_split(argv[n+5], &XX, &YY, &CC);

      /* check its shape */
      if( CC != 1 ) error("images must be single channel");
      if( n == 0 ) /* first image */
        {
          /* set image size */
          X = XX;
          Y = YY;

          /* get memory */
          images = (double *) xmalloc( X * Y * 2 * N * sizeof(double) );
        }
      else
        if( X != XX || Y != YY ) error("all images must be of equal size");

      /* copy the new image into the common buffer */
      for(i=0; i<X*Y; i++) images[i + n*X*Y] = new_img[i];

      /* free temporary buffer */
      free( (void *) new_img );
    }

  /* CHH read mask (255 for valid, 0 for invalid) */
  int XXX,YYY;
  unsigned char * input_mask = iio_read_image_uint8(argv[4], &XXX, &YYY);
  if(X != XXX || Y != YYY) error("input mask size does not match input bands");
  iio_write_image_uint8_vec((char*)"/tmp/in_mask_check.png",input_mask,XXX,YYY,1);
  fprintf(stderr, "X=%d Y=%d W=%d D=%d T=%f N=%d\n", X, Y, W, D, T, N);

  /* call cloud detector */
  cloud_mask = (double *) xmalloc( X * Y * sizeof(double) );
  parallax_cloud_detector(cloud_mask, images, X, Y, N, W, D, T, input_mask);

  /* write output */
  iio_write_image_double(argv[argc-1], cloud_mask, X, Y);

  /* free memory */
  free( (void *) images );
  free( (void *) cloud_mask );

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
