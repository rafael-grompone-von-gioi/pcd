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
#ifndef PARALLAX_CLOUD_DETECTOR_HEADER
#define PARALLAX_CLOUD_DETECTOR_HEADER

/*----------------------------------------------------------------------------*/
#define PARALLAX_CLOUD_DETECTOR_VERSION "0.7.dev0 (December 17, 2021)"

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
 */
void parallax_cloud_detector( double * cloud_mask, double * images,
                              int X, int Y, int N, int W, int D, double T,
                              unsigned char * input_mask);

#endif /* !PARALLAX_CLOUD_DETECTOR_HEADER */
/*----------------------------------------------------------------------------*/
