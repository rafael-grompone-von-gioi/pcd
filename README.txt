Cloud detector by parallax between the channels of push-broom satellite images
==============================================================================

Version 0.7 - April 15, 2020
by rafael grompone von gioi <grompone@gmail.com>


Files
-----

README.txt                 - This file
COPYING                    - GNU AFFERO GENERAL PUBLIC LICENSE Version 3
Makefile                   - Compilation instructions for 'make'
parallax_cloud_detector.c  - ANSI C89 code for the cloud detector
parallax_cloud_detector.h  - ANSI C89 header for the cloud detector
main.c                     - Command line interface for Devernay, ANSI C89
iio.c                      - Input/Output functions for command interface
iio.h                      - Input/Output functions for command interface
cR.tif                     - Red channel of the sample satellite image
cG.tif                     - Green channel of the sample satellite image
cB.tif                     - Blue channel of the sample satellite image


Requirements
------------

This program requires a C language compiler.  By default, the libraries
libtiff, libjpeg, and libpng are used for image input/output.

In GNU/Linux these libraries may be installed with a command similar to the
following:

  apt install libtiff-dev libpng-dev libjpeg-dev

The exact command may depend on the particular distribution used; the exact
package names may also vary depending on the distribution.

On MacOS systems using Homebrew (https://brew.sh), these libraries may be
installed with the following command:

  brew install libtiff libpng libjpeg

See below how to compile without installing these libraries.


Compiling
---------

The compiling instruction is just

  make

from the directory where the source codes and the Makefile are located.


Compiling without libtiff, libpng, and libjpeg
----------------------------------------------

This program can be compiled and used without installing any library.  In this
case, only image file formats built-in in IIO are available; this includes:
PGM, ASC, NPY, etc.

In this case, two files must be modified before compiling.  On the one hand,
the following lines must be commented in iio.c:

  #define I_CAN_HAS_LIBPNG
  #define I_CAN_HAS_LIBJPEG
  #define I_CAN_HAS_LIBTIFF

On the other hand, those three libraries (-lpng -ltiff -ljpeg) must be removed
from the compilation command, resulting in:

  cc -O3 -o parallax_cloud_detector main.c parallax_cloud_detector.c iio.c -lm


Test
----

To verify a correct compilation you can apply the algorithm to the test images
provided with:

  make test


Running
-------

A typical execution is:

  ./parallax_cloud_detector 10 20 cR.tif cG.tif cG.tif cB.tif output.png

where cR.tif, cG.tif and cG.tif are the three channels of a push-broom
satellite image, and 10 and 20 are two parameters: the size of the correlation
window and the maximal cloud displacement in pixels, respectively.  The output
cloud mask is written to output.png and it contains 0 values in the background
and 255 in the regions detected as clouds.

The input and output can be in any image format handled by Enric Meinhardt's
IIO library (https://github.com/mnhrdt/iio).


Copyright and License
---------------------

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
