from ctypes import c_int, c_double, c_ubyte, CDLL
from numpy.ctypeslib import ndpointer
import numpy as np
import os.path as op


def parallax_cloud_detector(
        images, win_radius=10, max_offset=20, grad_thres=0.1, data_mask=None):
    """
    Rafael Grompone von Gioi's parallax cloud detector

    Args:
        images: (ndarray, or list of ndarrays) list of N pairs of images, of
                size 2*N x h x w, with h and w the height and width of the
                images. Pairs are made using consecutive images in the list.
                For example use [R, G, G, B] to compute a mask based on the
                pairs (R, G) and (G, B).
        win_radius: (int) the correlation window size is 2W+1
        max_offset: (int) maximal offset length to be evaluated
        grad_thres: (float) threshold on minimal gradient norm allowed when
                    computing gradient angle
        data_mask: (ndarray, shape=(h,w), dtype=uint8) mask indicating where to
                   look for clouds (255 where data, 0 elsewhere)

    Returns:
        mask: (boolean ndarray) cloud mask, True where clouds are detected
    """

    images = np.ascontiguousarray(images, dtype=np.float64)

    M, Y, X = images.shape  #  number of images, height, width
    N = int(M / 2)  # number of pairs

    if M % 2:
        raise ValueError('Detector requires a list of *pairs* of images')

    if win_radius > min(X, Y) -2 -2 * max_offset:
        raise ValueError(f'Image is too small wrt parameters win_radius and '
                         'max_offset. Detector requires win_radius <= '
                         'min(width, height) -2 -2*max_offset, and '
                         f'width={X}, heigth={Y}, win_radius={win_radius}, '
                         f'max_offset={max_offset}')

    if data_mask is None: data_mask = np.full((Y, X), 255, dtype=np.uint8)
    assert data_mask.shape == (Y, X), (data_mask.shape, (Y, X))
    data_mask = np.ascontiguousarray(data_mask, dtype=np.uint8)

    # load library
    here = op.dirname(op.abspath(__file__))
    libparallax_path = op.join(here, 'libparallax.so')
    try:
        libparallax = CDLL(libparallax_path)
    except OSError:
        print(f'Error loading library {libparallax_path}.\nDid you compile it? '
              'This library is required for the method grompone-disparity-v2.')
        raise

    # C function parallax_cloud_detector prototype
    libparallax.parallax_cloud_detector.argtypes = (
            ndpointer(dtype=c_double, shape=(Y, X)),        # *cloud_mask
            ndpointer(dtype=c_double, shape=(M, Y, X)),     # *images
            c_int, c_int, c_int, c_int, c_int, c_double,    # X, Y, N, W, D, T
            ndpointer(dtype=c_ubyte, shape=(Y, X)))         # *input_mask

    # allocate memory for output mask
    mask = np.empty((Y, X), dtype=np.float64, order='C')

    # apply cloud detector
    libparallax.parallax_cloud_detector(
        mask, images, X, Y, N, win_radius, max_offset, grad_thres, data_mask)

    return mask != 0


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4 or len(sys.argv) % 2 == 1:
        print(f'Usage: {sys.argv[0]} <img1A> <img1B> [...imgNA imgNB] '
              f'<cloud_mask>', file=sys.stderr)
        sys.exit(1)
    import rasterio as rio
    ims = []
    for path in sys.argv[1:-1]:
        with rio.open(path, 'r') as src:
            ims.append(src.read(1))
            transform = src.transform
            crs = src.crs
    mask = parallax_cloud_detector(ims)
    p = {'driver': 'GTiff', 'tiled': True, 'blockxsize': 256,
         'blockysize': 256, 'compress': 'deflate', 'predictor': 2, 'zlevel': 9,
         'count': 1, 'height': mask.shape[0], 'width': mask.shape[1],
         'dtype': np.uint8, 'transform': transform, 'crs': crs}
    with rio.open(sys.argv[-1], 'w', **p) as dst:
        dst.write(mask.astype(np.uint8) * 255, 1)
