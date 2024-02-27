IIO=-DNDEBUG -std=c99 -lpng -ltiff -ljpeg
OPT= -O3

parallax_cloud_detector: main.c parallax_cloud_detector.c iio.c
	cc $(OPT) -o $@ $^ $(IIO) -lm

test: parallax_cloud_detector
	./parallax_cloud_detector  5 20 cR.tif cG.tif cG.tif cB.tif out_D05.png
	./parallax_cloud_detector 10 20 cR.tif cG.tif cG.tif cB.tif out_D10.png
	./parallax_cloud_detector 20 20 cR.tif cG.tif cG.tif cB.tif out_D20.png
	./parallax_cloud_detector 40 20 cR.tif cG.tif cG.tif cB.tif out_D40.png

clean:
	rm -f parallax_cloud_detector
	rm -f out_D05.png out_D10.png out_D20.png out_D40.png
