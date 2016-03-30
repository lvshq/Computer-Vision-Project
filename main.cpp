#include "edgeExtraction.h"

int main() {
	CImg<unsigned char> img1 = CImg<unsigned char>("Dataset/1.bmp");
	CImg<unsigned char> _img1 = egdeExtract(img1);

	_img1.display("Image 1");

	CImg<unsigned char> img2 = CImg<unsigned char>("Dataset/2.bmp");
	CImg<unsigned char> _img2 = egdeExtract(img2);

	_img2.display("Image 2");

	CImg<unsigned char> img3 = CImg<unsigned char>("Dataset/3.bmp");
	CImg<unsigned char> _img3 = egdeExtract(img3);

	_img3.display("Image 3");

	CImg<unsigned char> img4 = CImg<unsigned char>("Dataset/4.bmp");
	CImg<unsigned char> _img4 = egdeExtract(img4);

	_img4.display("Image 4");

	CImg<unsigned char> img5 = CImg<unsigned char>("Dataset/5.bmp");
	CImg<unsigned char> _img5 = egdeExtract(img5);

	_img5.display("Image 5");

	CImg<unsigned char> img6 = CImg<unsigned char>("Dataset/6.bmp");
	CImg<unsigned char> _img6 = egdeExtract(img6);

	_img6.display("Image 6");

	return 0;
}
