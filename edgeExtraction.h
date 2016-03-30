#ifndef EDGEEXTRACTION_H
#define EDGEEXTRACTION_H

#include "CImg.h"
#include<cmath>
#include<algorithm>
#include<iostream>
using namespace cimg_library;

/*
输入图像：
普通 A4 打印纸，上面可能有手写笔记或者打印内容，但是拍照时可能角度不正。
输出：
1. 图像的边缘
2. 计算 A4 纸边缘的各直线方程
3. 提取 A4 纸的四个角点

算法步骤：
1. 将图片缩小从而减少处理时间。
2. 把彩色图像灰度化，以便计算梯度。
3. 进行高斯模糊。
4. 计算灰度图的梯度。
5. 把梯度高的坐标映射到Hough图坐标上。
6. 提取Hough图峰值。
7. 计算直线方程和角点。
*/

const double PI = 3.1415926535898;

//用来表示极坐标参数空间里的数据
//其中角度angle指直线与x轴之间的夹角，r为原点到直线的垂直距离。
struct Para {
	int count;
	int angle;
	int r;
};

//极坐标空间里的点
struct Point {
	int angle;
	int r;
};

/*
用于缩小图片尺寸从而加快图片处理时间。在边缘提取结束后再恢复图片大小。
参数times表示改变的倍数，如0.5表示缩小成原来的一半。
*/
CImg<unsigned char> scaling(CImg<unsigned char> srcImg, float times) {
	//倍数非正数为非法
	if (times <= 0) return srcImg;

	//获得宽、高
	int width = srcImg.width()*times;
	int height = srcImg.height()*times;
	int srcWidth = srcImg.width();
	int srcHeight = srcImg.height();

	CImg<unsigned char> dstImg(width, height, 1, 3);  //目标图像dst：z值为1，3通道（RGB）

	//用双线性插值法缩放
	for (int i = 0; i < width; i++) {
		int x = i / times;
		double dx = i / times - x;
		for (int j = 0; j < height; j++) {
			int y = j / times;
			double dy = j / times - y;
			if (x < 0 || x + 1 >= srcWidth || y < 0 || y + 1 >= srcHeight) {  // 越界部分的像素值直接赋值为（0,0,0）
				dstImg(i, j, 0, 0) = 0;
				dstImg(i, j, 0, 1) = 0;
				dstImg(i, j, 0, 2) = 0;
			}
			else {
				//分别对三个通道使用双线性插值法
				dstImg(i, j, 0, 0) = (1 - dx) * (1 - dy) * srcImg(x, y, 0, 0) + (1 - dx) * dy * srcImg(x, y + 1, 0, 0) +
						dx * (1 - dy) * srcImg(x + 1, y, 0, 0) + dx * dy * srcImg(x + 1, y + 1, 0, 0);
				dstImg(i, j, 0, 1) = (1 - dx) * (1 - dy) * srcImg(x, y, 0, 1) + (1 - dx) * dy * srcImg(x, y + 1, 0, 1) +
						dx * (1 - dy) * srcImg(x + 1, y, 0, 1) + dx * dy * srcImg(x + 1, y + 1, 0, 1);
				dstImg(i, j, 0, 2) = (1 - dx) * (1 - dy) * srcImg(x, y, 0, 2) + (1 - dx) * dy * srcImg(x, y + 1, 0, 2) +
						dx * (1 - dy) * srcImg(x + 1, y, 0, 2) + dx * dy * srcImg(x + 1, y + 1, 0, 2);
			}
		}
	}
	return dstImg;
}

//用sobel算子计算梯度
CImg<unsigned char> sobel(CImg<unsigned char> src) {
	//调用CImg库中的函数进行高斯模糊
	src.blur(3);

	//两个方向的sobel filter
	float mask1[3][3] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };
	float mask2[3][3] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };

	//宽、高
	int width = src.width();
	int height = src.height();

	CImg<unsigned char> dst = CImg<unsigned char>(width, height, 1, 3);

	//x方向上的梯度，二维数组
	float **gradientX = new float*[width];
	//y方向上的梯度，二维数组
	float **gradientY = new float*[width];
	//整体梯度，二维数组
	int **gradient = new int*[width];
	//用于存放原来像素的灰度值
	int **origin = new int*[width];
	//初始化
	for (int i = 0; i < width; i++) {
		gradientX[i] = new float[height];
		gradientY[i] = new float[height];
		gradient[i] = new int[height];
		origin[i] = new int[height];
	}

	//计算灰度值
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			origin[i][j] = 0.299*src(i, j, 0, 0) + 0.587*src(i, j, 0, 1) + 0.114*src(i, j, 0, 2);
		}
	}

	// 用sobel算子滤波。注意下标不要越界。
	int sum = 0;
	int count = 0;
	for (int i = 1; i < width - 1; i++) {
		for (int j = 1; j < height - 1; j++) {
			double sum1 = 0;
			double sum2 = 0;
			for (int p = i - 1; p <= i + 1; p++) {
				for (int q = j - 1; q <= j + 1; q++) {
					sum1 += origin[p][q] * mask1[p - i + 1][q - j + 1];
					sum2 += origin[p][q] * mask2[p - i + 1][q - j + 1];
				}
			}
			gradientX[i][j] = abs(sum1);
			gradientY[i][j] = abs(sum2);
			gradient[i][j] = gradientX[i][j] + gradientY[i][j];  // 总梯度等于x与y方向梯度之和

			if (gradient[i][j] > 0) {  // 同时计算非零灰度值的平均值，用于计算阈值
				sum += gradient[i][j];
				count++;
			}
		}
	}

	int mean = sum / count;
	//手工试验不同阈值效果后，将阈值设为平均值的2倍
	int thresh = 2 * mean;

	//用阈值将梯度二值化
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			if (gradient[i][j] >= thresh) {
				gradient[i][j] = 255;
			}
			else {
				gradient[i][j] = 0;
			}
			//用梯度填充图片
			dst(i, j, 0, 0) = gradient[i][j];
			dst(i, j, 0, 1) = gradient[i][j];
			dst(i, j, 0, 2) = gradient[i][j];
		}
	}
	dst.display("Sobel Effect");
	return dst;
}


/*
进行hough变换：参数src是缩小后的图，source是缩小前的图，scale是放缩比例
*/
CImg<unsigned char> hough(CImg<unsigned char> src, CImg<unsigned char> source, int scale) {
	//srcWidth、srcHeight表示src的宽和高
	int srcWidth = src.width();
	int srcHeight = src.height();

	//width、height表示参数空间的宽和高
	//将角度分为500份（手工选定，在计算速度与角度的精度之间取折中），距离的取值范围是[0, sqrt(srcWidth * srcWidth + srcHeight * srcHeight)]
	int width = 500;
	int parallelBottomValue = 50;
	int parallelTopValue = width - 50;
	int height = sqrt(srcWidth * srcWidth + srcHeight * srcHeight);

	//用于生成hough空间（极坐标形式）的图片
	CImg<unsigned char> dst(width, height, 1, 3);
	dst.fill(0);

	//用于存放hough空间的点，并统计极大值。每个点包含参数r、angle（极坐标），和经过的直线个数count
	Para *houghSpace = new Para[width*height];

	//以原图中心为原点
	int centerX = srcWidth / 2;
	int centerY = srcHeight / 2;

	//data存放hough空间的值，进行初始化
	int **data = new int*[height];
	for (int i = 0; i < height; i++) {
		data[i] = new int[width];
		for (int j = 0; j < width; j++)
			data[i][j] = 0;
	}

	/*
	遍历原图src（已经过缩小）中的每个像素，对于每一个点，对应到hough空间有一条直线，该直线在hough空间经过的点对应的data值加一。
	方法：对每个取值范围内的角度angle（对应极坐标的横轴，即遍历width的每个取值），都计算相应的r，然后 data[r*width + angle]++。
	坐标转换时以原图中心为原点。
	*/
	double Theta;
	int newX, newY;
	for (int i = 0; i < srcWidth; i++) {
		for (int j = 0; j < srcHeight; j++) {
			if (src(i, j) == 255) { // 二值化后，只取大于阈值的点
				for (int k = 0; k < width; k++) {
					// 直角坐标与极坐标对应转换
					// x被ρcosθ代替，y被ρsinθ代替。x = i-centerX, y = j - centerY. 
					Theta = k / float(width) * PI; // 转换为对应角度
					newX = i - centerX;
					newY = j - centerY;
					int r = newX * cos(Theta) + newY * sin(Theta) + height / 2; // 加height/2为了保证r大于0
					data[r][k]++;
				}
			}
		}
	}

	//记录极大值的个数
	int total = 0;

	//data中的极大值即可能取的直线。将他们放入houghSpace，并统计个数total
	for (int r = 1; r < height - 1; r++) {
		for (int k = 1; k < width - 1; k++) {
			if (data[r][k] >= data[r][k+1]	&& data[r][k] >= data[r][k-1]
				&& data[r][k] > data[r+1][k] && data[r][k] >= data[r-1][k]) { // 比周围四个点大
				houghSpace[total].count = data[r][k];
				houghSpace[total].angle = k;
				houghSpace[total].r = r;
				total++;
			}
		}
	}

	//find表示当前找到的平行线对数
	int find = 0;
	Point p[10];

	//临时存放四个满足条件的角点
	Point temp[4];
	//表示p中的点是否与p中的其他线平行，初始化为false
	bool parallel[10];
	for (int i = 0; i < 10; i++)
		parallel[i] = false;

	//找到两对平行线即可
	for (int k = 0; find < 2 && k < 10; k++) {
		//最大投票值
		int m = 0;
		//最大值对应的坐标
		int n = 0;
		//找到最大值对应的坐标
		for (int i = 0; i < total; i++) {
			if (m < houghSpace[i].count) {
				m = houghSpace[i].count;
				n = i;
			}
		}

		//去除与已记录线接近重合的线
		bool overlap = true;
		for (int i = 0; i < k; i++) {
			if ((houghSpace[n].angle - p[i].angle) * (houghSpace[n].angle - p[i].angle) +
				(houghSpace[n].r - p[i].r) * (houghSpace[n].r - p[i].r) / 10 <= 100) { // 此处100为自主选定的阈值
				overlap = false;
				break;
			}
		}
		if (!overlap) {
			k--;
			houghSpace[n].count = 0;
			continue;
		}

		//记下非重合的线
		p[k].angle = houghSpace[n].angle;
		p[k].r = houghSpace[n].r;

		//判断与之前的线是不是平行线
		
		for (int i = 0; i < k; i++) {
			//平行条件。
			if (abs(p[i].angle - houghSpace[n].angle) < parallelBottomValue || abs(p[i].angle - houghSpace[n].angle) > parallelTopValue) {
				if (parallel[i] == false) {
					//若平行，将这些线放入temp中
					temp[2 * find].angle = p[i].angle;
					temp[2 * find].r = p[i].r;

					temp[2 * find + 1].angle = p[k].angle;
					temp[2 * find + 1].r = p[k].r;

					parallel[i] = true;
					parallel[k] = true;
					find++;
					break;
				}
			}
		}

		//将记录的线置为0，便于找到下一个最大值
		houghSpace[n].count = 0;
	}

	//将临时变量，即四个满足条件的点赋给p
	for (int i = 0; i < 4; i++) {
		p[i].angle = temp[i].angle;
		p[i].r = temp[i].r;
	}

	//记录转换成直角坐标后的四个顶点， xy[i][0]、xy[i][1]分别表示x、y坐标
	int xy[4][2];
	//记录点（xy[i][0]，xy[i][1]）是直线P[ flag[4][0] ]、P[ flag[4][1] ]的交点
	int flag[4][2];

	int cur = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			//只计算不平行的线的交点
			if (abs(p[i].angle - p[j].angle) < parallelBottomValue || abs(p[i].angle - p[j].angle) > parallelTopValue) continue;
			double x, y;

			if (p[i].angle == 0) { // 考虑角度为0的特殊情况
				x = (p[i].r - height / 2);
				y = ((p[j].r - height / 2) - x*cos(p[j].angle / float(width) * PI))
					/ sin(p[i].angle / float(width) * PI);
			}
			else if (p[j].angle == 0) {
				x = (p[j].r - height / 2);
				y = ((p[i].r - height / 2) - x*cos(p[i].angle / float(width) * PI))
					/ sin(p[i].angle / float(width) * PI);
			}
			else {
				//计算x值
				x = ((p[i].r - height / 2)*sin(p[j].angle / float(width) * PI) -
					(p[j].r - height / 2)*sin(p[i].angle / float(width) * PI)) /
					(cos(p[i].angle / float(width) * PI)*sin(p[j].angle / float(width) * PI) -
					cos(p[j].angle / float(width) * PI)*sin(p[i].angle / float(width) * PI));

				//计算y值
				y = ((p[i].r - height / 2) - x*cos(p[i].angle / float(width) * PI))
					/ sin(p[i].angle / float(width) * PI);
			}

			//x、y是以中心为原点的，需加上偏移值
			xy[cur][0] = x + centerX;
			xy[cur][1] = y + centerY;

			//记录点所在的两条直线
			flag[cur][0] = i;
			flag[cur][1] = j;

			//画点
			unsigned char c[3] = { 255, 0, 0 };
			// 此处处理的时缩小scale倍后的图片，因而应乘以scale
			source.draw_circle(xy[cur][0] * scale, xy[cur][1] * scale, 30, c, 2);
			cur++;
		}
	}

	//画四条直线，并输出四条直线方程
	std::cout << "四条边的直线方程：" << std::endl;
	unsigned char color[3] = { 255, 0, 0 };
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			//选两个能形成直线的点（避免画出对角线）
			if (flag[i][0] == flag[j][0] || flag[i][1] == flag[j][1]) {
				//为了便于观察，将线宽设为5（即同时画五个像素）
				source.draw_line(xy[i][0] * scale, xy[i][1] * scale, xy[j][0] * scale, xy[j][1] * scale, color, 1);
				source.draw_line(xy[i][0] * scale - 1, xy[i][1] * scale, xy[j][0] * scale - 1, xy[j][1] * scale, color, 1);
				source.draw_line(xy[i][0] * scale - 2, xy[i][1] * scale, xy[j][0] * scale - 2, xy[j][1] * scale, color, 1);
				source.draw_line(xy[i][0] * scale + 1, xy[i][1] * scale, xy[j][0] * scale + 1, xy[j][1] * scale, color, 1);
				source.draw_line(xy[i][0] * scale + 2, xy[i][1] * scale, xy[j][0] * scale + 2, xy[j][1] * scale, color, 1);

				//输出方程
				if (xy[i][0] == xy[j][0]) std::cout << "x = " << xy[i][0] * scale << std::endl;
				else std::cout << "y = " << (xy[i][1] - xy[j][1]) * 1.0 / (xy[i][0] - xy[j][0]) << " * x + "
					<< scale * (xy[i][1] - (xy[i][1] - xy[j][1]) * 1.0 / (xy[i][0] - xy[j][0]) * xy[i][0]) << std::endl;

			}
		}
	}

	//输出四个点
	std::cout << std::endl << "四个角点：" << std::endl;
	for (int i = 0; i < 4; i++) {
		std::cout << "角点" << i+1 << "：" << "( " << xy[i][0] * scale << " , " << xy[i][1] * scale << " )" << std::endl;
	}
	std::cout << std::endl;

	delete[]houghSpace;
	delete[]data;

	return source;

}

//提取边缘函数入口
CImg<unsigned char> egdeExtract(CImg<unsigned char> src) {
	//scale表示对原图缩小的倍数
	int scale = 2;
	//依次对图片进行缩小、sobel滤波、hough变换并还原回原来大小
	return hough(sobel(scaling(src, 1.0 / scale)), src, scale);
}

#endif
