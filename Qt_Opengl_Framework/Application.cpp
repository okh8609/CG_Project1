#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <omp.h>
#define SET0TO255(X)  ((X) > (255) ? (255) : ((X) < (0) ? (0) : (X)))

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene(void)
{

	ui_instance = Qt_Opengl_Framework::getInstance();

}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage(QString filePath)
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath)
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB(void)
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (!img_data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0; j < img_width; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0; i < 3; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			uint8_t gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			img_data[offset_rgba + rr] = rgb[offset_rgb + rr] / 32 * 32;
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg] / 32 * 32;
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb] / 32 * 32;
			img_data[offset_rgba + aa] = WHITE;
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
struct popuCell
{
	long long data; //最大可能到總像素
	uint8_t ir; //index red (0-31)
	uint8_t ig; //index green (0-31)
	uint8_t ib; //index blue (0-31)
};

typedef struct popuCell popuCell;

int popuCellCmp(const void *a, const void *b)
{
	return ((popuCell *)a)->data < ((popuCell *)b)->data ? 1 : -1;
}

void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	//宣告並初始化各個位置資訊
	popuCell ccc[32][32][32] = { 0 };
	for (size_t r = 0; r < 32; ++r)
	{
		for (size_t g = 0; g < 32; ++g)
		{
			for (size_t b = 0; b < 32; ++b)
			{
				popuCell* c = &ccc[r][g][b];
				c->ir = r, c->ig = g, c->ib = b;
			}
		}
	}

	//計算出現的數量
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;

			int r = rgb[offset_rgb + rr] / 8;
			int g = rgb[offset_rgb + gg] / 8;
			int b = rgb[offset_rgb + bb] / 8;
			++(ccc[r][g][b].data);
		}
	}

	//排序
	qsort(ccc, 32 * 32 * 32, sizeof(ccc[0][0][0]), popuCellCmp);

	//迭代整張圖
#pragma omp parallel for
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//找最近的顏色
			int indexOfColor = 0;
			int minColorError = INT_MAX;
			for (size_t k = 0; k < 256; ++k)
			{
				popuCell c = ((popuCell *)ccc)[k];
				int nowColorError = (pow(rgb[offset_rgb + rr] - c.ir * 8, 2) + pow(rgb[offset_rgb + gg] - c.ig * 8, 2) + pow(rgb[offset_rgb + bb] - c.ib * 8, 2));
				if (nowColorError < minColorError)
				{
					minColorError = nowColorError;
					indexOfColor = k;
				}
			}
			popuCell c = ((popuCell *)ccc)[indexOfColor];
			img_data[offset_rgba + rr] = c.ir * 8;
			img_data[offset_rgba + gg] = c.ig * 8;
			img_data[offset_rgba + bb] = c.ib * 8;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
// 二值化
void Application::Dither_Threshold()
{
	const uint8_t threshold = 127;

	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//先轉灰階
			uint8_t gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			//二值化
			gray = (gray > threshold) ? WHITE : BLACK;

			img_data[offset_rgba + rr] = gray;
			img_data[offset_rgba + gg] = gray;
			img_data[offset_rgba + bb] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();

	long long hstg[256] = { 0 }; //Histogram

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;

			//先轉灰階
			uint8_t gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			++hstg[gray]; //統計++

			rgb[offset_rgb + 0] = gray; //暫放灰階值在第0個通道!!
		}
	}

	//統計出threshold
	uint8_t threshold = 127;
	double acc = 0;
	for (int i = 0; i < 255; ++i)
	{
		acc += hstg[i];
	}
	threshold = 255 - (acc / img_height * img_width);

#pragma omp parallel for
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//二值化
			uint8_t gray = rgb[offset_rgb + 0];  //取出暫放的灰階值
			gray = (gray > threshold) ? WHITE : BLACK;

			img_data[offset_rgba + rr] = gray;
			img_data[offset_rgba + gg] = gray;
			img_data[offset_rgba + bb] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	//亂數種子
	srand(time(NULL));
	//亂數範圍/
	int ranMin = -0.2 * 256;
	int ranMax = 0.2 * 256;
	// 產生 [min , max) 亂數
	// int x = (max - min) * rand() / (RAND_MAX + 1.0) + min;

	const uint8_t threshold = 127; // xxxxxxxxx
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//先轉灰階 加上亂數
			int gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb]
				+ ((ranMax - ranMin) * rand() / (RAND_MAX + 1.0) + ranMin);
			//二值化
			gray = (gray > threshold) ? WHITE : BLACK;

			img_data[offset_rgba + rr] = gray;
			img_data[offset_rgba + gg] = gray;
			img_data[offset_rgba + bb] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	UINT8 I[4][4] =
	{ { 0.7059 * 255, 0.3529 * 255, 0.5882 * 255, 0.2353 * 255 },
	  { 0.0588 * 255, 0.9412 * 255, 0.8235 * 255, 0.4118 * 255 },
	  { 0.4706 * 255, 0.7647 * 255, 0.8824 * 255, 0.1176 * 255 },
	  { 0.1765 * 255, 0.5294 * 255, 0.2941 * 255, 0.6471 * 255 } };

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//先轉灰階
			uint8_t gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			//二值化
			gray = (gray > I[i % 4][j % 4]) ? WHITE : BLACK;

			img_data[offset_rgba + rr] = gray;
			img_data[offset_rgba + gg] = gray;
			img_data[offset_rgba + bb] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	const uint8_t threshold = 127;
	unsigned char *rgb = To_RGB();

	//先轉灰階
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			uint8_t gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			img_data[offset_rgba + rr] = gray;
			img_data[offset_rgba + gg] = gray;
			img_data[offset_rgba + bb] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	for (int i = 0; i < img_height - 1; i++)
	{
		for (int j = 1; j < img_width - 1; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			uint8_t oldpixel = img_data[offset_rgba + 0];
			uint8_t newpixel = (oldpixel > threshold) ? WHITE : BLACK;
			// img_data[offset_rgba + 0] = newpixel;
			double quant_error = oldpixel - newpixel;
			img_data[offset_rgba + 0 + 4] = SET0TO255(img_data[offset_rgba + 0 + 4] + quant_error * 7.0 / 16.0);
			img_data[offset_rgba + 0 + img_width * 4 - 4] = SET0TO255(img_data[offset_rgba + 0 + img_width * 4 - 4] + quant_error * 3.0 / 16.0);
			img_data[offset_rgba + 0 + img_width * 4] = SET0TO255(img_data[offset_rgba + 0 + img_width * 4] + quant_error * 5.0 / 16.0);
			img_data[offset_rgba + 0 + img_width * 4 + 4] = SET0TO255(img_data[offset_rgba + 0 + img_width * 4 + 4] + quant_error * 1.0 / 16.0);

			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = newpixel;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	const uint8_t threshold = 127;
	unsigned char *rgb = To_RGB();

#pragma omp parallel for
	for (int i = 0; i < img_height - 1; i++)
	{
		for (int j = 1; j < img_width - 1; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//要先去取樣建表?? 我是覺得結果不會差很多啦...

			img_data[offset_rgba + rr] = rgb[offset_rgb + rr];
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg];
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb];
			img_data[offset_rgba + aa] = WHITE;
		}
	}

#pragma omp parallel for
	for (int kk = 0; kk < 3; ++kk)
	{
		for (int i = 0; i < img_height - 1; i++)
		{
			for (int j = 1; j < img_width - 1; j++)
			{
				int offset_rgb = i * img_width * 3 + j * 3;
				int offset_rgba = i * img_width * 4 + j * 4;

				uint8_t oldpixel = img_data[offset_rgba + kk];
				uint8_t newpixel = (oldpixel > threshold) ? 255 : 0;
				img_data[offset_rgba + kk] = newpixel;
				double quant_error = oldpixel - newpixel;
				img_data[offset_rgba + kk + 4] = SET0TO255(img_data[offset_rgba + kk + 4] + quant_error * 7.0 / 16.0);
				img_data[offset_rgba + kk + img_width * 4 - 4] = SET0TO255(img_data[offset_rgba + kk + img_width * 4 - 4] + quant_error * 3.0 / 16.0);
				img_data[offset_rgba + kk + img_width * 4] = SET0TO255(img_data[offset_rgba + kk + img_width * 4] + quant_error * 5.0 / 16.0);
				img_data[offset_rgba + kk + img_width * 4 + 4] = SET0TO255(img_data[offset_rgba + kk + img_width * 4 + 4] + quant_error * 1.0 / 16.0);
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	INT8 ff[25] = {
			1,1,1,1,1,
			1,1,1,1,1,
			1,1,1,1,1,
			1,1,1,1,1,
			1,1,1,1,1 };
	int fff = 25;

	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 y = -2; y < 3; y++)
			{
				for (INT8 x = -2; x < 3; x++)
				{
					if (x + j >= 0 && y + i >= 0 && x + j < img_width && y + i < img_height)
					{
						rrr += rgb[offset_rgb + rr + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						ggg += rgb[offset_rgb + gg + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						bbb += rgb[offset_rgb + bb + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
					}
				}
			}

			img_data[offset_rgba + rr] = rrr / fff;
			img_data[offset_rgba + gg] = ggg / fff;
			img_data[offset_rgba + bb] = bbb / fff;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	INT8 ff[25] = {
			1,2,3,2,1,
			2,4,6,4,2,
			3,6,9,6,3,
			2,4,6,4,2,
			1,2,3,2,1 };
	int fff = 81;

	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 y = -2; y < 3; y++)
			{
				for (INT8 x = -2; x < 3; x++)
				{
					if (x + j >= 0 && y + i >= 0 && x + j < img_width && y + i < img_height)
					{
						rrr += rgb[offset_rgb + rr + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						ggg += rgb[offset_rgb + gg + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						bbb += rgb[offset_rgb + bb + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
					}
				}
			}

			img_data[offset_rgba + rr] = rrr / fff;
			img_data[offset_rgba + gg] = ggg / fff;
			img_data[offset_rgba + bb] = bbb / fff;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	INT8 ff[25] = {
		1,  4 ,  6,  4 , 1,
		4, 16 , 24, 16 , 4,
		6, 24 , 36, 24 , 6,
		4, 16 , 24, 16 , 4,
		1,  4 ,  6,  4 , 1 };
	int fff = 256;

	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 y = -2; y < 3; y++)
			{
				for (INT8 x = -2; x < 3; x++)
				{
					if (x + j >= 0 && y + i >= 0 && x + j < img_width && y + i < img_height)
					{
						rrr += rgb[offset_rgb + rr + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						ggg += rgb[offset_rgb + gg + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						bbb += rgb[offset_rgb + bb + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
					}
				}
			}

			img_data[offset_rgba + rr] = rrr / fff;
			img_data[offset_rgba + gg] = ggg / fff;
			img_data[offset_rgba + bb] = bbb / fff;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N(unsigned int N)
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	INT8 ff[25] = {
			1,  4 ,  6,  4 , 1,
			4, 16 , 24, 16 , 4,
			6, 24 , 36, 24 , 6,
			4, 16 , 24, 16 , 4,
			1,  4 ,  6,  4 , 1 };
	int fff = 256;

	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 y = -2; y < 3; y++)
			{
				for (INT8 x = -2; x < 3; x++)
				{
					if (x + j >= 0 && y + i >= 0 && x + j < img_width && y + i < img_height)
					{
						rrr += rgb[offset_rgb + rr + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						ggg += rgb[offset_rgb + gg + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						bbb += rgb[offset_rgb + bb + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
					}
				}
			}

			img_data[offset_rgba + rr] = SET0TO255(rgb[offset_rgb + rr] - (rrr / fff));
			img_data[offset_rgba + gg] = SET0TO255(rgb[offset_rgb + gg] - (ggg / fff));
			img_data[offset_rgba + bb] = SET0TO255(rgb[offset_rgb + bb] - (bbb / fff));
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	INT8 ff[25] = {
			1,  4 ,  6,  4 , 1,
			4, 16 , 24, 16 , 4,
			6, 24 , 36, 24 , 6,
			4, 16 , 24, 16 , 4,
			1,  4 ,  6,  4 , 1 };
	int fff = 256;

	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 y = -2; y < 3; y++)
			{
				for (INT8 x = -2; x < 3; x++)
				{
					if (x + j >= 0 && y + i >= 0 && x + j < img_width && y + i < img_height)
					{
						rrr += rgb[offset_rgb + rr + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						ggg += rgb[offset_rgb + gg + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
						bbb += rgb[offset_rgb + bb + y * img_width * 3 + x * 3] * ff[(y + 2) * 5 + (x + 2)];
					}
				}
			}

			img_data[offset_rgba + rr] = SET0TO255(2 * rgb[offset_rgb + rr] - (rrr / fff));
			img_data[offset_rgba + gg] = SET0TO255(2 * rgb[offset_rgb + gg] - (ggg / fff));
			img_data[offset_rgba + bb] = SET0TO255(2 * rgb[offset_rgb + bb] - (bbb / fff));

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	int new_width = img_width / 2;
	int new_height = img_height / 2;
	uint8_t *new_data = new uint8_t[new_width * new_height * 4];

	for (int i = 0; i < new_height; ++i)
	{
		for (int j = 0; j < new_width; ++j)
		{
			int offset_new = i * new_width * 4 + j * 4;
			int offset_rgba = (i * 2) * img_width * 4 + (j * 2) * 4;

			new_data[offset_new + rr] = img_data[offset_rgba + rr];
			new_data[offset_new + gg] = img_data[offset_rgba + gg];
			new_data[offset_new + bb] = img_data[offset_rgba + bb];
			new_data[offset_new + aa] = img_data[offset_rgba + aa];
		}
	}

	//delete[] img_data;

	img_data = new_data;
	img_width = new_width;
	img_height = new_height;

	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);

	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	int new_width = img_width * 2;
	int new_height = img_height * 2;
	uint8_t *new_data = new uint8_t[new_width * new_height * 4];

	//效能可能比較差 多了蠻多次重複的記憶體存取的 不過速度尚可
	for (int i = 0; i < new_height; ++i)
	{
		for (int j = 0; j < new_width; ++j)
		{
			int offset_new = i * new_width * 4 + j * 4;
			int offset_rgba = (i / 2) * img_width * 4 + (j / 2) * 4;

			new_data[offset_new + rr] = img_data[offset_rgba + rr];
			new_data[offset_new + gg] = img_data[offset_rgba + gg];
			new_data[offset_new + bb] = img_data[offset_rgba + bb];
			new_data[offset_new + aa] = img_data[offset_rgba + aa];
		}
	}

	//delete[] img_data;

	img_data = new_data;
	img_width = new_width;
	img_height = new_height;

	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);

	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize(float scale)
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate(float angleDegrees)
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge(QString filePath)
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image(int tMethod)
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::NPR_Paint_Layer(unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke(const Stroke& s)
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++)
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++)
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height))
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared)
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				}
				else if (dist_squared == radius_squared + 1)
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
	radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}





