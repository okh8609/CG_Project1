#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <omp.h>

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * ��l�e���A�����Ntust.png����
// 
//============================================================================
void Application::createScene(void)
{

	ui_instance = Qt_Opengl_Framework::getInstance();

}

//****************************************************************************
//
// * ���}���w����
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
// * ��s�e��
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
// * �e����l��
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
// * �x�s����
// 
//============================================================================
void Application::saveImage(QString filePath)
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * �N���ɸ���ഫ��RGB��m���
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
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

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
	long long data; //�̤j�i����`����
	uint8_t ir; //index red
	uint8_t ig; //index green
	uint8_t ib; //index blue
};

typedef struct popuCell popuCell;

int popuCellCmp(const void *a, const void *b)
{
	return ((popuCell *)a)->data < ((popuCell *)b)->data ? 1 : -1;
}

void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	//�ŧi�ê�l�ƦU�Ӧ�m��T
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

	//�p��X�{���ƶq
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			// int offset_rgba = i * img_width * 4 + j * 4;

			int r = rgb[offset_rgb + rr] / 8;
			int g = rgb[offset_rgb + gg] / 8;
			int b = rgb[offset_rgb + bb] / 8;
			++(ccc[r][g][b].data);
		}
	}

	//�Ƨ�
	qsort(ccc, 32 * 32 * 32, sizeof(ccc[0][0][0]), popuCellCmp);

	//���N��i��
#pragma omp parallel for
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//��̪��C��
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
// �G�Ȥ�
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

			//����Ƕ�
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			//�G�Ȥ�
			gray = (gray > threshold) ? WHITE : BLACK;

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
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	uint8_t threshold = 127;
	double sum = 0;

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//����Ƕ�
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			sum += gray;

			img_data[offset_rgba + 0] = gray; //�ȩ�Ƕ��Ȧb�Ĥ@�ӳq�D
		}
	}

	threshold = sum / img_height / img_width;


#pragma omp parallel for
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//�G�Ȥ�
			unsigned char gray = img_data[offset_rgba + 0];  //���X�ȩ񪺦Ƕ���
			gray = (gray > threshold) ? WHITE : BLACK;

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
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	//�üƺؤl
	srand(time(NULL));
	//�üƽd��/
	int ranMin = -0.2 * 256;
	int ranMax = 0.2 * 256;
	// ���� [min , max) �ü�
	// int x = (max - min) * rand() / (RAND_MAX + 1.0) + min;

	const uint8_t threshold = 127; // xxxxxxxxx
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			//����Ƕ� �[�W�ü�
			int gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb]
				+ ((ranMax - ranMin) * rand() / (RAND_MAX + 1.0) + ranMin);
			//�G�Ȥ�
			gray = (gray > threshold) ? WHITE : BLACK;

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
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();



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
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////


//void Application::filtering(double filter[][5])
//{
//	unsigned char *rgb = this->To_RGB();
//
//
//
//	delete[] rgb;
//	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
//	renew();
//}
//
//void Application::filtering(double **filter, int n)
//{
//	unsigned char *rgb = this->To_RGB();
//
//
//
//	delete[] rgb;
//	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
//	renew();
//}


///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	unsigned char *rgb = To_RGB();

	for (int i = 2; i < img_height - 2; i++)
	{
		for (int j = 2; j < img_width - 2; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0;
			double ggg = 0;
			double bbb = 0;

			for (INT8 x = -2; x < 3; x++)
			{
				for (INT8 y = -2; y < 3; y++)
				{
					rrr += rgb[offset_rgb + rr + x * img_width * 3 + y * 3];
					ggg += rgb[offset_rgb + gg + x * img_width * 3 + y * 3];
					bbb += rgb[offset_rgb + bb + x * img_width * 3 + y * 3];
				}
			}

			img_data[offset_rgba + rr] = rrr / 25;
			img_data[offset_rgba + gg] = ggg / 25;
			img_data[offset_rgba + bb] = bbb / 25;
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

	for (int i = 2; i < img_height - 2; i++)
	{
		for (int j = 2; j < img_width - 2; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 x = -2; x < 3; x++)
			{
				for (INT8 y = -2; y < 3; y++)
				{
					rrr += rgb[offset_rgb + rr + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					ggg += rgb[offset_rgb + gg + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					bbb += rgb[offset_rgb + bb + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
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

	for (int i = 2; i < img_height - 2; i++)
	{
		for (int j = 2; j < img_width - 2; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 x = -2; x < 3; x++)
			{
				for (INT8 y = -2; y < 3; y++)
				{
					rrr += rgb[offset_rgb + rr + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					ggg += rgb[offset_rgb + gg + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					bbb += rgb[offset_rgb + bb + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
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

	for (int i = 2; i < img_height - 2; i++)
	{
		for (int j = 2; j < img_width - 2; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 x = -2; x < 3; x++)
			{
				for (INT8 y = -2; y < 3; y++)
				{
					rrr += rgb[offset_rgb + rr + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					ggg += rgb[offset_rgb + gg + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					bbb += rgb[offset_rgb + bb + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
				}
			}

			img_data[offset_rgba + rr] = rgb[offset_rgb + rr] - (rrr / fff);
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg] - (ggg / fff);
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb] - (bbb / fff);
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

	for (int i = 2; i < img_height - 2; i++)
	{
		for (int j = 2; j < img_width - 2; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			double rrr = 0, ggg = 0, bbb = 0;

			for (INT8 x = -2; x < 3; x++)
			{
				for (INT8 y = -2; y < 3; y++)
				{
					rrr += rgb[offset_rgb + rr + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					ggg += rgb[offset_rgb + gg + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
					bbb += rgb[offset_rgb + bb + x * img_width * 3 + y * 3] * ff[(x + 2) * 5 + (y + 2)];
				}
			}

			rrr = 2 * rgb[offset_rgb + rr] - (rrr / fff);
			ggg = 2 * rgb[offset_rgb + gg] - (ggg / fff);
			bbb = 2 * rgb[offset_rgb + bb] - (bbb / fff);

			img_data[offset_rgba + rr] = (rrr > 255) ? 255 : rrr;
			img_data[offset_rgba + gg] = (ggg > 255) ? 255 : ggg;
			img_data[offset_rgba + bb] = (bbb > 255) ? 255 : bbb;

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
	uchar *new_data = new uchar[new_width * new_height * 4];

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
	uchar *new_data = new uchar[new_width * new_height * 4];

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



