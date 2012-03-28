// InPaintPort.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#if defined WIN32
#if defined _DEBUG
#pragma comment(lib,"opencv_core231d.lib")
#pragma comment(lib,"opencv_imgproc231d.lib")
#pragma comment(lib,"opencv_highgui231d.lib")
#pragma comment(lib, "opencv_features2d231d.lib")
#pragma comment(lib, "opencv_flann231d.lib")
#pragma comment(lib, "opencv_calib3d231d.lib")
#else
#pragma comment(lib,"opencv_core231.lib")
#pragma comment(lib,"opencv_imgproc231.lib")
#pragma comment(lib,"opencv_highgui231.lib")
#pragma comment(lib, "opencv_features2d231.lib")
#pragma comment(lib, "opencv_flann231.lib")
#pragma comment(lib, "opencv_calib3d231.lib")
#endif

#pragma warning(disable: 4251)
#pragma warning(disable: 4996)
#endif


#include <opencv2/opencv.hpp>

#include "Inpaint.h"
using namespace cv;

using namespace std;

void help()
{
	cout << "\nCool inpainging demo. Inpainting repairs damage to images by floodfilling the damage \n"
		<< "with surrounding image areas.\n"
		"Using OpenCV version %s\n" << CV_VERSION << "\n"
		"Usage:\n"
		"./inpaint [image_name -- Default fruits.jpg]\n" << endl;

	cout << "Hot keys: \n"
		"\tESC - quit the program\n"
		"\tr - restore the original image\n"
		"\ti or SPACE - run inpainting algorithm\n"
		"\t\t(before running it, paint something on the image)\n" << endl;
}

Mat img, inpaintMask;
Point prevPt(-1,-1);

void onMouse( int event, int x, int y, int flags, void* )
{
	if( event == CV_EVENT_LBUTTONUP || !(flags & CV_EVENT_FLAG_LBUTTON) )
		prevPt = Point(-1,-1);
	else if( event == CV_EVENT_LBUTTONDOWN )
		prevPt = Point(x,y);
	else if( event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON) )
	{
		Point pt(x,y);
		if( prevPt.x < 0 )
			prevPt = pt;

		int lineWidth = 20;
		line( inpaintMask, prevPt, pt, Scalar::all(255), lineWidth, 8, 0 );
		line( img, prevPt, pt, Scalar::all(255), lineWidth, 8, 0 );
		prevPt = pt;
		imshow("image", img);
	}
}


int main( int argc, char** argv )
{
	char* filename = argc >= 2 ? argv[1] : (char*)"tarja.bmp";
	Mat img0 = imread(filename, -1);
	if(img0.empty())
	{
		cout << "Couldn't open the image " << filename << ". Usage: inpaint <image_name>\n" << endl;
		return 0;
	}

	help();

	namedWindow( "image", 1 );

	img = img0.clone();
	inpaintMask = Mat::zeros(img.size(), CV_8U);

	imshow("image", img);
	setMouseCallback( "image", onMouse, 0 );

	for(;;)
	{
		char c = (char)waitKey();

		if( c == 27 )
			break;

		if( c == 'r' )
		{
			inpaintMask = Scalar::all(0);
			img0.copyTo(img);
			imshow("image", img);
		}

		if( c == 'i' || c == ' ' )
		{
			Mat src = img;
			Mat inpainted;

			Inpaint(src.data, inpaintMask.data, src.rows, src.cols, 3);
			imshow("inpainted image", src);
		}
	}

	return 0;
}


