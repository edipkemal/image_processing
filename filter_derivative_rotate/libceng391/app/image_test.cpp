//Edip Kemal Sardoğan
//240201026
#include <cstdlib>
#include <iostream>

#include "ceng391/image.hpp"

using namespace std;
using ceng391::Image;

int main(int argc, char** argv)
{	
	int test=9;
        
        if (argc != 2) {
                cerr << "Usage: " << argv[0] << " <filename>" << endl;
                return EXIT_FAILURE;
        }
        
        Image *img =  Image::from_png(argv[1], true);
        int filter_size;
	float sigma_x, sigma_y, theta;
	cout <<endl<< "*********"<<endl<<"Test cases will be overwritten the image."<<endl<<"Please reset the image(9) to see each test seperately."<<endl<<"*********"<<endl; 
        while(test!=0){

		switch(test){
			case(1):// Ex 1-a
				cout << "Enter gaussian filter size: ";
				cin >>filter_size;
				cout << "Enter sigma x value: ";
				cin >>sigma_x;
				
				img->smooth_x(filter_size,sigma_x);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(2):// Ex 1-b
				cout << "Enter gaussian filter size: ";
				cin >>filter_size;
				cout << "Enter sigma y value: ";
				cin >>sigma_y;
				
				img->smooth_y(filter_size,sigma_y);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(3):// Ex 1-c
				cout << "Enter gaussian filter size: ";
				cin >>filter_size;
				cout << "Enter sigma x value: ";
				cin >>sigma_x;
				cout << "Enter sigma y value: ";
				cin >>sigma_y;
				
				img->smooth(filter_size,sigma_x,sigma_y);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(4):// Ex 2-a
				img->derive_x(); //returns short type array
				img->xsave_png("/tmp/test_img.png");
				break;
			case(5):// Ex 2-b
				img->derive_y(); //returns short type array
				img->xsave_png("/tmp/test_img.png");
				break;
			case(6):// Ex 3-a
				cout << "Please enter rotation angle: ";
				cin >>theta;
				img->rotate(theta);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(7):// Ex 3-b
				break;
			case(8):// Ex 3-c
				cout << "Please enter rotation angle: ";
				cin >>theta;
				img->rotate_center(theta);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(9):// reset image
				img =  Image::from_png(argv[1], true);
				img->xsave_png("/tmp/test_img.png");
				break;
			case(0):
				break;
			default:
				cout << "Invalid input !"<< endl;
				break;
		}
	
		cout << endl << "User Manual For Test Case:" <<endl<<endl;
		cout << " 1 - Ex.1-a (smoot_x)" << endl;
		cout << " 2 - Ex.1-b (smoot_y)" << endl;
		cout << " 3 - Ex.1-c (smooth)" << endl;
		cout << " 4 - Ex.2-a (deriv_x)" << endl;
		cout << " 5 - Ex.2-b (deriv_y)" << endl;
		cout << " 6 - Ex.3-a (rotate)" << endl;
		cout << " 7 - Ex.3-b (bilinear_sampling)" << endl;
		cout << " 8 - Ex.3-c (rotate_center)" << endl;
		cout << " 9 - Reset the image" << endl;
		cout << " 0 - Exit" << endl;
		cout << endl << "Please enter the number of test: " << endl;
		
		cin >>test;	
	}
	
        delete img;
        return EXIT_SUCCESS;
}

/// Local Variables:
/// mode: c++
/// compile-command: "make -C ../build image-test"
/// End:
