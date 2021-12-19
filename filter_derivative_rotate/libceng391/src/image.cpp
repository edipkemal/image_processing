//Edip Kemal Sardoğan
//240201026
#include "ceng391/image.hpp"

#include <cstdlib>
#include <cstdio>
#include<fstream>
#include<iostream>
#include<cstring>
#include <stdexcept>
#include <math.h> 

#include <png.h>

#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif

#ifndef png_infopp_NULL
#  define png_infopp_NULL (png_infopp)NULL
#endif

#ifndef int_p_NULL
# define int_p_NULL (int*)NULL
#endif

using namespace std;

namespace ceng391 {

Image::Image(int width, int height, int n_channels, int step)
{
        m_width = width;
        m_height = height;
        m_n_channels = n_channels;
        int row_len = m_width * m_n_channels;
        if (step < row_len)
                step = row_len;
        m_step = step;

        m_data = new uchar[m_step * m_height];
}

Image::~Image()
{
        delete [] m_data;
}

void Image::reallocate(int width, int height, int n_channels)
{
        if (width  != this->m_width ||
            height != this->m_height ||
            n_channels != this->m_n_channels) {
                delete [] m_data;
                int step = width * n_channels;
                m_step = step;
                m_data = new uchar[height * m_step];
                m_width = width;
                m_height = height;
                m_n_channels = n_channels;
        }
}

void Image::set_rect_rgba(int x_tl, int y_tl, int width, int height,
                          uchar red, uchar green, uchar blue, uchar alpha)
{
        if (x_tl < 0) {
                width += x_tl;
                x_tl = 0;
        }

        if (y_tl < 0) {
                height += y_tl;
                y_tl = 0;
        }

        size_t length = min(width, m_width - x_tl);
        int y_max = min(y_tl + height, m_height);
        if (m_n_channels == 4) {
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        for (int x = x_tl; x < x_tl + (int)length; ++x) {
                                row_y[4*x]   = red;
                                row_y[4*x+1] = green;
                                row_y[4*x+2] = blue;
                                row_y[4*x+3] = alpha;
                        }
                }
        } else if (m_n_channels == 3) {
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        for (int x = x_tl; x < x_tl + (int)length; ++x) {
                                row_y[3*x]   = red;
                                row_y[3*x+1] = green;
                                row_y[3*x+2] = blue;
                        }
                }
        } else if (m_n_channels == 1) {
                int value = 0.3f * red + 0.59f * green + 0.11f * blue;
                if (value > 255)
                        value = 255;
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        memset(row_y, value, length);
                }
        } else {
                cerr << "Unexpected number of channels ("
                     << m_n_channels << ") in call to set_rect_xxx()" << endl;
                exit(EXIT_FAILURE);
        }
}

void Image::xsave_pnm(const  std::string& filename) const
{
        if (m_n_channels != 1) {
                cerr << "xsave_pnm() works only with grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        const string magic_head = "P5";
        ofstream fout;
        string extended_name = filename + ".pgm";
        fout.open(extended_name.c_str(), ios::out | ios::binary);
        fout << magic_head << "\n";
        fout << m_width << " " << m_height << " 255\n";
        for (int y = 0; y < m_height; ++y) {
                const uchar *row_data = this->data(y);
                fout.write(reinterpret_cast<const char*>(row_data),
                           m_width*sizeof(uchar));
        }
        fout.close();
}

void Image::xsave_png(const std::string &filename) const
{
        // We open the output file with C style IO since we are using libpng
        // C-API
        FILE *fout = fopen(filename.c_str(), "w+b");
        if (!fout) {
                cerr << "[ERROR][CENG391::Image] Failed open file for writing: " << filename << endl;
                exit(EXIT_FAILURE);
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                      0, 0, 0);
        if (!png_ptr) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG write structure!" << endl;
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG info structure!" << endl;
                png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG jump buffer!" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        int color_type;
        switch (this->m_n_channels) {
        case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
        case 3: color_type = PNG_COLOR_TYPE_RGB; break;
        case 4: color_type = PNG_COLOR_TYPE_RGBA; break;
        default:
                cerr << "[ERROR][CENG391::Image] Unsupported image type for saving as PNG!" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        png_init_io(png_ptr, fout);
        png_set_IHDR(png_ptr, info_ptr, this->m_width, this->m_height, 8,
                     color_type, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);

        png_bytepp row_pointers = (png_bytepp)malloc(this->m_height * sizeof(png_bytep));
        if (!row_pointers) {
                cerr << "[ERROR][CENG391::Image]Error creating PNG row pointers" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        for (png_int_32 k = 0; k < this->m_height; k++) {
                row_pointers[k] = (png_bytep)(this->data(k));
        }

        png_write_image(png_ptr, row_pointers);
        png_write_end(png_ptr, info_ptr);

        png_destroy_write_struct(&png_ptr, &info_ptr);
        free(row_pointers);
        fclose(fout);
}


void Image::xload_png(const std::string &filename, bool load_as_grayscale)
{
        // We open the output file with C style IO since we are using libpng
        // C-API
        FILE *fin = fopen(filename.c_str(), "r+b");
        if (!fin) {
                cerr << "[ERROR][CENG391::Image] Failed to open file for reading: " << filename << endl;
                exit(EXIT_FAILURE);
        }

        png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                                     NULL, NULL, NULL);
        if (!png_ptr) {
                cerr << "[ERROR][CENG391::Image] Could not create PNG read structure" << endl;
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
                cerr << "[ERROR][CENG391::Image] Could not create PNG info pointer" << endl;
                png_destroy_read_struct(&png_ptr, png_infopp_NULL,
                                        png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
                cerr << "[ERROR][CENG391::Image] Could not set jump point for reading PNG file" << endl;
                png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        png_init_io(png_ptr, fin);
        png_read_info(png_ptr, info_ptr);

        png_uint_32 stream_width, stream_height;
        int bit_depth, color_type, interlace_type;
        png_get_IHDR(png_ptr, info_ptr, &stream_width, &stream_height, &bit_depth, &color_type,
                     &interlace_type, int_p_NULL, int_p_NULL);

        png_set_strip_16(png_ptr);
        if (color_type == PNG_COLOR_TYPE_GA) {
                png_set_strip_alpha(png_ptr); /*(not recommended). */
        }

        png_set_packing(png_ptr);
        if (color_type == PNG_COLOR_TYPE_PALETTE) {
                png_set_palette_to_rgb(png_ptr);
        }

        png_set_expand(png_ptr);

        // Depending on the type of image in the file and the load_as_grayscale
        // flag, we determine the desired number of channels of the output
        // image.
        int desired_n_channels = 4;
        if (load_as_grayscale) {
                desired_n_channels = 1;
                png_set_rgb_to_gray_fixed(png_ptr, 1, 30000, 59000);
                png_set_strip_alpha(png_ptr); /*(not recommended). */
        } else {
                if (color_type == PNG_COLOR_TYPE_GRAY ||
                    color_type == PNG_COLOR_TYPE_GA) {
                        desired_n_channels = 1;
                }

                if(color_type == PNG_COLOR_TYPE_RGB ||
                   color_type == PNG_COLOR_TYPE_PALETTE) {
                        png_set_add_alpha(png_ptr, 255, PNG_FILLER_AFTER);
                }
        }

        // If the current image dimensions do not match the image to be loaded,
        // then reallocate with the desired dimensions.
        reallocate(stream_width, stream_height, desired_n_channels);

        png_bytepp row_pointers = (png_bytepp)malloc(this->m_height * sizeof(png_bytep));
        if (!row_pointers) {
                cerr << "[ERROR][CENG391::Image]Error creating PNG row pointers" << endl;
                png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }
        for (int k = 0; k < this->m_height; k++) {
                row_pointers[k] = (png_bytep)(this->data(k));
        }

        png_read_image(png_ptr, row_pointers);
        png_read_end(png_ptr, info_ptr);

        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

        free(row_pointers);

        fclose(fin);
}

void copy_row_to_buffer(uchar* buffer, const uchar* row, int n, int border_size)
{
        memcpy(reinterpret_cast<void *>(buffer + border_size), reinterpret_cast<const void *>(row), n*sizeof(*buffer));

        for (int i = 0; i < border_size; ++i) {
                buffer[i] = buffer[border_size];
                buffer[n-i+(2*border_size)-1] = buffer[n+border_size-1];
        }
}

void copy_column_to_buffer(uchar* buffer, const uchar* column, int n, int step, int border_size)
{
        for (int i = 0; i < n; ++i) {
                buffer[border_size + i] = column[i * step];
        }

        for (int i = 0; i < border_size; ++i) {
                buffer[i] = buffer[border_size];
                buffer[n-i+(2*border_size)-1] = buffer[n+border_size-1];
        }
}

void box_filter_buffer(int n, uchar* buffer, int filter_size)
{
        for(int i = 0; i < n; ++i) {
                float sum = 0;
                for (int j = 0; j < filter_size; ++j)
                        sum += buffer[i + j];
                sum /= filter_size;
                buffer[i] = sum;
        }
}

void Image::box_filter_x(int filter_size)
{
        if (m_n_channels != 1) {
                cerr << "Can only box filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int border_size = filter_size / 2;
        
        filter_size = 2*border_size + 1;

        int lbuffer = 2*border_size + m_width;
        uchar *buffer = new uchar[lbuffer];

        for (int y = 0; y < m_height; ++y) {
                copy_row_to_buffer(buffer, data(y), m_width, border_size);
                box_filter_buffer(m_width, buffer, filter_size);
                memcpy(reinterpret_cast<void *>(data(y)), reinterpret_cast<const void *>(buffer), m_width*sizeof(*buffer));
        }

        delete [] buffer;
}

void Image::box_filter_y(int filter_size)
{
        if (m_n_channels != 1) {
                cerr << "Can only box filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int border_size = filter_size / 2;
        filter_size = 2*border_size + 1;

        int lbuffer = 2*border_size + m_height;
        uchar *buffer = new uchar[lbuffer];

        for (int x = 0; x < m_width; ++x) {
                copy_column_to_buffer(buffer, m_data + x, m_height, m_step, border_size);
                box_filter_buffer(m_height, buffer, filter_size);
                for (int y = 0; y < m_height; ++y)
                        m_data[x + y*m_step] = buffer[y];
        }

        delete [] buffer;
}

void Image::box_filter(int filter_size)
{
        box_filter_x(filter_size);
        box_filter_y(filter_size);
}




void Image::smooth_buffer(int n, uchar* buffer, int filter_size, float* gaussian_filter)
{	
	
	int border_size=filter_size/2;
	
	uchar* temp_buffer=new uchar[n+(2*border_size)];
	
        for(int i = border_size; i < n+(2*border_size); i++) {
                int sum = 0;
                for (int j = -border_size; j < border_size+1; ++j)
                        sum += buffer[i + j]*gaussian_filter[j+border_size];
                temp_buffer[i-border_size] = sum;
        }
        for(int i=0;i<n+(2*border_size);i++){
        	buffer[i]=temp_buffer[i];
        }
}

void Image::smooth_x(int filter_size,float sigma_x)
{
        if (m_n_channels != 1) {
                cerr << "Can only gaussian filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        
        int border_size = filter_size / 2;

        float *gaussian_filter = new float[filter_size];
        float gauss_filt_sum=0;
        for (int i = -border_size; i<border_size+1;i++){
        	float a = exp((-0.5)*(pow(i,2)/pow(sigma_x,2)));
        	gaussian_filter[i+border_size] = a;
        	gauss_filt_sum+=a;
        }
	for (int i=0;i<filter_size;i++){
		gaussian_filter[i]/=gauss_filt_sum;
	}

        int lbuffer = 2*border_size + m_width;
        uchar *buffer = new uchar[lbuffer];
        
        for (int y = 0; y < m_height; ++y) {
                copy_row_to_buffer(buffer, data(y), m_width, border_size);
                smooth_buffer(m_width, buffer, filter_size, gaussian_filter);
                memcpy(reinterpret_cast<void *>(data(y)), reinterpret_cast<const void *>(buffer), m_width*sizeof(*buffer));
        }

        delete [] buffer;
}


void Image::smooth_y(int filter_size, float sigma_y)
{
        if (m_n_channels != 1) {
                cerr << "Can only gaussian filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int border_size = filter_size / 2;

        float *gaussian_filter = new float[filter_size];
        float gauss_filt_sum=0;
        for (int i = -border_size; i<border_size+1;i++){
        	float a = exp((-0.5)*(pow(i,2)/pow(sigma_y,2)));
        	gaussian_filter[i+border_size] = a;
        	gauss_filt_sum+=a;
        }

	for (int i=0;i<filter_size;i++){
		gaussian_filter[i]/=gauss_filt_sum;
	}
	
        int lbuffer = 2*border_size + m_height;
        uchar *buffer = new uchar[lbuffer];

        for (int x = 0; x < m_width; ++x) {
                copy_column_to_buffer(buffer, m_data + x, m_height, m_step, border_size);
                smooth_buffer(m_height, buffer, filter_size, gaussian_filter);
                for (int y = 0; y < m_height; ++y)
                        m_data[x + y*m_step] = buffer[y];
        }

        delete [] buffer;
}


void Image::smooth(int filter_size, float sigma_x, float sigma_y)
{
        smooth_x(filter_size,sigma_x);
        smooth_y(filter_size,sigma_y);
} 


void derive_buffer(int n, uchar* buffer, int filter_size, short* filter)
{	

	int border_size=filter_size/2;
	
	short* temp_buffer=new short[n+(2*border_size)];
	
        for(int i = 0; i < n+(border_size); i++) {
                float sum = 0;
                for (int j = 0; j < filter_size; ++j)
                        sum += buffer[i + j]*filter[j];
                if (sum>255) sum=255;
                if (sum<0) sum=0;
                temp_buffer[i-border_size] = sum;
        }
        for(int i=0;i<n+(2*border_size);i++){
        	buffer[i]=(uchar)temp_buffer[i];
        }
}

short* Image::derive_x(){

        if (m_n_channels != 1) {
                cerr << "Can only derive grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int filter_size=3;
        int border_size = filter_size / 2;
        
        int lbuffer = 2*border_size + m_width;
        uchar *buffer = new uchar[lbuffer];
        short *result = new short[lbuffer];
        short filter_x[] = {1, 2, 1};
        short filter_y[] = {-1, 0, 1};
        
        
        
        for (int y = 0; y < m_height; ++y) {
                copy_row_to_buffer(buffer, data(y), m_width, border_size);
                derive_buffer(m_width, buffer, filter_size,filter_x);
                memcpy(reinterpret_cast<void *>(data(y)), reinterpret_cast<const void *>(buffer), m_width*sizeof(*buffer));
        }
        for (int x = 0; x < m_width; ++x) {
                copy_column_to_buffer(buffer, m_data + x, m_height, m_step, border_size);
                derive_buffer(m_height, buffer, filter_size,filter_y);
		for (int y = 0; y < m_height; ++y)
                        m_data[x + y*m_step] = buffer[y];
        }
	for(int i=0;i<lbuffer;i++){
		result[i]=m_data[i];
	}

        delete [] buffer;
        return result;
}

short* Image::derive_y(){
        if (m_n_channels != 1) {
                cerr << "Can only derive grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int filter_size=3;
        int border_size = filter_size / 2;
        
        int lbuffer = 2*border_size + m_width;
        uchar *buffer = new uchar[lbuffer];
        short *result = new short[lbuffer];
        short filter_y[] = {1, 2, 1};
        short filter_x[] = {-1, 0, 1};
        
        for (int y = 0; y < m_height; ++y) {
                copy_row_to_buffer(buffer, data(y), m_width, border_size);
                derive_buffer(m_width, buffer, filter_size,filter_x);
                memcpy(reinterpret_cast<void *>(data(y)), reinterpret_cast<const void *>(buffer), m_width*sizeof(*buffer));
        }
        for (int x = 0; x < m_width; ++x) {
                copy_column_to_buffer(buffer, m_data + x, m_height, m_step, border_size);
                derive_buffer(m_height, buffer, filter_size,filter_y);
		for (int y = 0; y < m_height; ++y)
                        m_data[x + y*m_step] = buffer[y];
        }
        
	
	for(int i=0;i<lbuffer;i++){
		result[i]=m_data[i];
	}

	
        delete [] buffer;
        return result;
}


void Image::rotate(float theta){
	if (m_n_channels != 1) {
                cerr << "Can only rotate grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        
	float c=cos((theta / 180.0) * ((double) M_PI));
	float s=sin((theta / 180.0) * ((double) M_PI));
	
	uchar *new_data= new uchar[m_height*m_width];
	for(int i=0; i<m_height;i++){
		for(int j=0; j<m_width; j++){
			float new_x=c*j+s*i; // x' = cos(t)*x + sin(t)*y
			float new_y=(-s)*j+c*i; // y' = -sin(t)*x + cos(t)*y 
			
			int xx=round(new_x);
			int yy=round(new_y);
			
			if ((xx>m_width) | (yy>m_height) | (xx<0) | (yy<0))
				continue;
			
			new_data[i*m_width+j]=m_data[yy*m_width+xx];
		}	
	}
	for(int i=0;i<m_height*m_width;i++){
		m_data[i]=new_data[i];
	}
}
void Image::rotate_center(float theta){
	if (m_n_channels != 1) {
                cerr << "Can only rotate grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        
	float c=cos((theta / 180.0) * ((double) M_PI));
	float s=sin((theta / 180.0) * ((double) M_PI));
	
	uchar *new_data= new uchar[m_height*m_width];
	for(int i=0; i<m_height;i++){
		float ii=i-m_height/2;
		for(int j=0; j<m_width; j++){
			float jj=j-m_width/2;
			float new_x=c*(jj)+s*(ii); // x' = cos(t)*x + sin(t)*y
			float new_y=(-s)*(jj)+c*(ii); // y' = -sin(t)*x + cos(t)*y 
			
			new_x += m_width/2;
			new_y += m_height/2;
			
			int xx=round(new_x);
			int yy=round(new_y);
			
			uchar new_pixel;
			if ((xx>m_width) | (yy>m_height) | (xx<0) | (yy<0)){
				new_pixel=0;
			}else{
				new_pixel=m_data[yy*m_step+xx];
			}
			new_data[i*m_step+j]=new_pixel;
		}	
	}
	for(int i=0;i<m_height*m_width;i++){
		m_data[i]=new_data[i];
	}
}

}
