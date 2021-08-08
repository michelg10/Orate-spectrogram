//
//  writeOutImage.hpp
//  writeOutImage
//
//  Created by LegitMichel777 on 2021/8/8.
//

#ifndef writeOutImage_hpp
#define writeOutImage_hpp

#include <stdio.h>
#include <png.h>

using namespace std;
typedef long long ll;

bool writeOutImage(png_bytepp pngdt,string file,ll height,ll width) {
    //this is mostly copied and formatted from someone else's code lmao
    png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    if (png_ptr==NULL) { return 0; }
    png_infop info_ptr=png_create_info_struct(png_ptr);
    if (info_ptr==NULL) {png_destroy_write_struct(&png_ptr,(png_infopp)NULL); return 0;};
    
    if (setjmp(png_jmpbuf(png_ptr))) {png_destroy_write_struct(&png_ptr,&info_ptr); return 0;}
    
    FILE *ofs=fopen(file.c_str(),"wb");
    if (ofs==NULL) {
        return 0;
    }
    png_init_io(png_ptr,ofs);
    
    png_set_IHDR(png_ptr, info_ptr,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
    
    png_write_info(png_ptr,info_ptr);
    
    png_write_image(png_ptr,pngdt);
    
    png_write_end(png_ptr,NULL);
    
    png_destroy_write_struct(&png_ptr,&info_ptr);
    
    fclose(ofs);
    return true;
}

#endif /* writeOutImage_hpp */
