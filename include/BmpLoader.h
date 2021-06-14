#ifndef BMPLOADER_H
#define BMPLOADER_H
#include<windows.h>

class BmpLoader
{
    public:
        unsigned char* textureData;
        int iHeight, iWidth;

        BmpLoader(const char*);
        ~BmpLoader();

    private:
        BITMAPFILEHEADER bfh;
        BITMAPINFOHEADER bih;
};

#endif // BMPLOADER_H
