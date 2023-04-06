#pragma once
#include <fstream>
#include <cassert>
#include <vector>

#define USE_MULTITHREAD

#ifdef USE_MULTITHREAD
    #include <thread>
    #include <mutex>
    #include <condition_variable>
    #include <chrono>
#endif

const struct BmpHeader {
    char bitmapSignatureBytes[2] = { 'B', 'M' };
    uint32_t sizeOfBitmapFile = 54 + (3*1920*1080);
    uint32_t reservedBytes = 0;
    uint32_t pixelDataOffset = 54;

    bool write_to_file(std::ofstream& file) const
    {
        assert(file.is_open());

        if (file.is_open())
        {
            file.write(this->bitmapSignatureBytes, 2);
            file.write((const char*)&this->sizeOfBitmapFile, sizeof(uint32_t));
            file.write((const char*)&this->reservedBytes, sizeof(uint32_t));
            file.write((const char*)&this->pixelDataOffset, sizeof(uint32_t));
            return true;
        }

        return false;
    }
} bmpHeader;

const struct BmpInfoHeader {
    uint32_t sizeOfThisHeader = 40;
    int32_t width = 1920; // in pixels
    int32_t height = 1080; // in pixels, negative so that we start drawing at the top left
    uint16_t numberOfColorPlanes = 1; // must be 1
    uint16_t colorDepth = 24;
    uint32_t compressionMethod = 0;
    uint32_t rawBitmapDataSize = 3*1920*1080; // generally ignored
    int32_t horizontalResolution = 3780; // in pixel per meter
    int32_t verticalResolution = 3780; // in pixel per meter
    uint32_t colorTableEntries = 0;
    uint32_t importantColors = 0;

    bool write_to_file(std::ofstream& file) const
    {
        assert(file.is_open());

        if (file.is_open())
        {
            file.write((const char*)&this->sizeOfThisHeader, sizeof(uint32_t));
            file.write((const char*)&this->width, sizeof(int32_t));
            file.write((const char*)&this->height, sizeof(int32_t));
            file.write((const char*)&this->numberOfColorPlanes, sizeof(uint16_t));
            file.write((const char*)&this->colorDepth, sizeof(uint16_t));
            file.write((const char*)&this->compressionMethod, sizeof(uint32_t));
            file.write((const char*)&this->rawBitmapDataSize, sizeof(uint32_t));
            file.write((const char*)&this->horizontalResolution, sizeof(int32_t));
            file.write((const char*)&this->verticalResolution, sizeof(int32_t));
            file.write((const char*)&this->colorTableEntries, sizeof(uint32_t));
            file.write((const char*)&this->importantColors, sizeof(uint32_t));
            return true;
        }

        return false;
    }
} bmpInfoHeader;

struct Pixel {
    uint8_t blue = 0x0;
    uint8_t green = 0x0;
    uint8_t red = 0x0;

    //TVECTOR color;

    Pixel() = default;
    Pixel(uint8_t _red, uint8_t _green, uint8_t _blue) 
        : red(_red), green(_green), blue(_blue) {}

    bool write_to_file(std::ofstream& file) const
    {
        assert(file.is_open());

        if (file.is_open())
        {
            file.write((const char*)&this->blue, sizeof(uint8_t));
            file.write((const char*)&this->green, sizeof(uint8_t));
            file.write((const char*)&this->red, sizeof(uint8_t));
            return true;
        }
        return false;
    }

    bool read_file(std::ifstream& file) const
    {
        assert(file.is_open());

        if (file.is_open())
        {
            file.read((char*)&this->blue, sizeof(uint8_t));
            file.read((char*)&this->green, sizeof(uint8_t));
            file.read((char*)&this->red, sizeof(uint8_t));

            return true;
        }
        return false;
    }
};
std::vector<Pixel> pixels;

const int32_t canvas_width = bmpInfoHeader.width;
const int32_t canvas_height = bmpInfoHeader.height;
const int32_t viewport_width = 1;
const int32_t viewport_height = 1;
const int32_t projection_plane_z = 1;

#define BACKGROUND_COLOR {0.75f, 0.75f, 0.75f, 1.0f};

#ifdef USE_MULTITHREAD
struct ThreadData
{
    std::mutex* mutex;
    std::condition_variable* cnd;
    int32_t min[2]{ 0, };
    int32_t max[2]{ 0, };
};
#endif