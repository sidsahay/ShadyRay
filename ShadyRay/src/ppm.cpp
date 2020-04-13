//
// Created by walksbynight on 26/3/18.
//

#include "ppm.hpp"

namespace ppm {

    PPMImage* ReadPPM(const char *source) {
        FILE *fp;
        if ((fp = fopen(source, "rb")) == nullptr)
        {
            printf("Fopen failed.\n");
            fflush(stdin);
            return nullptr;
        }

        char buffer[1024] = { 0 };
        char ch;
        int i = 0;
        while (!isspace((ch = fgetc(fp))))
            buffer[i++] = ch;

        if (buffer[0] != 'P' || buffer[1] != '6')
        {
            printf("Invalid file.\n");
            fclose(fp);
            return nullptr;
        }

        while (isspace((ch = fgetc(fp))));
        memset(buffer, 0, sizeof(buffer));
        i = 0;
        buffer[i++] = ch;
        while (!isspace(ch = fgetc(fp)))
            buffer[i++] = ch;

        int width = atoi(buffer);


        while (isspace((ch = fgetc(fp))));
        memset(buffer, 0, sizeof(buffer));
        i = 0;
        buffer[i++] = ch;
        while (!isspace(ch = fgetc(fp)))
            buffer[i++] = ch;

        int height = atoi(buffer);


        while (isspace((ch = fgetc(fp))));
        memset(buffer, 0, sizeof(buffer));
        i = 0;
        buffer[i++] = ch;
        while (!isspace(ch = fgetc(fp)))
            buffer[i++] = ch;

        int maxVal = atoi(buffer);

        printf("I: H,W,V = %d, %d, %d\n", height, width, maxVal);
        auto *data = new Pixel[width * height];

        int pos1 = ftell(fp);
        fseek(fp, 0, SEEK_END);
        int pos2 = ftell(fp);
        printf("Length: %d\n", pos2 - pos1);
        fseek(fp, pos1, SEEK_SET);
        fread(data, sizeof(Pixel), height * width, fp);
        fclose(fp);

        auto *image = (PPMImage*)calloc(1, sizeof(PPMImage));
        image->height = height;
        image->width = width;
        image->maxVal = maxVal;
        image->data = data;
        return image;
    }

    bool WritePPM(PPMImage* image, const char *dest) {
        FILE *fp;
        if ((fp = fopen(dest, "wb")) == nullptr)
        {
            printf("Fopen failed.\n");
            fflush(stdin);
            return false;
        }

        char buffer[1024] = { 0 };
        printf("O: H,W,V = %d, %d, %d\n", image->height, image->width, image->maxVal);
        fprintf(fp, "P6\n%d %d\n%d\n", image->width, image->height, image->maxVal);
        fwrite(image->data, sizeof(Pixel), image->height * image->width, fp);
        fflush(fp);
        fclose(fp);
        return true;
    }

}