#include "InputManager.h"
#include "game.h"
#include "vector"
#include "../res/includes/glm/glm.hpp"
#include "cmath"
#include "fstream"
#include <iostream>
#include "stb_image.h"
#include "algorithm"
#define PI 3.14159265359

unsigned char* applyGaussianFilter(unsigned char* input, unsigned char* output, int width, int height) {
    int Gaussian_kernel[3][3] = {{1, 2, 1},
                        {2, 4, 2},
                        {1, 2, 1}};

    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 1; j < width - 1; j++)
        {
            float sum = 0;
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                    sum = sum + Gaussian_kernel[x][y] * (input)[((i+ x - 1) * width)*4 + (j+ y - 1)*4 ];

            sum = sum/16;
            unsigned char new_pixel = (unsigned char) std::min(std::max(sum, (float)0), (float)255);
            //set_val(output,to_index(i,j),new_pixel);
            output[i * width * 4 + j * 4] = new_pixel;
            output[i * width * 4 + j * 4 + 1] = new_pixel;
            output[i * width * 4 + j * 4 + 2] = new_pixel;
            output[i * width * 4 + j * 4 + 3] = 255;
        }
    }
return output;
}

unsigned char* applySoble(unsigned char* input, int width, int height) {

    auto* Xoutput = new unsigned char[256*256*4];
    auto* Youtput = new unsigned char[256*256*4];
    auto* Grad_output = new unsigned char[256*256*4];

    int Xkernel[3][3] = {{-1, 0, 1},
                         {-2, 0, 2},
                         {-1, 0, 1}};

    int Ykernel[3][3] = {{-1, -2, -1},
                         {0, 0, 0},
                         {1, 2, 1}};

    for (int i = 1; i < height - 1; i = i+1)
    {
        for (int j = 1; j < width - 1; j = j+1)
        {
            float Xsum = 0;
            float Ysum = 0;
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                {
                    Xsum = Xsum + Xkernel[x][y] * (input)[((i+ x - 1) * width)*4 + (j+ y - 1)*4];
                    Ysum = Ysum + Ykernel[x][y] * (input)[((i+ x - 1) * width)*4 + (j+ y - 1)*4];
                }
            unsigned char Xnew_pixel = (unsigned char) std::min(std::max(Xsum, (float)0), (float)255);
            unsigned char Ynew_pixel = (unsigned char) std::min(std::max(Ysum, (float)0), (float)255);
            //set_val(output,to_index(i,j),new_pixel);
            Xoutput[i * width * 4 + j * 4] = Xnew_pixel;
            Xoutput[i * width * 4 + j * 4 + 1] = Xnew_pixel;
            Xoutput[i * width * 4 + j * 4 + 2] = Xnew_pixel;
            Xoutput[i * width * 4 + j * 4 + 3] = Xnew_pixel;

            Youtput[i * width * 4 + j * 4] = Ynew_pixel;
            Youtput[i * width * 4 + j * 4 + 1] = Ynew_pixel;
            Youtput[i * width * 4 + j * 4 + 2] = Ynew_pixel;
            Youtput[i * width * 4 + j * 4 + 3] = Ynew_pixel;
        }
    }

    for (int i = 1; i < height; i++)
    {
        for (int j = 1; j < width; j++)
        {
            unsigned char Xpixel = (Xoutput)[(i * width)*4 + j*4];
            unsigned char Ypixel = (Youtput)[(i * width)*4 + j*4];
            unsigned char new_pixel= (unsigned char)((int)sqrt((Xpixel * Xpixel) + (Ypixel * Ypixel)));
            Grad_output[i * width * 4 + j * 4] = new_pixel;
            Grad_output[i * width * 4 + j * 4 + 1] = new_pixel;
            Grad_output[i * width * 4 + j * 4 + 2] = new_pixel;
            Grad_output[i * width * 4 + j * 4 + 3] = new_pixel;
        }
    }
    return Grad_output;
}
double customRound(double value, int factor) {
    if (factor <= 0)
        return 0;

    double roundedValue = std::round(value / factor) * factor;
    return std::round(value / factor) * factor;
}
void writePixel(unsigned char* input, int i, int j, unsigned char value){
    input[i*256*4 + 4*j] = value;
    input[i*256*4 + 4*j + 1] = value;
    input[i*256*4 + 4*j + 2] = value;
    input[i*256*4 + 4*j + 3] = 255;

}
unsigned char* nonMaxSuppression(unsigned char* input, unsigned char* Gaussinput, int width, int height) {
    //input is after sobel and Gaussinput i after Gaussian filter
    // I grab the gaussian char to use it in the gradiant calculation
    //we can do the sobel part and the non max supprition in the same func but I want it separately
    unsigned char* output = new unsigned char[width * height * 4];
    int Xkernel[3][3] = {{-1, 0, 1},
                         {-2, 0, 2},
                         {-1, 0, 1}};

    int Ykernel[3][3] = {{-1, -2, -1},
                         {0, 0, 0},
                         {1, 2, 1}};

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            float Xsum = 0;
            float Ysum = 0;

            for (int x = 0; x < 3; x++) {
                for (int y = 0; y < 3; y++) {
                    Xsum = Xsum + Xkernel[x][y] * (Gaussinput[((i + x - 1) * width) * 4 + (j + y - 1) * 4]);
                    Ysum = Ysum + Ykernel[x][y] * (Gaussinput[((i + x - 1) * width) * 4 + (j + y - 1) * 4]);
                }
            }

            float dx = Xsum;
            float dy = Ysum;
            float angle = atan2(dy, dx) * 180.0 / PI;
            int quantized = ((int(customRound(angle, 45)) + 180) % 180) / 45;
            float gradient = sqrt(dx * dx + dy * dy);

            if (quantized == 0) {
                if (gradient > input[(i * width + j - 1) * 4] &&
                    gradient > input[(i * width + j + 1) * 4]) {
                    writePixel(output, i, j, gradient);
                }
            } else if (quantized == 1) {
                if (gradient > input[((i - 1) * width + j + 1) * 4] &&
                    gradient > input[((i + 1) * width + j - 1) * 4]) {
                    writePixel(output, i, j, gradient);
                }
            } else if (quantized == 2) {
                if (gradient > input[((i - 1) * width + j) * 4] &&
                    gradient > input[((i + 1) * width + j) * 4]) {
                    writePixel(output, i, j, gradient);
                }
            } else if (quantized == 3) {
                if (gradient > input[((i - 1) * width + j - 1) * 4] &&
                    gradient > input[((i + 1) * width + j + 1) * 4]) {
                    writePixel(output, i, j, gradient);
                }
            }
        }
    }
    return output;
}

unsigned char* Thresholding(unsigned char* input, int width, int height) {

    unsigned char T1 = (unsigned char)( (255 * 0.1));
    unsigned char T2 = (unsigned char)( (255 * 0.5));
    auto* output = new unsigned char[256*256*4];

    for (int i = 0; i < height - 1; i++) {
        for (int j = 0; j < width - 1; j++) {
            if ((input)[(i * width) * 4 + j * 4] >= T2) {
                writePixel(output,i,j, 255);
            } else if ((input)[(i * width) * 4 + j * 4] < T1) {
                writePixel(output,i,j, 0);
            }
        }
    }
    return output;
}

unsigned char* halftone(unsigned char* input, int width, int height) {
    auto* output = new unsigned char[width * height * 4 * 4 ];
    for(int i = 0; i < height; i = i + 2){
        for(int j = 0; j < width; j = j + 2){
            int average = (input[i*width*4 + 4*j] + input[i*width*4 + 4*j + 1] + input[i*width*4 + 4*j + 2])/3;
            if(average < 0.2 * 255){
                writePixel(output, i, j, 0);
                writePixel(output, i, j+1, 0);
                writePixel(output, i+1, j, 0);
                writePixel(output, i+1, j+1, 0);
            }
            else if(average < 0.4 * 255){
                writePixel(output, i, j, 255);
                writePixel(output, i, j+1, 0);
                writePixel(output, i+1, j, 0);
                writePixel(output, i+1, j+1, 0);
            }
            else if(average < 0.6 * 255){
                writePixel(output, i, j, 255);
                writePixel(output, i, j+1, 0);
                writePixel(output, i+1, j, 0);
                writePixel(output, i+1, j+1, 255);
            }
            else if(average < 0.8 * 255){
                writePixel(output, i, j, 255);
                writePixel(output, i, j+1, 255);
                writePixel(output, i+1, j, 0);
                writePixel(output, i+1, j+1, 255);
            }
            else{
                writePixel(output, i, j, 255);
                writePixel(output, i, j+1, 255);
                writePixel(output, i+1, j, 255);
                writePixel(output, i+1, j+1, 255);
            }
        }
    }
    return output;
}


unsigned char* Floyd_Steinberg(unsigned char* input, int width, int height) {
    double a = 7.0 / 16, b = 3.0 / 16, g = 5.0 / 16, d = 1.0 / 16;

    unsigned char* output = new unsigned char[width * height * 4];
    unsigned char* CopyInput = new unsigned char [width * height * 4];
    for(int i = 0; i < width * height * 4; i++){
        CopyInput[i] = input[i];
    }

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {

            int index = (j * width + i) * 4;


            int prev_Pixel = CopyInput[index];
            int new_Pixel = customRound(prev_Pixel, 16);
            int err = prev_Pixel - new_Pixel;


            output[index] = new_Pixel;
            output[index + 1] = new_Pixel;
            output[index + 2] = new_Pixel;
            output[index + 3] = 255;


            if (i < width - 1) {

                int rightIndex = index + 4;
                CopyInput[rightIndex] += (unsigned char)(err * a);
            }
            if (i > 0 && j < height - 1) {

                int bottomLeftIndex = index + (width - 1) * 4;
                CopyInput[bottomLeftIndex] += (unsigned char)(err * b);
            }
            if (j < height - 1) {

                int bottomIndex = index + width * 4;
                CopyInput[bottomIndex] += (unsigned char)(err * g);
            }
            if (i < width - 1 && j < height - 1) {

                int bottomRightIndex = index + (width + 1) * 4;
                CopyInput[bottomRightIndex] += (unsigned char)(err * d);
            }
        }
    }
    return output;
}

void savePixelValuesIntoFile(const char* f, const unsigned char* pixels, int width, int height) {
    std::ofstream file(f);
    if (!file.is_open()) {
        std::cerr << "cannt open the file " << f << std::endl;
        return;
    }

    // Saving pixel values separated by commas
    int s=width * height;
    for (int i = 0; i < s; ++i) {
        file << static_cast<int>(pixels[i]);
        if (i < width * height - 1) {
            file << ",";
        }
    }
    // closing the file
    file.close();
}


int main(int argc,char *argv[])
{
    const int DISPLAY_WIDTH = 512;
    const int DISPLAY_HEIGHT = 512;
    const float CAMERA_ANGLE = 0.0f;
    const float NEAR = 1.0f;
    const float FAR = 100.0f;


    Game *scn = new Game(CAMERA_ANGLE,(float)DISPLAY_WIDTH/DISPLAY_HEIGHT,NEAR,FAR);

    Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");

    Init(display);

    scn->Init();

    display.SetScene(scn);

    int width, height, numComponents;
    unsigned char* input = stbi_load("../res/textures/lena256.jpg",&width, &height, &numComponents, 4);

    while(!display.CloseWindow())
    {
        /// GRAY_SCALE PART
        scn->AddTexture(width,height,input);
        scn->SetShapeTex(0,0);
        scn->Draw(1, 0, scn->BACK, true, false , 0);

        /// CANNY_EDGE_DETECTION ALGO
        auto* output = new unsigned char[256*256*4];
        auto* Gauss = new unsigned char[256*256*4];
        Gauss= applyGaussianFilter(input, Gauss, width, height);
        output = applySoble(Gauss,  width, height);
        output = nonMaxSuppression(output,Gauss, width, height);
        output = Thresholding(output, width, height);

        scn->AddTexture(256, 256, output);
        scn->SetShapeTex(0,1);
        scn->Draw(1, 0, scn->BACK, false, false , 1);
        savePixelValuesIntoFile("img4.txt", output, 256, 256);

        /// HALFTONE ALGO
        auto* output2 = halftone(input, width, height);
        scn->AddTexture(256, 256, output2);
        scn->SetShapeTex(0,2);
        scn->Draw(1, 0, scn->BACK, false, false , 2);

        //Halftone Saving pixel values for
        savePixelValuesIntoFile("img5.txt", output2, 256, 256);

        /// FLOYD_STEINBERG ALGO
        auto* output3 = Floyd_Steinberg(input, width, height);
        scn->AddTexture(256, 256, output3);
        scn->SetShapeTex(0,3);
        scn->Draw(1, 0, scn->BACK, false, false , 3);

        //Floyd_Steinberg Saving pixel values for
        savePixelValuesIntoFile("img6.txt", output3, 256, 256);


        scn->Motion();
        display.SwapBuffers();

//		scn->Draw(1,0,scn->BACK,true,false);
//		scn->Motion();
//		display.SwapBuffers();
        display.PollEvents();
        delete Gauss;
        delete output;
        delete output2;
        delete output3;

    }
    delete scn;
    return 0;
}
