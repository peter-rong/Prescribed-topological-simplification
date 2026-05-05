#ifndef MRC_READER_H
#define MRC_READER_H

// Suppress security warnings on Windows
#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdio>

class MRCReader {
public:
    /* Initializer */
    MRCReader(const char* fname) {
        sprintf(mrcfile, "%s", fname);
        
        FILE* fin = fopen(fname, "rb");
        if (fin == NULL) {
            printf("Cannot open file %s.\n", fname);
            throw std::runtime_error("Cannot open MRC file: " + std::string(fname));
        }

        // Parse header
        fread(&totx, sizeof(int), 1, fin);
        fread(&toty, sizeof(int), 1, fin);
        fread(&totz, sizeof(int), 1, fin);

        fread(&mode, sizeof(int), 1, fin);

        fread(&offx, sizeof(int), 1, fin);
        fread(&offy, sizeof(int), 1, fin);
        fread(&offz, sizeof(int), 1, fin);

        fread(&dimx, sizeof(int), 1, fin);
        fread(&dimy, sizeof(int), 1, fin);
        fread(&dimz, sizeof(int), 1, fin);
        dimx++;
        dimy++;
        dimz++;

        fread(&angsx, sizeof(float), 1, fin);
        fread(&angsy, sizeof(float), 1, fin);
        fread(&angsz, sizeof(float), 1, fin);

        fread(&anglex, sizeof(float), 1, fin);
        fread(&angley, sizeof(float), 1, fin);
        fread(&anglez, sizeof(float), 1, fin);

        fseek(fin, 12, SEEK_CUR);

        fread(&dmin, sizeof(float), 1, fin);
        fread(&dmax, sizeof(float), 1, fin);
        fread(&dmean, sizeof(float), 1, fin);

        fseek(fin, 4 * 32, SEEK_CUR);

        fread(&drms, sizeof(float), 1, fin);
        fclose(fin);

        dimx = totx;
        dimy = toty;
        dimz = totz;

        printf("Dimension: %d %d %d\n", dimx, dimy, dimz);
        printf("Mode: %d\n", mode);
        printf("Density: from %f to %f, mean at %f, rms at %f\n", dmin, dmax, dmean, drms);
        printf("Cell size: %f %f %f\n", angsx / (dimx - 1), angsy / (dimy - 1), angsz / (dimz - 1));
        printf("Cell angles: %f %f %f\n", anglex, angley, anglez);

        if (mode > 2) {
            printf("Complex mode not supported.\n");
            throw std::runtime_error("Complex mode not supported");
        }
    }

    /* Read volume */
    std::vector<std::vector<std::vector<double>>> getVolume(double adjustment = 0.0) {
        FILE* fin = fopen(mrcfile, "rb");
        fseek(fin, 1024, SEEK_SET);

        char chard;
        short shortd;
        float floatd;
        double d;

        // Create 3D array
        std::vector<std::vector<std::vector<double>>> data(dimz);
        for (int i = 0; i < dimz; ++i) {
            data[i].resize(dimy);
            for (int j = 0; j < dimy; ++j) {
                data[i][j].resize(dimx, 0.0);
            }
        }

        for (int i = 0; i < dimz; i++)
            for (int j = 0; j < dimy; j++)
                for (int k = 0; k < dimx; k++) {
                    switch (mode) {
                        case 0:
                            fread(&chard, sizeof(char), 1, fin);
                            d = (double)chard;
                            break;
                        case 1:
                            fread(&shortd, sizeof(short), 1, fin);
                            d = (double)shortd;
                            break;
                        case 2:
                            fread(&floatd, sizeof(float), 1, fin);
                            d = (double)floatd;
                            break;
                    }
                    

                    data[i][j][k] = -(d - adjustment);

                    if (data[i][j][k] == -0.0) {
                        data[i][j][k] = 0.0;
                    }
                }
        fclose(fin);

        printf("Done reading.\n");
        return data;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = angsx / (dimx - 1);
        ay = angsy / (dimy - 1);
        az = angsz / (dimz - 1);
    }

    // Static method for compatibility with existing code
    static std::vector<std::vector<std::vector<double>>> readMRCFile(const std::string& filename, double adjustment) {
        MRCReader reader(filename.c_str());
        return reader.getVolume(adjustment);
    }

private:
    int totx, toty, totz;
    int offx, offy, offz;
    int dimx, dimy, dimz;

    float angsx, angsy, angsz;
    float anglex, angley, anglez;
    float dmin, dmax, dmean, drms;

    int mode;

    char mrcfile[1024];
};

#endif // MRC_READER_H