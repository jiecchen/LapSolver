#pragma once

void grid2(int * restrict src, int * restrict dst, double * restrict weight,
           const int height, const int width, const double verticalWeight);

void grid3(int * restrict src, int * restrict dst, double * restrict weight,
           const int x, const int y, const int z,
           const double yweight, const double zweight);

void triangleGrid2(int * restrict src, int * restrict dst, double * restrict weight,
                   const int height, const int width);
