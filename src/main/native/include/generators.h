#pragma once

void grid2(int * restrict src, int * restrict dst, double * restrict weight,
           const int height, const int width, const double verticalWeight);

void triangleGrid2(int * restrict src, int * restrict dst, double * restrict weight,
                   const int height, const int width);
