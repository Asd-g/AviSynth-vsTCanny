## Description

Builds an edge map using canny edge detection.

This is [a port of the VapourSynth plugin TCanny](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-TCanny).

### Requirements:

- AviSynth 2.60 / AviSynth+ 3.4 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases)) (Windows only)

### Usage:

```
vsTCanny (clip, float "sigmaY", float "sigmaU", float "sigmaV", float sigma_vY", float "sigma_vU", float "sigma_vV", float "t_h", float "t_l", int "mode", int "op", float "gmmax", int "y", int "u", int "v", int "opt")
```

### Parameters:

- clip\
    A clip to process. All planar formats are supported.
    
- sigmaY, sigmaU, sigmaV\
    Standard deviation of horizontal gaussian blur.\
    Must be greater than 0.0.\
    Default: sigmaY = 1.5; sigmaU = sigmaV = 0.75 for half resolution chroma and sigmaU = sigmaV = 1.5 for full resolution chroma.

- sigma_vY, sigma_vU, sigma_vV\
    Standard deviation of vertical gaussian blur.\
    Must be greater than 0.0.\
    Default: sigma_vY = 1.5; sigma_vU = sigma_vV = 0.75 for half resolution chroma and sigma_vU = sigma_vV = 1.5 for full resolution chroma.
    
- t_h\
    High gradient magnitude threshold for hysteresis.\
    Default: 8.0.
    
- t_l\
    Low gradient magnitude threshold for hysteresis.\
    Must be lower than t_h.\
    Default: 1.0.
    
- mode\
    Sets output format.\
    -1: Gaussian blur only.\
    0: Thresholded edge map (2^bitdepth-1 for edge, 0 for non-edge).\
    1: Gradient magnitude map.\
    Default: 0.
    
- op\
    Sets the operator for edge detection.\
    0: The operator used in tritical's original filter.\
    1:  the Prewitt operator whose use is proposed by P. Zhou et al. [1]
    2: The Sobel operator.\
    3: The Scharr operator.\
    Default: 1.

- gmmax\
    Used for scaling gradient magnitude into [0, 2^bitdepth-1] for mode=1.\
    Musbe greater than or equal to 1.0.\
    Default: 50.0.
    
- y, u, v\
    Planes to process.\
    1: Return garbage.\
    2: Copy plane.\
    3: Process plane. Always process planes when the clip is RGB.\
    Default: y = u = v = 3.
    
- opt\
    Sets which cpu optimizations to use.\
    -1: Auto-detect.\
    0: Use C++ code.\
    1: Use SSE2 code.\
    2: Use AVX2 code.\
    Default: -1.

[1]: Zhou, P., Ye, W., & Wang, Q. (2011). An Improved Canny Algorithm for Edge Detection. Journal of Computational Information Systems, 7(5), 1516-1523.

### Building:

- Windows\
    Use solution files.

- Linux
    ```
    Requirements:
        - Git
        - C++17 compiler
        - CMake >= 3.16
    ```
    ```
    git clone https://github.com/Asd-g/AviSynth-vsTCanny && \
    cd AviSynth-vsTCanny && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    sudo make install
    ```
