#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "avisynth.h"

static constexpr float M_PIF{ 3.14159265358979323846f };
static constexpr float M_1_PIF{ 0.318309886183790671538f };
static constexpr float fltMax{ std::numeric_limits<float>::max() };
static constexpr float fltLowest{ std::numeric_limits<float>::lowest() };

enum Operator
{
    TRITICAL,
    PREWITT,
    SOBEL,
    SCHARR,
    KROON,
    KIRSCH,
    FDOG
};

class vsTCanny : public GenericVideoFilter
{
    float t_h_;
    float t_l_;
    int mode_;
    float op_;
    float scale;
    bool process[3];
    int radiusH[3];
    int radiusV[3];
    std::unique_ptr<float[]> weightsH[3];
    std::unique_ptr<float[]> weightsV[3];
    int peak;
    int alignment;
    int radiusAlign;
    float* blur;
    float* gradient;
    int* direction;
    std::unique_ptr<bool[]> found;
    bool has_at_least_v8;

    template<typename T>
    void filter_c(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
    template<typename T>
    void filter_sse2(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
    template<typename T>
    void filter_avx2(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
    template<typename T>
    void filter_avx512(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;

    void (vsTCanny::* filter)(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;

public:
    vsTCanny(PClip _child, float sigmaY, float sigmaU, float sigmaV, float sigma_vY, float sigma_vU, float sigma_vV, float t_h, float t_l, int mode, int op, float gmmax, int y, int u, int v, int opt, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
    }
    ~vsTCanny();
};

static void hysteresis(float* __restrict srcp, bool* __restrict found, const int width, const int height, const int stride, const float t_h, const float t_l) noexcept
{
    std::fill_n(found, width * height, false);
    std::vector<std::pair<int, int>> coordinates;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            if (!found[width * y + x] && srcp[stride * y + x] >= t_h)
            {
                srcp[stride * y + x] = fltMax;
                found[width * y + x] = true;

                coordinates.emplace_back(std::make_pair(x, y));

                while (!coordinates.empty())
                {
                    const auto pos = coordinates.back();
                    coordinates.pop_back();

                    const int xxStart{ std::max(pos.first - 1, 0) };
                    const int xxStop{ std::min(pos.first + 1, width - 1) };
                    const int yyStart{ std::max(pos.second - 1, 0) };
                    const int yyStop{ std::min(pos.second + 1, height - 1) };

                    for (int yy{ yyStart }; yy <= yyStop; ++yy)
                    {
                        for (int xx{ xxStart }; xx <= xxStop; ++xx)
                        {
                            if (!found[width * yy + xx] && srcp[stride * yy + xx] >= t_l)
                            {
                                srcp[stride * yy + xx] = fltMax;
                                found[width * yy + xx] = true;

                                coordinates.emplace_back(std::make_pair(xx, yy));
                            }
                        }
                    }
                }
            }
        }
    }
}
