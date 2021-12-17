#include <cmath>
#include <string>

#include "vsTCanny.h"

AVS_FORCEINLINE void* aligned_malloc(size_t size, size_t align)
{
    void* result = [&]() {
#ifdef _MSC_VER 
        return _aligned_malloc(size, align);
#else 
        if (posix_memalign(&result, align, size))
            return result = nullptr;
        else
            return result;
#endif
    }();

    return result;
}

AVS_FORCEINLINE void aligned_free(void* ptr)
{
#ifdef _MSC_VER 
    _aligned_free(ptr);
#else 
    free(ptr);
#endif
}

template<typename T>
static void copyPlane(const T* srcp, float* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
            dstp[x] = srcp[x];

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlur(const T* _srcp, float* __restrict temp, float* __restrict dstp, const float* weightsH, const float* weightsV, const int width, const int height, const int srcStride, const int dstStride, const int radiusH, const int radiusV) noexcept
{
    const int diameter{ radiusV * 2 + 1 };
    std::unique_ptr<const T* []> srcp{ std::make_unique<const T* []>(diameter) };

    srcp[radiusV] = _srcp;
    for (int i{ 1 }; i <= radiusV; ++i)
        srcp[radiusV - i] = srcp[radiusV + i] = srcp[radiusV] + srcStride * i;

    weightsH += radiusH;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            float sum{ 0.0f };

            for (int i{ 0 }; i < diameter; ++i)
                sum += srcp[i][x] * weightsV[i];

            temp[x] = sum;
        }

        for (int i{ 1 }; i <= radiusH; ++i)
        {
            temp[-i] = temp[i];
            temp[width - 1 + i] = temp[width - 1 - i];
        }

        for (int x{ 0 }; x < width; ++x)
        {
            float sum{ 0.0f };

            for (int i = -radiusH; i <= radiusH; ++i)
                sum += temp[x + i] * weightsH[i];

            dstp[x] = sum;
        }

        for (int i{ 0 }; i < diameter - 1; ++i)
            srcp[i] = srcp[i + 1];

        srcp[diameter - 1] += (y < height - 1 - radiusV) ? srcStride : -srcStride;;

        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlurV(const T* _srcp, float* __restrict dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius) noexcept
{
    const int diameter{ radius * 2 + 1 };
    std::unique_ptr<const T* []> srcp{ std::make_unique<const T* []>(diameter) };

    srcp[radius] = _srcp;
    for (int i{ 1 }; i <= radius; ++i)
        srcp[radius - i] = srcp[radius + i] = srcp[radius] + srcStride * i;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            float sum{ 0.0f };

            for (int i{ 0 }; i < diameter; ++i)
                sum += srcp[i][x] * weights[i];

            dstp[x] = sum;
        }

        for (int i{ 0 }; i < diameter - 1; ++i)
            srcp[i] = srcp[i + 1];

        srcp[diameter - 1] += (y < height - 1 - radius) ? srcStride : -srcStride;

        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlurH(const T* srcp, float* __restrict temp, float* __restrict dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius) noexcept
{
    weights += radius;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
            temp[x] = srcp[x];

        for (int i{ 1 }; i <= radius; ++i)
        {
            temp[-i] = temp[i];
            temp[width - 1 + i] = temp[width - 1 - i];
        }

        for (int x{ 0 }; x < width; ++x)
        {
            float sum{ 0.0f };

            for (int i{ -radius }; i <= radius; ++i)
                sum += temp[x + i] * weights[i];

            dstp[x] = sum;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

static void detectEdge(float* __restrict blur, float* __restrict gradient, int* __restrict direction, const int width, const int height, const int stride, const int bgStride, const int mode, const int op, const float scale) noexcept
{
    float* __restrict cur{ blur };
    float* __restrict next{ blur + bgStride };
    float* __restrict next2{ blur + bgStride * 2 };
    float* __restrict prev{ next };
    float* __restrict prev2{ next2 };

    cur[-1] = cur[1];
    cur[width] = cur[width - 2];

    if (op == FDOG)
    {
        cur[-2] = cur[2];
        cur[width + 1] = cur[width - 3];
    }

    for (int y{ 0 }; y < height; ++y)
    {
        next[-1] = next[1];
        next[width] = next[width - 2];

        if (op == FDOG)
        {
            next[-2] = next[2];
            next[width + 1] = next[width - 3];

            next2[-1] = next2[1];
            next2[-2] = next2[2];
            next2[width] = next2[width - 2];
            next2[width + 1] = next2[width - 3];
        }

        for (int x{ 0 }; x < width; ++x)
        {
            float gx, gy;

            if (op != FDOG)
            {
                const float c1{ prev[x - 1] };
                const float c2{ prev[x] };
                const float c3{ prev[x + 1] };
                const float c4{ cur[x - 1] };
                const float c6{ cur[x + 1] };
                const float c7{ next[x - 1] };
                const float c8{ next[x] };
                const float c9{ next[x + 1] };

                switch (op)
                {
                    case TRITICAL:
                    {
                        gx = c6 - c4;
                        gy = c2 - c8;
                        break;
                    }
                    case PREWITT:
                    {
                        gx = (c3 + c6 + c9 - c1 - c4 - c7) / 2.0f;
                        gy = (c1 + c2 + c3 - c7 - c8 - c9) / 2.0f;
                        break;
                    }
                    case SOBEL:
                    {
                        gx = c3 + 2.0f * c6 + c9 - c1 - 2.0f * c4 - c7;
                        gy = c1 + 2.0f * c2 + c3 - c7 - 2.0f * c8 - c9;
                        break;
                    }
                    case SCHARR:
                    {
                        gx = 3.0f * c3 + 10.0f * c6 + 3.0f * c9 - 3.0f * c1 - 10.0f * c4 - 3.0f * c7;
                        gy = 3.0f * c1 + 10.0f * c2 + 3.0f * c3 - 3.0f * c7 - 10.0f * c8 - 3.0f * c9;
                        break;
                    }
                    case KROON:
                    {
                        gx = 17.0f * c3 + 61.0f * c6 + 17.0f * c9 - 17.0f * c1 - 61.0f * c4 - 17.0f * c7;
                        gy = 17.0f * c1 + 61.0f * c2 + 17.0f * c3 - 17.0f * c7 - 61.0f * c8 - 17.0f * c9;
                        break;
                    }
                    case KIRSCH:
                    {
                        const float g1{ 5.0f * c1 + 5.0f * c2 + 5.0f * c3 - 3.0f * c4 - 3.0f * c6 - 3.0f * c7 - 3.0f * c8 - 3.0f * c9 };
                        const float g2{ 5.0f * c1 + 5.0f * c2 - 3.0f * c3 + 5.0f * c4 - 3.0f * c6 - 3.0f * c7 - 3.0f * c8 - 3.0f * c9 };
                        const float g3{ 5.0f * c1 - 3.0f * c2 - 3.0f * c3 + 5.0f * c4 - 3.0f * c6 + 5.0f * c7 - 3.0f * c8 - 3.0f * c9 };
                        const float g4{ -3.0f * c1 - 3.0f * c2 - 3.0f * c3 + 5.0f * c4 - 3.0f * c6 + 5.0f * c7 + 5.0f * c8 - 3.0f * c9 };
                        const float g5{ -3.0f * c1 - 3.0f * c2 - 3.0f * c3 - 3.0f * c4 - 3.0f * c6 + 5.0f * c7 + 5.0f * c8 + 5.0f * c9 };
                        const float g6{ -3.0f * c1 - 3.0f * c2 - 3.0f * c3 - 3.0f * c4 + 5.0f * c6 - 3.0f * c7 + 5.0f * c8 + 5.0f * c9 };
                        const float g7{ -3.0f * c1 - 3.0f * c2 + 5.0f * c3 - 3.0f * c4 + 5.0f * c6 - 3.0f * c7 - 3.0f * c8 + 5.0f * c9 };
                        const float g8{ -3.0f * c1 + 5.0f * c2 + 5.0f * c3 - 3.0f * c4 + 5.0f * c6 - 3.0f * c7 - 3.0f * c8 - 3.0f * c9 };
                        const float g{ std::max({ std::abs(g1), std::abs(g2), std::abs(g3), std::abs(g4), std::abs(g5), std::abs(g6), std::abs(g7), std::abs(g8) }) };
                        gradient[x] = g * scale;
                        break;
                    }
                }
            }
            else
            {
                const float c1{ prev2[x - 2] };
                const float c2{ prev2[x - 1] };
                const float c3{ prev2[x] };
                const float c4{ prev2[x + 1] };
                const float c5{ prev2[x + 2] };
                const float c6{ prev[x - 2] };
                const float c7{ prev[x - 1] };
                const float c8{ prev[x] };
                const float c9{ prev[x + 1] };
                const float c10{ prev[x + 2] };
                const float c11{ cur[x - 2] };
                const float c12{ cur[x - 1] };
                const float c14{ cur[x + 1] };
                const float c15{ cur[x + 2] };
                const float c16{ next[x - 2] };
                const float c17{ next[x - 1] };
                const float c18{ next[x] };
                const float c19{ next[x + 1] };
                const float c20{ next[x + 2] };
                const float c21{ next2[x - 2] };
                const float c22{ next2[x - 1] };
                const float c23{ next2[x] };
                const float c24{ next2[x + 1] };
                const float c25{ next2[x + 2] };

                gx = c5 + 2.0f * c10 + 3.0f * c15 + 2.0f * c20 + c25 + c4 + 2.0f * c9 + 3.0f * c14 + 2.0f * c19 + c24
                    - c2 - 2.0f * c7 - 3.0f * c12 - 2.0f * c17 - c22 - c1 - 2.0f * c6 - 3.0f * c11 - 2.0f * c16 - c21;
                gy = c1 + 2.0f * c2 + 3.0f * c3 + 2.0f * c4 + c5 + c6 + 2.0f * c7 + 3.0f * c8 + 2.0f * c9 + c10
                    - c16 - 2.0f * c17 - 3.0f * c18 - 2.0f * c19 - c20 - c21 - 2.0f * c22 - 3.0f * c23 - 2.0f * c24 - c25;
            }

            if (op != KIRSCH)
            {
                gx *= scale;
                gy *= scale;
                gradient[x] = std::sqrt(gx * gx + gy * gy);
            }

            if (mode == 0)
            {
                float dr{ std::atan2(gy, gx) };

                if (dr < 0.0f)
                    dr += M_PIF;

                const int bin{ static_cast<int>(dr * 4.0f * M_1_PIF + 0.5f) };
                direction[x] = (bin >= 4) ? 0 : bin;
            }
        }

        prev2 = prev;
        prev = cur;
        cur = next;

        if (op != FDOG)
            next += (y < height - 2) ? bgStride : -bgStride;
        else
        {
            next = next2;
            next2 += (y < height - 3) ? bgStride : -bgStride;
        }

        gradient += bgStride;
        direction += stride;
    }
}

static void nonMaximumSuppression(const int* direction, float* __restrict gradient, float* __restrict blur, const int width, const int height, const int stride, const int bgStride, const int radiusAlign) noexcept
{
    const int offsets[]{ 1, -bgStride + 1, -bgStride, -bgStride - 1 };

    gradient[-1] = gradient[1];
    gradient[-1 + bgStride * (height - 1)] = gradient[1 + bgStride * (height - 1)];
    gradient[width] = gradient[width - 2];
    gradient[width + bgStride * (height - 1)] = gradient[width - 2 + bgStride * (height - 1)];
    std::copy_n(gradient - radiusAlign + bgStride, width + radiusAlign * 2, gradient - radiusAlign - bgStride);
    std::copy_n(gradient - radiusAlign + bgStride * (height - 2), width + radiusAlign * 2, gradient - radiusAlign + bgStride * height);

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            const int offset{ offsets[direction[x]] };
            blur[x] = (gradient[x] >= std::max(gradient[x + offset], gradient[x - offset])) ? gradient[x] : fltLowest;
        }

        direction += stride;
        gradient += bgStride;
        blur += bgStride;
    }
}

template<typename T>
static void binarizeCE(const float* srcp, T* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const int peak) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = (srcp[x] == fltMax) ? static_cast<T>(peak) : 0;
            else
                dstp[x] = (srcp[x] == fltMax) ? 1.0f : 0.0f;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void discretizeGM(const float* srcp, T* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const int peak) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = static_cast<T>(std::min(static_cast<int>(srcp[x] + 0.5f), peak));
            else
                dstp[x] = srcp[x];
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void vsTCanny::filter_c(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept
{
    const int planes_y[3]{ PLANAR_Y, PLANAR_U, PLANAR_V };
    const int planes_r[3]{ PLANAR_G, PLANAR_B, PLANAR_R };
    const int* current_planes{ (vi.IsRGB()) ? planes_r : planes_y };
    const int planecount{ std::min(vi.NumComponents(), 3) };

    for (int i{ 0 }; i < planecount; ++i)
    {
        const int height{ src->GetHeight(current_planes[i]) };

        if (process[i])
        {
            const size_t stride{ src->GetPitch(current_planes[i]) / sizeof(T) };
            const size_t bgStride{ stride + radiusAlign * 2 };
            const size_t dst_stride{ dst->GetPitch(current_planes[i]) / sizeof(T) };
            const size_t width{ src->GetRowSize(current_planes[i]) / sizeof(T) };
            const T* srcp{ reinterpret_cast<const T*>(src->GetReadPtr(current_planes[i])) };
            T* dstp{ reinterpret_cast<T*>(dst->GetWritePtr(current_planes[i])) };

            float* blur{ vsTCanny::blur + radiusAlign };
            float* gradient{ vsTCanny::gradient + bgStride + radiusAlign };

            if (radiusV[i] && radiusH[i])
                gaussianBlur(srcp, gradient, blur, weightsH[i].get(), weightsV[i].get(), width, height, stride, bgStride, radiusH[i], radiusV[i]);
            else if (radiusV[i])
                gaussianBlurV(srcp, blur, weightsV[i].get(), width, height, stride, bgStride, radiusV[i]);
            else if (radiusH[i])
                gaussianBlurH(srcp, gradient, blur, weightsH[i].get(), width, height, stride, bgStride, radiusH[i]);
            else
                copyPlane(srcp, blur, width, height, stride, bgStride);

            if (mode_ != -1)
            {
                detectEdge(blur, gradient, direction, width, height, stride, bgStride, mode_, op_, scale);

                if (mode_ == 0)
                {
                    nonMaximumSuppression(direction, gradient, blur, width, height, stride, bgStride, radiusAlign);
                    hysteresis(blur, found.get(), width, height, bgStride, t_h_, t_l_);
                }
            }

            if (mode_ == 0)
                binarizeCE(blur, dstp, width, height, bgStride, dst_stride, peak);
            else
                discretizeGM((mode_ == 1) ? gradient : blur, dstp, width, height, bgStride, dst_stride, peak);
        }
        else
            env->BitBlt(dst->GetWritePtr(current_planes[i]), dst->GetPitch(current_planes[i]), src->GetReadPtr(current_planes[i]), src->GetPitch(current_planes[i]), src->GetRowSize(current_planes[i]), height);
    }
}

static float* gaussianWeights(const float sigma, int& radius) noexcept
{
    const int diameter{ std::max(static_cast<int>(sigma * 3.0f + 0.5f), 1) * 2 + 1 };
    radius = diameter / 2;

    float* weights{ new float[diameter]() };
    float sum{ 0.0f };

    for (int k = -radius; k <= radius; ++k)
    {
        const float w{ std::exp(-(k * k) / (2.0f * sigma * sigma)) };
        weights[k + radius] = w;
        sum += w;
    }

    for (int k{ 0 }; k < diameter; ++k)
        weights[k] /= sum;

    return weights;
}

vsTCanny::vsTCanny(PClip _child, float sigmaY, float sigmaU, float sigmaV, float sigma_vY, float sigma_vU, float sigma_vV, float t_h, float t_l, int mode, int op, float scale_, int y, int u, int v, int opt, IScriptEnvironment* env)
    : GenericVideoFilter(_child), t_h_(t_h), t_l_(t_l), mode_(mode), op_(op), scale(scale_)
{
    if (!vi.IsPlanar())
        env->ThrowError("vsTCanny: the clip is not in planar format.");

    const int height{ vi.height };
    const int width{ vi.width };

    if (height < 3)
        env->ThrowError("vsTCanny: the clip's height must be at least 3.");
    if (t_l_ >= t_h_)
        env->ThrowError("vsTCanny: t_h must be greater than t_l.");
    if (mode_ < -1 || mode_ > 1)
        env->ThrowError("vsTCanny: mode must be -1, 0, or 1.");
    if (op_ < 0 || op_ > 6)
        env->ThrowError("vsTCanny: op must be 0, 1, 2, 3, 4, 5 or 6.");
    if (op_ == 5 && mode == 0)
        env->ThrowError("vsTCanny: op=5 cannot be used when mode=0.");
    if (scale <= 0.0f)
        env->ThrowError("vsTCanny: scale must be greater than 0.0.");
    if (opt < -1 || opt > 3)
        env->ThrowError("vsTCanny: opt must be between -1..3.");

    const bool avx512{ !!(env->GetCPUFlags() & CPUF_AVX512F) };
    const bool avx2{ !!(env->GetCPUFlags() & CPUF_AVX2) };
    const bool sse2{ !!(env->GetCPUFlags() & CPUF_SSE2) };

    if (!avx512 && opt == 3)
        env->ThrowError("vsTCanny: opt=2 requires AVX512.");
    if (!avx2 && opt == 2)
        env->ThrowError("vsTCanny: opt=2 requires AVX2.");
    if (!sse2 && opt == 1)
        env->ThrowError("vsTCanny: opt=1 requires SSE2.");

    const bool rgb{ vi.IsRGB() };
    const int sw{ vi.GetPlaneWidthSubsampling(PLANAR_U) };
    const int sh{ vi.GetPlaneHeightSubsampling(PLANAR_U) };
    const int planecount{ std::min(vi.NumComponents(), 3) };

    if (planecount > 1)
    {
        if (sigmaU == -1354.4f)
            sigmaU = (rgb) ? sigmaY : (sigmaY / (1 << sw));
        if (sigmaV == -1354.4f)
            sigmaV = sigmaU;
        if (sigma_vU == -1354.4f)
            sigma_vU = (rgb) ? sigma_vY : (sigma_vY / (1 << sh));
        if (sigma_vV == -1354.4f)
            sigma_vV = sigma_vU;
    }

    const float sigmaH[3]{ sigmaY, sigmaU, sigmaV };
    const float sigmaV_[3]{ sigma_vY, sigma_vU, sigma_vV };
    const int planes[3]{ y, u, v };

    for (int i{ 0 }; i < planecount; ++i)
    {
        if (rgb)
            process[i] = true;
        else
        {
            switch (planes[i])
            {
                case 3: process[i] = true; break;
                case 2: process[i] = false; break;
            }
        }

        if (sigmaH[i] < 0.0f)
        {
            const std::string sigmaOrder[3]{ "sigmaY", "sigmaU", "sigmaV" };
            env->ThrowError(std::string{ "vsTCanny: " + sigmaOrder[i] + " must be greater than or equal to 0.0." }.c_str());
        }
        if (sigmaV_[i] < 0.0f)
        {
            const std::string sigmaVOrder[3]{ "sigma_vY", "sigma_vU", "sigma_vV" };
            env->ThrowError(std::string{ "vsTCanny: " + sigmaVOrder[i] + " must be greater than or equal to 0.0." }.c_str());
        }
        if (planes[i] < 1 || planes[i] > 3)
            env->ThrowError("vsTCanny: y, u, v must be between 1..3.");

        if (process[i])
        {
            if (sigmaH[i])
            {
                weightsH[i].reset(gaussianWeights(sigmaH[i], radiusH[i]));

                const int width_{ (i && !rgb) ? (width >> sw) : width };
                if (width_ < radiusH[i] + 1)
                {
                    const std::string planeOrder[3]{ "first", "second", "third" };
                    env->ThrowError(std::string{ "vsTCanny: the " + planeOrder[i] + " plane's width must be greater than or equal to " + std::to_string(radiusH[i] + 1) + " for specified sigma." }.c_str());
                }
            }
            else
                radiusH[i] = 0;

            if (sigmaV_[i])
            {
                weightsV[i].reset(gaussianWeights(sigmaV_[i], radiusV[i]));

                const int height_{ (i && !rgb) ? (height >> sh) : height };
                if (height_ < radiusV[i] + 1)
                {
                    const std::string planeOrder[3]{ "first", "second", "third" };
                    env->ThrowError(std::string{ "vsTCanny: the " + planeOrder[i] + " plane's height must be greater than or equal to " + std::to_string(radiusV[i] + 1) + " for specified sigma_v." }.c_str());
                }
            }
            else
                radiusV[i] = 0;
        }
    }

    const int comp_size{ vi.ComponentSize() };
    if (comp_size < 4)
    {
        peak = (1 << vi.BitsPerComponent()) - 1;
        const float scale_{ peak / 255.0f };
        t_h_ *= scale_;
        t_l_ *= scale_;
    }
    else
    {
        t_h_ /= 255.0f;
        t_l_ /= 255.0f;
    }

    int vectorSize;

    if ((avx512 && opt < 0) || opt == 3)
    {
        vectorSize = 16;
        alignment = 64;

        switch (comp_size)
        {
            case 1: filter = &vsTCanny::filter_avx512<uint8_t>; break;
            case 2: filter = &vsTCanny::filter_avx512<uint16_t>; break;
            default: filter = &vsTCanny::filter_avx512<float>; break;
        }
    }
    else if ((avx2 && opt < 0) || opt == 2)
    {
        vectorSize = 8;
        alignment = 32;

        switch (comp_size)
        {
            case 1: filter = &vsTCanny::filter_avx2<uint8_t>; break;
            case 2: filter = &vsTCanny::filter_avx2<uint16_t>; break;
            default: filter = &vsTCanny::filter_avx2<float>; break;
        }
    }
    else if ((sse2 && opt < 0) || opt == 1)
    {
        vectorSize = 4;
        alignment = 16;

        switch (comp_size)
        {
            case 1: filter = &vsTCanny::filter_sse2<uint8_t>; break;
            case 2: filter = &vsTCanny::filter_sse2<uint16_t>; break;
            default: filter = &vsTCanny::filter_sse2<float>; break;
        }
    }
    else
    {
        vectorSize = 1;
        alignment = 4;

        switch (comp_size)
        {
            case 1: filter = &vsTCanny::filter_c<uint8_t>; break;
            case 2: filter = &vsTCanny::filter_c<uint16_t>; break;
            default: filter = &vsTCanny::filter_c<float>; break;
        }
    }

    radiusAlign = (std::max({ radiusH[0], radiusH[1], radiusH[2], (op == FDOG) ? 2 : 1 }) + vectorSize - 1) & ~(vectorSize - 1);

    PVideoFrame clip{ child->GetFrame(0, env) };
    const int pitch{ clip->GetPitch() };

    blur = reinterpret_cast<float*>(aligned_malloc((pitch / comp_size + radiusAlign * 2) * height * sizeof(float), alignment));
    if (!blur)
        env->ThrowError("vsTCanny: malloc failure (blur).");

    gradient = reinterpret_cast<float*>(aligned_malloc((pitch / comp_size + radiusAlign * 2) * (height + 2) * sizeof(float), alignment));
    if (!gradient)
        env->ThrowError("vsTCanny: malloc failure (gradient).");

    if (mode_ == 0)
    {
        direction = reinterpret_cast<int*>(aligned_malloc(pitch / comp_size * height * sizeof(int), alignment));
        if (!direction)
            env->ThrowError("vsTCanny: malloc failure (direction).");

        found = std::make_unique<bool[]>(width * height);
    }
    else
        direction = nullptr;

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; }
}

vsTCanny::~vsTCanny()
{
    aligned_free(blur);
    aligned_free(gradient);

    if (direction)
        aligned_free(direction);
}

PVideoFrame __stdcall vsTCanny::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src{ child->GetFrame(n, env) };
    PVideoFrame dst{ (has_at_least_v8) ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi) };

    (this->*filter)(src, dst, env);

    return dst;
}

AVSValue __cdecl Create_vsTCanny(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    const float sigmaY{ args[1].AsFloatf(1.5f) };
    if (sigmaY < 0.0f)
        env->ThrowError("vsTCanny: sigmaY must be greater than or equal to 0.0.");

    const float sigmaU{ args[2].AsFloatf(-1354.4f) };
    if (sigmaU < 0.0f && sigmaU != -1354.4f)
        env->ThrowError("vsTCanny: sigmaU must be greater than or equal to 0.0.");

    const float sigmaV{ args[3].AsFloatf(sigmaU) };
    if (sigmaV < 0.0f && sigmaV != -1354.4f)
        env->ThrowError("vsTCanny: sigmaV must be greater than or equal to 0.0.");

    const float sigma_vY{ args[4].AsFloatf(sigmaY) };
    if (sigma_vY < 0.0f)
        env->ThrowError("vsTCanny: sigma_vY must be greater than or equal to 0.0.");

    const float sigma_vU{ args[5].AsFloatf(-1354.4f) };
    if (sigma_vU < 0.0f && sigma_vU != -1354.4f)
        env->ThrowError("vsTCanny: sigma_vU must be greater than or equal to 0.0.");

    const float sigma_vV{ args[6].AsFloatf(sigma_vU) };
    if (sigma_vV < 0.0f && sigma_vV != -1354.4f)
        env->ThrowError("vsTCanny: sigma_vV must be greater than or equal to 0.0.");

    return new vsTCanny(
        args[0].AsClip(),
        sigmaY,
        sigmaU,
        sigmaV,
        sigma_vY,
        sigma_vU,
        sigma_vV,
        args[7].AsFloatf(8.0f),
        args[8].AsFloatf(1.0f),
        args[9].AsInt(0),
        args[10].AsInt(1),
        args[11].AsFloatf(1.0f),
        args[12].AsInt(3),
        args[13].AsInt(3),
        args[14].AsInt(3),
        args[15].AsInt(-1),
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("vsTCanny", "c[sigmaY]f[sigmaU]f[sigmaV]f[sigma_vY]f[sigma_vU]f[sigma_vV]f[t_h]f[t_l]f[mode]i[op]i[scale]f[y]i[u]i[v]i[opt]i", Create_vsTCanny, 0);

    return "vsTCanny";
}
