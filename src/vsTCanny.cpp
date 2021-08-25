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
static void copyPlane(const T* srcp, float* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = srcp[x];
            else
                dstp[x] = srcp[x] + offset;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlur(const T* _srcp, float* __restrict temp, float* __restrict dstp, const float* weightsH, const float* weightsV, const int width, const int height, const int srcStride, const int dstStride, const int radiusH, const int radiusV, const float offset) noexcept
{
    const int diameter = radiusV * 2 + 1;
    const T** srcp = new const T * [diameter];

    srcp[radiusV] = _srcp;
    for (int i = 1; i <= radiusV; ++i)
    {
        srcp[radiusV - i] = srcp[radiusV - 1 + i];
        srcp[radiusV + i] = srcp[radiusV] + srcStride * static_cast<int64_t>(i);
    }

    weightsH += radiusH;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float sum = 0.f;

            for (int i = 0; i < diameter; ++i)
            {
                if constexpr (std::is_integral_v<T>)
                    sum += srcp[i][x] * weightsV[i];
                else
                    sum += (srcp[i][x] + offset) * weightsV[i];
            }

            temp[x] = sum;
        }

        for (int i = 1; i <= radiusH; ++i)
        {
            temp[-i] = temp[-1 + i];
            temp[width - 1 + i] = temp[width - i];
        }

        for (int x = 0; x < width; ++x)
        {
            float sum = 0.f;

            for (int i = -radiusH; i <= radiusH; ++i)
                sum += temp[x + i] * weightsH[i];

            dstp[x] = sum;
        }

        for (int i = 0; i < diameter - 1; ++i)
            srcp[i] = srcp[i + 1];
        if (y < height - 1 - radiusV)
            srcp[diameter - 1] += srcStride;
        else if (y > height - 1 - radiusV)
            srcp[diameter - 1] -= srcStride;
        dstp += dstStride;
    }

    delete[] srcp;
}

template<typename T>
static void gaussianBlurV(const T* _srcp, float* __restrict dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius, const float offset) noexcept
{
    const int diameter = radius * 2 + 1;
    const T** srcp = new const T * [diameter];

    srcp[radius] = _srcp;
    for (int i = 1; i <= radius; ++i)
    {
        srcp[radius - i] = srcp[radius - 1 + i];
        srcp[radius + i] = srcp[radius] + srcStride * i;
    }

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float sum = 0.f;

            for (int i = 0; i < diameter; ++i)
            {
                if constexpr (std::is_integral_v<T>)
                    sum += srcp[i][x] * weights[i];
                else
                    sum += (srcp[i][x] + offset) * weights[i];
            }

            dstp[x] = sum;
        }

        for (int i = 0; i < diameter - 1; ++i)
            srcp[i] = srcp[i + 1];
        if (y < height - 1 - radius)
            srcp[diameter - 1] += srcStride;
        else if (y > height - 1 - radius)
            srcp[diameter - 1] -= srcStride;
        dstp += dstStride;
    }

    delete[] srcp;
}

template<typename T>
static void gaussianBlurH(const T* srcp, float* __restrict temp, float* __restrict dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius, const float offset) noexcept
{
    weights += radius;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                temp[x] = srcp[x];
            else
                temp[x] = srcp[x] + offset;
        }

        for (int i = 1; i <= radius; ++i)
        {
            temp[-i] = temp[-1 + i];
            temp[width - 1 + i] = temp[width - i];
        }

        for (int x = 0; x < width; ++x)
        {
            float sum = 0.f;

            for (int i = -radius; i <= radius; ++i)
                sum += temp[x + i] * weights[i];

            dstp[x] = sum;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

static void detectEdge(float* blur, float* __restrict gradient, unsigned* __restrict direction, const int width, const int height, const int stride, const int bgStride, const int mode, const int op) noexcept
{
    float* __restrict srcpp = blur;
    float* __restrict srcp = blur;
    float* __restrict srcpn = blur + bgStride;

    srcp[-1] = srcp[0];
    srcp[width] = srcp[width - 1];

    for (int y = 0; y < height; ++y)
    {
        srcpn[-1] = srcpn[0];
        srcpn[width] = srcpn[width - 1];

        for (int x = 0; x < width; ++x)
        {
            float gx, gy;

            if (op == 0)
            {
                gx = srcp[x + 1] - srcp[x - 1];
                gy = srcpp[x] - srcpn[x];
            }
            else if (op == 1)
            {
                gx = (srcpp[x + 1] + srcp[x + 1] + srcpn[x + 1] - srcpp[x - 1] - srcp[x - 1] - srcpn[x - 1]) / 2.f;
                gy = (srcpp[x - 1] + srcpp[x] + srcpp[x + 1] - srcpn[x - 1] - srcpn[x] - srcpn[x + 1]) / 2.f;
            }
            else if (op == 2)
            {
                gx = srcpp[x + 1] + 2.f * srcp[x + 1] + srcpn[x + 1] - srcpp[x - 1] - 2.f * srcp[x - 1] - srcpn[x - 1];
                gy = srcpp[x - 1] + 2.f * srcpp[x] + srcpp[x + 1] - srcpn[x - 1] - 2.f * srcpn[x] - srcpn[x + 1];
            }
            else
            {
                gx = 3.f * srcpp[x + 1] + 10.f * srcp[x + 1] + 3.f * srcpn[x + 1] - 3.f * srcpp[x - 1] - 10.f * srcp[x - 1] - 3.f * srcpn[x - 1];
                gy = 3.f * srcpp[x - 1] + 10.f * srcpp[x] + 3.f * srcpp[x + 1] - 3.f * srcpn[x - 1] - 10.f * srcpn[x] - 3.f * srcpn[x + 1];
            }

            gradient[x] = std::sqrt(gx * gx + gy * gy);

            if (mode == 0) {
                float dr = std::atan2(gy, gx);
                if (dr < 0.f)
                    dr += M_PIF;

                const unsigned bin = static_cast<unsigned>(dr * 4.f * M_1_PIF + 0.5f);
                direction[x] = (bin >= 4) ? 0 : bin;
            }
        }

        srcpp = srcp;
        srcp = srcpn;
        if (y < height - 2)
            srcpn += bgStride;
        gradient += bgStride;
        direction += stride;
    }
}

static void nonMaximumSuppression(const unsigned* direction, float* __restrict gradient, float* __restrict blur, const int width, const int height, const int stride, const int bgStride, const int radiusAlign) noexcept
{
    const int offsets[] = { 1, -bgStride + 1, -bgStride, -bgStride - 1 };

    gradient[-1] = gradient[0];
    gradient[-1 + bgStride * (height - 1)] = gradient[bgStride * (height - 1)];
    gradient[width] = gradient[width - 1];
    gradient[width + bgStride * (height - 1)] = gradient[width - 1 + bgStride * (height - 1)];
    std::copy_n(gradient - radiusAlign, width + radiusAlign * 2, gradient - radiusAlign - bgStride);
    std::copy_n(gradient - radiusAlign + bgStride * (height - static_cast<int64_t>(1)), width + radiusAlign * 2, gradient - radiusAlign + bgStride * static_cast<int64_t>(height));

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            const int offset = offsets[direction[x]];
            blur[x] = (gradient[x] >= std::max(gradient[x + offset], gradient[x - offset])) ? gradient[x] : fltLowest;
        }

        direction += stride;
        gradient += bgStride;
        blur += bgStride;
    }
}

template<typename T>
static void outputGB(const float* srcp, T* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = static_cast<T>(std::min(static_cast<unsigned>(lrintf(srcp[x])), static_cast<unsigned>(peak)));
            else
                dstp[x] = srcp[x] - offset;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void binarizeCE(const float* srcp, T* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const uint16_t peak, const float lower, const float upper) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = static_cast<T>((srcp[x] == fltMax) ? peak : 0);
            else
                dstp[x] = (srcp[x] == fltMax) ? upper : lower;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void discretizeGM(const float* srcp, T* __restrict dstp, const int width, const int height, const int srcStride, const int dstStride, const float magnitude, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if constexpr (std::is_integral_v<T>)
                dstp[x] = static_cast<T>(std::min(static_cast<unsigned>(lrintf(srcp[x] * magnitude)), static_cast<unsigned>(peak)));
            else
                dstp[x] = srcp[x] * magnitude - offset;
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void vsTCanny::filter_c(PVideoFrame& src, PVideoFrame& dst, const vsTCanny* const __restrict, IScriptEnvironment* env) noexcept
{
    int planes_y[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    int planes_r[3] = { PLANAR_G, PLANAR_B, PLANAR_R };
    const int* current_planes = vi.IsRGB() ? planes_r : planes_y;
    for (int i = 0; i < planecount; ++i)
    {
        const int plane = current_planes[i];
        const int height = src->GetHeight(plane);

        if (process[i] == 3)
        {
            const int stride = src->GetPitch(plane) / sizeof(T);
            const int bgStride = stride + radiusAlign * 2;
            const int dst_stride = dst->GetPitch(plane) / sizeof(T);
            const int width = src->GetRowSize(plane) / sizeof(T);
            const T* srcp = reinterpret_cast<const T*>(src->GetReadPtr(plane));
            T* dstp = reinterpret_cast<T*>(dst->GetWritePtr(plane));

            float* blur = vsTCanny::blur + radiusAlign;
            float* gradient = vsTCanny::gradient + bgStride + radiusAlign;

            if (radiusV[i] && radiusH[i])
                gaussianBlur(srcp, gradient, blur, weightsH[i], weightsV[i], width, height, stride, bgStride, radiusH[i], radiusV[i], offset[i]);
            else if (radiusV[i])
                gaussianBlurV(srcp, blur, weightsV[i], width, height, stride, bgStride, radiusV[i], offset[i]);
            else if (radiusH[i])
                gaussianBlurH(srcp, gradient, blur, weightsH[i], width, height, stride, bgStride, radiusH[i], offset[i]);
            else
                copyPlane(srcp, blur, width, height, stride, bgStride, offset[i]);

            if (mode_ != -1)
            {
                detectEdge(blur, gradient, direction, width, height, stride, bgStride, mode_, op_);

                if (mode_ == 0)
                {
                    nonMaximumSuppression(direction, gradient, blur, width, height, stride, bgStride, radiusAlign);
                    hysteresis(blur, found, width, height, bgStride, t_h_, t_l_);
                }
            }

            if (mode_ == -1)
                outputGB(blur, dstp, width, height, bgStride, dst_stride, peak, offset[i]);
            else if (mode_ == 0)
                binarizeCE(blur, dstp, width, height, bgStride, dst_stride, peak, lower[i], upper[i]);
            else
                discretizeGM(gradient, dstp, width, height, bgStride, dst_stride, magnitude, peak, offset[i]);
        }
        else if (process[i] == 2)
            env->BitBlt(dst->GetWritePtr(plane), dst->GetPitch(plane), src->GetReadPtr(plane), src->GetPitch(plane), src->GetRowSize(plane), height);
    }
}

static inline float* gaussianWeights(const float sigma, int& radius) noexcept
{
    const int diameter = std::max(static_cast<int>(lrintf(sigma * 3.f)), 1) * 2 + 1;
    radius = diameter / 2;

    float* __restrict weights = new (std::nothrow) float[diameter];

    if (!weights)
        return nullptr;

    float sum = 0.f;

    for (int k = -radius; k <= radius; ++k)
    {
        const float w = std::exp(-(k * k) / (2.f * sigma * sigma));
        weights[k + radius] = w;
        sum += w;
    }
    
    for (int k = 0; k < diameter; ++k)
        weights[k] /= sum;

    return weights;
}

vsTCanny::vsTCanny(PClip _child, float sigmaY, float sigmaU, float sigmaV, float sigma_vY, float sigma_vU, float sigma_vV, float t_h, float t_l, int mode, int op, float gmmax, int y, int u, int v, int opt, IScriptEnvironment* env)
    : GenericVideoFilter(_child), t_h_(t_h), t_l_(t_l), mode_(mode), op_(op), opt_(opt)
{
    if (!vi.IsPlanar())
        env->ThrowError("vsTCanny: the clip is not in planar format.");

    const int height = vi.height;

    if (height < 2)
        env->ThrowError("vsTCanny: the clip's height must be greater than or equal to 2.");
    if (t_l_ >= t_h_)
        env->ThrowError("vsTCanny: t_h must be greater than t_l.");
    if (mode_ < -1 || mode_ > 1)
        env->ThrowError("vsTCanny: mode must be -1, 0, or 1.");
    if (op_ < 0 || op_ > 3)
        env->ThrowError("vsTCanny: op must be 0, 1, 2, or 3.");
    if (gmmax < 1.f)
        env->ThrowError("vsTCanny: gmmax must be greater than or equal to 1.0.");
    if (opt_ < -1 || opt_ > 2)
        env->ThrowError("vsTCanny: opt must be between -1..2.");
    if (!(env->GetCPUFlags() & CPUF_AVX2) && opt_ == 2)
        env->ThrowError("vsTCanny: opt=2 requires AVX2.");
    if (!(env->GetCPUFlags() & CPUF_SSE2) && opt_ == 1)
        env->ThrowError("vsTCanny: opt=1 requires SSE2.");

    const float sigmaH[3] = { sigmaY, sigmaU, sigmaV };
    const float sigmaV_[3] = { sigma_vY, sigma_vU, sigma_vV };
    const std::string planeOrder[3] = { "first", "second", "third" };
    const std::string sigmaOrder[3] = { "sigmaY", "sigmaU", "sigmaV" };
    const std::string sigmaVOrder[3] = { "sigma_vY", "sigma_vU", "sigma_vV" };

    const int planes[3] = { y, u, v };
    planecount = std::min(vi.NumComponents(), 3);
    for (int i = 0; i < planecount; ++i)
    {
        if (vi.IsRGB())
            process[i] = 3;
        else
        {
            switch (planes[i])
            {
                case 3: process[i] = 3; break;
                case 2: process[i] = 2; break;
                default: process[i] = 1; break;
            }
        }

        if (sigmaH[i] < 0.f)
            env->ThrowError(std::string("vsTCanny: " + sigmaOrder[i] + " must be greater than or equal to 0.0.").c_str());
        if (sigmaV_[i] < 0.f)
            env->ThrowError(std::string("vsTCanny: " + sigmaVOrder[i] + " must be greater than or equal to 0.0.").c_str());
        if (planes[i] < 1 || planes[i] > 3)
            env->ThrowError("vsTCanny: y, u, v must be between 1..3.");

        if (i == 0 || vi.IsRGB())
        {
            offset[i] = 0.f;
            lower[i] = 0.f;
            upper[i] = 1.f;
        }
        else
        {
            offset[i] = 0.5f;
            lower[i] = -0.5f;
            upper[i] = 0.5f;
        }

        if (process[i] == 3)
        {
            if (sigmaH[i])
            {
                weightsH[i] = gaussianWeights(sigmaH[i], radiusH[i]);
                if (!weightsH[i])
                    env->ThrowError("vsTCanny: malloc failure (weightsH).");

                const int width = vi.width >> ((i && !vi.IsRGB()) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0);
                
                if (width < radiusH[i] + 1)
                    env->ThrowError(std::string("vsTCanny: the " + planeOrder[i] + " plane's width must be greater than or equal to " + std::to_string(radiusH[i] + 1) + " for specified sigma.").c_str());
            }

            if (sigmaV_[i])
            {
                weightsV[i] = gaussianWeights(sigmaV_[i], radiusV[i]);
                if (!weightsV[i])
                    env->ThrowError("vsTCanny: malloc failure (weightsV).");

                const int height_ = height >> ((i && !vi.IsRGB()) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0);

                if (height_ < radiusV[i] + 1)
                    env->ThrowError(std::string("vsTCanny: the " + planeOrder[i] + " plane's height must be greater than or equal to " + std::to_string(radiusV[i] + 1) + " for specified sigma_v.").c_str());
            }
        }
    }

    const int comp_size = vi.ComponentSize();

    if (comp_size != 4)
    {
        peak = (1 << vi.BitsPerComponent()) - 1;
        const float scale = peak / 255.f;
        t_h_ *= scale;
        t_l_ *= scale;
    }
    else
    {
        peak = 0;
        t_h_ /= 255.f;
        t_l_ /= 255.f;
    }

    avx2 = ((!!(env->GetCPUFlags() & CPUF_AVX2) && opt_ < 0) || opt_ == 2);
    sse2 = ((!!(env->GetCPUFlags() & CPUF_SSE2) && opt_ < 0) || opt_ == 1);

    if (avx2)
    {
        vectorSize = 8;
        alignment = 32;
    }
    else if (sse2)
    {
        vectorSize = 4;
        alignment = 16;
    }
    else
    {
        vectorSize = 1;
        alignment = 4;
    }

    radiusAlign = (std::max({ radiusH[0], radiusH[1], radiusH[2], 1 }) + vectorSize - 1) & -vectorSize;
    magnitude = 255.f / gmmax;

    PVideoFrame clip = child->GetFrame(0, env);
    const int pitch = clip->GetPitch(PLANAR_Y);

    blur = reinterpret_cast<float*>(aligned_malloc((pitch / comp_size + radiusAlign * static_cast<int64_t>(2)) * height * sizeof(float), alignment));
    if (!blur)
        env->ThrowError("vsTCanny: malloc failure (blur).");

    gradient = reinterpret_cast<float*>(aligned_malloc((pitch / comp_size + radiusAlign * static_cast<int64_t>(2)) * (height + static_cast<int64_t>(2)) * sizeof(float), alignment));
    if (!gradient)
        env->ThrowError("vsTCanny: malloc failure (gradient).");

    if (mode_ == 0)
    {
        direction = reinterpret_cast<unsigned*>(aligned_malloc(pitch / comp_size * static_cast<int64_t>(height) * sizeof(unsigned), alignment));
        if (!direction)
            env->ThrowError("vsTCanny: malloc failure (direction).");

        found = new (std::nothrow) bool[vi.width * static_cast<int64_t>(height)];
        if (!found)
            env->ThrowError("vsTCanny: malloc failure (found).");
    }
    else
    {
        direction = nullptr;
        found = nullptr;
    }

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; }
}

vsTCanny::~vsTCanny()
{
    for (int i = 0; i < planecount; ++i)
    {
        if (process[i] == 3)
        {
            delete[] weightsH[i];
            delete[] weightsV[i];
        }
    }

    aligned_free(blur);
    aligned_free(gradient);
    aligned_free(direction);
    delete[] found;
}

PVideoFrame __stdcall vsTCanny::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src = child->GetFrame(n, env);
    PVideoFrame dst = has_at_least_v8 ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi);

    if (avx2)
    {
        switch (vi.ComponentSize())
        {
            case 1: filter_avx2<uint8_t>(src, dst, 0, env); break;
            case 2: filter_avx2<uint16_t>(src, dst, 0, env); break;
            default: filter_avx2<float>(src, dst, 0, env); break;
        }
    }
    else if (sse2)
    {
        switch (vi.ComponentSize())
        {
            case 1: filter_sse2<uint8_t>(src, dst, 0, env); break;
            case 2: filter_sse2<uint16_t>(src, dst, 0, env); break;
            default: filter_sse2<float>(src, dst, 0, env); break;
        }
    }
    else
    {
        switch (vi.ComponentSize())
        {
            case 1: filter_c<uint8_t>(src, dst, 0, env); break;
            case 2: filter_c<uint16_t>(src, dst, 0, env); break;
            default: filter_c<float>(src, dst, 0, env); break;
        }
    }

    return dst;
}

AVSValue __cdecl Create_vsTCanny(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    PClip clip = args[0].AsClip();
    const VideoInfo& vi = clip->GetVideoInfo();
    const float sigmaY = args[1].AsFloatf(1.5f);
    const float sigma_vY = args[4].AsFloatf(1.5f);

    const float sigmaU = [&]() {
        if (vi.NumComponents() > 1)
            return (vi.IsRGB() || vi.Is444()) ? args[2].AsFloatf(sigmaY) : args[2].AsFloatf(sigmaY / (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)));
        else
            return 0.0f;
    }();

    const float sigma_vU = [&]() {
        if (vi.NumComponents() > 1)
            return (vi.IsRGB() || vi.Is444()) ? args[5].AsFloatf(sigma_vY) : args[5].AsFloatf(sigma_vY / (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)));
        else
            return 0.0f;
    }();

    return new vsTCanny(
        clip,
        sigmaY,
        sigmaU,
        args[3].AsFloatf(sigmaU),
        sigma_vY,
        sigma_vU,
        args[6].AsFloatf(sigma_vU),
        args[7].AsFloatf(8.0),
        args[8].AsFloatf(1.0),
        args[9].AsInt(0),
        args[10].AsInt(1),
        args[11].AsFloatf(50.0),
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

    env->AddFunction("vsTCanny", "c[sigmaY]f[sigmaU]f[sigmaV]f[sigma_vY]f[sigma_vU]f[sigma_vV]f[t_h]f[t_l]f[mode]i[op]i[gmmax]f[y]i[u]i[v]i[opt]i", Create_vsTCanny, 0);

    return "vsTCanny";
}
