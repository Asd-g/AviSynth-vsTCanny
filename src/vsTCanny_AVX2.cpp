#include "vsTCanny.h"
#include "VCL2/vectorclass.h"
#include "VCL2/vectormath_trig.h"

template<typename T>
static void copyPlane(const T* srcp, float* dstp, const int width, const int height, const int srcStride, const int dstStride, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            if constexpr (std::is_same_v<T, uint8_t>)
                to_float(Vec8i().load_8uc(srcp + x)).store_nt(dstp + x);
            else if constexpr (std::is_same_v<T, uint16_t>)
                to_float(Vec8i().load_8us(srcp + x)).store_nt(dstp + x);
            else
                (Vec8f().load_a(reinterpret_cast<const float*>(srcp + x)) + offset).store_nt(dstp + x);
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlur(const T* __srcp, float* temp, float* dstp, const float* weightsH, const float* weightsV, const int width, const int height,
    const int srcStride, const int dstStride, const int radiusH, const int radiusV, const float offset) noexcept
{
    const int diameter = radiusV * 2 + 1;
    const T** _srcp = new const T * [diameter];

    _srcp[radiusV] = __srcp;
    for (int i = 1; i <= radiusV; ++i)
    {
        _srcp[radiusV - i] = _srcp[radiusV - 1 + i];
        _srcp[radiusV + i] = _srcp[radiusV] + srcStride * static_cast<int64_t>(i);
    }

    weightsH += radiusH;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            Vec8f sum = zero_8f();

            for (int i = 0; i < diameter; ++i)
            {
                if constexpr (std::is_same_v<T, uint8_t>)
                {
                    const Vec8f srcp = to_float(Vec8i().load_8uc(_srcp[i] + x));
                    sum = mul_add(srcp, weightsV[i], sum);
                }
                else if constexpr (std::is_same_v<T, uint16_t>)
                {
                    const Vec8f srcp = to_float(Vec8i().load_8us(_srcp[i] + x));
                    sum = mul_add(srcp, weightsV[i], sum);
                }
                else
                {
                    const Vec8f srcp = Vec8f().load_a(reinterpret_cast<const float*>(_srcp[i] + x));
                    sum = mul_add(srcp + offset, weightsV[i], sum);
                }
            }

            sum.store_a(temp + x);
        }

        for (int i = 1; i <= radiusH; ++i)
        {
            temp[-i] = temp[-1 + i];
            temp[width - 1 + i] = temp[width - i];
        }

        for (int x = 0; x < width; x += 8)
        {
            Vec8f sum = zero_8f();

            for (int i = -radiusH; i <= radiusH; ++i)
            {
                const Vec8f srcp = Vec8f().load(temp + x + i);
                sum = mul_add(srcp, weightsH[i], sum);
            }

            sum.store_nt(dstp + x);
        }

        for (int i = 0; i < diameter - 1; ++i)
            _srcp[i] = _srcp[i + 1];
        if (y < height - 1 - radiusV)
            _srcp[diameter - 1] += srcStride;
        else if (y > height - 1 - radiusV)
            _srcp[diameter - 1] -= srcStride;
        dstp += dstStride;
    }

    delete[] _srcp;
}

template<typename T>
static void gaussianBlurV(const T* __srcp, float* dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride,
    const int radius, const float offset) noexcept
{
    const int diameter = radius * 2 + 1;
    const T** _srcp = new const T * [diameter];

    _srcp[radius] = __srcp;
    for (int i = 1; i <= radius; ++i)
    {
        _srcp[radius - i] = _srcp[radius - 1 + i];
        _srcp[radius + i] = _srcp[radius] + srcStride * i;
    }

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            Vec8f sum = zero_8f();

            for (int i = 0; i < diameter; ++i)
            {
                if constexpr (std::is_same_v<T, uint8_t>)
                {
                    const Vec8f srcp = to_float(Vec8i().load_8uc(_srcp[i] + x));
                    sum = mul_add(srcp, weights[i], sum);
                }
                else if constexpr (std::is_same_v<T, uint16_t>)
                {
                    const Vec8f srcp = to_float(Vec8i().load_8us(_srcp[i] + x));
                    sum = mul_add(srcp, weights[i], sum);
                }
                else
                {
                    const Vec8f srcp = Vec8f().load_a(reinterpret_cast<const float*>(_srcp[i] + x));
                    sum = mul_add(srcp + offset, weights[i], sum);
                }
            }

            sum.store_nt(dstp + x);
        }

        for (int i = 0; i < diameter - 1; ++i)
            _srcp[i] = _srcp[i + 1];
        if (y < height - 1 - radius)
            _srcp[diameter - 1] += srcStride;
        else if (y > height - 1 - radius)
            _srcp[diameter - 1] -= srcStride;
        dstp += dstStride;
    }

    delete[] _srcp;
}

template<typename T>
static void gaussianBlurH(const T* _srcp, float* temp, float* dstp, const float* weights, const int width, const int height,
    const int srcStride, const int dstStride, const int radius, const float offset) noexcept
{
    weights += radius;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            if constexpr (std::is_same_v<T, uint8_t>)
                to_float(Vec8i().load_8uc(_srcp + x)).store_a(temp + x);
            else if (std::is_same_v<T, uint16_t>)
                to_float(Vec8i().load_8us(_srcp + x)).store_a(temp + x);
            else
                (Vec8f().load_a(reinterpret_cast<const float*>(_srcp + x)) + offset).store_a(temp + x);
        }

        for (int i = 1; i <= radius; ++i)
        {
            temp[-i] = temp[-1 + i];
            temp[width - 1 + i] = temp[width - i];
        }

        for (int x = 0; x < width; x += 8)
        {
            Vec8f sum = zero_8f();

            for (int i = -radius; i <= radius; ++i)
            {
                const Vec8f srcp = Vec8f().load(temp + x + i);
                sum = mul_add(srcp, weights[i], sum);
            }

            sum.store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

static void detectEdge(float* blur, float* gradient, unsigned* direction, const int width, const int height, const int stride, const int bgStride,
    const int mode, const int op) noexcept
{
    float* srcpp = blur;
    float* srcp = blur;
    float* srcpn = blur + bgStride;

    srcp[-1] = srcp[0];
    srcp[width] = srcp[width - 1];

    for (int y = 0; y < height; ++y)
    {
        srcpn[-1] = srcpn[0];
        srcpn[width] = srcpn[width - 1];

        for (int x = 0; x < width; x += 8)
        {
            const Vec8f topLeft = Vec8f().load(srcpp + x - 1);
            const Vec8f top = Vec8f().load_a(srcpp + x);
            const Vec8f topRight = Vec8f().load(srcpp + x + 1);
            const Vec8f left = Vec8f().load(srcp + x - 1);
            const Vec8f right = Vec8f().load(srcp + x + 1);
            const Vec8f bottomLeft = Vec8f().load(srcpn + x - 1);
            const Vec8f bottom = Vec8f().load_a(srcpn + x);
            const Vec8f bottomRight = Vec8f().load(srcpn + x + 1);

            Vec8f gx, gy;

            if (op == 0)
            {
                gx = right - left;
                gy = top - bottom;
            }
            else if (op == 1)
            {
                gx = (topRight + right + bottomRight - topLeft - left - bottomLeft) * 0.5f;
                gy = (topLeft + top + topRight - bottomLeft - bottom - bottomRight) * 0.5f;
            }
            else if (op == 2)
            {
                gx = topRight + mul_add(2.f, right, bottomRight) - topLeft - mul_add(2.f, left, bottomLeft);
                gy = topLeft + mul_add(2.f, top, topRight) - bottomLeft - mul_add(2.f, bottom, bottomRight);
            }
            else
            {
                gx = mul_add(3.f, topRight, mul_add(10.f, right, 3.f * bottomRight)) - mul_add(3.f, topLeft, mul_add(10.f, left, 3.f * bottomLeft));
                gy = mul_add(3.f, topLeft, mul_add(10.f, top, 3.f * topRight)) - mul_add(3.f, bottomLeft, mul_add(10.f, bottom, 3.f * bottomRight));
            }

            sqrt(mul_add(gx, gx, gy * gy)).store_nt(gradient + x);

            if (mode == 0)
            {
                Vec8f dr = atan2(gy, gx);
                dr = if_add(dr < 0.f, dr, M_PIF);

                const Vec8ui bin = Vec8ui(truncatei(mul_add(dr, 4.f * M_1_PIF, 0.5f)));
                select(bin >= 4, zero_si256(), bin).store_nt(direction + x);
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

static void nonMaximumSuppression(const unsigned* _direction, float* _gradient, float* blur, const int width, const int height,
    const int stride, const int bgStride, const int radiusAlign) noexcept
{
    _gradient[-1] = _gradient[0];
    _gradient[-1 + bgStride * (height - 1)] = _gradient[bgStride * (height - 1)];
    _gradient[width] = _gradient[width - 1];
    _gradient[width + bgStride * (height - 1)] = _gradient[width - 1 + bgStride * (height - 1)];
    std::copy_n(_gradient - radiusAlign, width + radiusAlign * 2, _gradient - radiusAlign - bgStride);
    std::copy_n(_gradient - radiusAlign + bgStride * (height - static_cast<int64_t>(1)), width + radiusAlign * 2, _gradient - radiusAlign + bgStride * static_cast<int64_t>(height));

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            const Vec8ui direction = Vec8ui().load_a(_direction + x);

            Vec8fb mask = Vec8fb(direction == 0);
            Vec8f gradient = max(Vec8f().load(_gradient + x + 1), Vec8f().load(_gradient + x - 1));
            Vec8f result = gradient & mask;

            mask = Vec8fb(direction == 1);
            gradient = max(Vec8f().load(_gradient + x - bgStride + 1), Vec8f().load(_gradient + x + bgStride - 1));
            result |= gradient & mask;

            mask = Vec8fb(direction == 2);
            gradient = max(Vec8f().load_a(_gradient + x - bgStride), Vec8f().load_a(_gradient + x + bgStride));
            result |= gradient & mask;

            mask = Vec8fb(direction == 3);
            gradient = max(Vec8f().load(_gradient + x - bgStride - 1), Vec8f().load(_gradient + x + bgStride + 1));
            result |= gradient & mask;

            gradient = Vec8f().load_a(_gradient + x);
            select(gradient >= result, gradient, fltLowest).store_nt(blur + x);
        }

        _direction += stride;
        _gradient += bgStride;
        blur += bgStride;
    }
}

template<typename T> static void outputGB(const float*, T*, const int, const int, const int, const int, const uint16_t, const float) noexcept;

template<>
void outputGB(const float* _srcp, uint8_t* dstp, const int width, const int height, const int srcStride, const int dstStride, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 32)
        {
            const Vec8i srcp_8i_0 = truncatei(Vec8f().load_a(_srcp + x) + 0.5f);
            const Vec8i srcp_8i_1 = truncatei(Vec8f().load_a(_srcp + x + 8) + 0.5f);
            const Vec8i srcp_8i_2 = truncatei(Vec8f().load_a(_srcp + x + 16) + 0.5f);
            const Vec8i srcp_8i_3 = truncatei(Vec8f().load_a(_srcp + x + 24) + 0.5f);
            const Vec16s srcp_16s_0 = compress_saturated(srcp_8i_0, srcp_8i_1);
            const Vec16s srcp_16s_1 = compress_saturated(srcp_8i_2, srcp_8i_3);
            const Vec32uc srcp = compress_saturated_s2u(srcp_16s_0, srcp_16s_1);
            srcp.store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void outputGB(const float* _srcp, uint16_t* dstp, const int width, const int height, const int srcStride, const int dstStride, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 16)
        {
            const Vec8i srcp_8i_0 = truncatei(Vec8f().load_a(_srcp + x) + 0.5f);
            const Vec8i srcp_8i_1 = truncatei(Vec8f().load_a(_srcp + x + 8) + 0.5f);
            const Vec16us srcp = compress_saturated_s2u(srcp_8i_0, srcp_8i_1);
            min(srcp, peak).store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void outputGB(const float* _srcp, float* dstp, const int width, const int height, const int srcStride, const int dstStride, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            const Vec8f srcp = Vec8f().load_a(_srcp + x);
            (srcp - offset).store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T> static void binarizeCE(const float*, T*, const int, const int, const int, const int, const uint16_t, const float, const float) noexcept;

template<>
void binarizeCE(const float* srcp, uint8_t* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const uint16_t peak, const float lower, const float upper) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 32)
        {
            const Vec8ib mask_8ib_0 = Vec8ib(Vec8f().load_a(srcp + x) == fltMax);
            const Vec8ib mask_8ib_1 = Vec8ib(Vec8f().load_a(srcp + x + 8) == fltMax);
            const Vec8ib mask_8ib_2 = Vec8ib(Vec8f().load_a(srcp + x + 16) == fltMax);
            const Vec8ib mask_8ib_3 = Vec8ib(Vec8f().load_a(srcp + x + 24) == fltMax);
            const Vec16sb mask_16sb_0 = Vec16sb(compress_saturated(mask_8ib_0, mask_8ib_1));
            const Vec16sb mask_16sb_1 = Vec16sb(compress_saturated(mask_8ib_2, mask_8ib_3));
            const Vec32cb mask = Vec32cb(compress_saturated(mask_16sb_0, mask_16sb_1));
            select(mask, Vec32uc(255), zero_si256()).store_nt(dstp + x);
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void binarizeCE(const float* srcp, uint16_t* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const uint16_t peak, const float lower, const float upper) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 16)
        {
            const Vec8ib mask_8ib_0 = Vec8ib(Vec8f().load_a(srcp + x) == fltMax);
            const Vec8ib mask_8ib_1 = Vec8ib(Vec8f().load_a(srcp + x + 8) == fltMax);
            const Vec16sb mask = Vec16sb(compress_saturated(mask_8ib_0, mask_8ib_1));
            select(mask, Vec16us(peak), zero_si256()).store_nt(dstp + x);
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void binarizeCE(const float* srcp, float* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const uint16_t peak, const float lower, const float upper) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            const Vec8fb mask = (Vec8f().load_a(srcp + x) == fltMax);
            select(mask, Vec8f(upper), Vec8f(lower)).store_nt(dstp + x);
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T> static void discretizeGM(const float*, T*, const int, const int, const int, const int, const float, const uint16_t, const float) noexcept;

template<>
void discretizeGM(const float* _srcp, uint8_t* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const float magnitude, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 32)
        {
            const Vec8f srcp_8f_0 = Vec8f().load_a(_srcp + x);
            const Vec8f srcp_8f_1 = Vec8f().load_a(_srcp + x + 8);
            const Vec8f srcp_8f_2 = Vec8f().load_a(_srcp + x + 16);
            const Vec8f srcp_8f_3 = Vec8f().load_a(_srcp + x + 24);
            const Vec8i srcp_8i_0 = truncatei(mul_add(srcp_8f_0, magnitude, 0.5f));
            const Vec8i srcp_8i_1 = truncatei(mul_add(srcp_8f_1, magnitude, 0.5f));
            const Vec8i srcp_8i_2 = truncatei(mul_add(srcp_8f_2, magnitude, 0.5f));
            const Vec8i srcp_8i_3 = truncatei(mul_add(srcp_8f_3, magnitude, 0.5f));
            const Vec16s srcp_16s_0 = compress_saturated(srcp_8i_0, srcp_8i_1);
            const Vec16s srcp_16s_1 = compress_saturated(srcp_8i_2, srcp_8i_3);
            const Vec32uc srcp = compress_saturated_s2u(srcp_16s_0, srcp_16s_1);
            srcp.store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void discretizeGM(const float* _srcp, uint16_t* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const float magnitude, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 16)
        {
            const Vec8f srcp_8f_0 = Vec8f().load_a(_srcp + x);
            const Vec8f srcp_8f_1 = Vec8f().load_a(_srcp + x + 8);
            const Vec8i srcp_8i_0 = truncatei(mul_add(srcp_8f_0, magnitude, 0.5f));
            const Vec8i srcp_8i_1 = truncatei(mul_add(srcp_8f_1, magnitude, 0.5f));
            const Vec16us srcp = compress_saturated_s2u(srcp_8i_0, srcp_8i_1);
            min(srcp, peak).store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<>
void discretizeGM(const float* _srcp, float* dstp, const int width, const int height, const int srcStride, const int dstStride,
    const float magnitude, const uint16_t peak, const float offset) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x += 8)
        {
            const Vec8f srcp = Vec8f().load_a(_srcp + x);
            mul_sub(srcp, magnitude, offset).store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void vsTCanny::filter_avx2(PVideoFrame& src, PVideoFrame& dst, const vsTCanny* const __restrict, IScriptEnvironment* env) noexcept
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

template void vsTCanny::filter_avx2<uint8_t>(PVideoFrame& src, PVideoFrame& dst, const vsTCanny* const __restrict, IScriptEnvironment* env) noexcept;
template void vsTCanny::filter_avx2<uint16_t>(PVideoFrame& src, PVideoFrame& dst, const vsTCanny* const __restrict, IScriptEnvironment* env) noexcept;
template void vsTCanny::filter_avx2<float>(PVideoFrame& src, PVideoFrame& dst, const vsTCanny* const __restrict, IScriptEnvironment* env) noexcept;
