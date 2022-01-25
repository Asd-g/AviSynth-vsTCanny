#include "VCL2/vectormath_trig.h"
#include "vsTCanny.h"

template<typename T>
static void copyPlane(const T* srcp, float* dstp, const int width, const int height, const int srcStride, const int dstStride) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            if constexpr (std::is_same_v<T, uint8_t>)
                to_float(Vec16i().load_16uc(srcp + x)).store_nt(dstp + x);
            else if constexpr (std::is_same_v<T, uint16_t>)
                to_float(Vec16i().load_16us(srcp + x)).store_nt(dstp + x);
            else
                Vec16f().load_a(srcp + x).store_nt(dstp + x);
        }

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlur(const T* __srcp, float* temp, float* dstp, const float* weightsH, const float* weightsV, const int width, const int height, const int srcStride, const int dstStride, const int radiusH, const int radiusV) noexcept
{
    const int diameter{ radiusV * 2 + 1 };
    std::unique_ptr<const T* []> _srcp{ std::make_unique<const T* []>(diameter) };

    _srcp[radiusV] = __srcp;
    for (int i{ 1 }; i <= radiusV; ++i)
        _srcp[radiusV - i] = _srcp[radiusV + i] = _srcp[radiusV] + srcStride * i;

    weightsH += radiusH;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            Vec16f sum{ zero_16f() };

            for (int i{ 0 }; i < diameter; ++i)
            {
                if constexpr (std::is_same_v<T, uint8_t>)
                {
                    const Vec16f srcp{ to_float(Vec16i().load_16uc(_srcp[i] + x)) };
                    sum = mul_add(srcp, weightsV[i], sum);
                }
                else if constexpr (std::is_same_v<T, uint16_t>)
                {
                    const Vec16f srcp{ to_float(Vec16i().load_16us(_srcp[i] + x)) };
                    sum = mul_add(srcp, weightsV[i], sum);
                }
                else
                {
                    const Vec16f srcp{ Vec16f().load_a(_srcp[i] + x) };
                    sum = mul_add(srcp, weightsV[i], sum);
                }
            }

            sum.store_a(temp + x);
        }

        for (int i{ 1 }; i <= radiusH; ++i)
        {
            temp[-i] = temp[i];
            temp[width - 1 + i] = temp[width - 1 - i];
        }

        for (int x{ 0 }; x < width; x += 16)
        {
            Vec16f sum{ zero_16f() };

            for (int i{ -radiusH }; i <= radiusH; ++i)
            {
                const Vec16f srcp{ Vec16f().load(temp + x + i) };
                sum = mul_add(srcp, weightsH[i], sum);
            }

            sum.store_nt(dstp + x);
        }

        for (int i{ 0 }; i < diameter - 1; ++i)
            _srcp[i] = _srcp[i + 1];

        _srcp[diameter - 1] += (y < height - 1 - radiusV) ? srcStride : -srcStride;

        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlurV(const T* __srcp, float* dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius) noexcept
{
    const int diameter{ radius * 2 + 1 };
    std::unique_ptr<const T* []> _srcp{ std::make_unique<const T* []>(diameter) };

    _srcp[radius] = __srcp;
    for (int i{ 1 }; i <= radius; ++i)
        _srcp[radius - i] = _srcp[radius + i] = _srcp[radius] + srcStride * i;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            Vec16f sum{ zero_16f() };

            for (int i{ 0 }; i < diameter; ++i)
            {
                if constexpr (std::is_same_v<T, uint8_t>)
                {
                    const Vec16f srcp{ to_float(Vec16i().load_16uc(_srcp[i] + x)) };
                    sum = mul_add(srcp, weights[i], sum);
                }
                else if constexpr (std::is_same_v<T, uint16_t>)
                {
                    const Vec16f srcp{ to_float(Vec16i().load_16us(_srcp[i] + x)) };
                    sum = mul_add(srcp, weights[i], sum);
                }
                else
                {
                    const Vec16f srcp{ Vec16f().load_a(_srcp[i] + x) };
                    sum = mul_add(srcp, weights[i], sum);
                }
            }

            sum.store_nt(dstp + x);
        }

        for (int i{ 0 }; i < diameter - 1; ++i)
            _srcp[i] = _srcp[i + 1];

        _srcp[diameter - 1] += (y < height - 1 - radius) ? srcStride : -srcStride;

        dstp += dstStride;
    }
}

template<typename T>
static void gaussianBlurH(const T* _srcp, float* temp, float* dstp, const float* weights, const int width, const int height, const int srcStride, const int dstStride, const int radius) noexcept
{
    weights += radius;

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            if constexpr (std::is_same_v<T, uint8_t>)
                to_float(Vec16i().load_16uc(_srcp + x)).store_a(temp + x);
            else if constexpr (std::is_same_v<T, uint16_t>)
                to_float(Vec16i().load_16us(_srcp + x)).store_a(temp + x);
            else
                Vec16f().load_a(_srcp + x).store_a(temp + x);
        }

        for (int i{ 1 }; i <= radius; ++i)
        {
            temp[-i] = temp[i];
            temp[width - 1 + i] = temp[width - 1 - i];
        }

        for (int x{ 0 }; x < width; x += 16)
        {
            Vec16f sum{ zero_16f() };

            for (int i{ -radius }; i <= radius; ++i)
            {
                const Vec16f srcp{ Vec16f().load(temp + x + i) };
                sum = mul_add(srcp, weights[i], sum);
            }

            sum.store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

static void detectEdge(float* blur, float* gradient, int* direction, const int width, const int height, const int stride, const int bgStride, const int mode, const int op, const float scale) noexcept
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

        for (int x{ 0 }; x < width; x += 16)
        {
            Vec16f gx, gy;

            if (op != FDOG)
            {
                const Vec16f c1{ Vec16f().load(prev + x - 1) };
                const Vec16f c2{ Vec16f().load_a(prev + x) };
                const Vec16f c3{ Vec16f().load(prev + x + 1) };
                const Vec16f c4{ Vec16f().load(cur + x - 1) };
                const Vec16f c6{ Vec16f().load(cur + x + 1) };
                const Vec16f c7{ Vec16f().load(next + x - 1) };
                const Vec16f c8{ Vec16f().load_a(next + x) };
                const Vec16f c9{ Vec16f().load(next + x + 1) };

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
                        gx = (c3 + c6 + c9 - c1 - c4 - c7) * 0.5f;
                        gy = (c1 + c2 + c3 - c7 - c8 - c9) * 0.5f;
                        break;
                    }
                    case SOBEL:
                    {
                        gx = c3 + mul_add(2.0f, c6, c9) - c1 - mul_add(2.0f, c4, c7);
                        gy = c1 + mul_add(2.0f, c2, c3) - c7 - mul_add(2.0f, c8, c9);
                        break;
                    }
                    case SCHARR:
                    {
                        gx = mul_add(3.0f, c3 + c9, 10.0f * c6) - mul_add(3.0f, c1 + c7, 10.0f * c4);
                        gy = mul_add(3.0f, c1 + c3, 10.0f * c2) - mul_add(3.0f, c7 + c9, 10.0f * c8);
                        break;
                    }
                    case KROON:
                    {
                        gx = mul_add(17.0f, c3 + c9, 61.0f * c6) - mul_add(17.0f, c1 + c7, 61.0f * c4);
                        gy = mul_add(17.0f, c1 + c3, 61.0f * c2) - mul_add(17.0f, c7 + c9, 61.0f * c8);
                        break;
                    }
                    case KIRSCH:
                    {
                        const Vec16f g1{ mul_sub(5.0f, c1 + c2 + c3, 3.0f * (c4 + c6 + c7 + c8 + c9)) };
                        const Vec16f g2{ mul_sub(5.0f, c1 + c2 + c4, 3.0f * (c3 + c6 + c7 + c8 + c9)) };
                        const Vec16f g3{ mul_sub(5.0f, c1 + c4 + c7, 3.0f * (c2 + c3 + c6 + c8 + c9)) };
                        const Vec16f g4{ mul_sub(5.0f, c4 + c7 + c8, 3.0f * (c1 + c2 + c3 + c6 + c9)) };
                        const Vec16f g5{ mul_sub(5.0f, c7 + c8 + c9, 3.0f * (c1 + c2 + c3 + c4 + c6)) };
                        const Vec16f g6{ mul_sub(5.0f, c6 + c8 + c9, 3.0f * (c1 + c2 + c3 + c4 + c7)) };
                        const Vec16f g7{ mul_sub(5.0f, c3 + c6 + c9, 3.0f * (c1 + c2 + c4 + c7 + c8)) };
                        const Vec16f g8{ mul_sub(5.0f, c2 + c3 + c6, 3.0f * (c1 + c4 + c7 + c8 + c9)) };
                        const Vec16f g{ max(max(max(abs(g1), abs(g2)), max(abs(g3), abs(g4))), max(max(abs(g5), abs(g6)), max(abs(g7), abs(g8)))) };
                        (g * scale).store_nt(gradient + x);
                        break;
                    }
                }
            }
            else
            {
                const Vec16f c1{ Vec16f().load(prev2 + x - 2) };
                const Vec16f c2{ Vec16f().load(prev2 + x - 1) };
                const Vec16f c3{ Vec16f().load(prev2 + x) };
                const Vec16f c4{ Vec16f().load(prev2 + x + 1) };
                const Vec16f c5{ Vec16f().load(prev2 + x + 2) };
                const Vec16f c6{ Vec16f().load(prev + x - 2) };
                const Vec16f c7{ Vec16f().load(prev + x - 1) };
                const Vec16f c8{ Vec16f().load(prev + x) };
                const Vec16f c9{ Vec16f().load(prev + x + 1) };
                const Vec16f c10{ Vec16f().load(prev + x + 2) };
                const Vec16f c11{ Vec16f().load(cur + x - 2) };
                const Vec16f c12{ Vec16f().load(cur + x - 1) };
                const Vec16f c14{ Vec16f().load(cur + x + 1) };
                const Vec16f c15{ Vec16f().load(cur + x + 2) };
                const Vec16f c16{ Vec16f().load(next + x - 2) };
                const Vec16f c17{ Vec16f().load(next + x - 1) };
                const Vec16f c18{ Vec16f().load(next + x) };
                const Vec16f c19{ Vec16f().load(next + x + 1) };
                const Vec16f c20{ Vec16f().load(next + x + 2) };
                const Vec16f c21{ Vec16f().load(next2 + x - 2) };
                const Vec16f c22{ Vec16f().load(next2 + x - 1) };
                const Vec16f c23{ Vec16f().load(next2 + x) };
                const Vec16f c24{ Vec16f().load(next2 + x + 1) };
                const Vec16f c25{ Vec16f().load(next2 + x + 2) };

                gx = c5 + c25 + c4 + c24 + mul_add(2.0f, c10 + c20 + c9 + c19, 3.0f * (c15 + c14))
                    - c2 - c22 - c1 - c21 - mul_add(2.0f, c7 + c17 + c6 + c16, 3.0f * (c12 + c11));
                gy = c1 + c5 + c6 + c10 + mul_add(2.0f, c2 + c4 + c7 + c9, 3.0f * (c3 + c8))
                    - c16 - c20 - c21 - c25 - mul_add(2.0f, c17 + c19 + c22 + c24, 3.0f * (c18 + c23));
            }

            if (op != KIRSCH)
            {
                gx *= scale;
                gy *= scale;
                sqrt(mul_add(gx, gx, gy * gy)).store_nt(gradient + x);
            }

            if (mode == 0)
            {
                Vec16f dr{ atan2(gy, gx) };
                dr = if_add(dr < 0.0f, dr, M_PIF);

                const Vec16i bin{ truncatei(mul_add(dr, 4.0f * M_1_PIF, 0.5f)) };
                select(bin >= 4, zero_si512(), bin).store_nt(direction + x);
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

static void nonMaximumSuppression(const int* _direction, float* _gradient, float* blur, const int width, const int height, const int stride, const int bgStride, const int radiusAlign) noexcept
{
    _gradient[-1] = _gradient[1];
    _gradient[-1 + bgStride * (height - 1)] = _gradient[1 + bgStride * (height - 1)];
    _gradient[width] = _gradient[width - 2];
    _gradient[width + bgStride * (height - 1)] = _gradient[width - 2 + bgStride * (height - 1)];
    std::copy_n(_gradient - radiusAlign + bgStride, width + radiusAlign * 2, _gradient - radiusAlign - bgStride);
    std::copy_n(_gradient - radiusAlign + bgStride * (height - 2), width + radiusAlign * 2, _gradient - radiusAlign + bgStride * static_cast<int64_t>(height));

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            const Vec16ui direction{ Vec16ui().load_a(_direction + x) };

            Vec16fb mask{ Vec16fb(direction == 0) };
            Vec16f gradient{ max(Vec16f().load(_gradient + x + 1), Vec16f().load(_gradient + x - 1)) };
            Vec16f result{ gradient & mask };

            mask = Vec16fb(direction == 1);
            gradient = max(Vec16f().load(_gradient + x - bgStride + 1), Vec16f().load(_gradient + x + bgStride - 1));
            result |= gradient & mask;

            mask = Vec16fb(direction == 2);
            gradient = max(Vec16f().load_a(_gradient + x - bgStride), Vec16f().load_a(_gradient + x + bgStride));
            result |= gradient & mask;

            mask = Vec16fb(direction == 3);
            gradient = max(Vec16f().load(_gradient + x - bgStride - 1), Vec16f().load(_gradient + x + bgStride + 1));
            result |= gradient & mask;

            gradient = Vec16f().load_a(_gradient + x);
            select(gradient >= result, gradient, fltLowest).store_nt(blur + x);
        }

        _direction += stride;
        _gradient += bgStride;
        blur += bgStride;
    }
}

template<typename T>
static void binarizeCE(const float* _srcp, T* dstp, const int width, const int height, const int srcStride, const int dstStride, const int peak) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            const Vec16f srcp{ Vec16f().load_a(_srcp + x) };

            if constexpr (std::is_same_v<T, uint8_t>)
            {
                const Vec16cb mask{ Vec16cb(srcp == fltMax) };
                select(mask, Vec16uc(255), zero_si128()).store_nt(dstp + x);
            }
            else if constexpr (std::is_same_v<T, uint16_t>)
            {
                const Vec16sb mask{ Vec16sb(srcp == fltMax) };
                select(mask, Vec16us(peak), zero_si256()).store_nt(dstp + x);
            }
            else
            {
                const Vec16fb mask{ srcp == fltMax };
                select(mask, Vec16f(1.0f), Vec16f(0.0f)).store_nt(dstp + x);
            }
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T, bool clampFP = true>
static void discretizeGM(const float* _srcp, T* dstp, const int width, const int height, const int srcStride, const int dstStride, const int peak) noexcept
{
    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; x += 16)
        {
            const Vec16f srcp{ Vec16f().load_a(_srcp + x) };

            if constexpr (std::is_same_v<T, uint8_t>)
            {
                const Vec16uc result{ compress_saturated_s2u(compress_saturated(truncatei(srcp + 0.5f), zero_si512()), zero_si512()).get_low().get_low() };
                result.store_nt(dstp + x);
            }
            else if constexpr (std::is_same_v<T, uint16_t>)
            {
                const Vec16us result{ compress_saturated_s2u(truncatei(srcp + 0.5f), zero_si512()).get_low() };
                min(result, peak).store_nt(dstp + x);
            }
            else if constexpr (clampFP)
                min(max(srcp, 0.0f), 1.0f).store_nt(dstp + x);
            else
                srcp.store_nt(dstp + x);
        }

        _srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void vsTCanny::filter_avx512(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept
{
    const int planes_y[3]{ PLANAR_Y, PLANAR_U, PLANAR_V };
    const int planes_r[3]{ PLANAR_G, PLANAR_B, PLANAR_R };
    const int* current_planes{ (vi.IsRGB()) ? planes_r : planes_y };
    const int planecount{ std::min(vi.NumComponents(), 3) };

    for (int i{ 0 }; i < planecount; ++i)
    {
        const int height{ src->GetHeight(current_planes[i]) };

        if (process[i] == 3)
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

            switch (mode_)
            {
                case 0: binarizeCE(blur, dstp, width, height, bgStride, dst_stride, peak); break;
                case 1: discretizeGM(gradient, dstp, width, height, bgStride, stride, peak); break;
                default: discretizeGM<T, false>(blur, dstp, width, height, bgStride, dst_stride, peak); break;
            }
        }
        else if (process[i] == 2)
            env->BitBlt(dst->GetWritePtr(current_planes[i]), dst->GetPitch(current_planes[i]), src->GetReadPtr(current_planes[i]), src->GetPitch(current_planes[i]), src->GetRowSize(current_planes[i]), height);
    }
}

template void vsTCanny::filter_avx512<uint8_t>(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
template void vsTCanny::filter_avx512<uint16_t>(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
template void vsTCanny::filter_avx512<float>(PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) noexcept;
