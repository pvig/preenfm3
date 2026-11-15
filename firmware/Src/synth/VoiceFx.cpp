/*
 * Copyright 2013 Xavier Hosxe
 *
 * Author: Xavier Hosxe (xavier . hosxe (at) gmail . com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "Voice.h"
#include "Timbre.h"

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC push_options
#pragma GCC optimize("Ofast", "fast-math")
#endif

// for bitwise manipulations
#define FLOAT2SHORT 32768.0f
#define SHORT2FLOAT (1.0f / 32768.0f)

#define RATIOINV (1.0f / 131072.0f)

#define SVFRANGE 1.23f
#define SVFOFFSET 0.151f
#define SVFGAINOFFSET 0.3f

#define LP2OFFSET -0.045f
#define min(a, b) ((a) < (b) ? (a) : (b))

extern float noise[32];

#define filterWindowMin 0.01f
#define filterWindowMax 0.99f

inline float expf_fast(float a)
{
    // https://github.com/ekmett/approximate/blob/master/cbits/fast.c
    union
    {
        float f;
        int x;
    } u;
    u.x = (int)(12102203 * a + 1064866805);
    return u.f;
}
inline float exp_pole(float p)
{
    // High-accuracy Padé (2,2) scaling & squaring for z-plane pole mapping exp(p):
    // Exact 1.0000000 at p = 0 (critical for high-Q TB-303 poles near the unit circle)
    const float h = p * 0.25f;
    const float h2 = h * h * 0.083333336f;
    const float hp = 0.5f * h;
    float e = (1.0f + hp + h2) / (1.0f - hp + h2);
    e *= e;
    return e * e;
}
inline float sqrt3(const float x)
{
    union
    {
        int i;
        float x;
    } u;

    u.x = x;
    u.i = (1 << 29) + (u.i >> 1) - (1 << 22);
    return u.x;
}
inline float fastInvSqrtVoice(float x)
{
    float xhalf = 0.5f * x;
    union
    {
        float f;
        uint32_t i;
    } u;
    u.f = x;
    u.i = 0x5f3759dfu - (u.i >> 1);
    u.f = u.f * (1.5f - xhalf * u.f * u.f);
    return u.f;
}
inline float softLimitVoice(float x)
{
    return x * fastInvSqrtVoice(1.0f + x * x);
}
inline float tanh3(float x)
{
    // return x / (0.7f + fabsf(0.8f * x));
    return 1.5f * x / (1.7f + fabsf(0.34f * x * x));
}
inline float tanh4(float x)
{
    return x / sqrt3(x * x + 1);
}
inline float tanh4ADAA1_F(float x)
{
    return sqrtf(x * x + 1.0f);
}
inline float tanh4ADAA(float x, float &prevX, float &prevF)
{
    const float f = tanh4ADAA1_F(x);
    const float diff = x - prevX;
    float out;
    if (likely(fabsf(diff) > 1e-5f))
    {
        out = (f - prevF) / diff;
    }
    else
    {
        out = softLimitVoice(0.5f * (x + prevX));
    }
    prevX = x;
    prevF = f;
    return out;
}
inline float sat25(float x)
{
    if (unlikely(fabsf(x) > 4))
        return 0;
    return x * (1.f - fabsf(x * 0.25f));
}
inline float sat25ADAA1_F(float x)
{
    const float ax = fabsf(x);
    if (unlikely(ax > 4.0f))
    {
        return 2.6666667f;
    }
    return x * x * (0.5f - 0.083333336f * ax);
}
inline float sat25ADAA(float x, float &prevX, float &prevF)
{
    const float f = sat25ADAA1_F(x);
    const float diff = x - prevX;
    float out;
    if (likely(fabsf(diff) > 1e-5f))
    {
        out = (f - prevF) / diff;
    }
    else
    {
        out = sat25(0.5f * (x + prevX));
    }
    prevX = x;
    prevF = f;
    return out;
}
inline float sat66(float x)
{
    return x * (1 - (x * x * 0.055f));
}
inline float sat33(float x)
{
    if (unlikely(fabsf(x) > 2.58f))
        return 0;
    return x * (1 - x * x * 0.15f);
}
inline float clamp(float d, float min, float max)
{
    const float t = unlikely(d < min) ? min : d;
    return unlikely(t > max) ? max : t;
}
inline float satSin(float x, float drive, float posF)
{
    // https://www.desmos.com/calculator/pdsqpi5lp6
    // Linearly interpolated lookup on sinTable to avoid table-quantization & modulation stepping aliasing
    const float phase = fabsf(x) * 360.0f + posF;
    const uint_fast16_t idx = (uint_fast16_t)phase;
    const float frac = phase - (float)idx;
    const float s0 = sinTable[idx & 0x3FF];
    const float s1 = sinTable[(idx + 1) & 0x3FF];
    return x + clamp(x, -0.5f * drive, drive) * (s0 + frac * (s1 - s0));
}
// https://www.musicdsp.org/en/latest/Other/120-saturation.html
inline float sigmoid(float x)
{
    const float xc = clamp(x, -1.0f, 1.0f);
    return xc * (1.5f - 0.5f * xc * xc);
}
inline float sigmoidADAA1_F(float x)
{
    const float xc = clamp(x, -1.0f, 1.0f);
    const float x2 = xc * xc;
    return x2 * (0.75f - 0.125f * x2) + (fabsf(x) - fabsf(xc));
}
inline float sigmoidADAA(float x, float &prevX, float &prevF)
{
    const float f = sigmoidADAA1_F(x);
    const float diff = x - prevX;
    float out;
    if (likely(fabsf(diff) > 1e-4f))
    {
        out = (f - prevF) / diff;
    }
    else
    {
        out = sigmoid(0.5f * (x + prevX));
    }
    prevX = x;
    prevF = f;
    return out;
}
inline float sigmoidPos(float x)
{
    // x : 0 -> 1
    return (sigmoid((x * 2) - 1) + 1) * 0.5f;
}
inline float fold(float x4)
{
    // https://www.desmos.com/calculator/ge2wvg2wgj
    // x4 = x / 4
    return (fabsf(x4 + 0.25f - roundf(x4 + 0.25f)) - 0.25f);
}
inline float foldADAA1_F(float x4)
{
    const float u = x4 + 0.25f;
    const float d = u - roundf(u);
    return d * (0.5f * fabsf(d) - 0.25f);
}
inline float foldADAA(float x4, float &prevX4, float &prevF)
{
    const float f = foldADAA1_F(x4);
    const float diff = x4 - prevX4;
    float out;
    if (likely(fabsf(diff) > 1e-5f))
    {
        out = (f - prevF) / diff;
    }
    else
    {
        out = fold(0.5f * (x4 + prevX4));
    }
    prevX4 = x4;
    prevF = f;
    return out;
}
inline float wrap(float x)
{
    // Centered symmetric phase wrapper: wrap(0) = 0, wraps at +/-1 into [-1, 1)
    const float u = 0.5f * x + 0.5f;
    return (u - floorf(u)) * 2.0f - 1.0f;
}
inline float wrapADAA1_F(float x)
{
    const float u = 0.5f * x + 0.5f;
    const float d = u - floorf(u);
    return 2.0f * d * (d - 1.0f);
}
inline float wrapADAA(float x, float &prevX, float &prevF)
{
    const float f = wrapADAA1_F(x);
    const float diff = x - prevX;
    float out;
    if (likely(fabsf(diff) > 1e-5f))
    {
        out = (f - prevF) / diff;
    }
    else
    {
        out = wrap(0.5f * (x + prevX));
    }
    prevX = x;
    prevF = f;
    return out;
}

// TB303 filter emulation ------
const float filterpoles[10][64] = {
    {-1.42475857e-02, -1.10558351e-02, -9.58097367e-03,
     -8.63568249e-03, -7.94942757e-03, -7.41570560e-03,
     -6.98187179e-03, -6.61819537e-03, -6.30631927e-03,
     -6.03415378e-03, -5.79333654e-03, -5.57785533e-03,
     -5.38325013e-03, -5.20612558e-03, -5.04383985e-03,
     -4.89429884e-03, -4.75581571e-03, -4.62701254e-03,
     -4.50674977e-03, -4.39407460e-03, -4.28818259e-03,
     -4.18838855e-03, -4.09410427e-03, -4.00482112e-03,
     -3.92009643e-03, -3.83954259e-03, -3.76281836e-03,
     -3.68962181e-03, -3.61968451e-03, -3.55276681e-03,
     -3.48865386e-03, -3.42715236e-03, -3.36808777e-03,
     -3.31130196e-03, -3.25665127e-03, -3.20400476e-03,
     -3.15324279e-03, -3.10425577e-03, -3.05694308e-03,
     -3.01121207e-03, -2.96697733e-03, -2.92415989e-03,
     -2.88268665e-03, -2.84248977e-03, -2.80350622e-03,
     -2.76567732e-03, -2.72894836e-03, -2.69326825e-03,
     -2.65858922e-03, -2.62486654e-03, -2.59205824e-03,
     -2.56012496e-03, -2.52902967e-03, -2.49873752e-03,
     -2.46921570e-03, -2.44043324e-03, -2.41236091e-03,
     -2.38497108e-03, -2.35823762e-03, -2.33213577e-03,
     -2.30664208e-03, -2.28173430e-03, -2.25739130e-03,
     -2.23359302e-03},
    {1.63323670e-16, -1.61447133e-02, -1.99932070e-02,
     -2.09872000e-02, -2.09377795e-02, -2.04470150e-02,
     -1.97637613e-02, -1.90036975e-02, -1.82242987e-02,
     -1.74550383e-02, -1.67110053e-02, -1.59995606e-02,
     -1.53237941e-02, -1.46844019e-02, -1.40807436e-02,
     -1.35114504e-02, -1.29747831e-02, -1.24688429e-02,
     -1.19916965e-02, -1.15414484e-02, -1.11162818e-02,
     -1.07144801e-02, -1.03344362e-02, -9.97465446e-03,
     -9.63374867e-03, -9.31043725e-03, -9.00353710e-03,
     -8.71195702e-03, -8.43469084e-03, -8.17081077e-03,
     -7.91946102e-03, -7.67985179e-03, -7.45125367e-03,
     -7.23299254e-03, -7.02444481e-03, -6.82503313e-03,
     -6.63422244e-03, -6.45151640e-03, -6.27645413e-03,
     -6.10860728e-03, -5.94757730e-03, -5.79299303e-03,
     -5.64450848e-03, -5.50180082e-03, -5.36456851e-03,
     -5.23252970e-03, -5.10542063e-03, -4.98299431e-03,
     -4.86501921e-03, -4.75127814e-03, -4.64156716e-03,
     -4.53569463e-03, -4.43348032e-03, -4.33475462e-03,
     -4.23935774e-03, -4.14713908e-03, -4.05795659e-03,
     -3.97167614e-03, -3.88817107e-03, -3.80732162e-03,
     -3.72901453e-03, -3.65314257e-03, -3.57960420e-03,
     -3.50830319e-03},
    {-1.83545593e-06, -1.35008051e-03, -1.51527847e-03,
     -1.61437715e-03, -1.68536679e-03, -1.74064961e-03,
     -1.78587681e-03, -1.82410854e-03, -1.85719118e-03,
     -1.88632533e-03, -1.91233586e-03, -1.93581405e-03,
     -1.95719818e-03, -1.97682215e-03, -1.99494618e-03,
     -2.01177700e-03, -2.02748155e-03, -2.04219657e-03,
     -2.05603546e-03, -2.06909331e-03, -2.08145062e-03,
     -2.09317612e-03, -2.10432901e-03, -2.11496056e-03,
     -2.12511553e-03, -2.13483321e-03, -2.14414822e-03,
     -2.15309131e-03, -2.16168985e-03, -2.16996830e-03,
     -2.17794867e-03, -2.18565078e-03, -2.19309254e-03,
     -2.20029023e-03, -2.20725864e-03, -2.21401130e-03,
     -2.22056055e-03, -2.22691775e-03, -2.23309332e-03,
     -2.23909688e-03, -2.24493730e-03, -2.25062280e-03,
     -2.25616099e-03, -2.26155896e-03, -2.26682328e-03,
     -2.27196010e-03, -2.27697514e-03, -2.28187376e-03,
     -2.28666097e-03, -2.29134148e-03, -2.29591970e-03,
     -2.30039977e-03, -2.30478562e-03, -2.30908091e-03,
     -2.31328911e-03, -2.31741351e-03, -2.32145721e-03,
     -2.32542313e-03, -2.32931406e-03, -2.33313263e-03,
     -2.33688133e-03, -2.34056255e-03, -2.34417854e-03,
     -2.34773145e-03},
    {-2.96292613e-06, 6.75138822e-04, 6.96581050e-04,
     7.04457808e-04, 7.07837502e-04, 7.09169651e-04,
     7.09415480e-04, 7.09031433e-04, 7.08261454e-04,
     7.07246872e-04, 7.06074484e-04, 7.04799978e-04,
     7.03460301e-04, 7.02080606e-04, 7.00678368e-04,
     6.99265907e-04, 6.97852005e-04, 6.96442963e-04,
     6.95043317e-04, 6.93656323e-04, 6.92284301e-04,
     6.90928882e-04, 6.89591181e-04, 6.88271928e-04,
     6.86971561e-04, 6.85690300e-04, 6.84428197e-04,
     6.83185182e-04, 6.81961088e-04, 6.80755680e-04,
     6.79568668e-04, 6.78399727e-04, 6.77248505e-04,
     6.76114631e-04, 6.74997722e-04, 6.73897392e-04,
     6.72813249e-04, 6.71744904e-04, 6.70691972e-04,
     6.69654071e-04, 6.68630828e-04, 6.67621875e-04,
     6.66626854e-04, 6.65645417e-04, 6.64677222e-04,
     6.63721940e-04, 6.62779248e-04, 6.61848835e-04,
     6.60930398e-04, 6.60023644e-04, 6.59128290e-04,
     6.58244058e-04, 6.57370684e-04, 6.56507909e-04,
     6.55655483e-04, 6.54813164e-04, 6.53980718e-04,
     6.53157918e-04, 6.52344545e-04, 6.51540387e-04,
     6.50745236e-04, 6.49958895e-04, 6.49181169e-04,
     6.48411873e-04},
    {-1.00014774e+00, -1.35336624e+00, -1.42048887e+00,
     -1.46551548e+00, -1.50035433e+00, -1.52916086e+00,
     -1.55392254e+00, -1.57575858e+00, -1.59536715e+00,
     -1.61321568e+00, -1.62963377e+00, -1.64486333e+00,
     -1.65908760e+00, -1.67244897e+00, -1.68506052e+00,
     -1.69701363e+00, -1.70838333e+00, -1.71923202e+00,
     -1.72961221e+00, -1.73956855e+00, -1.74913935e+00,
     -1.75835773e+00, -1.76725258e+00, -1.77584919e+00,
     -1.78416990e+00, -1.79223453e+00, -1.80006075e+00,
     -1.80766437e+00, -1.81505964e+00, -1.82225940e+00,
     -1.82927530e+00, -1.83611794e+00, -1.84279698e+00,
     -1.84932127e+00, -1.85569892e+00, -1.86193740e+00,
     -1.86804360e+00, -1.87402388e+00, -1.87988413e+00,
     -1.88562983e+00, -1.89126607e+00, -1.89679760e+00,
     -1.90222885e+00, -1.90756395e+00, -1.91280679e+00,
     -1.91796101e+00, -1.92303002e+00, -1.92801704e+00,
     -1.93292509e+00, -1.93775705e+00, -1.94251559e+00,
     -1.94720328e+00, -1.95182252e+00, -1.95637561e+00,
     -1.96086471e+00, -1.96529188e+00, -1.96965908e+00,
     -1.97396817e+00, -1.97822093e+00, -1.98241904e+00,
     -1.98656411e+00, -1.99065768e+00, -1.99470122e+00,
     -1.99869613e+00},
    {1.30592376e-04, 3.54780202e-01, 4.22050344e-01,
     4.67149412e-01, 5.02032084e-01, 5.30867858e-01,
     5.55650170e-01, 5.77501296e-01, 5.97121154e-01,
     6.14978238e-01, 6.31402872e-01, 6.46637440e-01,
     6.60865515e-01, 6.74229755e-01, 6.86843408e-01,
     6.98798009e-01, 7.10168688e-01, 7.21017938e-01,
     7.31398341e-01, 7.41354603e-01, 7.50925074e-01,
     7.60142923e-01, 7.69037045e-01, 7.77632782e-01,
     7.85952492e-01, 7.94016007e-01, 8.01841009e-01,
     8.09443333e-01, 8.16837226e-01, 8.24035549e-01,
     8.31049962e-01, 8.37891065e-01, 8.44568531e-01,
     8.51091211e-01, 8.57467223e-01, 8.63704040e-01,
     8.69808551e-01, 8.75787123e-01, 8.81645657e-01,
     8.87389629e-01, 8.93024133e-01, 8.98553916e-01,
     9.03983409e-01, 9.09316756e-01, 9.14557836e-01,
     9.19710291e-01, 9.24777540e-01, 9.29762800e-01,
     9.34669099e-01, 9.39499296e-01, 9.44256090e-01,
     9.48942030e-01, 9.53559531e-01, 9.58110882e-01,
     9.62598250e-01, 9.67023698e-01, 9.71389181e-01,
     9.75696562e-01, 9.79947614e-01, 9.84144025e-01,
     9.88287408e-01, 9.92379299e-01, 9.96421168e-01,
     1.00041442e+00},
    {-2.96209812e-06, -2.45794824e-04, -8.18027564e-04,
     -1.19157447e-03, -1.46371229e-03, -1.67529045e-03,
     -1.84698016e-03, -1.99058664e-03, -2.11344205e-03,
     -2.22039065e-03, -2.31478873e-03, -2.39905115e-03,
     -2.47496962e-03, -2.54390793e-03, -2.60692676e-03,
     -2.66486645e-03, -2.71840346e-03, -2.76809003e-03,
     -2.81438252e-03, -2.85766225e-03, -2.89825096e-03,
     -2.93642247e-03, -2.97241172e-03, -3.00642174e-03,
     -3.03862912e-03, -3.06918837e-03, -3.09823546e-03,
     -3.12589065e-03, -3.15226077e-03, -3.17744116e-03,
     -3.20151726e-03, -3.22456591e-03, -3.24665644e-03,
     -3.26785166e-03, -3.28820859e-03, -3.30777919e-03,
     -3.32661092e-03, -3.34474723e-03, -3.36222800e-03,
     -3.37908995e-03, -3.39536690e-03, -3.41109012e-03,
     -3.42628855e-03, -3.44098902e-03, -3.45521647e-03,
     -3.46899410e-03, -3.48234354e-03, -3.49528498e-03,
     -3.50783728e-03, -3.52001812e-03, -3.53184405e-03,
     -3.54333061e-03, -3.55449241e-03, -3.56534320e-03,
     -3.57589590e-03, -3.58616273e-03, -3.59615520e-03,
     -3.60588419e-03, -3.61536000e-03, -3.62459235e-03,
     -3.63359049e-03, -3.64236316e-03, -3.65091867e-03,
     -3.65926491e-03},
    {-7.75894750e-06, 3.11294169e-03, 3.41779455e-03,
     3.52160375e-03, 3.55957019e-03, 3.56903631e-03,
     3.56431495e-03, 3.55194570e-03, 3.53526954e-03,
     3.51613008e-03, 3.49560287e-03, 3.47434152e-03,
     3.45275527e-03, 3.43110577e-03, 3.40956242e-03,
     3.38823540e-03, 3.36719598e-03, 3.34648945e-03,
     3.32614343e-03, 3.30617351e-03, 3.28658692e-03,
     3.26738515e-03, 3.24856568e-03, 3.23012330e-03,
     3.21205091e-03, 3.19434023e-03, 3.17698219e-03,
     3.15996727e-03, 3.14328577e-03, 3.12692791e-03,
     3.11088400e-03, 3.09514449e-03, 3.07970007e-03,
     3.06454165e-03, 3.04966043e-03, 3.03504790e-03,
     3.02069585e-03, 3.00659636e-03, 2.99274180e-03,
     2.97912486e-03, 2.96573849e-03, 2.95257590e-03,
     2.93963061e-03, 2.92689635e-03, 2.91436713e-03,
     2.90203718e-03, 2.88990095e-03, 2.87795312e-03,
     2.86618855e-03, 2.85460234e-03, 2.84318974e-03,
     2.83194618e-03, 2.82086729e-03, 2.80994883e-03,
     2.79918673e-03, 2.78857707e-03, 2.77811607e-03,
     2.76780009e-03, 2.75762559e-03, 2.74758919e-03,
     2.73768761e-03, 2.72791768e-03, 2.71827634e-03,
     2.70876064e-03},
    {-9.99869423e-01, -6.38561407e-01, -5.69514530e-01,
     -5.23990915e-01, -4.89176780e-01, -4.60615628e-01,
     -4.36195579e-01, -4.14739573e-01, -3.95520699e-01,
     -3.78056805e-01, -3.62010728e-01, -3.47136887e-01,
     -3.33250504e-01, -3.20208824e-01, -3.07899106e-01,
     -2.96230641e-01, -2.85129278e-01, -2.74533563e-01,
     -2.64391946e-01, -2.54660728e-01, -2.45302512e-01,
     -2.36285026e-01, -2.27580207e-01, -2.19163487e-01,
     -2.11013226e-01, -2.03110249e-01, -1.95437482e-01,
     -1.87979648e-01, -1.80723016e-01, -1.73655197e-01,
     -1.66764971e-01, -1.60042136e-01, -1.53477393e-01,
     -1.47062234e-01, -1.40788856e-01, -1.34650080e-01,
     -1.28639289e-01, -1.22750366e-01, -1.16977645e-01,
     -1.11315866e-01, -1.05760138e-01, -1.00305900e-01,
     -9.49488960e-02, -8.96851464e-02, -8.45109223e-02,
     -7.94227260e-02, -7.44172709e-02, -6.94914651e-02,
     -6.46423954e-02, -5.98673139e-02, -5.51636250e-02,
     -5.05288741e-02, -4.59607376e-02, -4.14570134e-02,
     -3.70156122e-02, -3.26345497e-02, -2.83119399e-02,
     -2.40459880e-02, -1.98349851e-02, -1.56773019e-02,
     -1.15713843e-02, -7.51574873e-03, -3.50897732e-03,
     4.50285508e-04},
    {1.13389002e-04, 3.50509549e-01, 4.19971782e-01,
     4.66835760e-01, 5.03053790e-01, 5.32907131e-01,
     5.58475931e-01, 5.80942937e-01, 6.01050219e-01,
     6.19296203e-01, 6.36032925e-01, 6.51518847e-01,
     6.65949666e-01, 6.79477330e-01, 6.92222311e-01,
     7.04281836e-01, 7.15735567e-01, 7.26649641e-01,
     7.37079603e-01, 7.47072578e-01, 7.56668915e-01,
     7.65903438e-01, 7.74806427e-01, 7.83404383e-01,
     7.91720644e-01, 7.99775871e-01, 8.07588450e-01,
     8.15174821e-01, 8.22549745e-01, 8.29726527e-01,
     8.36717208e-01, 8.43532720e-01, 8.50183021e-01,
     8.56677208e-01, 8.63023619e-01, 8.69229911e-01,
     8.75303138e-01, 8.81249811e-01, 8.87075954e-01,
     8.92787154e-01, 8.98388600e-01, 9.03885123e-01,
     9.09281227e-01, 9.14581119e-01, 9.19788738e-01,
     9.24907772e-01, 9.29941684e-01, 9.34893728e-01,
     9.39766966e-01, 9.44564285e-01, 9.49288407e-01,
     9.53941905e-01, 9.58527211e-01, 9.63046630e-01,
     9.67502344e-01, 9.71896424e-01, 9.76230838e-01,
     9.80507456e-01, 9.84728057e-01, 9.88894335e-01,
     9.93007906e-01, 9.97070310e-01, 1.00108302e+00,
     1.00504744e+00}};

void Voice::fxAfterBlock()
{
    float ratioTimbres = 1.0f;
    float matrixFilterFrequency = matrix.getDestination(FILTER1_PARAM1);
    float matrixFilterParam2 = matrix.getDestination(FILTER1_PARAM2);
    float matrixFilterAmp = matrix.getDestination(FILTER1_AMP);

    // LP Algo
    int effectType = currentTimbre->params_.effect1.type;
    if (unlikely(fxParamB1 != (float)effectType))
    {
        v0L = v1L = v2L = v3L = v4L = v5L = v6L = v7L = v8L = v0R = v1R = v2R = v3R = v4R = v5R = v6R = v7R = v8R = 0.0f;
        fxParamA1 = fxParamA2 = fxParamB2 = 0.0f;
        fxPhase = 0.0f;
        fxParam1 = clamp(currentTimbre->params_.effect1.param1, 0.0f, 1.0f);
        fxParam2 = clamp(currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParamB1 = (float)effectType;
    }
    float gainTmp = clamp(currentTimbre->params_.effect1.param3 + matrixFilterAmp, 0, 16);
    mixerGain = 0.02f * gainTmp + .98f * mixerGain;

    switch (effectType)
    {
    case FILTER_LP:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 4.0f * prevParam1) * 0.2f, 0.0f, 1.0f);

        float res = clamp(currentTimbre->params_.effect1.param2 + matrixFilterParam2, 0.0f, 1.0f);
        const float g0 = 0.0025f + prevParam1 * (0.42f + 0.73f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.42f + 0.73f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 2.0f - 1.97f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float satK = 0.04f + 0.14f * res;
        const float outGain = (1.0f + 0.22f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = *sp;
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            *sp++ = clamp(v2L_ * outGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            *sp++ = clamp(v2R_ * outGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_HP:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(currentTimbre->params_.effect1.param2 + matrixFilterParam2, 0.0f, 1.0f);
        const float g0 = 0.002f + prevParam1 * (0.38f + 0.72f * prevParam1);
        const float g1 = 0.002f + fxParam1 * (0.38f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 2.0f - 1.96f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float satK = 0.04f + 0.15f * res;
        const float outGain = (1.0f + 0.18f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = *sp;
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            float hpL = inL - k * v1L_ - v2L_;
            *sp++ = clamp(hpL * outGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            float hpR = inR - k * v1R_ - v2R_;
            *sp++ = clamp(hpR * outGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_BASS:
    {
        float freqParam = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = (freqParam + 9.0f * fxParam1) * 0.1f;
        float boostParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // Tunable resonant sub-bass contour (approx 45 Hz to 320 Hz) with warm saturation
        const float g = 0.003f + fxParam1 * (0.012f + 0.028f * fxParam1);
        const float k = 1.35f - 0.55f * boostParam;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;
        const float boostGain = boostParam * (1.8f + 1.4f * boostParam);
        const float dryComp = 1.0f - 0.18f * boostParam;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            float inL = *sp;
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * ic1eqL + a2 * v3L_;
            float v2L_ = ic2eqL + a2 * ic1eqL + a3 * v3L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            float lowSatL = softLimitVoice(v2L_ * 1.4f);
            *sp++ = clamp((inL * dryComp + lowSatL * boostGain) * mixerGain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * ic1eqR + a2 * v3R_;
            float v2R_ = ic2eqR + a2 * ic1eqR + a3 * v3R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            float lowSatR = softLimitVoice(v2R_ * 1.4f);
            *sp++ = clamp((inR * dryComp + lowSatR * boostGain) * mixerGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_MIXER:
    {
        float panParam = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        float pan = panParam * 2.0f - 1.0f;
        float *sp = this->sampleBlock;
        float sampleR, sampleL;
        if (pan <= 0)
        {
            float onePlusPan = 1 + pan;
            float minusPan = -pan;
            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                sampleL = *(sp);
                sampleR = *(sp + 1);

                *sp = (sampleL + sampleR * minusPan) * mixerGain;
                sp++;
                *sp = sampleR * onePlusPan * mixerGain;
                sp++;
            }
        }
        else if (pan > 0)
        {
            float oneMinusPan = 1 - pan;
            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                sampleL = *(sp);
                sampleR = *(sp + 1);

                *sp = sampleL * oneMinusPan * mixerGain;
                sp++;
                *sp = (sampleR + sampleL * pan) * mixerGain;
                sp++;
            }
        }
    }
    break;
    case FILTER_CRUSHER:
    {
        // Lo-Fi Crusher: symmetric mid-tread quantizer & sample-and-hold decimator
        // with adaptive reconstruction smoothing to tame harsh zero-crossing clicks
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency + 0.005f, 0.005f, 1.0f);
        fxParamA1 = clamp((fxParamTmp + 9.0f * fxParamA1) * 0.1f, 0.0f, 1.0f);
        fxParam1 = fxParamA1;
        const float fxFreq = fxParamA1 * fxParamA1 * (0.35f + 0.65f * fxParamA1);

        float *sp = this->sampleBlock;
        float localPhase = fxPhase;

        float param2Sum = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = param2Sum;
        // Continuous bit-depth scaling from 2 levels (1-bit) to 32768 levels (15-bit)
        int bits = (int)(1.0f + 14.0f * param2Sum);
        float subBit = 1.0f + ((1.0f + 14.0f * param2Sum) - (float)bits);
        const float localPower = (float)(1 << bits) * subBit;
        const float localStep = 1.0f / localPower;
        // Adaptive vintage DAC reconstruction filter
        const float smooth = clamp(0.28f + 0.72f * fxParamA1, 0.28f, 1.0f);

        float holdL = v0L, holdR = v0R;
        float outL = v1L, outR = v1R;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            localPhase += fxFreq;
            if (unlikely(localPhase >= 1.0f))
            {
                localPhase -= 1.0f;
                holdL = roundf((*sp) * localPower) * localStep;
                holdR = roundf((*(sp + 1)) * localPower) * localStep;
            }
            outL += smooth * (holdL - outL);
            outR += smooth * (holdR - outR);

            *sp++ = clamp(outL * mixerGain, -ratioTimbres, ratioTimbres);
            *sp++ = clamp(outR * mixerGain, -ratioTimbres, ratioTimbres);
        }
        v0L = holdL;
        v0R = holdR;
        v1L = outL;
        v1R = outR;
        fxPhase = localPhase;
    }
    break;
    case FILTER_BP:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float qParam = clamp(currentTimbre->params_.effect1.param2 + matrixFilterParam2, 0.0f, 1.0f);
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 2.0f - 1.97f * qParam;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float satK = 0.03f + 0.12f * qParam;
        // Normalized bandpass gain (k + min floor) so high Q stays full and musical
        const float bpNorm = (0.38f + 0.78f * k) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            float inL = *sp;
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            *sp++ = clamp(v1L_ * bpNorm, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            *sp++ = clamp(v1R_ * bpNorm, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        break;
    }

    case FILTER_LP2:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(LP2OFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float g0 = 0.0025f + prevParam1 * (0.42f + 0.73f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.42f + 0.73f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.95f - 1.93f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float drive = 1.15f + 0.85f * res;
        const float satK = 0.10f + 0.32f * res;
        const float outGain = (0.95f + 0.20f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice: MS-20 / Polivoks style asymmetric driven ZDF SVF
            float inL = softLimitVoice((*sp) * drive + 0.04f * res) - 0.038f * res;
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            *sp++ = clamp(v2L_ * outGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = softLimitVoice((*sp) * drive + 0.04f * res) - 0.038f * res;
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            *sp++ = clamp(v2R_ * outGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_HP2:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(LP2OFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float g0 = 0.002f + prevParam1 * (0.38f + 0.72f * prevParam1);
        const float g1 = 0.002f + fxParam1 * (0.38f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.95f - 1.93f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float drive = 1.15f + 0.75f * res;
        const float satK = 0.10f + 0.30f * res;
        const float outGain = (0.95f + 0.18f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = softLimitVoice((*sp) * drive);
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            float hpL = inL - k * v1L_ - v2L_;
            *sp++ = clamp(hpL * outGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = softLimitVoice((*sp) * drive);
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            float hpR = inR - k * v1R_ - v2R_;
            *sp++ = clamp(hpR * outGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_BP2:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float qParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.95f - 1.93f * qParam;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float drive = 1.25f + 0.90f * qParam;
        const float satK = 0.10f + 0.28f * qParam;
        const float bpNorm = (0.42f + 0.72f * k) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left
            float inL = softLimitVoice((*sp) * drive);
            float v3L_ = inL - ic2eqL;
            float v1L_ = a1 * (ic1eqL + g * v3L_);
            float v2L_ = ic2eqL + g * v1L_;
            ic1eqL = 2.0f * v1L_ - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2L_ - ic2eqL;
            *sp++ = clamp(v1L_ * bpNorm, -ratioTimbres, ratioTimbres);

            // Right
            float inR = softLimitVoice((*sp) * drive);
            float v3R_ = inR - ic2eqR;
            float v1R_ = a1 * (ic1eqR + g * v3R_);
            float v2R_ = ic2eqR + g * v1R_;
            ic1eqR = 2.0f * v1R_ - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2R_ - ic2eqR;
            *sp++ = clamp(v1R_ * bpNorm, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        break;
    }
    case FILTER_LP3:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // True 24 dB/oct 2-stage cascaded ZDF SVF with staggered resonance
        const float g0 = 0.0025f + prevParam1 * (0.42f + 0.70f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.42f + 0.70f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k1 = 1.45f - 0.85f * res;
        const float k2 = 1.85f - 1.81f * res;
        const float a1_1_0 = 1.0f / (1.0f + g0 * (g0 + k1));
        const float a1_1_1 = 1.0f / (1.0f + g1 * (g1 + k1));
        const float da1_1 = (a1_1_1 - a1_1_0) * 0.03125f;
        float a1_1_lin = a1_1_0;
        const float a1_2_0 = 1.0f / (1.0f + g0 * (g0 + k2));
        const float a1_2_1 = 1.0f / (1.0f + g1 * (g1 + k2));
        const float da1_2 = (a1_2_1 - a1_2_0) * 0.03125f;
        float a1_2_lin = a1_2_0;
        const float satK = 0.08f + 0.24f * res;
        const float svfGain = (1.0f + 0.28f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float s1_ic1L = v0L, s1_ic2L = v1L, s2_ic1L = v2L, s2_ic2L = v3L;
        float s1_ic1R = v0R, s1_ic2R = v1R, s2_ic1R = v2R, s2_ic2R = v3R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_1_lin += da1_1;
            a1_2_lin += da1_2;
            const float a1_1 = a1_1_lin * (2.0f - (1.0f + g * (g + k1)) * a1_1_lin);
            const float a1_2 = a1_2_lin * (2.0f - (1.0f + g * (g + k2)) * a1_2_lin);

            // Left voice - Stage 1
            float inL = *sp;
            float v3 = inL - s1_ic2L;
            float v1 = a1_1 * (s1_ic1L + g * v3);
            float v2 = s1_ic2L + g * v1;
            s1_ic1L = 2.0f * v1 - s1_ic1L;
            s1_ic2L = 2.0f * v2 - s1_ic2L;
            // Left voice - Stage 2
            v3 = v2 - s2_ic2L;
            v1 = a1_2 * (s2_ic1L + g * v3);
            v2 = s2_ic2L + g * v1;
            s2_ic1L = 2.0f * v1 - s2_ic1L;
            s2_ic1L -= satK * s2_ic1L * (s2_ic1L * s2_ic1L) * fastInvSqrtVoice(1.0f + s2_ic1L * s2_ic1L);
            s2_ic2L = 2.0f * v2 - s2_ic2L;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice - Stage 1
            float inR = *sp;
            v3 = inR - s1_ic2R;
            v1 = a1_1 * (s1_ic1R + g * v3);
            v2 = s1_ic2R + g * v1;
            s1_ic1R = 2.0f * v1 - s1_ic1R;
            s1_ic2R = 2.0f * v2 - s1_ic2R;
            // Right voice - Stage 2
            v3 = v2 - s2_ic2R;
            v1 = a1_2 * (s2_ic1R + g * v3);
            v2 = s2_ic2R + g * v1;
            s2_ic1R = 2.0f * v1 - s2_ic1R;
            s2_ic1R -= satK * s2_ic1R * (s2_ic1R * s2_ic1R) * fastInvSqrtVoice(1.0f + s2_ic1R * s2_ic1R);
            s2_ic2R = 2.0f * v2 - s2_ic2R;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = s1_ic1L;
        v1L = s1_ic2L;
        v2L = s2_ic1L;
        v3L = s2_ic2L;
        v0R = s1_ic1R;
        v1R = s1_ic2R;
        v2R = s2_ic1R;
        v3R = s2_ic2R;
    }
    break;
    case FILTER_HP3:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // True 24 dB/oct 2-stage cascaded ZDF Highpass
        const float g0 = 0.002f + prevParam1 * (0.38f + 0.70f * prevParam1);
        const float g1 = 0.002f + fxParam1 * (0.38f + 0.70f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k1 = 1.45f - 0.85f * res;
        const float k2 = 1.85f - 1.81f * res;
        const float a1_1_0 = 1.0f / (1.0f + g0 * (g0 + k1));
        const float a1_1_1 = 1.0f / (1.0f + g1 * (g1 + k1));
        const float da1_1 = (a1_1_1 - a1_1_0) * 0.03125f;
        float a1_1_lin = a1_1_0;
        const float a1_2_0 = 1.0f / (1.0f + g0 * (g0 + k2));
        const float a1_2_1 = 1.0f / (1.0f + g1 * (g1 + k2));
        const float da1_2 = (a1_2_1 - a1_2_0) * 0.03125f;
        float a1_2_lin = a1_2_0;
        const float satK = 0.08f + 0.24f * res;
        const float svfGain = (1.0f + 0.22f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float s1_ic1L = v0L, s1_ic2L = v1L, s2_ic1L = v2L, s2_ic2L = v3L;
        float s1_ic1R = v0R, s1_ic2R = v1R, s2_ic1R = v2R, s2_ic2R = v3R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_1_lin += da1_1;
            a1_2_lin += da1_2;
            const float a1_1 = a1_1_lin * (2.0f - (1.0f + g * (g + k1)) * a1_1_lin);
            const float a1_2 = a1_2_lin * (2.0f - (1.0f + g * (g + k2)) * a1_2_lin);

            // Left voice - Stage 1
            float inL = *sp;
            float v3 = inL - s1_ic2L;
            float v1 = a1_1 * (s1_ic1L + g * v3);
            float v2 = s1_ic2L + g * v1;
            s1_ic1L = 2.0f * v1 - s1_ic1L;
            s1_ic2L = 2.0f * v2 - s1_ic2L;
            float hp1L = inL - k1 * v1 - v2;
            // Left voice - Stage 2
            v3 = hp1L - s2_ic2L;
            v1 = a1_2 * (s2_ic1L + g * v3);
            v2 = s2_ic2L + g * v1;
            s2_ic1L = 2.0f * v1 - s2_ic1L;
            s2_ic1L -= satK * s2_ic1L * (s2_ic1L * s2_ic1L) * fastInvSqrtVoice(1.0f + s2_ic1L * s2_ic1L);
            s2_ic2L = 2.0f * v2 - s2_ic2L;
            float hp2L = hp1L - k2 * v1 - v2;
            *sp++ = clamp(hp2L * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice - Stage 1
            float inR = *sp;
            v3 = inR - s1_ic2R;
            v1 = a1_1 * (s1_ic1R + g * v3);
            v2 = s1_ic2R + g * v1;
            s1_ic1R = 2.0f * v1 - s1_ic1R;
            s1_ic2R = 2.0f * v2 - s1_ic2R;
            float hp1R = inR - k1 * v1 - v2;
            // Right voice - Stage 2
            v3 = hp1R - s2_ic2R;
            v1 = a1_2 * (s2_ic1R + g * v3);
            v2 = s2_ic2R + g * v1;
            s2_ic1R = 2.0f * v1 - s2_ic1R;
            s2_ic1R -= satK * s2_ic1R * (s2_ic1R * s2_ic1R) * fastInvSqrtVoice(1.0f + s2_ic1R * s2_ic1R);
            s2_ic2R = 2.0f * v2 - s2_ic2R;
            float hp2R = hp1R - k2 * v1 - v2;
            *sp++ = clamp(hp2R * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = s1_ic1L;
        v1L = s1_ic2L;
        v2L = s2_ic1L;
        v3L = s2_ic2L;
        v0R = s1_ic1R;
        v1R = s1_ic2R;
        v2R = s2_ic1R;
        v3R = s2_ic2R;
    }
    break;
    case FILTER_BP3:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // 4-pole twin-peak vocal / surgical bandpass (2 slightly spread ZDF SVF stages)
        const float gBase0 = 0.0025f + prevParam1 * (0.38f + 0.68f * prevParam1);
        const float gBase1 = 0.0025f + fxParam1 * (0.38f + 0.68f * fxParam1);
        const float dgBase = (gBase1 - gBase0) * 0.03125f;
        const float spread = 1.0f + 0.14f * (1.0f - 0.5f * res);
        float g1 = gBase0;
        float g2 = gBase0 * spread;
        const float g2_1 = gBase1 * spread;
        const float dg1 = dgBase;
        const float dg2 = dgBase * spread;
        const float k = 1.65f - 1.58f * res;
        const float a1_1_0 = 1.0f / (1.0f + g1 * (g1 + k));
        const float a1_1_1 = 1.0f / (1.0f + gBase1 * (gBase1 + k));
        const float da1_1 = (a1_1_1 - a1_1_0) * 0.03125f;
        float a1_1_lin = a1_1_0;
        const float a1_2_0 = 1.0f / (1.0f + g2 * (g2 + k));
        const float a1_2_1 = 1.0f / (1.0f + g2_1 * (g2_1 + k));
        const float da1_2 = (a1_2_1 - a1_2_0) * 0.03125f;
        float a1_2_lin = a1_2_0;
        const float satK = 0.08f + 0.22f * res;
        const float bpNorm = (0.55f + 0.85f * k) * (0.75f + 0.55f * k) * mixerGain;

        float *sp = this->sampleBlock;
        float s1_ic1L = v0L, s1_ic2L = v1L, s2_ic1L = v2L, s2_ic2L = v3L;
        float s1_ic1R = v0R, s1_ic2R = v1R, s2_ic1R = v2R, s2_ic2R = v3R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g1 += dg1;
            g2 += dg2;
            a1_1_lin += da1_1;
            a1_2_lin += da1_2;
            const float a1_1 = a1_1_lin * (2.0f - (1.0f + g1 * (g1 + k)) * a1_1_lin);
            const float a1_2 = a1_2_lin * (2.0f - (1.0f + g2 * (g2 + k)) * a1_2_lin);

            // Left voice
            float inL = *sp;
            float v3 = inL - s1_ic2L;
            float v1 = a1_1 * (s1_ic1L + g1 * v3);
            float v2 = s1_ic2L + g1 * v1;
            s1_ic1L = 2.0f * v1 - s1_ic1L;
            s1_ic2L = 2.0f * v2 - s1_ic2L;

            v3 = v1 - s2_ic2L;
            v1 = a1_2 * (s2_ic1L + g2 * v3);
            v2 = s2_ic2L + g2 * v1;
            s2_ic1L = 2.0f * v1 - s2_ic1L;
            s2_ic1L -= satK * s2_ic1L * (s2_ic1L * s2_ic1L) * fastInvSqrtVoice(1.0f + s2_ic1L * s2_ic1L);
            s2_ic2L = 2.0f * v2 - s2_ic2L;
            *sp++ = clamp(v1 * bpNorm, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            v3 = inR - s1_ic2R;
            v1 = a1_1 * (s1_ic1R + g1 * v3);
            v2 = s1_ic2R + g1 * v1;
            s1_ic1R = 2.0f * v1 - s1_ic1R;
            s1_ic2R = 2.0f * v2 - s1_ic2R;

            v3 = v1 - s2_ic2R;
            v1 = a1_2 * (s2_ic1R + g2 * v3);
            v2 = s2_ic2R + g2 * v1;
            s2_ic1R = 2.0f * v1 - s2_ic1R;
            s2_ic1R -= satK * s2_ic1R * (s2_ic1R * s2_ic1R) * fastInvSqrtVoice(1.0f + s2_ic1R * s2_ic1R);
            s2_ic2R = 2.0f * v2 - s2_ic2R;
            *sp++ = clamp(v1 * bpNorm, -ratioTimbres, ratioTimbres);
        }

        v0L = s1_ic1L;
        v1L = s1_ic2L;
        v2L = s2_ic1L;
        v3L = s2_ic2L;
        v0R = s1_ic1R;
        v1R = s1_ic2R;
        v2R = s2_ic1R;
        v3R = s2_ic2R;
    }
    break;
    case FILTER_PEAK:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // True ZDF Peak / Resonant Accent filter (LP - HP + resonant BP peak)
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.90f - 1.86f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float peakGain = 0.65f + 1.85f * res;
        const float satK = 0.06f + 0.20f * res;
        const float svfGain = (0.85f + 0.15f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2 - ic2eqL;
            // Peak = in - k*v1 - 2*v2 + resonant emphasis
            float outL = (inL - k * v1 - 2.0f * v2) * 0.55f + v1 * peakGain * 0.45f;
            *sp++ = clamp(outL * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2 - ic2eqR;
            float outR = (inR - k * v1 - 2.0f * v2) * 0.55f + v1 * peakGain * 0.45f;
            *sp++ = clamp(outR * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_NOTCH:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // True ZDF Notch (Band-Reject): out = in - k * v1
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.95f - 1.88f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = 0; kSample < BLOCK_SIZE; kSample++)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            float notchL = inL - k * v1;
            *sp++ = clamp(notchL * mixerGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            float notchR = inR - k * v1;
            *sp++ = clamp(notchR * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_BELL:
    {
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // Musical Bell EQ: -12 dB (A=0.5) at 0, flat (A=1.0) at 0.25, up to +14 dB (A=2.25) at 1.0
        const float A = 0.5f + filterParam2 * (1.25f + 0.5f * filterParam2);
        const float g = 0.0025f + fxParam1 * (0.38f + 0.68f * fxParam1);
        const float k = 1.0f / (0.55f * A);
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;
        const float m1 = k * (A * A - 1.0f);
        const float eqComp = (1.0f / (0.65f + 0.35f * A)) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp((inL + m1 * v1) * eqComp, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp((inR + m1 * v1) * eqComp, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_LOWSHELF:
    {
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float A = 0.45f + filterParam2 * (1.35f + 0.55f * filterParam2);
        const float k = 1.4142f;
        const float g = (0.0025f + fxParam1 * (0.34f + 0.58f * fxParam1)) / sqrt3(A);
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;
        const float m1 = k * (A - 1.0f);
        const float m2 = (A * A - 1.0f);
        const float svfGain = (1.0f / (0.65f + 0.35f * A)) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp((inL + m1 * v1 + m2 * v2) * svfGain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp((inR + m1 * v1 + m2 * v2) * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_HIGHSHELF:
    {
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float A = 0.45f + filterParam2 * (1.35f + 0.55f * filterParam2);
        const float k = 1.4142f;
        const float g = (0.0025f + fxParam1 * (0.34f + 0.58f * fxParam1)) * sqrt3(A);
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;
        const float m0 = A * A;
        const float m1 = k * (1.0f - A) * A;
        const float m2 = (1.0f - A * A);
        const float svfGain = (1.0f / (0.65f + 0.35f * A)) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp((m0 * inL + m1 * v1 + m2 * v2) * svfGain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp((m0 * inR + m1 * v1 + m2 * v2) * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_LPHP:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // Oberheim SEM style continuous morphing ZDF SVF: LP -> Notch/BP -> HP across Frequency sweep
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.95f - 1.91f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;

        // Continuous unquantized equal-power crossfade (eliminates 7-bit panTable stepping)
        const float invP = 1.0f - fxParam1;
        const float mixLP = 0.96f * (1.0f - fxParam1 * fxParam1) + 0.04f;
        const float mixHP = 0.96f * (1.0f - invP * invP) + 0.04f;
        const float satK = 0.06f + 0.22f * res;
        const float gain = (1.0f + 0.25f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2 - ic2eqL;
            float hpL = inL - k * v1 - v2;
            float outL = softLimitVoice(v2 * mixLP + hpL * mixHP);
            *sp++ = clamp(outL * gain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2 - ic2eqR;
            float hpR = inR - k * v1 - v2;
            float outR = softLimitVoice(v2 * mixLP + hpR * mixHP);
            *sp++ = clamp(outR * gain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_BPds:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // Saturated dirty bandpass ZDF SVF without input feedback corruption
        const float g0 = 0.0025f + prevParam1 * (0.40f + 0.72f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.40f + 0.72f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.90f - 1.86f * res;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float satDrive = 1.35f + 1.85f * res;
        const float satK = 0.12f + 0.30f * res;
        const float svfGain = (0.65f + 0.55f * k) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic1eqL -= satK * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2 - ic2eqL;
            float outL = softLimitVoice(v1 * satDrive + 0.06f * res) - 0.055f * res;
            *sp++ = clamp(outL * svfGain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic1eqR -= satK * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2 - ic2eqR;
            float outR = softLimitVoice(v1 * satDrive + 0.06f * res) - 0.055f * res;
            *sp++ = clamp(outR * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_LPWS:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // ZDF 2-pole resonant Lowpass followed by ADAA warm tube waveshaper with Dry/Wet Mix
        const float g0 = 0.0025f + prevParam1 * (0.18f + 1.35f * prevParam1);
        const float g1 = 0.0025f + fxParam1 * (0.18f + 1.35f * fxParam1);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.15f - 0.35f * filterParam2;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;

        // Continuous unquantized equal-power wet/dry crossfade
        const float invW = 1.0f - filterParam2;
        const float mixWetGain = (1.0f - invW * invW) * mixerGain;
        const float mixDryGain = (1.0f - filterParam2 * filterParam2) * mixerGain;
        const float wsDrive = 1.8f + 1.4f * filterParam2;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float prevXL = v2L, prevXR = v2R;
        float prevFL = sat25ADAA1_F(prevXL), prevFR = sat25ADAA1_F(prevXR);

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            float shapedL = softLimitVoice(sat25ADAA(v2 * wsDrive + 0.06f, prevXL, prevFL) * 1.35f - 0.04f);
            *sp++ = clamp(shapedL * mixWetGain + inL * mixDryGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            float shapedR = softLimitVoice(sat25ADAA(v2 * wsDrive + 0.06f, prevXR, prevFR) * 1.35f - 0.04f);
            *sp++ = clamp(shapedR * mixWetGain + inR * mixDryGain, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v2L = prevXL;
        v2R = prevXR;
    }
    break;
    case FILTER_TILT:
    {
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float modParam = clamp(currentTimbre->params_.effect1.param2 + matrixFilterParam2, 0.0f, 1.0f);
        // True bipolar spectral Tilt EQ around pivot frequency fxParam1:
        // modParam < 0.5 -> warm bass tilt (+lows / -highs)
        // modParam = 0.5 -> flat
        // modParam > 0.5 -> bright air tilt (-lows / +highs)
        const float tilt = (modParam - 0.5f) * 2.0f; // [-1..+1]
        const float g = 0.004f + fxParam1 * (0.08f + 0.45f * fxParam1);
        const float k = 1.35f;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        const float lowGain = (tilt < 0.0f) ? (1.0f - 1.65f * tilt) : (1.0f - 0.65f * tilt);
        const float highGain = (tilt > 0.0f) ? (1.0f + 1.65f * tilt) : (1.0f + 0.65f * tilt);
        const float midGain = 1.0f;
        const float normGain = (0.88f - 0.12f * fabsf(tilt)) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            float inL = *sp;
            float v3 = inL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            float hpL = inL - k * v1 - v2;
            float outL = v2 * lowGain + (k * v1) * midGain + hpL * highGain;
            *sp++ = clamp(outL * normGain, -ratioTimbres, ratioTimbres);

            float inR = *sp;
            v3 = inR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            float hpR = inR - k * v1 - v2;
            float outR = v2 * lowGain + (k * v1) * midGain + hpR * highGain;
            *sp++ = clamp(outR * normGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_STEREO:
    {
        // Binaural Stereo Widener & Complementary 3-Stage Allpass Phaser ("Pann": Pos, Sprd)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = (fxParamTmp + 19.0f * fxParam1) * 0.05f;

        float spreadTmp = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = (spreadTmp + 19.0f * fxParam2) * 0.05f;

        // Opposite L/R frequency offsets so Sprd > 0 creates genuine stereo decorrelation
        const float spread = fxParam2 * 0.28f;
        const float panTilt = (fxParam1 - 0.5f) * 0.35f * fxParam2;
        const float gainL = (1.0f - panTilt) * mixerGain;
        const float gainR = (1.0f + panTilt) * mixerGain;
        const float wetMix = 0.55f * fxParam2;

        const float f1L = clamp(fxParam1 - spread, 0.02f, 0.96f) * 1.75f - 0.875f;
        const float f2L = clamp(fxParam1 * 0.85f + spread * 0.6f, 0.02f, 0.96f) * 1.75f - 0.875f;
        const float f3L = clamp(fxParam1 * 1.15f - spread * 1.3f, 0.02f, 0.96f) * 1.75f - 0.875f;

        const float f1R = clamp(fxParam1 + spread, 0.02f, 0.96f) * 1.75f - 0.875f;
        const float f2R = clamp(fxParam1 * 0.85f - spread * 0.6f, 0.02f, 0.96f) * 1.75f - 0.875f;
        const float f3R = clamp(fxParam1 * 1.15f + spread * 1.3f, 0.02f, 0.96f) * 1.75f - 0.875f;

        float *sp = this->sampleBlock;
        float lowL = v0L, bandL = v1L;
        float lowR = v0R, bandR = v1R;
        float lowL2 = v2L, bandL2 = v3L;
        float lowR2 = v2R, bandR2 = v3R;
        float lowL3 = v4L, bandL3 = v5L;
        float lowR3 = v4R, bandR3 = v5R;

        for (int k = BLOCK_SIZE; k--;)
        {
            // Left voice
            float inL = *sp;
            bandL = lowL - f1L * inL;
            lowL = inL + f1L * bandL;
            bandL2 = lowL2 - f2L * bandL;
            lowL2 = bandL + f2L * bandL2;
            bandL3 = lowL3 - f3L * bandL2;
            lowL3 = bandL2 + f3L * bandL3;
            float outL = inL * (1.0f - 0.35f * wetMix) + bandL3 * wetMix;
            *sp++ = clamp(outL * gainL, -ratioTimbres, ratioTimbres);

            // Right voice (complementary polarity & opposite allpass poles)
            float inR = *sp;
            bandR = lowR - f1R * inR;
            lowR = inR + f1R * bandR;
            bandR2 = lowR2 - f2R * bandR;
            lowR2 = bandR + f2R * bandR2;
            bandR3 = lowR3 - f3R * bandR2;
            lowR3 = bandR2 + f3R * bandR3;
            float outR = inR * (1.0f - 0.35f * wetMix) - bandR3 * wetMix;
            *sp++ = clamp(outR * gainR, -ratioTimbres, ratioTimbres);
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;
        v2L = lowL2;
        v3L = bandL2;
        v2R = lowR2;
        v3R = bandR2;
        v4L = lowL3;
        v5L = bandL3;
        v4R = lowR3;
        v5R = bandR3;
    }
    break;
    case FILTER_SAT:
    {
        // Warm Asymmetric Tape / Triode Saturator ("Sat ": Thrs = param1, Tone = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;
        float toneParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // As Thrs (fxParam1) increases, saturation drive increases with automatic gain makeup
        const float drive = 1.0f + fxParam1 * (2.2f + 2.8f * fxParam1);
        const float asym = 0.08f * fxParam1;
        const float makeup = (1.0f / (0.85f + 0.38f * fxParam1)) * mixerGain;

        // Post-saturation ZDF 2-pole Tone filter (sweeps from warm tape rolloff to open presence)
        const float g = 0.015f + toneParam * (0.22f + 1.15f * toneParam);
        const float k = 1.35f;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float prevXL = v2L, prevXR = v2R;
        float prevFL = sigmoidADAA1_F(prevXL), prevFR = sigmoidADAA1_F(prevXR);

        for (int k = BLOCK_SIZE; k--;)
        {
            // Left
            float xL = clamp((*sp) * drive + asym, -1.45f, 1.45f);
            float satL = sigmoidADAA(xL, prevXL, prevFL) - 0.85f * asym;
            float v3 = satL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp(v2 * makeup, -ratioTimbres, ratioTimbres);

            // Right
            float xR = clamp((*sp) * drive + asym, -1.45f, 1.45f);
            float satR = sigmoidADAA(xR, prevXR, prevFR) - 0.85f * asym;
            v3 = satR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp(v2 * makeup, -ratioTimbres, ratioTimbres);
        }
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v2L = prevXL;
        v2R = prevXR;
    }
    break;
    case FILTER_SIGMOID:
    {
        // S-Curve Soft-Clipped Overdrive & Resonant Tone Shaper ("Sigm": Driv = param1, Tone = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float toneParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float drive = 1.1f + fxParam1 * (3.5f + 8.5f * fxParam1);
        const float gainComp = (0.95f / (1.0f + 0.22f * sqrt3(fxParam1))) * mixerGain;

        // ZDF SVF Tone filter with subtle midrange focus
        const float g = 0.008f + toneParam * (0.20f + 1.15f * toneParam);
        const float k = 1.10f - 0.35f * fxParam1;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float prevXL = v2L, prevXR = v2R;
        float prevFL = tanh4ADAA1_F(prevXL), prevFR = tanh4ADAA1_F(prevXR);

        for (int k = BLOCK_SIZE; k--;)
        {
            // Left voice
            float satL = tanh4ADAA((*sp) * drive, prevXL, prevFL);
            float v3 = satL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            float outL = v2 + 0.22f * v1;
            *sp++ = clamp(outL * gainComp, -ratioTimbres, ratioTimbres);

            // Right voice
            float satR = tanh4ADAA((*sp) * drive, prevXR, prevFR);
            v3 = satR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            float outR = v2 + 0.22f * v1;
            *sp++ = clamp(outR * gainComp, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v2L = prevXL;
        v2R = prevXR;
    }
    break;
    case FILTER_FOLD:
    {
        // West-Coast Anti-Aliased Wavefolder + Input Phase Allpass + ZDF Resonant Lowpass Gate ("Fold": Driv = param1, Tone = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float toneParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        const float drive = fxParam1 * (0.45f + 0.55f * fxParam1);
        const float gain = (1.0f + 14.0f * drive) * 0.25f;
        const float finalGain = (0.92f / (1.0f + 0.35f * sqrt3(drive))) * mixerGain;

        // Input allpass disperser shifting harmonic symmetry before folding
        const float f1 = clamp(0.18f + 0.52f * toneParam - 0.15f * fxParam1, filterWindowMin, filterWindowMax);
        const float coef1 = (1.0f - f1) / (1.0f + f1);

        // ZDF 2-pole Lowpass controlled by Tone
        const float g = 0.006f + toneParam * (0.38f + 0.78f * toneParam);
        const float k = 1.15f - 0.35f * fxParam1;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float _ly1L = v2L, _ly1R = v2R;
        float _lx1L = v3L, _lx1R = v3R;
        float prevX4L = v4L, prevX4R = v4R;
        float prevFL = foldADAA1_F(prevX4L), prevFR = foldADAA1_F(prevX4R);

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            // Left
            float inL = *sp;
            _ly1L = coef1 * (_ly1L + inL) - _lx1L;
            _lx1L = inL;
            float foldedL = foldADAA(_ly1L * gain, prevX4L, prevFL) * 4.0f;
            float v3 = foldedL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp(v2 * finalGain, -ratioTimbres, ratioTimbres);

            // Right
            float inR = *sp;
            _ly1R = coef1 * (_ly1R + inR) - _lx1R;
            _lx1R = inR;
            float foldedR = foldADAA(_ly1R * gain, prevX4R, prevFR) * 4.0f;
            v3 = foldedR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp(v2 * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v2L = _ly1L;
        v2R = _ly1R;
        v3L = _lx1L;
        v3R = _lx1R;
        v4L = prevX4L;
        v4R = prevX4R;
    }
    break;
    case FILTER_WRAP:
    {
        // Anti-Aliased Phase-Wrap Distortion + ZDF Tone Filter ("Wrap": Driv = param1, Tone = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = clamp((fxParamTmp + 9.0f * clamp(fxParam1, 0.0f, 1.0f)) * 0.1f, 0.0f, 1.0f);

        float toneParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        // Progressive convex drive: fine control around the +/-1 wrap threshold, deep multi-wrap at high values
        const float drive = fxParam1 * (0.35f + 0.65f * fxParam1);
        const float gain = 0.90f + 5.70f * drive;
        const float wetMix = clamp(fxParam1 * 2.0f, 0.0f, 1.0f);
        const float dryMix = 1.0f - wetMix;
        const float finalGain = (0.95f / (1.0f + 0.65f * sqrt3(drive))) * mixerGain;

        // ZDF 2-pole Lowpass controlled by Tone
        const float g = 0.006f + toneParam * (0.20f + 1.20f * toneParam);
        const float k = 1.25f - 0.30f * drive;
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float prevXL = v4L, prevXR = v4R;
        float prevFL = wrapADAA1_F(prevXL), prevFR = wrapADAA1_F(prevXR);

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            // Left
            float inL = *sp;
            float wrappedL = inL * dryMix + wrapADAA(inL * gain, prevXL, prevFL) * wetMix;
            float v3 = wrappedL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp(v2 * finalGain, -ratioTimbres, ratioTimbres);

            // Right
            float inR = *sp;
            float wrappedR = inR * dryMix + wrapADAA(inR * gain, prevXR, prevFR) * wetMix;
            v3 = wrappedR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp(v2 * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v4L = prevXL;
        v4R = prevXR;
    }
    break;
    case FILTER_XOR:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * .1f, 0, 1);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0, 1);

        float p1sq = sqrt3(fxParam1);
        float prevP1sq = sqrt3(prevParam1);
        // Continuous unquantized S-curve replacing 64-step integer panTable lookup
        float pos = sigmoid(prevP1sq * 2.0f - 1.0f);
        const float posTarget = sigmoid(p1sq * 2.0f - 1.0f);
        const float dPos = (posTarget - pos) * 0.03125f;

        float *sp = this->sampleBlock;
        float localv0L = v0L;
        float localv0R = v0R;

        const float a = 0.95f - filterParam2 * 0.5f;
        const float b = 1 - a;

        int digitsA, digitsB;
        const float threshold = (0.66f - p1sq * 0.66f) * 1.0f;

        float _ly2L = v1L, _ly2R = v1R;
        float _lx2L = v5L, _lx2R = v5R;

        const float f2 = clamp(0.33f + fxParam1 * 0.5f, filterWindowMin, filterWindowMax);
        float coef2 = (1.0f - f2) / (1.0f + f2);

        for (int k = BLOCK_SIZE; k--;)
        {
            pos += dPos;

            // Left
            if (fabsf(*sp) > threshold)
            {
                localv0L = (*sp) - pos * (localv0L + pos * (*sp));
                digitsA = FLOAT2SHORT * (*sp);
                digitsB = FLOAT2SHORT * localv0L;
                localv0L = SHORT2FLOAT * (float)(digitsA ^ (digitsB & 0xfff));
            }
            else
            {
                localv0L = *sp;
            }

            localv0L = (*sp * b) + (localv0L * a);

            _ly2L = coef2 * (_ly2L + localv0L) - _lx2L; // allpass
            _lx2L = localv0L;

            *sp++ = clamp(_ly2L * mixerGain, -ratioTimbres, ratioTimbres);

            // Right
            if (fabsf(*sp) > threshold)
            {
                localv0R = (*sp) - pos * (localv0R + pos * (*sp));
                digitsA = FLOAT2SHORT * (*sp);
                digitsB = FLOAT2SHORT * localv0R;
                localv0R = SHORT2FLOAT * (float)(digitsA ^ (digitsB & 0xfff));
            }
            else
            {
                localv0R = *sp;
            }

            localv0R = (*sp * b) + (localv0R * a);

            _ly2R = coef2 * (_ly2R + localv0R) - _lx2R; // allpass
            _lx2R = localv0R;

            *sp++ = clamp(_ly2R * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v0L = localv0L;
        v0R = localv0R;

        v1L = _ly2L;
        v1R = _ly2R;
        v5L = _lx2L;
        v5R = _lx2R;
    }
    break;
    case FILTER_TEXTURE1:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * .1f, 0, 1);

        float filterParam2 = sqrt3(1 - clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0, 1) * 0.99f);

        float *sp = this->sampleBlock;
        float lowL = v0L, highL = 0, bandL = v1L;
        float lowR = v0R, highR = 0, bandR = v1R;
        float notch;

        const float f0 = (prevParam1 * prevParam1 * prevParam1) * 0.5f;
        const float f1Target = (fxParam1 * fxParam1 * fxParam1) * 0.5f;
        const float df = (f1Target - f0) * 0.03125f;
        float f = f0;
        const float fb = filterParam2;
        const float scale = sqrt3(fb);

        const int highBits = 0xFFFFFD4F;
        const int lowBits = ~(highBits);

        const int ll = (int)(fxParam1 * lowBits);
        int digitsL, digitsR;
        int lowDigitsL, lowDigitsR;

        float _ly1L = v2L, _ly1R = v2R;
        float _lx1L = v3L, _lx1R = v3R;
        const float f1 = clamp(0.33f + f1Target * 0.47f, filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);

        for (int k = BLOCK_SIZE; k--;)
        {
            f += df;

            // LEFT
            lowL = lowL + f * bandL;
            highL = scale * (*sp) - lowL - fb * bandL;
            bandL = f * highL + bandL;
            notch = lowL + highL;

            digitsL = FLOAT2SHORT * bandL;
            lowDigitsL = (digitsL & lowBits);
            bandL = SHORT2FLOAT * (float)((digitsL & highBits) ^ ((lowDigitsL ^ ll) & 0x1FFF));

            lowL = lowL + f * bandL;
            highL = scale * (notch)-lowL - fb * bandL;
            bandL = f * highL + bandL;

            _ly1L = coef1 * (_ly1L + notch) - _lx1L; // allpass
            _lx1L = notch;

            *sp++ = clamp(_ly1L * mixerGain, -ratioTimbres, ratioTimbres);

            // RIGHT
            lowR = lowR + f * bandR;
            highR = scale * (*sp) - lowR - fb * bandR;
            bandR = f * highR + bandR;
            notch = lowR + highR;

            digitsR = FLOAT2SHORT * bandR;
            lowDigitsR = (digitsR & lowBits);
            bandR = SHORT2FLOAT * (float)((digitsR & highBits) ^ ((lowDigitsR ^ ll) & 0x1FFF));

            lowR = lowR + f * bandR;
            highR = scale * (notch)-lowR - fb * bandR;
            bandR = f * highR + bandR;

            _ly1R = coef1 * (_ly1R + notch) - _lx1R; // allpass
            _lx1R = notch;

            *sp++ = clamp(_ly1R * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;

        v2L = _ly1L;
        v2R = _ly1R;
        v3L = _lx1L;
        v3R = _lx1R;
    }
    break;
    case FILTER_TEXTURE2:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * .1f, 0, 1);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0, 1);

        float *sp = this->sampleBlock;
        float lowL = v0L, highL = 0, bandL = v1L;
        float lowR = v0R, highR = 0, bandR = v1R;
        float notch;

        const float f0 = (prevParam1 * prevParam1 * prevParam1) * 0.5f;
        const float f1Target = (fxParam1 * fxParam1 * fxParam1) * 0.5f;
        const float df = (f1Target - f0) * 0.03125f;
        float f = f0;
        const float fb = filterParam2;
        const float scale = sqrt3(fb);

        const int highBits = 0xFFFFFFAF;
        const int lowBits = ~(highBits);

        const int ll = fxParam1 * lowBits;
        int digitsL, digitsR;
        int lowDigitsL, lowDigitsR;

        float _ly1L = v2L, _ly1R = v2R;
        float _lx1L = v3L, _lx1R = v3R;
        float _ly2L = v4L, _ly2R = v4R;
        float _lx2L = v5L, _lx2R = v5R;
        const float f1 = clamp(0.57f - fxParam1 * 0.33f, filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);
        const float f2 = clamp(0.50f - fxParam1 * 0.39f, filterWindowMin, filterWindowMax);
        float coef2 = (1.0f - f2) / (1.0f + f2);

        const float r = 0.9940f;

        for (int k = BLOCK_SIZE; k--;)
        {
            f += df;

            // LEFT
            *sp = _ly2L + ((-_ly2L + *sp) * r);

            lowL = lowL + f * bandL;
            highL = scale * (*sp) - lowL - fb * bandL;
            bandL = f * highL + bandL;
            notch = lowL + highL;

            digitsL = FLOAT2SHORT * (bandL);
            lowDigitsL = (digitsL & lowBits);
            bandL = SHORT2FLOAT * (float)((digitsL & highBits) ^ ((lowDigitsL * ll) & 0x1FFF));

            _ly1L = coef1 * (_ly1L + notch) - _lx1L; // allpass
            _lx1L = notch;

            _ly2L = coef2 * (_ly2L + _ly1L) - _lx2L; // allpass 2
            _lx2L = _ly1L;

            *sp++ = clamp(_ly2L * mixerGain, -ratioTimbres, ratioTimbres);

            // RIGHT
            *sp = _ly2R + ((-_ly2R + *sp) * r);

            lowR = lowR + f * bandR;
            highR = scale * (*sp) - lowR - fb * bandR;
            bandR = f * highR + bandR;
            notch = lowR + highR;

            digitsR = FLOAT2SHORT * (bandR);
            lowDigitsR = (digitsR & lowBits);
            bandR = SHORT2FLOAT * (float)((digitsR & highBits) ^ ((lowDigitsR * ll) & 0x1FFF));

            _ly1R = coef1 * (_ly1R + notch) - _lx1R; // allpass
            _lx1R = notch;

            _ly2R = coef2 * (_ly2R + _ly1R) - _lx2R; // allpass 2
            _lx2R = _ly1R;

            *sp++ = clamp(_ly2R * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;

        v2L = _ly1L;
        v2R = _ly1R;
        v3L = _lx1L;
        v3R = _lx1R;
        v4L = _ly2L;
        v4R = _ly2R;
        v5L = _lx2L;
        v5R = _lx2R;
    }
    break;
    case FILTER_LPXOR:
    {
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        const uint_fast8_t random = (*(uint_fast8_t *)noise) & 0xff;
        if (random > 252)
        {
            fxParam1 += ((random & 1) * 0.007874015748031f);
        }

        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * .1f, 0, 1);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0, 1);

        const float a = 1.f - fxParam1;
        const float a2f = a * SHORT2FLOAT;
        const float b = 1.f - a;

        float *sp = this->sampleBlock;
        float localv0L = v0L;
        float localv0R = v0R;
        float digitized;

        float drive = (filterParam2 * filterParam2);
        float gain = (1 + 2 * drive);

        int digitsAL, digitsBL, digitsAR, digitsBR;
        digitsAL = digitsBL = digitsAR = digitsBR = 0;

        const short bitmask = 0xfac;

        float _ly1L = v2L, _ly1R = v2R;
        float _lx1L = v3L, _lx1R = v3R;
        float prevX4L = v4L, prevX4R = v4R;
        float prevFL = foldADAA1_F(prevX4L), prevFR = foldADAA1_F(prevX4R);
        const float f1 = clamp(0.27f + fxParam1 * 0.2f, filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);

        const float r = 0.9940f;

        for (int k = BLOCK_SIZE; k--;)
        {
            // Left voice
            *sp = _ly1L + ((-_ly1L + *sp) * r);

            digitsAL = FLOAT2SHORT * foldADAA(localv0L * gain, prevX4L, prevFL);
            digitsBL = FLOAT2SHORT * (*sp);
            digitized = (float)(digitsAL ^ (digitsBL & bitmask));
            localv0L = (*sp * b) + (digitized * a2f);

            _ly1L = coef1 * (_ly1L + localv0L) - _lx1L; // allpass
            _lx1L = localv0L;

            *sp++ = clamp(_ly1L * mixerGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = _ly1R + ((-_ly1R + *sp) * r);

            digitsAR = FLOAT2SHORT * foldADAA(localv0R * gain, prevX4R, prevFR);
            digitsBR = FLOAT2SHORT * (*sp);
            digitized = (float)(digitsAR ^ (digitsBR & bitmask));
            localv0R = (*sp * b) + (digitized * a2f);

            _ly1R = coef1 * (_ly1R + localv0R) - _lx1R; // allpass
            _lx1R = localv0R;

            *sp++ = clamp(_ly1R * mixerGain, -ratioTimbres, ratioTimbres);
        }
        v0L = localv0L;
        v0R = localv0R;

        v2L = _ly1L;
        v2R = _ly1R;
        v3L = _lx1L;
        v3R = _lx1R;
        v4L = prevX4L;
        v4R = prevX4R;
    }
    break;
    case FILTER_LPXOR2:
    {
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        const uint8_t random = (*(uint8_t *)noise) & 0xff;
        if (random > 252)
        {
            fxParam1 += ((random & 1) * 0.007874015748031f);
        }

        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * .1f, 0, 1);
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0, 1);

        const float a = 1.f - fxParam1;
        const float a2f = a * SHORT2FLOAT;
        const float b = 1.f - a;

        float *sp = this->sampleBlock;
        float localv0L = v0L;
        float localv0R = v0R;
        float digitized;

        float drive = (filterParam2 * filterParam2);
        float gain = (1 + 8 * (drive)) * 0.25f;

        int digitsAL, digitsBL, digitsAR, digitsBR;
        digitsAL = digitsBL = digitsAR = digitsBR = 0;

        const short bitmask = 0xfba - (int)(fxParam1 * 7);

        float _ly1L = v2L, _ly1R = v2R;
        float _lx1L = v3L, _lx1R = v3R;
        float prevX4L = v4L, prevX4R = v4R;
        float prevFL = foldADAA1_F(prevX4L), prevFR = foldADAA1_F(prevX4R);
        const float f1 = clamp(0.27f + fxParam1 * 0.2f, filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);

        const float r = 0.9940f;

        for (int k = BLOCK_SIZE; k--;)
        {
            // Left voice
            *sp = _ly1L + ((-_ly1L + *sp) * r);

            digitsAL = FLOAT2SHORT * localv0L;
            digitsBL = FLOAT2SHORT * foldADAA(*sp * gain, prevX4L, prevFL);
            digitized = ((digitsAL | (digitsBL & bitmask)));
            localv0L = (*sp * b) + (digitized * a2f);

            _ly1L = coef1 * (_ly1L + localv0L) - _lx1L; // allpass
            _lx1L = localv0L;

            *sp++ = clamp(_ly1L * mixerGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = _ly1R + ((-_ly1R + *sp) * r);

            digitsAR = FLOAT2SHORT * localv0R;
            digitsBR = FLOAT2SHORT * foldADAA(*sp * gain, prevX4R, prevFR);
            digitized = ((digitsAR | (digitsBR & bitmask)));
            localv0R = (*sp * b) + (digitized * a2f);

            _ly1R = coef1 * (_ly1R + localv0R) - _lx1R; // allpass
            _lx1R = localv0R;

            *sp++ = clamp(_ly1R * mixerGain, -ratioTimbres, ratioTimbres);
        }
        v0L = localv0L;
        v0R = localv0R;

        v2L = _ly1L;
        v2R = _ly1R;
        v3L = _lx1L;
        v3R = _lx1R;
        v4L = prevX4L;
        v4R = prevX4R;
    }
    break;
    case FILTER_LPSIN:
    {
        // Sine-Table Harmonic Waveshaper + ZDF Resonant Lowpass ("LpSn": Pos = param1, Freq = param2)
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        const float prevParam2 = clamp(fxParam2, 0.0f, 1.0f);
        float freqParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = clamp((freqParam + 9.0f * prevParam2) * 0.1f, 0.0f, 1.0f);

        float posF = prevParam1 * 2048.0f;
        const float dPosF = (fxParam1 - prevParam1) * (2048.0f * 0.03125f);
        const float drive = 0.12f + fxParam1 * fxParam1 * 0.55f;

        // ZDF 2-pole resonant Lowpass controlled by Freq (fxParam2) with per-sample interpolation
        const float g0 = 0.0025f + prevParam2 * (0.18f + 1.40f * prevParam2);
        const float g1 = 0.0025f + fxParam2 * (0.18f + 1.40f * fxParam2);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 0.85f;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float svfGain = (1.0f + 0.15f * fxParam2) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            posF += dPosF;
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = softLimitVoice(satSin(*sp, drive, posF));
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic1eqL -= 0.08f * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = softLimitVoice(satSin(*sp, drive, posF));
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic1eqR -= 0.08f * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_HPSIN:
    {
        // Sine-Table Harmonic Waveshaper + ZDF Resonant Highpass ("HpSn": Pos = param1, Freq = param2)
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);

        const float prevParam2 = clamp(fxParam2, 0.0f, 1.0f);
        float freqParam = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = clamp((freqParam + 9.0f * prevParam2) * 0.1f, 0.0f, 1.0f);

        float posF = prevParam1 * 2048.0f;
        const float dPosF = (fxParam1 - prevParam1) * (2048.0f * 0.03125f);
        const float drive = 0.12f + fxParam1 * fxParam1 * 0.55f;

        // ZDF 2-pole resonant Highpass controlled by Freq (fxParam2) with per-sample interpolation
        const float g0 = 0.002f + prevParam2 * (0.16f + 1.30f * prevParam2);
        const float g1 = 0.002f + fxParam2 * (0.16f + 1.30f * fxParam2);
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;
        const float k = 1.05f;
        const float a1_0 = 1.0f / (1.0f + g0 * (g0 + k));
        const float a1_1 = 1.0f / (1.0f + g1 * (g1 + k));
        const float da1 = (a1_1 - a1_0) * 0.03125f;
        float a1_lin = a1_0;
        const float svfGain = (1.0f + 0.15f * fxParam2) * mixerGain;

        float *sp = this->sampleBlock;
        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            posF += dPosF;
            g += dg;
            a1_lin += da1;
            const float a1 = a1_lin * (2.0f - (1.0f + g * (g + k)) * a1_lin);

            // Left voice
            float inL = softLimitVoice(satSin(*sp, drive, posF));
            float v3 = inL - ic2eqL;
            float v1 = a1 * (ic1eqL + g * v3);
            float v2 = ic2eqL + g * v1;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic1eqL -= 0.08f * ic1eqL * (ic1eqL * ic1eqL) * fastInvSqrtVoice(1.0f + ic1eqL * ic1eqL);
            ic2eqL = 2.0f * v2 - ic2eqL;
            float hpL = inL - k * v1 - v2;
            *sp++ = clamp(hpL * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = softLimitVoice(satSin(*sp, drive, posF));
            v3 = inR - ic2eqR;
            v1 = a1 * (ic1eqR + g * v3);
            v2 = ic2eqR + g * v1;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic1eqR -= 0.08f * ic1eqR * (ic1eqR * ic1eqR) * fastInvSqrtVoice(1.0f + ic1eqR * ic1eqR);
            ic2eqR = 2.0f * v2 - ic2eqR;
            float hpR = inR - k * v1 - v2;
            *sp++ = clamp(hpR * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
    }
    break;
    case FILTER_QUADNOTCH:
    {
        // Stereo 4-Stage Spread Notch Phaser ("Nch4": Freq = param1, Sprd = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = (filterParam2 + 9.0f * fxParam2) * 0.1f;

        const float bipolarf = (fxParam1 - 0.5f);
        const float folded = sigmoid(bipolarf * 2.5f * fxParam2);

        const float offset = fxParam2 * fxParam2 * 0.17f;
        const float lrDelta = 0.006f * folded;
        const float range = 0.88f;

        const float windowMin = 0.008f, windowMax = 0.82f;

        const float f1LTarget = clamp((fxParam1 - offset - lrDelta) * range, windowMin, windowMax);
        const float f2LTarget = clamp((fxParam1 + offset + lrDelta) * range, windowMin, windowMax);
        const float f3LTarget = clamp((fxParam1 - (offset * 2.0f) - lrDelta) * range, windowMin, windowMax);
        const float f4LTarget = clamp((fxParam1 + (offset * 2.0f) + lrDelta) * range, windowMin, windowMax);
        const float f1RTarget = clamp((f1LTarget - lrDelta) * range, windowMin, windowMax);
        const float f2RTarget = clamp((f2LTarget + lrDelta) * range, windowMin, windowMax);
        const float f3RTarget = clamp((f3LTarget - lrDelta) * range, windowMin, windowMax);
        const float f4RTarget = clamp((f4LTarget + lrDelta) * range, windowMin, windowMax);
        const float q = 0.55f;

        float f1L = (v8L > 0.0f) ? v8L : f1LTarget;
        const float df1L = (f1LTarget - f1L) * 0.03125f;
        const float d2L = f2LTarget - f1LTarget, d3L = f3LTarget - f1LTarget, d4L = f4LTarget - f1LTarget;
        const float d1R = f1RTarget - f1LTarget, d2R = f2RTarget - f1LTarget, d3R = f3RTarget - f1LTarget, d4R = f4RTarget - f1LTarget;

        float *sp = this->sampleBlock;
        float lowL = v0L, highL = 0, bandL = v1L;
        float lowR = v0R, highR = 0, bandR = v1R;
        float lowL2 = v2L, highL2 = 0, bandL2 = v3L;
        float lowR2 = v2R, highR2 = 0, bandR2 = v3R;
        float lowL3 = v4L, highL3 = 0, bandL3 = v5L;
        float lowR3 = v4R, highR3 = 0, bandR3 = v5R;
        float lowL4 = v6L, highL4 = 0, bandL4 = v7L;
        float lowR4 = v6R, highR4 = 0, bandR4 = v7R;

        for (int k = BLOCK_SIZE; k--;)
        {
            f1L += df1L;
            const float f2L = clamp(f1L + d2L, windowMin, windowMax);
            const float f3L = clamp(f1L + d3L, windowMin, windowMax);
            const float f4L = clamp(f1L + d4L, windowMin, windowMax);
            const float f1R = clamp(f1L + d1R, windowMin, windowMax);
            const float f2R = clamp(f1L + d2R, windowMin, windowMax);
            const float f3R = clamp(f1L + d3R, windowMin, windowMax);
            const float f4R = clamp(f1L + d4R, windowMin, windowMax);

            // Left voice
            *sp = (highL4 + lowL4) + ((-(highL4 + lowL4) + *sp) * 0.9840f);

            lowL = lowL + f1L * bandL;
            highL = (*sp) - lowL - q * bandL;
            bandL = f1L * highL + bandL;

            lowL2 = lowL2 + f2L * bandL2;
            highL2 = (highL + lowL) - lowL2 - q * bandL2;
            bandL2 = f2L * highL2 + bandL2;

            lowL3 = lowL3 + f3L * bandL3;
            highL3 = (highL2 + lowL2) - lowL3 - q * bandL3;
            bandL3 = f3L * highL3 + bandL3;

            lowL4 = lowL4 + f4L * bandL4;
            highL4 = (highL3 + lowL3) - lowL4 - q * bandL4;
            bandL4 = f4L * highL4 + bandL4;

            *sp++ = clamp((highL4 + lowL4) * mixerGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = (highR4 + lowR4) + ((-(highR4 + lowR4) + *sp) * 0.9840f);

            lowR = lowR + f1R * bandR;
            highR = (*sp) - lowR - q * bandR;
            bandR = f1R * highR + bandR;

            lowR2 = lowR2 + f2R * bandR2;
            highR2 = (highR + lowR) - lowR2 - q * bandR2;
            bandR2 = f2R * highR2 + bandR2;

            lowR3 = lowR3 + f3R * bandR3;
            highR3 = (highR2 + lowR2) - lowR3 - q * bandR3;
            bandR3 = f3R * highR3 + bandR3;

            lowR4 = lowR4 + f4R * bandR4;
            highR4 = (highR3 + lowR3) - lowR4 - q * bandR4;
            bandR4 = f4R * highR4 + bandR4;

            *sp++ = clamp((highR4 + lowR4) * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;

        v2L = lowL2;
        v3L = bandL2;
        v2R = lowR2;
        v3R = bandR2;

        v4L = lowL3;
        v5L = bandL3;
        v4R = lowR3;
        v5R = bandR3;

        v6L = lowL4;
        v7L = bandL4;
        v6R = lowR4;
        v7R = bandR4;
        v8L = f1LTarget;
    }
    break;
    case FILTER_AP4:
    {
        // http://denniscronin.net/dsp/vst.html
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;

        const float bipolarf = (fxParam1 - 0.5f);
        const float folded = sigmoid(bipolarf * 2.2f * fxParam2); // smooth bounded [-1, 1]

        const float offset = fxParam2 * fxParam2 * 0.17f;
        const float lrDelta = 0.005f * folded;
        const float range = 0.90f;

        const float f1L = clamp(0.02f + (fxParam1 - offset - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f2L = clamp(0.02f + (fxParam1 + offset + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f3L = clamp(0.02f + (fxParam1 - (offset * 2.0f) - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f4L = clamp(0.02f + (fxParam1 + (offset * 2.0f) + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float c1LTarget = (1.0f - f1L) / (1.0f + f1L);
        const float c2LTarget = (1.0f - f2L) / (1.0f + f2L);
        const float c3LTarget = (1.0f - f3L) / (1.0f + f3L);
        const float c4LTarget = (1.0f - f4L) / (1.0f + f4L);
        const float f1R = clamp(0.02f + (fxParam1 - offset + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f2R = clamp(0.02f + (fxParam1 + offset - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f3R = clamp(0.02f + (fxParam1 - (offset * 2.0f) + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f4R = clamp(0.02f + (fxParam1 + (offset * 2.0f) - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float c1RTarget = (1.0f - f1R) / (1.0f + f1R);
        const float c2RTarget = (1.0f - f2R) / (1.0f + f2R);
        const float c3RTarget = (1.0f - f3R) / (1.0f + f3R);
        const float c4RTarget = (1.0f - f4R) / (1.0f + f4R);

        // Per-sample coefficient interpolation to eliminate 1.5 kHz stepping on phaser sweeps
        float coef1L = (v8L != 0.0f) ? v8L : c1LTarget;
        const float dCoef = (c1LTarget - coef1L) * 0.03125f;
        float coef2L = c2LTarget - (c1LTarget - coef1L);
        float coef3L = c3LTarget - (c1LTarget - coef1L);
        float coef4L = c4LTarget - (c1LTarget - coef1L);
        float coef1R = c1RTarget - (c1LTarget - coef1L);
        float coef2R = c2RTarget - (c1LTarget - coef1L);
        float coef3R = c3RTarget - (c1LTarget - coef1L);
        float coef4R = c4RTarget - (c1LTarget - coef1L);

        float *sp = this->sampleBlock;

        float inmix;
        float _ly1L = v0L, _ly1R = v0R;
        float _ly2L = v1L, _ly2R = v1R;
        float _ly3L = v2L, _ly3R = v2R;
        float _ly4L = v3L, _ly4R = v3R;
        float _lx1L = v4L, _lx1R = v4R;
        float _lx2L = v5L, _lx2R = v5R;
        float _lx3L = v6L, _lx3R = v6R;
        float _lx4L = v7L, _lx4R = v7R;
        float _feedback = 0.22f + 0.35f * fxParam2;
        float _crossFeedback = 0.18f;
        float finalGain = mixerGain * 0.52f;

        for (int k = BLOCK_SIZE; k--;)
        {
            coef1L += dCoef;
            coef2L += dCoef;
            coef3L += dCoef;
            coef4L += dCoef;
            coef1R += dCoef;
            coef2R += dCoef;
            coef3R += dCoef;
            coef4R += dCoef;

            // Left voice
            inmix = (*sp) + softLimitVoice(_ly3L * _feedback - _ly2R * _crossFeedback);

            _ly1L = coef1L * (_ly1L + inmix) - _lx1L; // do 1st filter
            _lx1L = inmix;
            _ly2L = coef2L * (_ly2L + _ly1L) - _lx2L; // do 2nd filter
            _lx2L = _ly1L;
            _ly3L = coef3L * (_ly3L + _ly2L) - _lx3L; // do 3rd filter
            _lx3L = _ly2L;
            _ly4L = coef4L * (_ly4L + _ly3L) - _lx4L; // do 4th filter
            _lx4L = _ly3L;

            *sp++ = clamp((*sp + _ly4L) * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            inmix = (*sp) + softLimitVoice(_ly3R * _feedback - _ly2L * _crossFeedback);

            _ly1R = coef1R * (_ly1R + inmix) - _lx1R; // do 1st filter
            _lx1R = inmix;
            _ly2R = coef2R * (_ly2R + _ly1R) - _lx2R; // do 2nd filter
            _lx2R = _ly1R;
            _ly3R = coef3R * (_ly3R + _ly2R) - _lx3R; // do 3rd filter
            _lx3R = _ly2R;
            _ly4R = coef4R * (_ly4R + _ly3R) - _lx4R; // do 4th filter
            _lx4R = _ly3R;

            *sp++ = clamp((*sp + _ly4R) * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = _ly1L;
        v0R = _ly1R;
        v1L = _ly2L;
        v1R = _ly2R;
        v2L = _ly3L;
        v2R = _ly3R;
        v3L = _ly4L;
        v3R = _ly4R;
        v4L = _lx1L;
        v4R = _lx1R;
        v5L = _lx2L;
        v5R = _lx2R;
        v6L = _lx3L;
        v6R = _lx3R;
        v7L = _lx4L;
        v7R = _lx4R;
        v8L = c1LTarget;
    }
    break;
    case FILTER_AP4B:
    {
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;

        const float bipolarf = (fxParam1 - 0.5f);
        const float folded = sigmoid(bipolarf * 3.2f * fxParam2); // smooth bounded [-1, 1]

        const float offset = fxParam2 * fxParam2 * 0.17f;
        const float lrDelta = 0.005f * folded;
        const float range = 0.90f;

        const float f1L = clamp(0.02f + (fxParam1 - offset - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f2L = clamp(0.02f + (fxParam1 + offset + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f3L = clamp(0.02f + (fxParam1 - (offset * 2.0f) - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f4L = clamp(0.02f + (fxParam1 + (offset * 2.0f) + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float c1LTarget = (1.0f - f1L) / (1.0f + f1L);
        const float c2LTarget = (1.0f - f2L) / (1.0f + f2L);
        const float c3LTarget = (1.0f - f3L) / (1.0f + f3L);
        const float c4LTarget = (1.0f - f4L) / (1.0f + f4L);
        const float f1R = clamp(0.02f + (fxParam1 - offset + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f2R = clamp(0.02f + (fxParam1 + offset - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f3R = clamp(0.02f + (fxParam1 - (offset * 2.0f) + lrDelta) * range, filterWindowMin, filterWindowMax);
        const float f4R = clamp(0.02f + (fxParam1 + (offset * 2.0f) - lrDelta) * range, filterWindowMin, filterWindowMax);
        const float c1RTarget = (1.0f - f1R) / (1.0f + f1R);
        const float c2RTarget = (1.0f - f2R) / (1.0f + f2R);
        const float c3RTarget = (1.0f - f3R) / (1.0f + f3R);
        const float c4RTarget = (1.0f - f4R) / (1.0f + f4R);

        // Per-sample coefficient interpolation to eliminate 1.5 kHz stepping on phaser sweeps
        float coef1L = (v8L != 0.0f) ? v8L : c1LTarget;
        const float dCoef = (c1LTarget - coef1L) * 0.03125f;
        float coef2L = c2LTarget - (c1LTarget - coef1L);
        float coef3L = c3LTarget - (c1LTarget - coef1L);
        float coef4L = c4LTarget - (c1LTarget - coef1L);
        float coef1R = c1RTarget - (c1LTarget - coef1L);
        float coef2R = c2RTarget - (c1LTarget - coef1L);
        float coef3R = c3RTarget - (c1LTarget - coef1L);
        float coef4R = c4RTarget - (c1LTarget - coef1L);

        float *sp = this->sampleBlock;

        float _ly1L = v0L, _ly1R = v0R;
        float _ly2L = v1L, _ly2R = v1R;
        float _ly3L = v2L, _ly3R = v2R;
        float _ly4L = v3L, _ly4R = v3R;
        float _lx1L = v4L, _lx1R = v4R;
        float _lx2L = v5L, _lx2R = v5R;
        float _lx3L = v6L, _lx3R = v6R;
        float _lx4L = v7L, _lx4R = v7R;

        float _feedback = 0.68f;
        float _crossFeedback = 0.17f;
        float inmix;
        float finalGain = mixerGain * 0.5f;

        for (int k = BLOCK_SIZE; k--;)
        {
            coef1L += dCoef;
            coef2L += dCoef;
            coef3L += dCoef;
            coef4L += dCoef;
            coef1R += dCoef;
            coef2R += dCoef;
            coef3R += dCoef;
            coef4R += dCoef;

            // Left voice
            inmix = (*sp) + softLimitVoice(-_feedback * _ly3L + _crossFeedback * _ly3R);

            _ly1L = coef1L * (_ly1L + inmix) - _lx1L; // do 1st filter
            _lx1L = inmix;
            _ly2L = coef2L * (_ly2L + _ly1L) - _lx2L; // do 2nd filter
            _lx2L = _ly1L;
            _ly3L = coef3L * (_ly3L + _ly2L) - _lx3L; // do 3rd filter
            _lx3L = _ly2L;
            _ly4L = coef4L * (_ly4L + _ly3L) - _lx4L; // do 4nth filter
            _lx4L = _ly3L;

            *sp++ = clamp((*sp + _ly4L) * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            inmix = (*sp) + softLimitVoice(-_feedback * _ly3R + _crossFeedback * _ly3L);

            _ly1R = coef1R * (_ly1R + inmix) - _lx1R; // do 1st filter
            _lx1R = inmix;
            _ly2R = coef2R * (_ly2R + _ly1R) - _lx2R; // do 2nd filter
            _lx2R = _ly1R;
            _ly3R = coef3R * (_ly3R + _ly2R) - _lx3R; // do 3rd filter
            _lx3R = _ly2R;
            _ly4R = coef4R * (_ly4R + _ly3R) - _lx4R; // do 4nth filter
            _lx4R = _ly3R;

            *sp++ = clamp((*sp + _ly4R) * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = _ly1L;
        v0R = _ly1R;
        v1L = _ly2L;
        v1R = _ly2R;
        v2L = _ly3L;
        v2R = _ly3R;
        v3L = _ly4L;
        v3R = _ly4R;
        v4L = _lx1L;
        v4R = _lx1R;
        v5L = _lx2L;
        v5R = _lx2R;
        v6L = _lx3L;
        v6R = _lx3R;
        v7L = _lx4L;
        v7R = _lx4R;
        v8L = c1LTarget;
    }
    break;
    case FILTER_AP4D:
    {
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;

        const float offset = fxParam2 * fxParam2 * 0.17f;

        const float range1 = 0.33f;
        const float range2 = 0.14f;
        const float range3 = 0.05f;

        const float lowlimit = 0.027f;

        const float f1 = clamp(((lowlimit + fxParam1 + offset) * range1), filterWindowMin, filterWindowMax);
        const float f2 = clamp(((lowlimit + fxParam1 - offset) * range2), filterWindowMin, filterWindowMax);
        const float f3 = clamp(((lowlimit + fxParam1 + offset) * range3), filterWindowMin, filterWindowMax);
        const float c1Target = (1.0f - f1) / (1.0f + f1);
        const float c2Target = (1.0f - f2) / (1.0f + f2);
        float coef3 = (1.0f - f3) / (1.0f + f3);

        float coef1 = (v8L != 0.0f) ? v8L : c1Target;
        float coef2 = (v8R != 0.0f) ? v8R : c2Target;
        const float dCoef1 = (c1Target - coef1) * 0.03125f;
        const float dCoef2 = (c2Target - coef2) * 0.03125f;

        float *sp = this->sampleBlock;

        float _ly1L = v0L, _ly1R = v0R;
        float _ly2L = v1L, _ly2R = v1R;
        float _ly3L = v2L, _ly3R = v2R;
        float _ly4L = v3L, _ly4R = v3R;
        float _lx1L = v4L, _lx1R = v4R;
        float _lx2L = v5L, _lx2R = v5R;
        float _lx3L = v6L, _lx3R = v6R;
        float _lx4L = clamp(v7L, -0.96f, 0.96f), _lx4R = clamp(v7R, -0.96f, 0.96f);

        const float _feedback = 0.15f;
        float inmix;
        float finalGain = mixerGain * 0.5f;

        const float a = 0.85f;
        const float b = 1.f - a;

        for (int k = BLOCK_SIZE; k--;)
        {
            coef1 += dCoef1;
            coef2 += dCoef2;

            // Left voice
            float coef3L = clamp(coef3 + (*sp * 0.0075f), -0.96f, 0.96f); // nested f audio-rate xmod
            inmix = (*sp) + _feedback * _ly4L;
            _ly1L = coef1 * (_ly1L + inmix) - _lx1L; // do 1st filter
            _lx1L = inmix;
            _ly2L = coef2 * (_ly3L + _ly1L) - _lx2L; // do 2nd filter
            _lx2L = _ly1L;
            _ly3L = _lx4L * (_ly3L + _ly2L) - _lx3L; // nested filter
            _lx3L = _ly2L;

            _lx4L = (_lx4L * 31.0f + coef3L) * 0.03125f; // xmod damp

            _ly4L = (_ly3L * b) + (_ly4L * a); // lowpass

            *sp++ = clamp((_ly1L + _ly2L - *sp) * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float coef3R = clamp(coef3 - (*sp * 0.0075f), -0.96f, 0.96f); // stereo nested f audio-rate xmod
            inmix = (*sp) + _feedback * _ly4R;
            _ly1R = coef1 * (_ly1R + inmix) - _lx1R; // do 1st filter
            _lx1R = inmix;
            _ly2R = coef2 * (_ly3R + _ly1R) - _lx2R; // do 2nd filter
            _lx2R = _ly1R;
            _ly3R = _lx4R * (_ly3R + _ly2R) - _lx3R; // nested filter
            _lx3R = _ly2R;

            _lx4R = (_lx4R * 31.0f + coef3R) * 0.03125f; // xmod damp

            _ly4R = (_ly3R * b) + (_ly4R * a); // lowpass

            *sp++ = clamp((_ly1R + _ly2R - *sp) * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = _ly1L;
        v0R = _ly1R;
        v1L = _ly2L;
        v1R = _ly2R;
        v2L = _ly3L;
        v2R = _ly3R;
        v3L = _ly4L;
        v3R = _ly4R;
        v4L = _lx1L;
        v4R = _lx1R;
        v5L = _lx2L;
        v5R = _lx2R;
        v6L = _lx3L;
        v6R = _lx3R;
        v7L = _lx4L;
        v7R = _lx4R;
        v8L = c1Target;
        v8R = c2Target;
    }
    break;
    case FILTER_ORYX:
    {
        // https://www.musicdsp.org/en/latest/Filters/23-state-variable.html
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = (filterParam2 + 19.0f * fxParam2) * 0.05f;
        const float rootf = 0.011f + ((fxParam2 * fxParam2) - 0.005f) * 0.08f;
        const float frange = 0.14f;
        const float bipolarf = (fxParam1 - 0.5f);

        const float f1Target = rootf + sigmoidPos(0.5f - fabsf(bipolarf)) * frange;
        const float fb = 0.105f + fxParam1 * 0.015f;
        const float scale = 0.62f * fb + 0.38f * sqrt3(fb);

        const float spread = 1.0f + (fxParam2 * 0.13f);

        const float fold2 = 0.375f + sigmoid(bipolarf * 3.5f * fxParam2) * frange;
        const float f2Target = rootf + frange * 1.15f + sigmoidPos(clamp(1.0f - fxParam1 + fold2 * 0.25f, 0.0f, 1.0f)) * spread * frange * 1.25f;
        const float fb2 = clamp(0.22f - fold2 * 0.07f, 0.11f, 0.35f);
        const float scale2 = 0.62f * fb2 + 0.38f * sqrt3(fb2);

        float f1 = (v7L > 0.0f) ? v7L : f1Target;
        float f2 = (v7R > 0.0f) ? v7R : f2Target;
        const float df1 = (f1Target - f1) * 0.03125f;
        const float df2 = (f2Target - f2) * 0.03125f;

        float *sp = this->sampleBlock;
        float lowL = v0L, bandL = v1L;
        float lowR = v0R, bandR = v1R;

        float lowL2 = v2L, bandL2 = v3L;
        float lowR2 = v2R, bandR2 = v3R;

        const float svfGain = (1.15f + SVFGAINOFFSET) * mixerGain;

        float out;
        const float r = 0.990f;
        float drift = v6L;
        float nexDrift = noise[7] * 0.010f;
        float deltaD = (nexDrift - drift) * 0.03125f;

        const float sat = 1.20f;
        const float dampK = 0.065f;

        for (int k = BLOCK_SIZE; k--;)
        {
            f1 += df1;
            f2 += df2;
            float fb1L = clamp(fb - drift, 0.09f, 0.35f);
            float fb1R = clamp(fb + drift, 0.09f, 0.35f);

            // ----------- Left voice
            float inL = lowL + ((-lowL + *sp) * r);
            // bandpass 1 with self-limiting resonant damping
            lowL += f1 * bandL;
            bandL += f1 * (scale * sat66(inL * sat) - lowL - fb1L * bandL);
            bandL -= dampK * bandL * (bandL * bandL) * fastInvSqrtVoice(1.0f + bandL * bandL);

            // bandpass 2 with self-limiting resonant damping
            lowL2 += f2 * bandL2;
            bandL2 += f2 * (scale2 * inL - lowL2 - fb2 * bandL2);
            bandL2 -= dampK * bandL2 * (bandL2 * bandL2) * fastInvSqrtVoice(1.0f + bandL2 * bandL2);

            out = softLimitVoice(bandL + bandL2);

            *sp++ = clamp(out * svfGain, -ratioTimbres, ratioTimbres);

            // ----------- Right voice
            float inR = lowR + ((-lowR + *sp) * r);
            // bandpass 1 with self-limiting resonant damping
            lowR += f1 * bandR;
            bandR += f1 * (scale * sat66(inR * sat) - lowR - fb1R * bandR);
            bandR -= dampK * bandR * (bandR * bandR) * fastInvSqrtVoice(1.0f + bandR * bandR);

            // bandpass 2 with self-limiting resonant damping
            lowR2 += f2 * bandR2;
            bandR2 += f2 * (scale2 * inR - lowR2 - fb2 * bandR2);
            bandR2 -= dampK * bandR2 * (bandR2 * bandR2) * fastInvSqrtVoice(1.0f + bandR2 * bandR2);

            out = softLimitVoice(bandR + bandR2);

            *sp++ = clamp(out * svfGain, -ratioTimbres, ratioTimbres);

            drift += deltaD;
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;

        v2L = lowL2;
        v3L = bandL2;
        v2R = lowR2;
        v3R = bandR2;

        v6L = drift;
        v7L = f1Target;
        v7R = f2Target;
    }
    break;
    case FILTER_ORYX2:
    {
        // https://www.musicdsp.org/en/latest/Filters/23-state-variable.html
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = (filterParam2 + 19.0f * fxParam2) * 0.05f;
        const float rootf = 0.031f + (fxParam2) * 0.0465f;
        const float frange = 0.22f;
        const float bipolarf = sigmoid(fxParam1 - 0.5f);

        // Smooth C1 arch instead of discontinuous ternary at fxParam1 = 0.55f
        const float arch = fxParam1 * (1.0f - 0.72f * fxParam1);
        const float f1Target = rootf + sigmoidPos(arch) * frange;
        const float fb = 0.095f + fxParam1 * 0.015f;
        const float scale = 0.62f * fb + 0.38f * sqrt3(fb);

        const float spread = 1.0f + (fxParam2 * 0.13f);

        const float fold2 = 0.375f + sigmoid(bipolarf * 3.5f * fxParam2) * frange;
        const float f2Target = rootf + frange * 1.1f - (sigmoidPos(clamp(fxParam1 + fold2 * 0.25f, 0.0f, 1.0f)) - 0.5f) * spread * frange * 0.5f;
        const float fb2 = clamp(0.13f - f1Target * fold2 * 0.32f, 0.095f, 0.35f);
        const float scale2 = 0.62f * fb2 + 0.38f * sqrt3(fb2);

        const float fold3 = sigmoid(bipolarf * 2.8f * fxParam2) * frange;
        const float f3Target = rootf + frange * 1.35f - tanh3(((1.0f - fxParam1) * 2.0f - fold3 * 0.2f)) * spread * frange * 0.25f;
        const float fb3 = clamp(0.15f - f1Target * fold3 * 0.38f, 0.10f, 0.35f);
        const float scale3 = 0.62f * fb3 + 0.38f * sqrt3(fb3);

        float f1 = (v7L > 0.0f) ? v7L : f1Target;
        float f2 = (v7R > 0.0f) ? v7R : f2Target;
        float f3 = (v8L > 0.0f) ? v8L : f3Target;
        const float df1 = (f1Target - f1) * 0.03125f;
        const float df2 = (f2Target - f2) * 0.03125f;
        const float df3 = (f3Target - f3) * 0.03125f;

        float *sp = this->sampleBlock;
        float lowL = v0L, bandL = v1L;
        float lowR = v0R, bandR = v1R;

        float lowL2 = v2L, bandL2 = v3L;
        float lowR2 = v2R, bandR2 = v3R;

        float lowL3 = v4L, bandL3 = v5L;
        float lowR3 = v4R, bandR3 = v5R;

        const float svfGain = (1.10f + SVFGAINOFFSET) * mixerGain;

        float out;
        const float r = 0.990f;
        float drift = v6L;
        float nexDrift = noise[7] * 0.010f;
        float deltaD = (nexDrift - drift) * 0.03125f;
        const float dampK = 0.065f;

        for (int k = BLOCK_SIZE; k--;)
        {
            f1 += df1;
            f2 += df2;
            f3 += df3;
            float fb1L = clamp(fb - drift, 0.09f, 0.35f);
            float fb1R = clamp(fb + drift, 0.09f, 0.35f);
            float f2d = clamp(f2 + drift, 0.01f, 0.48f);

            // ----------- Left voice
            float inL = lowL + ((-lowL + *sp) * r);
            // bandpass 1
            lowL += f1 * bandL;
            bandL += f1 * (scale * inL - lowL - fb1L * bandL);
            bandL -= dampK * bandL * (bandL * bandL) * fastInvSqrtVoice(1.0f + bandL * bandL);

            // bandpass 2
            lowL2 += f2d * bandL2;
            bandL2 += f2d * (scale2 * inL - lowL2 - fb2 * bandL2);
            bandL2 -= dampK * bandL2 * (bandL2 * bandL2) * fastInvSqrtVoice(1.0f + bandL2 * bandL2);

            // bandpass 3
            lowL3 += f3 * bandL3;
            bandL3 += f3 * (scale3 * inL - lowL3 - fb3 * bandL3);
            bandL3 -= dampK * bandL3 * (bandL3 * bandL3) * fastInvSqrtVoice(1.0f + bandL3 * bandL3);

            out = softLimitVoice(bandL + bandL2 + bandL3);

            *sp++ = clamp(out * svfGain, -ratioTimbres, ratioTimbres);

            // ----------- Right voice
            float inR = lowR + ((-lowR + *sp) * r);
            // bandpass 1
            lowR += f1 * bandR;
            bandR += f1 * (scale * inR - lowR - fb1R * bandR);
            bandR -= dampK * bandR * (bandR * bandR) * fastInvSqrtVoice(1.0f + bandR * bandR);

            // bandpass 2
            lowR2 += f2d * bandR2;
            bandR2 += f2d * (scale2 * inR - lowR2 - fb2 * bandR2);
            bandR2 -= dampK * bandR2 * (bandR2 * bandR2) * fastInvSqrtVoice(1.0f + bandR2 * bandR2);

            // bandpass 3
            lowR3 += f3 * bandR3;
            bandR3 += f3 * (scale3 * inR - lowR3 - fb3 * bandR3);
            bandR3 -= dampK * bandR3 * (bandR3 * bandR3) * fastInvSqrtVoice(1.0f + bandR3 * bandR3);

            out = softLimitVoice(bandR + bandR2 + bandR3);

            *sp++ = clamp(out * svfGain, -ratioTimbres, ratioTimbres);

            drift += deltaD;
        }

        v0L = lowL;
        v1L = bandL;
        v0R = lowR;
        v1R = bandR;

        v2L = lowL2;
        v3L = bandL2;
        v2R = lowR2;
        v3R = bandR2;

        v4L = lowL3;
        v5L = bandL3;
        v4R = lowR3;
        v5R = bandR3;

        v6L = drift;
        v7L = f1Target;
        v7R = f2Target;
        v8L = f3Target;
    }
    break;
    case FILTER_ORYX3:
    {
        // Throat / Growl Binaural Talking Formant Filter:
        // 3 formants with smoothed laryngeal F1->F2/F3 cross-modulation, anti-formant phaser notches & binaural stereo spread
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParam2 = (filterParam2 + 19.0f * fxParam2) * 0.05f;

        // Distinct throat-to-nasal formant trajectory (wider and more aggressive than ORYX1/2)
        const float rootf = 0.014f + (fxParam2 * fxParam2) * 0.065f;
        const float frange = 0.165f;
        const float bipolarf = (fxParam1 - 0.5f);
        const float prevBipolarf = (prevParam1 - 0.5f);

        // F1: Pharynx / throat resonance
        const float f10 = clamp(rootf + sigmoidPos(0.5f - fabsf(prevBipolarf)) * frange * 1.15f, 0.012f, 0.32f);
        const float f1Target = clamp(rootf + sigmoidPos(0.5f - fabsf(bipolarf)) * frange * 1.15f, 0.012f, 0.32f);
        float f1 = f10;
        const float df1 = (f1Target - f10) * 0.03125f;
        const float fb = 0.105f + fxParam1 * 0.015f;
        const float scale = 0.62f * fb + 0.38f * sqrt3(fb);

        const float spread = 1.0f + (fxParam2 * 0.22f);

        // F2: Tongue / oral cavity contrapuntal sweep
        const float fold2 = 0.40f + sigmoid(bipolarf * 3.8f * (0.25f + 0.75f * fxParam2)) * frange;
        const float fold2_0 = 0.40f + sigmoid(prevBipolarf * 3.8f * (0.25f + 0.75f * fxParam2)) * frange;
        const float f20 = clamp(rootf + frange * 0.95f + sigmoidPos(clamp(1.0f - prevParam1 + fold2_0 * 0.35f, 0.0f, 1.0f)) * spread * frange * 1.45f, 0.035f, 0.44f);
        const float f2Target = clamp(rootf + frange * 0.95f + sigmoidPos(clamp(1.0f - fxParam1 + fold2 * 0.35f, 0.0f, 1.0f)) * spread * frange * 1.45f, 0.035f, 0.44f);
        float f2 = f20;
        const float df2 = (f2Target - f20) * 0.03125f;
        const float fb2 = clamp(0.125f - fold2 * 0.035f, 0.095f, 0.35f);
        const float scale2 = 0.62f * fb2 + 0.38f * sqrt3(fb2);

        // F3: High nasal / vocal-tract edge resonance
        const float fold3 = sigmoid(bipolarf * 3.0f * (0.25f + 0.75f * fxParam2)) * frange;
        const float fold3_0 = sigmoid(prevBipolarf * 3.0f * (0.25f + 0.75f * fxParam2)) * frange;
        const float f30 = clamp(rootf + frange * 2.10f + tanh3(1.0f - (prevParam1 * 2.0f - fold3_0 * 0.3f)) * spread * frange * 0.65f, 0.08f, 0.48f);
        const float f3Target = clamp(rootf + frange * 2.10f + tanh3(1.0f - (fxParam1 * 2.0f - fold3 * 0.3f)) * spread * frange * 0.65f, 0.08f, 0.48f);
        float f3 = f30;
        const float df3 = (f3Target - f30) * 0.03125f;
        const float fb3 = clamp(0.135f - fold3 * 0.04f, 0.10f, 0.35f);
        const float scale3 = 0.62f * fb3 + 0.38f * sqrt3(fb3);

        // Laryngeal F1 -> F2/F3 cross-modulation depth (band-limited via xmodEnv to prevent Nyquist FM foldover)
        const float xmodDepth = 0.024f + 0.052f * fxParam2;

        float *sp = this->sampleBlock;
        float low = v0L, band = v1L;
        float low2 = v2L, band2 = v3L;
        float low3 = v4L, band3 = v5L;
        float drift = v6L;
        float xmodEnv = v8L;

        float _ly1L = v0R, _lx1L = v1R;
        float _ly1R = v2R, _lx1R = v3R;

        float _ly2L = v4R, _lx2L = v5R;
        float _ly2R = v6R, _lx2R = v7R;

        const float svfGain = (1.08f + SVFGAINOFFSET) * mixerGain;

        float inL, inR, in, outCore;
        const float r = 0.990f;
        float nexDrift = noise[7] * 0.010f;
        float deltaD = (nexDrift - drift) * 0.03125f;

        // Asymmetric binaural allpass poles creating distinct L/R anti-formant notches
        const float apf1 = clamp(0.18f + fxParam1 * 0.32f + fxParam2 * 0.18f, filterWindowMin, filterWindowMax);
        const float coef1 = (1.0f - apf1) / (1.0f + apf1);
        const float apf2 = clamp(0.52f - fxParam1 * 0.28f - fxParam2 * 0.14f, filterWindowMin, filterWindowMax);
        const float coef2 = (1.0f - apf2) / (1.0f + apf2);

        const float sat = 1.45f + 0.45f * fxParam2;
        const float sat2 = 1.55f + 0.45f * fxParam2;
        const float dampK = 0.075f;
        float prevSatX = v7L, prevTanhX = v8R;
        float prevSatF = sat25ADAA1_F(prevSatX), prevTanhF = tanh4ADAA1_F(prevTanhX);

        for (int k = BLOCK_SIZE; k--;)
        {
            f1 += df1;
            f2 += df2;
            f3 += df3;
            float fb1 = clamp(fb - drift, 0.095f, 0.35f);

            inL = *sp;
            inR = *(sp + 1);
            in = (inL + inR) * 0.5f;

            in = low + ((-low + in) * r);
            // Formant 1 (Throat / Chest)
            low += f1 * band;
            band += f1 * (scale * in - low - fb1 * band);
            band -= dampK * band * (band * band) * fastInvSqrtVoice(1.0f + band * band);

            // Band-limited laryngeal cross-modulation from F1 into F2 & F3 (prevents ultrasonic FM foldover)
            xmodEnv += 0.12f * (band - xmodEnv);
            float f2m = clamp(f2 + xmodEnv * xmodDepth, 0.025f, 0.45f);
            float f3m = clamp(f3 - xmodEnv * (xmodDepth * 1.25f), 0.05f, 0.47f);

            // Formant 2 (Oral cavity, fed with subtle L/R difference)
            low2 += f2m * band2;
            band2 += f2m * (scale2 * (in + 0.20f * (inL - inR)) - low2 - fb2 * band2);
            band2 -= dampK * band2 * (band2 * band2) * fastInvSqrtVoice(1.0f + band2 * band2);

            // Formant 3 (Nasal / Sibilant edge)
            low3 += f3m * band3;
            band3 += f3m * (scale3 * (in - 0.20f * (inL - inR) + 0.18f * xmodEnv) - low3 - fb3 * band3);
            band3 -= dampK * band3 * (band3 * band3) * fastInvSqrtVoice(1.0f + band3 * band3);

            float throat = sat25ADAA(band * sat + low * 0.22f, prevSatX, prevSatF);
            float mouth = tanh4ADAA(band2 * sat2 + band3 * 0.85f, prevTanhX, prevTanhF);
            outCore = throat * 0.72f + mouth * 0.58f;

            // ----------- Left voice (binaural anti-formant notch + F2 emphasis)
            float inApL = outCore + 0.25f * band2;
            _ly1L = coef1 * (_ly1L + inApL) - _lx1L;
            _lx1L = inApL;
            _ly2L = coef1 * (_ly2L + _ly1L) - _lx2L;
            _lx2L = _ly1L;
            float outL = softLimitVoice(0.65f * _ly2L + 0.35f * (throat + band3 * 0.6f));

            *sp++ = clamp(outL * svfGain, -ratioTimbres, ratioTimbres);

            // ----------- Right voice (complementary binaural anti-formant notch + F3 emphasis)
            float inApR = outCore + 0.25f * band3;
            _ly1R = coef2 * (_ly1R + inApR) - _lx1R;
            _lx1R = inApR;
            _ly2R = coef2 * (_ly2R + _ly1R) - _lx2R;
            _lx2R = _ly1R;
            float outR = softLimitVoice(0.65f * _ly2R + 0.35f * (throat + band2 * 0.6f));

            *sp++ = clamp(outR * svfGain, -ratioTimbres, ratioTimbres);

            drift += deltaD;
        }

        v0L = low;
        v1L = band;
        v2L = low2;
        v3L = band2;
        v4L = low3;
        v5L = band3;
        v6L = drift;
        v7L = prevSatX;
        v8L = xmodEnv;

        v0R = _ly1L;
        v1R = _lx1L;
        v2R = _ly1R;
        v3R = _lx1R;

        v4R = _ly2L;
        v5R = _lx2L;
        v6R = _ly2R;
        v7R = _lx2R;
        v8R = prevTanhX;

        break;
    }
    case FILTER_18DB:
    {
        // https://www.musicdsp.org/en/latest/Filters/26-moog-vcf-variation-2.html
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);
        const float cut0 = clamp(sqrt3(prevParam1), filterWindowMin, 0.96f);
        const float cut1 = clamp(sqrt3(fxParam1), filterWindowMin, 0.96f);
        const float dCut = (cut1 - cut0) * 0.03125f;
        float cut = cut0;
        const float cut2 = cut1 * cut1;
        // True 3-pole DC normalization: 1 / (1 + inAtn)^3 = 1 / 1.3^3 = 0.455166 * cut^3
        const float attn0 = 0.455166f * cut0 * cut0 * cut0;
        const float attn1 = 0.455166f * cut2 * cut1;
        const float dAttn = (attn1 - attn0) * 0.03125f;
        float attn = attn0;
        const float inAtn = 0.3f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;
        float res = fxParam2 * (1.15f + 0.95f * fxParam2) + cut2 * 0.12f;
        float fb = res * (1.0f - 0.22f * cut2);

        float *sp = this->sampleBlock;

        float buf1L = v0L, buf2L = v1L, buf3L = v2L;
        float buf1R = v0R, buf2R = v1R, buf3R = v2R;

        float buf1inL = v4L, buf2inL = v5L, buf3inL = v6L;
        float buf1inR = v4R, buf2inR = v5R, buf3inR = v6R;

        float _ly1L = v3L, _ly1R = v3R;
        float _lx1L = v7L, _lx1R = v7R;
        float prevXL = v8L, prevXR = v8R;
        float prevFL = tanh4ADAA1_F(prevXL), prevFR = tanh4ADAA1_F(prevXR);

        const float f1 = clamp(0.37f + cut2 * 0.047f, filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);

        float finalGain = mixerGain * (1.0f + fxParam2 * 0.55f) * 1.15f;

        const float r = 0.987f;

        for (int k = BLOCK_SIZE; k--;)
        {
            cut += dCut;
            attn += dAttn;
            const float invCut = 1.0f - cut;

            // Left voice
            *sp = buf3L + ((-buf3L + *sp) * r);

            *sp = tanh4ADAA(*sp - buf3L * fb, prevXL, prevFL);
            *sp *= attn;
            buf1L = (*sp) + inAtn * buf1inL + invCut * buf1L; // Pole 1
            buf1inL = (*sp);
            buf2L = (buf1L) + inAtn * buf2inL + invCut * buf2L; // Pole 2
            buf2inL = buf1L;
            buf3L = (buf2L) + inAtn * buf3inL + invCut * buf3L; // Pole 3
            buf3inL = buf2L;

            _ly1L = (coef1) * (_ly1L + buf3L) - _lx1L; // allpass
            _lx1L = buf3L;

            *sp++ = clamp(_ly1L * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = buf3R + ((-buf3R + *sp) * r);

            *sp = tanh4ADAA(*sp - buf3R * fb, prevXR, prevFR);
            *sp *= attn;
            buf1R = (*sp) + inAtn * buf1inR + invCut * buf1R; // Pole 1
            buf1inR = (*sp);
            buf2R = (buf1R) + inAtn * buf2inR + invCut * buf2R; // Pole 2
            buf2inR = buf1R;
            buf3R = (buf2R) + inAtn * buf3inR + invCut * buf3R; // Pole 3
            buf3inR = buf2R;

            _ly1R = (coef1) * (_ly1R + buf3R) - _lx1R; // allpass
            _lx1R = buf3R;

            *sp++ = clamp(_ly1R * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = buf1L;
        v1L = buf2L;
        v2L = buf3L;

        v0R = buf1R;
        v1R = buf2R;
        v2R = buf3R;

        v4L = buf1inL;
        v5L = buf2inL;
        v6L = buf3inL;

        v4R = buf1inR;
        v5R = buf2inR;
        v6R = buf3inR;

        v3L = _ly1L;
        v3R = _ly1R;
        v7L = _lx1L;
        v7R = _lx1R;
        v8L = prevXL;
        v8R = prevXR;

        break;
    }
    break;
    case FILTER_LADDER:
    {
        // https://www.musicdsp.org/en/latest/Filters/26-moog-vcf-variation-2.html
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);
        const float cut0 = clamp(sqrt3(prevParam1), filterWindowMin, 0.96f);
        const float cut1 = clamp(sqrt3(fxParam1), filterWindowMin, 0.96f);
        const float dCut = (cut1 - cut0) * 0.03125f;
        float cut = cut0;
        const float cut0_2 = cut0 * cut0;
        const float cut2 = cut1 * cut1;
        const float attn0 = 0.35013f * cut0_2 * cut0_2;
        const float attn1 = 0.35013f * cut2 * cut2;
        const float dAttn = (attn1 - attn0) * 0.03125f;
        float attn = attn0;
        const float inAtn = 0.3f;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;
        float res = fxParam2 * (2.4f + 1.5f * fxParam2) * (0.8f + cut2 * 0.2f);
        float fb = res * (1.0f - 0.15f * cut2);

        const float r = 0.985f;

        float *sp = this->sampleBlock;

        float buf1L = v0L, buf2L = v1L, buf3L = v2L, buf4L = v3L;
        float buf1R = v0R, buf2R = v1R, buf3R = v2R, buf4R = v3R;

        float buf1inL = v4L, buf2inL = v5L, buf3inL = v6L, buf4inL = v7L;
        float buf1inR = v4R, buf2inR = v5R, buf3inR = v6R, buf4inR = v7R;

        const float inBoost = 1.0f + 0.25f * fxParam2;
        float finalGain = mixerGain * (1.0f + fxParam2 * 1.15f);

        for (int k = BLOCK_SIZE; k--;)
        {
            cut += dCut;
            attn += dAttn;
            const float invCut = 1.0f - cut;

            // Left voice
            *sp = buf4L + ((-buf4L + *sp) * r);

            *sp = softLimitVoice(*sp * inBoost - tanh3(buf4L) * (fb + buf4R * 0.05f));
            *sp *= attn;
            buf1L = *sp + inAtn * buf1inL + invCut * buf1L; // Pole 1
            buf1inL = *sp;
            buf2L = buf1L + inAtn * buf2inL + invCut * buf2L; // Pole 2
            buf2inL = buf1L;
            buf3L = buf2L + inAtn * buf3inL + invCut * buf3L; // Pole 3
            buf3inL = buf2L;
            buf4L = buf3L + inAtn * buf4inL + invCut * buf4L; // Pole 4
            buf4inL = buf3L;

            *sp++ = clamp(buf4L * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = buf4R + ((-buf4R + *sp) * r);

            *sp = softLimitVoice(*sp * inBoost - tanh3(buf4R) * (fb - buf4L * 0.05f));
            *sp *= attn;
            buf1R = *sp + inAtn * buf1inR + invCut * buf1R; // Pole 1
            buf1inR = *sp;
            buf2R = buf1R + inAtn * buf2inR + invCut * buf2R; // Pole 2
            buf2inR = buf1R;
            buf3R = buf2R + inAtn * buf3inR + invCut * buf3R; // Pole 3
            buf3inR = buf2R;
            buf4R = buf3R + inAtn * buf4inR + invCut * buf4R; // Pole 4
            buf4inR = buf3R;

            *sp++ = clamp(buf4R * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = buf1L;
        v1L = buf2L;
        v2L = buf3L;
        v3L = buf4L;

        v0R = buf1R;
        v1R = buf2R;
        v2R = buf3R;
        v3R = buf4R;

        v4L = buf1inL;
        v5L = buf2inL;
        v6L = buf3inL;
        v7L = buf4inL;

        v4R = buf1inR;
        v5R = buf2inR;
        v6R = buf3inR;
        v7R = buf4inR;

        break;
    }
    case FILTER_LADDER2:
    {
        // https://www.musicdsp.org/en/latest/Filters/240-karlsen-fast-ladder.html
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;

        fxParam1 = (fxParamTmp + 9.0f * fxParam1) * 0.1f;
        const float cut0 = clamp(sqrt3(prevParam1), filterWindowMin, 0.92f);
        const float cut1 = clamp(sqrt3(fxParam1), filterWindowMin, 0.92f);
        const float dCut = (cut1 - cut0) * 0.03125f;
        float cut = cut0;

        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        float OffsetTmp = fabsf(filterParam2);
        fxParam2 = (OffsetTmp + 9.0f * fxParam2) * 0.1f;
        float res = fxParam2 * (0.72f + 0.42f * fxParam2);

        float *sp = this->sampleBlock;

        float buf1L = v0L, buf2L = v1L, buf3L = v2L, buf4L = v3L;
        float buf1R = v0R, buf2R = v1R, buf3R = v2R, buf4R = v3R;

        float _ly1L = v4L, _lx1L = v5L;
        float _ly1R = v4R, _lx1R = v5R;

        float _ly2L = v6L, _ly2R = v6R;
        float _lx2L = v7L, _lx2R = v7R;
        float prevXL = v8L, prevXR = v8R;
        float prevFL = sat25ADAA1_F(prevXL), prevFR = sat25ADAA1_F(prevXR);

        const float bipolarf = fxParam1 - 0.5f;
        const float fold3 = (sigmoid(2.0f * bipolarf) * 0.25f + 0.125f) * 0.25f;
        const float f1 = clamp(0.33f + fold3 * (1.0f - cut1 * 0.5f), filterWindowMin, filterWindowMax);
        float coef1 = (1.0f - f1) / (1.0f + f1);
        const float f2 = clamp(f1 * f1, filterWindowMin, filterWindowMax);
        float coef2 = (1.0f - f2) / (1.0f + f2);

        float finalGain = mixerGain * (1.15f + 0.25f * fxParam2);

        for (int k = BLOCK_SIZE; k--;)
        {
            cut += dCut;

            // Left voice
            *sp = *sp - sat25ADAA(buf4L * (res + buf4R * 0.05f), prevXL, prevFL);
            buf1L = ((*sp - buf1L) * cut) + buf1L;
            buf2L = ((buf1L - buf2L) * cut) + buf2L;
            _ly1L = coef1 * (_ly1L + buf2L) - _lx1L; // allpass
            _lx1L = buf2L;
            buf3L = ((_ly1L - buf3L) * cut) + buf3L;
            _ly2L = coef2 * (_ly2L + buf3L) - _lx2L; // allpass
            _lx2L = buf3L;
            buf4L = softLimitVoice(((_ly2L - buf4L) * cut) + buf4L);
            *sp++ = clamp(buf4L * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = *sp - sat25ADAA(buf4R * (res - buf4L * 0.05f), prevXR, prevFR);
            buf1R = ((*sp - buf1R) * cut) + buf1R;
            buf2R = ((buf1R - buf2R) * cut) + buf2R;
            _ly1R = coef1 * (_ly1R + buf2R) - _lx1R; // allpass
            _lx1R = buf2R;
            buf3R = ((_ly1R - buf3R) * cut) + buf3R;
            _ly2R = coef2 * (_ly2R + buf3R) - _lx2R; // allpass
            _lx2R = buf3R;
            buf4R = softLimitVoice(((_ly2R - buf4R) * cut) + buf4R);
            *sp++ = clamp(buf4R * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = buf1L;
        v1L = buf2L;
        v2L = buf3L;
        v3L = buf4L;

        v0R = buf1R;
        v1R = buf2R;
        v2R = buf3R;
        v3R = buf4R;

        v4L = _ly1L;
        v5L = _lx1L;
        v4R = _ly1R;
        v5R = _lx1R;

        v6L = _ly2L;
        v6R = _ly2R;
        v7L = _lx2L;
        v7R = _lx2R;
        v8L = prevXL;
        v8R = prevXR;

        break;
    }
    case FILTER_DIOD:
    {
        // Coupled Asymmetric Diode Ladder / Acid Filter ("Diod": Freq = param1, Reso = param2)
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);
        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // Coupled 2-stage ZDF Diode topology with asymmetric diode clipping & vocal acid resonance
        const float resCurve = res * (0.45f + 0.55f * res);
        const float g1_0 = 0.002f + prevParam1 * (0.16f + 1.25f * prevParam1);
        const float g1_1 = 0.002f + fxParam1 * (0.16f + 1.25f * fxParam1);
        const float dg1 = (g1_1 - g1_0) * 0.03125f;
        float g1 = g1_0;
        float g2 = g1 * 1.18f;
        const float g2_1 = g1_1 * 1.18f;
        const float k1 = 1.65f - 0.85f * resCurve;
        const float k2 = 1.75f - 1.52f * resCurve;
        const float a1_1_0 = 1.0f / (1.0f + g1 * (g1 + k1));
        const float a1_1_1 = 1.0f / (1.0f + g1_1 * (g1_1 + k1));
        const float da1_1 = (a1_1_1 - a1_1_0) * 0.03125f;
        float a1_1_lin = a1_1_0;
        const float a1_2_0 = 1.0f / (1.0f + g2 * (g2 + k2));
        const float a1_2_1 = 1.0f / (1.0f + g2_1 * (g2_1 + k2));
        const float da1_2 = (a1_2_1 - a1_2_0) * 0.03125f;
        float a1_2_lin = a1_2_0;

        const float diodeBias = 0.12f + 0.10f * res;
        const float diodeComp = softLimitVoice(diodeBias);
        const float fbCouple = 0.14f * resCurve;
        const float drive = 1.05f + 0.45f * resCurve;
        const float svfGain = (1.0f + 0.32f * res) * mixerGain;
        const float dampK = 0.10f + 0.25f * res;

        float *sp = this->sampleBlock;
        float s1_ic1L = v0L, s1_ic2L = v1L, s2_ic1L = v2L, s2_ic2L = v3L;
        float s1_ic1R = v0R, s1_ic2R = v1R, s2_ic1R = v2R, s2_ic2R = v3R;

        for (int k = BLOCK_SIZE; k--;)
        {
            g1 += dg1;
            g2 = g1 * 1.18f;
            a1_1_lin += da1_1;
            a1_2_lin += da1_2;
            const float a1_1 = a1_1_lin * (2.0f - (1.0f + g1 * (g1 + k1)) * a1_1_lin);
            const float a1_2 = a1_2_lin * (2.0f - (1.0f + g2 * (g2 + k2)) * a1_2_lin);

            // Left voice
            float inL = softLimitVoice(((*sp) - fbCouple * s2_ic2L) * drive + diodeBias) - diodeComp;
            float v3 = inL - s1_ic2L;
            float v1 = a1_1 * (s1_ic1L + g1 * v3);
            float v2 = s1_ic2L + g1 * v1;
            s1_ic1L = 2.0f * v1 - s1_ic1L;
            s1_ic2L = 2.0f * v2 - s1_ic2L;

            float midL = softLimitVoice(v2 + 0.10f * v1 + diodeBias) - diodeComp;
            v3 = midL - s2_ic2L;
            v1 = a1_2 * (s2_ic1L + g2 * v3);
            v2 = s2_ic2L + g2 * v1;
            s2_ic1L = 2.0f * v1 - s2_ic1L;
            s2_ic1L -= dampK * s2_ic1L * (s2_ic1L * s2_ic1L) * fastInvSqrtVoice(1.0f + s2_ic1L * s2_ic1L);
            s2_ic2L = 2.0f * v2 - s2_ic2L;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = softLimitVoice(((*sp) - fbCouple * s2_ic2R) * drive + diodeBias) - diodeComp;
            v3 = inR - s1_ic2R;
            v1 = a1_1 * (s1_ic1R + g1 * v3);
            v2 = s1_ic2R + g1 * v1;
            s1_ic1R = 2.0f * v1 - s1_ic1R;
            s1_ic2R = 2.0f * v2 - s1_ic2R;

            float midR = softLimitVoice(v2 + 0.10f * v1 + diodeBias) - diodeComp;
            v3 = midR - s2_ic2R;
            v1 = a1_2 * (s2_ic1R + g2 * v3);
            v2 = s2_ic2R + g2 * v1;
            s2_ic1R = 2.0f * v1 - s2_ic1R;
            s2_ic1R -= dampK * s2_ic1R * (s2_ic1R * s2_ic1R) * fastInvSqrtVoice(1.0f + s2_ic1R * s2_ic1R);
            s2_ic2R = 2.0f * v2 - s2_ic2R;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = s1_ic1L;
        v1L = s1_ic2L;
        v2L = s2_ic1L;
        v3L = s2_ic2L;
        v0R = s1_ic1R;
        v1R = s1_ic2R;
        v2R = s2_ic1R;
        v3R = s2_ic2R;
    }
    break;
    case FILTER_KRMG:
    {
        // https://github.com/ddiakopoulos/MoogLadders/blob/master/src/KrajeskiModel.h
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(SVFOFFSET + currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        // Low pass... on the Frequency
        fxParam1 = clamp((fxParamTmp + 9.0f * fxParam1) * 0.1f, 0.0f, 1.0f);
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        const float wc0 = prevParam1 * (0.35f + 0.65f * sqrt3(prevParam1));
        const float wc0_2 = wc0 * wc0;
        const float wc0_3 = wc0_2 * wc0;
        const float g0 = 0.9892f * wc0 - 0.4342f * wc0_2 + 0.1381f * wc0_3 - 0.0202f * wc0_2 * wc0_2;

        const float wc = fxParam1 * (0.35f + 0.65f * sqrt3(fxParam1));
        const float wc2 = wc * wc;
        const float wc3 = wc2 * wc;
        const float g1 = 0.9892f * wc - 0.4342f * wc2 + 0.1381f * wc3 - 0.0202f * wc2 * wc2;
        const float dg = (g1 - g0) * 0.03125f;
        float g = g0;

        const float drive = 0.95f + 0.25f * filterParam2;
        const float gComp = 0.5f;
        const float resonance = clamp(filterParam2 * 0.86f - fxParam1 * 0.04f, 0.0f, 0.92f);
        const float gRes0 = 2.85f * resonance * (1.0029f + 0.0526f * wc0 - 0.926f * wc0_2 + 0.0218f * wc0_3);
        const float gRes1 = 2.85f * resonance * (1.0029f + 0.0526f * wc - 0.926f * wc2 + 0.0218f * wc3);
        const float dgRes = (gRes1 - gRes0) * 0.03125f;
        float gRes = gRes0;

        float *sp = this->sampleBlock;
        float state0L = v0L, state1L = v1L, state2L = v2L, state3L = v3L, state4L = v4L;
        float delay0L = v5L, delay1L = v6L, delay2L = v7L, delay3L = v8L;

        float state0R = v0R, state1R = v1R, state2R = v2R, state3R = v3R, state4R = v4R;
        float delay0R = v5R, delay1R = v6R, delay2R = v7R, delay3R = v8R;

        const float va1 = 0.2307692308f; // 0.3f / 1.3f;
        const float va2 = 0.7692307692f; // 1 / 1.3f;
        float r = 0.989f;

        float drift = fxParamB2;
        float nexDrift = noise[7] * 0.012f;
        float deltaD = (nexDrift - drift) * 0.03125f;

        float finalGain = mixerGain * (0.85f + 0.20f * filterParam2);

        for (int k = BLOCK_SIZE; k--;)
        {
            g += dg;
            gRes += dgRes;

            // Left voice
            *sp = (state4L + ((-state4L + *sp) * r));
            state0L = tanh4((drive - drift) * (*sp - gRes * (state4L - gComp * *sp)));

            state1L = g * (va1 * state0L + va2 * delay0L - state1L) + state1L;
            delay0L = state0L;

            state2L = g * (va1 * state1L + va2 * delay1L - state2L) + state2L;
            delay1L = state1L;

            state3L = g * (va1 * state2L + va2 * delay2L - state3L) + state3L;
            delay2L = state2L;

            state4L = g * (va1 * state3L + va2 * delay3L - state4L) + state4L;
            delay3L = state3L;

            *sp++ = clamp(state4L * finalGain, -ratioTimbres, ratioTimbres);

            // Right voice
            *sp = (state4R + ((-state4R + *sp) * r));
            state0R = tanh4((drive + drift) * (*sp - gRes * (state4R - gComp * *sp)));

            state1R = g * (va1 * state0R + va2 * delay0R - state1R) + state1R;
            delay0R = state0R;

            state2R = g * (va1 * state1R + va2 * delay1R - state2R) + state2R;
            delay1R = state1R;

            state3R = g * (va1 * state2R + va2 * delay2R - state3R) + state3R;
            delay2R = state2R;

            state4R = g * (va1 * state3R + va2 * delay3R - state4R) + state4R;
            delay3R = state3R;

            *sp++ = clamp(state4R * finalGain, -ratioTimbres, ratioTimbres);

            drift += deltaD;
        }

        v0L = state0L;
        v1L = state1L;
        v2L = state2L;
        v3L = state3L;
        v4L = state4L;

        v5L = delay0L;
        v6L = delay1L;
        v7L = delay2L;
        v8L = delay3L;

        v0R = state0R;
        v1R = state1R;
        v2R = state2R;
        v3R = state3R;
        v4R = state4R;

        v5R = delay0R;
        v6R = delay1R;
        v7R = delay2R;
        v8R = delay3R;

        fxParamB2 = drift;
    }
    break;
    case FILTER_TEEBEE:
    {
        // https://github.com/a1k0n/303
        //  Authentic TB-303 4-pole diode ladder + C13 "gummy capacitor" Accent circuit
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParam1 = clamp((fxParamTmp + 4.0f * fxParam1) * 0.2f, 0.0f, 1.0f);
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // Detect raw MIDI velocity >= 78 for authentic TB-303 Accent (independent of PreenFM3's engine1.velocity scaling)
        const int rawVelo = (uint8_t)this->midiVelocity;
        const bool isAccented = (rawVelo >= 78);
        const float accentStrength = isAccented ? clamp((float)(rawVelo - 72) * (1.0f / 55.0f), 0.15f, 1.0f) : 0.0f;

        // State variables for the 2-stage TB-303 Accent circuit:
        // v6L = short 200ms Accent decay envelope (R_accent)
        // fxParamA2 = C13 capacitor charge (accumulates across consecutive accented 16th notes)
        // fxParamA1 = smoothed C13 "wow" sweep output
        float accentDecayEnv = v6L;
        if (this->fxNoteTrigger_)
        {
            this->fxNoteTrigger_ = 0;
            if (isAccented)
            {
                // Trigger fast 200ms decay envelope and stack charge onto C13 capacitor
                accentDecayEnv = accentStrength;
                fxParamA2 = clamp(fxParamA2 * 0.55f + accentStrength * (0.45f + 0.35f * filterParam2), 0.0f, 1.65f);
            }
            else
            {
                // Unaccented note: dump accent envelope & C13 charge for sharp accented/unaccented contrast
                accentDecayEnv *= 0.15f;
                fxParamA2 *= 0.25f;
            }
        }

        // 1) Short 190ms TB-303 Accent decay (decays immediately after note-on, even while key is held!)
        accentDecayEnv *= 0.9965f;
        v6L = accentDecayEnv;

        // 2) C13 "Gummy Capacitor" RC charge/discharge network (~14ms rounded "wow" attack + resonant bleed)
        const float c13Target = accentDecayEnv * (0.65f + 0.60f * filterParam2);
        if (c13Target > fxParamA2)
        {
            fxParamA2 += 0.26f * (c13Target - fxParamA2); // ~14ms rounded rise ("wow" attack)
        }
        else
        {
            fxParamA2 = fxParamA2 * 0.9945f + 0.0055f * c13Target; // capacitor bleed
        }
        fxParamA1 = (fxParamA1 * 7.0f + fxParamA2) * 0.125f;

        float *sp = this->sampleBlock;
        float state0L = v0L, state1L = v1L, state2L = v2L, state3L = v3L, state4L = v4L;
        float state0R = v0R, state1R = v1R, state2R = v2R, state3R = v3R, state4R = v4R;

        // Continuous resonance index [0.0 .. 61.5] with linear interpolation across filterpoles table
        const float resoScaled = filterParam2 * (0.65f + 0.35f * filterParam2) * 61.5f;
        const int resoIdx = __USAT((int)resoScaled, 6);
        const int resoIdxNext = (resoIdx < 62) ? (resoIdx + 1) : 62;
        const float resoFrac = resoScaled - (float)resoIdx;
        const float resoInv = 1.0f - resoFrac;

        // Smooth C-infinity exponential cutoff via exp_pole (replaces piecewise-linear expf_fast octave corners)
        const float accentSweep = fxParamA1 * (0.42f + 0.65f * filterParam2);
        const float vcf_cutoff = clamp(fxParam1 + accentSweep, 0.0f, 1.55f);
        const float cutoffHz = 141.17496f * exp_pole(3.05f * vcf_cutoff) + 40.0f;
        const float wTarget = clamp(cutoffHz * (2.0f * 3.14159265358979f / PREENFM_FREQUENCY), 0.015f, 1.38f);
        const float wStart = (v7L > 0.0f) ? v7L : wTarget;

        // Interpolate closed-loop s-plane poles from a1k0n/303 table
        const float fp0 = filterpoles[0][resoIdx] * resoInv + filterpoles[0][resoIdxNext] * resoFrac;
        const float fp1 = filterpoles[1][resoIdx] * resoInv + filterpoles[1][resoIdxNext] * resoFrac;
        const float fp2 = filterpoles[2][resoIdx] * resoInv + filterpoles[2][resoIdxNext] * resoFrac;
        const float fp3 = filterpoles[3][resoIdx] * resoInv + filterpoles[3][resoIdxNext] * resoFrac;
        const float fp4 = filterpoles[4][resoIdx] * resoInv + filterpoles[4][resoIdxNext] * resoFrac;
        const float fp5 = filterpoles[5][resoIdx] * resoInv + filterpoles[5][resoIdxNext] * resoFrac;
        const float fp6 = filterpoles[6][resoIdx] * resoInv + filterpoles[6][resoIdxNext] * resoFrac;
        const float fp7 = filterpoles[7][resoIdx] * resoInv + filterpoles[7][resoIdxNext] * resoFrac;
        const float fp8 = filterpoles[8][resoIdx] * resoInv + filterpoles[8][resoIdxNext] * resoFrac;
        const float fp9 = filterpoles[9][resoIdx] * resoInv + filterpoles[9][resoIdxNext] * resoFrac;

        // Start-of-block and end-of-block biquad coefficients for 48 kHz sample-rate interpolation
        const float p0_0 = fp0 + wStart * fp1;
        const float p1r_0 = fp2 + wStart * fp4;
        const float p1i_0 = fp3 + wStart * fp5;
        const float p2r_0 = fp6 + wStart * fp8;
        const float p2i_0 = fp7 + wStart * fp9;

        const float p0_1 = fp0 + wTarget * fp1;
        const float p1r_1 = fp2 + wTarget * fp4;
        const float p1i_1 = fp3 + wTarget * fp5;
        const float p2r_1 = fp6 + wTarget * fp8;
        const float p2i_1 = fp7 + wTarget * fp9;

        const float reso_k = filterParam2 * 3.8f;
        const float targetgain = (1.85f / (1.0f + reso_k * 0.55f)) + 0.55f * filterParam2 + 0.50f * accentDecayEnv + 0.25f * fxParamA1;

        const float p0f_0 = exp_pole(p0_0);
        const float p0f_1 = exp_pole(p0_1);
        // f0gain = halfGain * (1 - f0a1) is linear in f0a1: derive it per sample instead of interpolating it separately
        const float halfGain = targetgain * 0.5f;
        float f0a1 = -p0f_0;
        const float df0a1 = (-p0f_1 - f0a1) * 0.03125f;

        // Filter section 1: non-resonant 2-pole diode skirt (+ subtle mid-bark boost on Accent)
        const float accentBark = 1.0f + 0.25f * fxParamA1;
        const float exp_p1r_0 = exp_pole(p1r_0);
        float f1a1 = -2.0f * exp_p1r_0 * cosf(p1i_0);
        float f1a2 = exp_p1r_0 * exp_p1r_0;

        const float exp_p1r_1 = exp_pole(p1r_1);
        const float f1a1_1 = -2.0f * exp_p1r_1 * cosf(p1i_1);
        const float f1a2_1 = exp_p1r_1 * exp_p1r_1;
        const float df1a1 = (f1a1_1 - f1a1) * 0.03125f;
        const float df1a2 = (f1a2_1 - f1a2) * 0.03125f;

        // Filter section 2: high-Q resonant 2-pole acid peak
        const float exp_p2r_0 = exp_pole(p2r_0);
        float f2a1 = -2.0f * exp_p2r_0 * cosf(p2i_0);
        float f2a2 = exp_p2r_0 * exp_p2r_0;

        const float exp_p2r_1 = exp_pole(p2r_1);
        const float f2a1_1 = -2.0f * exp_p2r_1 * cosf(p2i_1);
        const float f2a2_1 = exp_p2r_1 * exp_p2r_1;
        const float df2a1 = (f2a1_1 - f2a1) * 0.03125f;
        const float df2a2 = (f2a2_1 - f2a2) * 0.03125f;

        // Asymmetric diode-ladder inter-stage saturation & BA662 VCA output drive
        const float drive2 = 1.15f + 0.55f * filterParam2 + 0.65f * fxParamA1;
        const float asymBias = 0.06f * (0.4f + 0.6f * filterParam2 + 0.5f * fxParamA1);
        const float asymOut = asymBias * 1.45f;
        const float finalGain = mixerGain * (0.82f / (1.0f + 0.22f * filterParam2));
        float prevXL = v5L, prevXR = v5R;
        float prevFL = tanh4ADAA1_F(prevXL), prevFR = tanh4ADAA1_F(prevXR);

        for (int k = BLOCK_SIZE; k--;)
        {
            f0a1 += df0a1;
            f1a1 += df1a1;
            f1a2 += df1a2;
            f2a1 += df2a1;
            f2a2 += df2a2;
            const float f0gain = halfGain * (1.0f - f0a1);
            const float f1b = (1.0f + f1a1 + f1a2) * accentBark;
            const float f2b = 1.0f + f2a1 + f2a2;

            // -------- -------- -------- Left
            float inL = *sp;
            float y0L = f0gain * inL + state0L;
            state0L = -f0gain * inL - f0a1 * y0L;

            // First two-pole stage (diode skirt)
            float y1L = f1b * y0L + state1L;
            state1L = state2L - f1a1 * y1L;
            state2L = -f1a2 * y1L;

            // Subtle diode ladder inter-stage soft compression (tames extreme resonant whistling into liquid acid squelch)
            float midL = softLimitVoice(y1L * 0.65f + asymBias) * 1.5f - asymOut;

            // Second two-pole stage (acid resonant peak)
            float y2L = f2b * midL + state3L;
            state3L = state4L - f2a1 * y2L;
            state4L = -f2a2 * y2L;

            *sp++ = clamp(tanh4ADAA(y2L * drive2, prevXL, prevFL) * finalGain, -ratioTimbres, ratioTimbres);

            // -------- -------- -------- Right
            float inR = *sp;
            float y0R = f0gain * inR + state0R;
            state0R = -f0gain * inR - f0a1 * y0R;

            // First two-pole stage (diode skirt)
            float y1R = f1b * y0R + state1R;
            state1R = state2R - f1a1 * y1R;
            state2R = -f1a2 * y1R;

            float midR = softLimitVoice(y1R * 0.65f + asymBias) * 1.5f - asymOut;

            // Second two-pole stage (acid resonant peak)
            float y2R = f2b * midR + state3R;
            state3R = state4R - f2a1 * y2R;
            state4R = -f2a2 * y2R;

            *sp++ = clamp(tanh4ADAA(y2R * drive2, prevXR, prevFR) * finalGain, -ratioTimbres, ratioTimbres);
        }

        v0L = state0L;
        v1L = state1L;
        v2L = state2L;
        v3L = state3L;
        v4L = state4L;
        v5L = prevXL;
        v7L = wTarget;

        v0R = state0R;
        v1R = state1R;
        v2R = state2R;
        v3R = state3R;
        v4R = state4R;
        v5R = prevXR;
    }
    break;
    case FILTER_SVFLH:
    {
        const float prevParam1 = clamp(fxParam1, 0.0f, 1.0f);
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * prevParam1) * 0.1f, 0.0f, 1.0f);
        float res = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        // MS-20 style dual ZDF SVF + Allpass phase disperser: Resonant Highpass -> Allpass -> Driven Resonant Lowpass
        const float gLP0 = 0.0025f + prevParam1 * (0.42f + 0.72f * prevParam1);
        const float gLP1 = 0.0025f + fxParam1 * (0.42f + 0.72f * fxParam1);
        const float dgLP = (gLP1 - gLP0) * 0.03125f;
        float gLP = gLP0;

        const float gHP0 = 0.0015f + prevParam1 * (0.015f + 0.045f * res);
        const float gHP1 = 0.0015f + fxParam1 * (0.015f + 0.045f * res);
        const float dgHP = (gHP1 - gHP0) * 0.03125f;
        float gHP = gHP0;

        const float kLP = 1.92f - 1.88f * res;
        const float kHP = 1.45f - 1.15f * res;

        const float f1 = clamp(0.28f + 0.42f * fxParam1, filterWindowMin, filterWindowMax);
        const float coef1 = (1.0f - f1) / (1.0f + f1);

        const float a1_lp_0 = 1.0f / (1.0f + gLP0 * (gLP0 + kLP));
        const float a1_lp_1 = 1.0f / (1.0f + gLP1 * (gLP1 + kLP));
        const float da1_lp = (a1_lp_1 - a1_lp_0) * 0.03125f;
        float a1_lp_lin = a1_lp_0;

        const float a1_hp_0 = 1.0f / (1.0f + gHP0 * (gHP0 + kHP));
        const float a1_hp_1 = 1.0f / (1.0f + gHP1 * (gHP1 + kHP));
        const float da1_hp = (a1_hp_1 - a1_hp_0) * 0.03125f;
        float a1_hp_lin = a1_hp_0;

        const float sat = 1.10f + res * 0.65f;
        const float satK = 0.10f + 0.28f * res;
        const float svfGain = (0.95f + 0.25f * res) * mixerGain;

        float *sp = this->sampleBlock;
        float hp_ic1L = v0L, hp_ic2L = v1L, lp_ic1L = v2L, lp_ic2L = v3L;
        float hp_ic1R = v0R, hp_ic2R = v1R, lp_ic1R = v2R, lp_ic2R = v3R;
        float _ly1L = v4L, _ly1R = v4R;
        float _lx1L = v5L, _lx1R = v5R;
        float prevXL = v8L, prevXR = v8R;
        float prevFL = sat25ADAA1_F(prevXL), prevFR = sat25ADAA1_F(prevXR);

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            gLP += dgLP;
            gHP += dgHP;
            a1_lp_lin += da1_lp;
            a1_hp_lin += da1_hp;
            const float a1_lp = a1_lp_lin * (2.0f - (1.0f + gLP * (gLP + kLP)) * a1_lp_lin);
            const float a1_hp = a1_hp_lin * (2.0f - (1.0f + gHP * (gHP + kHP)) * a1_hp_lin);

            // Left voice: ADAA input drive -> Resonant HP -> Allpass -> Resonant LP
            float inL = sat25ADAA((*sp) * sat, prevXL, prevFL);
            float v3 = inL - hp_ic2L;
            float v1 = a1_hp * (hp_ic1L + gHP * v3);
            float v2 = hp_ic2L + gHP * v1;
            hp_ic1L = 2.0f * v1 - hp_ic1L;
            hp_ic2L = 2.0f * v2 - hp_ic2L;
            float hpOutL = inL - kHP * v1 - v2;

            _ly1L = coef1 * (_ly1L + hpOutL) - _lx1L;
            _lx1L = hpOutL;

            v3 = _ly1L - lp_ic2L;
            v1 = a1_lp * (lp_ic1L + gLP * v3);
            v2 = lp_ic2L + gLP * v1;
            lp_ic1L = 2.0f * v1 - lp_ic1L;
            lp_ic1L -= satK * lp_ic1L * (lp_ic1L * lp_ic1L) * fastInvSqrtVoice(1.0f + lp_ic1L * lp_ic1L);
            lp_ic2L = 2.0f * v2 - lp_ic2L;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);

            // Right voice
            float inR = sat25ADAA((*sp) * sat, prevXR, prevFR);
            v3 = inR - hp_ic2R;
            v1 = a1_hp * (hp_ic1R + gHP * v3);
            v2 = hp_ic2R + gHP * v1;
            hp_ic1R = 2.0f * v1 - hp_ic1R;
            hp_ic2R = 2.0f * v2 - hp_ic2R;
            float hpOutR = inR - kHP * v1 - v2;

            _ly1R = coef1 * (_ly1R + hpOutR) - _lx1R;
            _lx1R = hpOutR;

            v3 = _ly1R - lp_ic2R;
            v1 = a1_lp * (lp_ic1R + gLP * v3);
            v2 = lp_ic2R + gLP * v1;
            lp_ic1R = 2.0f * v1 - lp_ic1R;
            lp_ic1R -= satK * lp_ic1R * (lp_ic1R * lp_ic1R) * fastInvSqrtVoice(1.0f + lp_ic1R * lp_ic1R);
            lp_ic2R = 2.0f * v2 - lp_ic2R;
            *sp++ = clamp(v2 * svfGain, -ratioTimbres, ratioTimbres);
        }

        v0L = hp_ic1L;
        v1L = hp_ic2L;
        v2L = lp_ic1L;
        v3L = lp_ic2L;
        v0R = hp_ic1R;
        v1R = hp_ic2R;
        v2R = lp_ic1R;
        v3R = lp_ic2R;
        v4L = _ly1L;
        v4R = _ly1R;
        v5L = _lx1L;
        v5R = _lx1R;
        v8L = prevXL;
        v8R = prevXR;
    }
    break;
    case FILTER_CRUSH2:
    {
        // Dual-Clock Polyrhythmic S&H Decimator + Resonant ZDF Reconstruction Filter ("Alkx": Smp1 = param1, Smp2 = param2)
        float fxParamTmp = clamp(currentTimbre->params_.effect1.param1 + matrixFilterFrequency, 0.0f, 1.0f);
        fxParamTmp *= fxParamTmp;
        fxParam1 = clamp((fxParamTmp + 9.0f * clamp(fxParam1, 0.0f, 1.0f)) * 0.1f, 0.0f, 1.0f);
        float filterParam2 = clamp(matrixFilterParam2 + currentTimbre->params_.effect1.param2, 0.0f, 1.0f);

        const float b1inc = 0.008f + fxParam1 * fxParam1 * 0.992f;
        const float b2inc = 0.008f + filterParam2 * filterParam2 * 0.94f + b1inc * 0.05f;
        const float sat = 1.25f + fxParam1 * 0.35f;

        // Resonant ZDF 2-pole DAC reconstruction filter tracking Smp1 & Smp2
        const float cutoff = 0.30f + 0.45f * fxParam1 + 0.25f * filterParam2;
        const float g = 0.025f + cutoff * (0.22f + 0.85f * cutoff);
        const float k = 1.05f - 0.35f * (1.0f - fxParam1);
        const float a1 = 1.0f / (1.0f + g * (g + k));
        const float a2 = g * a1;
        const float a3 = g * a2;

        float *sp = this->sampleBlock;
        float phase1 = v8L;
        float phase2 = v8R;

        float ic1eqL = v0L, ic2eqL = v1L;
        float ic1eqR = v0R, ic2eqR = v1R;
        float hold1L = v2L, hold1R = v2R;
        float hold2L = v3L, hold2R = v3R;
        float smoothL = v6L, smoothR = v6R;

        for (int kSample = BLOCK_SIZE; kSample--;)
        {
            float rawL = *sp;
            float rawR = *(sp + 1);

            // Branchless sample & hold (phase stays in [0, 2) since increments are < 1):
            // avoids conditional VFP moves that trigger a GCC ARM ICE with -mslow-flash-data
            phase1 += b1inc;
            const float wrap1 = (float)(int)phase1;
            phase1 -= wrap1;
            hold1L += wrap1 * (rawL - hold1L);
            hold1R += wrap1 * (rawR - hold1R);

            phase2 += b2inc;
            const float wrap2 = (float)(int)phase2;
            phase2 -= wrap2;
            const float satInL = clamp(rawL * sat, -2.0f, 2.0f);
            const float satInR = clamp(rawR * sat, -2.0f, 2.0f);
            hold2L += wrap2 * (satInL * (1.0f - 0.25f * fabsf(satInL)) - hold2L);
            hold2R += wrap2 * (satInR * (1.0f - 0.25f * fabsf(satInR)) - hold2R);

            float destL = 0.52f * hold1L + 0.48f * hold2L;
            float destR = 0.52f * hold1R + 0.48f * hold2R;
            smoothL = (smoothL * 3.0f + destL) * 0.25f;
            smoothR = (smoothR * 3.0f + destR) * 0.25f;

            // Left ZDF reconstruction filter
            float v3 = smoothL - ic2eqL;
            float v1 = a1 * ic1eqL + a2 * v3;
            float v2 = ic2eqL + a2 * ic1eqL + a3 * v3;
            ic1eqL = 2.0f * v1 - ic1eqL;
            ic2eqL = 2.0f * v2 - ic2eqL;
            *sp++ = clamp(softLimitVoice(v2) * mixerGain, -ratioTimbres, ratioTimbres);

            // Right ZDF reconstruction filter
            v3 = smoothR - ic2eqR;
            v1 = a1 * ic1eqR + a2 * v3;
            v2 = ic2eqR + a2 * ic1eqR + a3 * v3;
            ic1eqR = 2.0f * v1 - ic1eqR;
            ic2eqR = 2.0f * v2 - ic2eqR;
            *sp++ = clamp(softLimitVoice(v2) * mixerGain, -ratioTimbres, ratioTimbres);
        }

        v8L = phase1;
        v8R = phase2;
        v0L = ic1eqL;
        v1L = ic2eqL;
        v0R = ic1eqR;
        v1R = ic2eqR;
        v2L = hold1L;
        v2R = hold1R;
        v3L = hold2L;
        v3R = hold2R;
        v6L = smoothL;
        v6R = smoothR;
    }
    break;

    case FILTER_OFF:
    {
        // Filter off has gain...
        if (unlikely(mixerGain != 1.0f))
        {
            float *sp = this->sampleBlock;
            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                *sp = (*sp) * mixerGain;
                sp++;
                *sp = (*sp) * mixerGain;
                sp++;
            }
        }
    }
    break;

    default:
        // NO EFFECT
        break;
    }
}

void Voice::setNewEffectParam(int encoder)
{
    if (encoder == 0)
    {
        v0L = v1L = v2L = v3L = v4L = v5L = v6L = v7L = v8L = v0R = v1R = v2R = v3R = v4R = v5R = v6R = v7R = v8R = 0.0f;
        fxParamA1 = fxParamA2 = fxParamB2 = 0.0f;
        fxPhase = 0.0f;
        fxParam1 = clamp(currentTimbre->params_.effect1.param1, 0.0f, 1.0f);
        fxParam2 = clamp(currentTimbre->params_.effect1.param2, 0.0f, 1.0f);
        fxParamB1 = (float)currentTimbre->params_.effect1.type;

        for (int k = 1; k < NUMBER_OF_ENCODERS_PFM2; k++)
        {
            setNewEffectParam(k);
        }
    }
    switch ((int)currentTimbre->params_.effect1.type)
    {
    case FILTER_BASS:
        fxParam1 = clamp(currentTimbre->params_.effect1.param1, 0.0f, 1.0f);
        break;
    default:
        switch (encoder)
        {
        case ENCODER_EFFECT_TYPE:
            break;
        case ENCODER_EFFECT_PARAM1:
            break;
        case ENCODER_EFFECT_PARAM2:
            break;
        }
    }
}

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC pop_options
#endif
