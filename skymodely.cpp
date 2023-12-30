#include <array>
#include <memory>
#include <numbers>
#include <print>
#include <string>
#include <vector>

#include "Arhosek12/ArHosekSkyModel.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


constexpr float k_pi = std::numbers::pi;

struct color_t
{
    union
    {
        struct { float r, g, b; };
        struct { float x, y, z; };
    };

    float& operator[](int i) { return ((float*)this)[i]; }

    color_t xyz_to_srgb()
    {
        // https://github.com/mmp/pbrt-v3/blob/master/src/core/spectrum.h#L56-L60

        color_t rgb{};
        color_t& xyz = *this;

        rgb[0] =  3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
        rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
        rgb[2] =  0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];

        return rgb;
    }
};

template<int wavelength_start_, int samples_, int delta_>
class base_spectrum_t
{
public:
    static constexpr int k_wavelength_start{ wavelength_start_ };
    static constexpr int k_wavelength_end{ wavelength_start_ + samples_ * delta_ };

    static constexpr int k_samples{ samples_ };
    static constexpr int k_delta{ delta_ };

public:
    float& operator[](int i) { return a_[i]; }
    float wavelength(int i) { return k_wavelength_start + i * delta_; }

    color_t to_srgb()
    {
        color_t color{};

        for (int i = 0; i < k_samples; ++i)
        {
            auto& spectrum = a_;
            float wavelen = wavelength(i);

            color.x += spectrum[i] * xFit_1931(wavelen) * k_delta;
            color.y += spectrum[i] * yFit_1931(wavelen) * k_delta;
            color.z += spectrum[i] * zFit_1931(wavelen) * k_delta;
        }

        return color.xyz_to_srgb();
    }

private:
    std::array<float, k_samples> a_{};

private:
    // Simple Analytic Approximationsto the CIE XYZ Color Matching Functions (JCGT) https://jcgt.org/published/0002/02/01/
    float xFit_1931(float wavelen)
    {
        float t1 = (wavelen - 442.0f) * ((wavelen < 442.0f) ? 0.0624f : 0.0374f);
        float t2 = (wavelen - 599.8f) * ((wavelen < 599.8f) ? 0.0264f : 0.0323f);
        float t3 = (wavelen - 501.1f) * ((wavelen < 501.1f) ? 0.0490f : 0.0382f);
        return 0.362f * expf(-0.5f * t1 * t1) + 1.056f * expf(-0.5f * t2 * t2) - 0.065f * expf(-0.5f * t3 * t3);
    }

    float yFit_1931(float wavelen)
    {
        float t1 = (wavelen - 568.8f) * ((wavelen < 568.8f) ? 0.0213f : 0.0247f);
        float t2 = (wavelen - 530.9f) * ((wavelen < 530.9f) ? 0.0613f : 0.0322f);
        return 0.821f * exp(-0.5f * t1 * t1) + 0.286f * expf(-0.5f * t2 * t2);
    }

    float zFit_1931(float wavelen)
    {
        float t1 = (wavelen - 437.0f) * ((wavelen < 437.0f) ? 0.0845f : 0.0278f);
        float t2 = (wavelen - 459.0f) * ((wavelen < 459.0f) ? 0.0385f : 0.0725f);
        return 1.217f * exp(-0.5f * t1 * t1) + 0.681f * expf(-0.5f * t2 * t2);
    }
};

using spectrum300_t = base_spectrum_t<400, 300,  1>;
using spectrum60_t  = base_spectrum_t<400,  60,  5>;
using spectrum30_t  = base_spectrum_t<400,  30, 10>;
using spectrum10_t  = base_spectrum_t<320,  10, 40>;

enum class spectrum_samples_enum_t
{
    n300, // 400nm~700nm, 300 sample,  1nm per sample
    n60,  // 400nm~700nm,  60 sample,  5nm per sample
    n30,  // 400nm~700nm,  30 sample, 10nm per sample
    n10,  // 320nm~720nm,  10 sample, 40nm per sample
    count
};


enum class skymodel_enum_t
{
    none,
    arhosek12_rgb = 1,
    arhosek12_xyz = 2,
    arhosek12_spectrum_sky = 4,
    arhosek12_spectrum_sun = 8,
    arhosek12_spectrum_skysun = arhosek12_spectrum_sky | arhosek12_spectrum_sun,
};
bool enum_have(skymodel_enum_t skymodel_enum, skymodel_enum_t target) { return (int)skymodel_enum & (int)target; }

class skymodel_t
{
public:
    virtual ~skymodel_t() {}

    virtual color_t radiance(double theta, double gamma) = 0;
};

class arhosek12_tristim_t : public skymodel_t
{
public:
    arhosek12_tristim_t(
        skymodel_enum_t skymodel_enum,
        double turbidity,
        color_t albedo,
        double solar_elevation)
    {
        skymodel_enum_ = skymodel_enum;

        auto skymodelstate_alloc_init = 
            skymodel_enum == skymodel_enum_t::arhosek12_rgb ?
            arhosek_rgb_skymodelstate_alloc_init:
            arhosek_xyz_skymodelstate_alloc_init;

        for (int i = 0; i < k_channels; ++i)
        {
            skymodel_state_[i] = skymodelstate_alloc_init(
                    turbidity,
                    albedo[i],
                    solar_elevation);
        }
    }

    ~arhosek12_tristim_t()
    {
        for (int i = 0; i < k_channels; ++i)
        {
            arhosekskymodelstate_free(skymodel_state_[i]);
        }
    }

    color_t radiance(double theta, double gamma) override
    {
        color_t color{};

        for (int i = 0; i < k_channels; ++i)
        {
            color[i] = arhosek_tristim_skymodel_radiance(
                skymodel_state_[i],
                theta,
                gamma,
                i);
        }

        if (skymodel_enum_ == skymodel_enum_t::arhosek12_xyz)
        {
            color = color.xyz_to_srgb();
        }

        return color;
    }

private:
    static constexpr int k_channels = 3;
    skymodel_enum_t skymodel_enum_{};
    std::array<ArHosekSkyModelState*, k_channels> skymodel_state_{};
};

template<typename spectrum_t>
//    requires(  )
class arhosek12_spectrum_t : public skymodel_t
{
public:
    arhosek12_spectrum_t(
        skymodel_enum_t skymodel_enum,
        double turbidity,
        color_t albedo,
        double solar_elevation)
    {
        skymodel_enum_ = skymodel_enum;

        // TODO: rgb to spectrum
        spectrum_t spectrum_albedo{};

        for (int i = 0; i < k_channels; ++i)
        {
            skymodel_state_[i] = arhosekskymodelstate_alloc_init(
                solar_elevation,
                turbidity,
                spectrum_albedo[i]);
        }
    }

    ~arhosek12_spectrum_t()
    {
        for (int i = 0; i < k_channels; ++i)
        {
            arhosekskymodelstate_free(skymodel_state_[i]);
        }
    }

    color_t radiance(double theta, double gamma) override
    {
        spectrum_t spectrum{};

        for (int i = 0; i < k_channels; ++i)
        {
            double wavelength = spectrum.wavelength(i);

            if (enum_have(skymodel_enum_, skymodel_enum_t::arhosek12_spectrum_sky))
                spectrum[i] += arhosekskymodel_radiance(skymodel_state_[i], theta, gamma, wavelength);

            if (enum_have(skymodel_enum_, skymodel_enum_t::arhosek12_spectrum_sun))
                spectrum[i] += arhosekskymodel_solar_radiance(skymodel_state_[i], theta, gamma, wavelength);
        }

        return spectrum.to_srgb();
    }

private:
    static constexpr int k_channels{ spectrum_t::k_samples };
    skymodel_enum_t skymodel_enum_{};
    std::array<ArHosekSkyModelState*, k_channels> skymodel_state_{};
};

std::unique_ptr<skymodel_t> create_skymodel(
    skymodel_enum_t skymodel_enum,
    spectrum_samples_enum_t spectrum_samples_eum,
    double turbidity,
    color_t albedo,
    double solar_elevation)
{
    switch (skymodel_enum)
    {
    case skymodel_enum_t::arhosek12_rgb:
    case skymodel_enum_t::arhosek12_xyz:
        return std::make_unique<arhosek12_tristim_t>(skymodel_enum, turbidity, albedo, solar_elevation);
    case skymodel_enum_t::arhosek12_spectrum_sky:
    case skymodel_enum_t::arhosek12_spectrum_sun:
    case skymodel_enum_t::arhosek12_spectrum_skysun:
    {
        switch (spectrum_samples_eum)
        {
        case spectrum_samples_enum_t::n300:
            return std::make_unique<arhosek12_spectrum_t<spectrum300_t>>(skymodel_enum, turbidity, albedo, solar_elevation);
        case spectrum_samples_enum_t::n60:
            return std::make_unique<arhosek12_spectrum_t<spectrum60_t>>(skymodel_enum, turbidity, albedo, solar_elevation);
        case spectrum_samples_enum_t::n30:
            return std::make_unique<arhosek12_spectrum_t<spectrum30_t>>(skymodel_enum, turbidity, albedo, solar_elevation);
        }
    }
    default:
        std::print("ERROR! no skymodel");
        return nullptr;
    }
}


// e.g. skymodely filename width height turbidity(1~10) albedo(RGB) solar_elevation(degree, 0~180) skymodel_enum spectrum_samples_enum
int main(int argc, char **argv)
{
    using namespace std;

    string filename = "skymodely.hdr";
    int width = 512;
    int height = 512;
    double turbidity = 1; // 1~10
    color_t albedo{ 0, 0, 0 };
    double solar_elevation = 90.0 / 180.0 * k_pi; // 0~180
    skymodel_enum_t skymodel_enum = skymodel_enum_t::arhosek12_spectrum_sky;
    spectrum_samples_enum_t spectrum_samples_eum = spectrum_samples_enum_t::n30;

    if (argc == 11)
    {
        filename = argv[1];
        width = atof(argv[2]);
        height = atof(argv[3]);
        turbidity = atof(argv[4]);
        albedo = color_t{ stof(argv[5]), stof(argv[6]), stof(argv[7]) };
        solar_elevation = atof(argv[8]) / 180.0 * k_pi;
        skymodel_enum = skymodel_enum_t(atoi(argv[9]));
        spectrum_samples_eum = spectrum_samples_enum_t(atoi(argv[10]));
    }

    auto skymodel = create_skymodel(skymodel_enum, spectrum_samples_eum, turbidity, albedo, solar_elevation);

    unique_ptr<color_t[]> img = make_unique<color_t[]>(width * height);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int cx = width / 2;
            int cy = height / 2;
            int rx = (x - cx);
            int ry = (y - cy);

            double nr = sqrt(rx * rx + ry * ry) / (width / 2.0);
            double th = nr * 0.5 * k_pi;
            double ph = atan2(rx, ry);

            double gamma = acos(cos(solar_elevation) * sin(th) * sin(ph) + sin(solar_elevation) * cos(th));
            double theta = th;

            if (nr < width / 2.0)
            {
                img[y * width + x] = skymodel->radiance(theta, gamma);
            }
        }
    }

    // TODO: Tone Mapping

    // TODO: output jpeg
    stbi_write_hdr(filename.c_str(), width, height, 3, (float*)img.get());
    //system("tex skymodely.hdr");

    return 0;
}