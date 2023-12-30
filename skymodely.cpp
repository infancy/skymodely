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

struct spectrum_t
{
    // TODO

    void to_rgb()
    {
    }
};


enum class skymodel_enum_t
{
    none,
    arhosek12_rgb,
    arhosek12_xyz,
    arhosek12_spectrum_sky,
    arhosek12_spectrum_sun,
    arhosek12_spectrum_skysun,
};

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
        for (int i = 0; i < k_channels; i++)
        {
            arhosekskymodelstate_free(skymodel_state_[i]);
        }
    }

    color_t radiance(double theta, double gamma) override
    {
        color_t color;

        for (int i = 0; i < k_channels; i++)
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
    ArHosekSkyModelState* skymodel_state_[k_channels]{};
};

std::unique_ptr<skymodel_t> create_skymodel(
    skymodel_enum_t skymodel_enum,
    double turbidity,
    color_t albedo,
    double solar_elevation)
{
    switch (skymodel_enum)
    {
    case skymodel_enum_t::arhosek12_rgb:
    case skymodel_enum_t::arhosek12_xyz:
        return std::make_unique<arhosek12_tristim_t>(skymodel_enum, turbidity, albedo, solar_elevation);
    default:
        std::print("ERROR! no skymodel");
        return nullptr;
    }
}


// e.g. skymodely filename width height turbidity(1~10) albedo(RGB) solar_elevation(degree, 0~180) skymodel_enum
int main(int argc, char **argv)
{
    using namespace std;

    string filename = "skymodely.hdr";
    int width = 512;
    int height = 512;
    double turbidity = 1; // 1~10
    color_t albedo{ 0, 0, 0 };
    double solar_elevation = 90.0 / 180.0 * k_pi; // 0~180
    skymodel_enum_t skymodel_enum = skymodel_enum_t::arhosek12_xyz;

    if (argc == 10)
    {
        filename = argv[1];
        width = atof(argv[2]);
        height = atof(argv[3]);
        turbidity = atof(argv[4]);
        albedo = color_t{ stof(argv[5]), stof(argv[6]), stof(argv[7]) };
        solar_elevation = atof(argv[8]) / 180.0 * k_pi;
        skymodel_enum = skymodel_enum_t(atoi(argv[9]));
    }

    auto skymodel = create_skymodel(skymodel_enum, turbidity, albedo, solar_elevation);

    unique_ptr<color_t[]> img = make_unique<color_t[]>(sizeof(color_t) * width * height);

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