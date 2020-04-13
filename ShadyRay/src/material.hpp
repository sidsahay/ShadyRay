//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_MATERIAL_HPP
#define SHADYRAY_MATERIAL_HPP

#pragma once

#include "linalg.hpp"
#include "spectrum.hpp"
#include "interaction.hpp"
#include "texture.hpp"

#include <vector>

Vector3r UnifromSampleHemisphere(const Point2r &u);
Vector3r CosineSampleHemisphere(const Point2r &u);

real MatCosTheta(const Vector3r &w);

real MatCos2Theta(const Vector3r &w);

real MatabsCosTheta(const Vector3r &w);

real MatSin2Theta(const Vector3r &w);

real MatSinTheta(const Vector3r &w);

real MatTanTheta(const Vector3r &w);

real MatTan2Theta(const Vector3r &w);

real MatCosPhi(const Vector3r &w);

real MatSinPhi(const Vector3r &w);

real MatCos2Phi(const Vector3r &w);

real MatSin2Phi(const Vector3r &w);

real MatCosDPhi(const Vector3r &wa, const Vector3r &wb);

real FrDielectric(real cos_theta_i, real eta_i, real eta_t);

bool MatSameHemisphere(const Vector3r &w, const Vector3r &wp);

Vector3r MatReflect(const Vector3r &wo, const Vector3r &n);

bool MatRefract(const Vector3r &wi, const Vector3r &n, real eta, Vector3r *wt);


class Fresnel {
public:
    virtual Spectrum Evaluate(real cos_i);
};

class FresnelDielectric : public Fresnel {
public:
    FresnelDielectric(real eta_i, real eta_t);

    Spectrum Evaluate(real cos_i);

protected:
    real eta_i;
    real eta_t;
};

enum BxDFType {
    kBsdfReflection = 1 << 0,
    kBsdfTransmission = 1 << 1,
    kBsdfDiffuse = 1 << 2,
    kBsdfGlossy = 1 << 3,
    kBsdfSpecular = 1 << 4,
    kBsdfAll = kBsdfReflection | kBsdfTransmission | kBsdfDiffuse | kBsdfGlossy | kBsdfSpecular
};


class MicrofacetDistribution {
public:
    virtual real D(const Vector3r &wh) = 0;

    virtual real Lambda(const Vector3r &w) = 0;

    virtual Vector3r SampleWh(const Vector3r &wo, const Point2r &u) = 0;

    real Pdf(const Vector3r &wo, const Vector3r &wh);

    real G1(const Vector3r &w);

    real G(const Vector3r &w0, const Vector3r &wi);

    virtual ~MicrofacetDistribution() = default;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
    explicit BeckmannDistribution(real alpha);

    real D(const Vector3r &wh) override;

    real Lambda(const Vector3r &w) override;

    Vector3r SampleWh(const Vector3r &wo, const Point2r &u) override;

private:
    real alpha;
};

// class TrowbridgeReitzDistribution : public MicrofacetDistribution {
// public:
//     TrowbridgeReitzDistribution(real alpha);
//     real D(const Vector3r& wh);
//     real Lambda(const Vector3r& w);

// private:
//     real alpha;
// };


class BxDF {
public:
    explicit BxDF(BxDFType type);

    bool MatchesFlags(BxDFType t);

    virtual Spectrum F(const Vector3r &wo, const Vector3r &wi) = 0;

    virtual Spectrum
    SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) = 0;

    virtual real Pdf(const Vector3r &wo, const Vector3r &wi) = 0;
    //virtual Spectrum Rho(const Vector3r& wo, int n_samples, const Point2r* samples) = 0;

    virtual ~BxDF() = default;

    BxDFType type;
};

class LambertianReflection : public BxDF {
public:
    LambertianReflection(const Spectrum &R);

    Spectrum F(const Vector3r &wo, const Vector3r &wi) override;

    Spectrum
    SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) override;

    real Pdf(const Vector3r &wo, const Vector3r &wi) override;

    ~LambertianReflection() override;

    Spectrum R;
};

class MicrofacetReflection : public BxDF {
public:
    MicrofacetReflection(const Spectrum &R, MicrofacetDistribution *distribution, Fresnel *fresnel);

    Spectrum F(const Vector3r &wo, const Vector3r &wi) override;

    Spectrum
    SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) override;

    real Pdf(const Vector3r &wo, const Vector3r &wi) override;

    ~MicrofacetReflection() override;

    Spectrum R;
    MicrofacetDistribution *distribution;
    Fresnel *fresnel;
};

class MicrofacetTransmission : public BxDF {
public:
    MicrofacetTransmission(real eta_a, real eta_b, const Spectrum &R, MicrofacetDistribution *distribution,
                           Fresnel *fresnel);

    Spectrum F(const Vector3r &wo, const Vector3r &wi) override;

    Spectrum
    SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) override;

    real Pdf(const Vector3r &wo, const Vector3r &wi) override;

    ~MicrofacetTransmission() override;

    //Spectrum Rho(const Vector3r& wo, int n_samples, const Point2r* samples) = 0;

protected:
    real eta_a, eta_b;
    Spectrum R;
    MicrofacetDistribution *distribution;
    Fresnel *fresnel;
};

class AshikhminShirleyReflection : public BxDF {
public:
    AshikhminShirleyReflection(const Spectrum &r_d, const Spectrum &r_s, MicrofacetDistribution *distribution);

    Spectrum SchlickFresnel(real cos_theta);

    Spectrum F(const Vector3r &wo, const Vector3r &wi);

protected:
    Spectrum r_d, r_s;
    MicrofacetDistribution *distribution;
};

static constexpr unsigned int kMaxBxDFs = 8;

class BSDF {
public:
    explicit BSDF(const Vector3r &normal, real eta = 1);

    Vector3r WorldToLocal(const Vector3r &v);

    Vector3r LocalToWorld(const Vector3r &v);

    void AddBxDF(BxDF *bxdf);

    void Choose(const Point2r &u);

    virtual Spectrum F(const Vector3r &wo_world, const Vector3r &wi_world);

    virtual Spectrum SampleF(const Vector3r &wo_world, Vector3r *wi_world, const Point2r &u, real *pdf);

    real eta;

    BxDF *bxdfs[kMaxBxDFs];

    int num_bxdfs = 0;
    int chosen_index = 0;

    Vector3r ns, ss, ts;

    ~BSDF();

protected:
};

typedef unsigned int MaterialType;

static constexpr MaterialType kMaterialNone = 1 << 0;
static constexpr MaterialType kMaterialDiffuse = 1 << 1;
static constexpr MaterialType kMaterialSpecular = 1 << 2;
static constexpr MaterialType kMaterialEmission = 1 << 3;
static constexpr MaterialType kMaterialTransparent = 1 << 4;

bool TestFlag(unsigned int expression, unsigned int flag);

class Material {
public:
    Material() = default;

    virtual BSDF *GetBSDF(const Interaction &interaction, real eta);

    ~Material() = default;

    MaterialType material_type = 0;

    Texture *color_map = nullptr; // for diffuse color
    Texture *normal_map = nullptr;

    Spectrum diffuse_color;
    Spectrum specular_color;
    Spectrum emission_color;
    Spectrum transmission_color;

    real alpha = (real) 1.;
    real eta = (real) 1.;
    real ior = (real)1.;

    bool is_lamp = false;
};


#endif //SHADYRAY_MATERIAL_HPP