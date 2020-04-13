//
// Created by walksbynight on 23/3/18.
//

#include "material.hpp"

real MatCosTheta(const Vector3r &w) {
    return w.values.z;
}

real MatCos2Theta(const Vector3r &w) {
    return w.values.z * w.values.z;
}

real MatAbsCosTheta(const Vector3r &w) {
    return std::abs(w.values.z);
}

real MatSin2Theta(const Vector3r &w) {
    return std::max((real) 0., (real) 1. - MatCos2Theta(w));
}

real MatSinTheta(const Vector3r &w) {
    return std::sqrt(MatSin2Theta(w));
}

real MatTanTheta(const Vector3r &w) {
    return MatSinTheta(w) / MatCosTheta(w);
}

real MatTan2Theta(const Vector3r &w) {
    return MatSin2Theta(w) / MatCos2Theta(w);
}

real MatCosPhi(const Vector3r &w) {
    real sin_theta = MatSinTheta(w);
    return (sin_theta == 0) ? (real) 1. : Clamp(w.values.x / sin_theta, (real) -1., (real) 1.);
}

real MatSinPhi(const Vector3r &w) {
    real sin_theta = MatSinTheta(w);
    return (sin_theta == 0) ? (real) 0. : Clamp(w.values.y / sin_theta, (real) -1., (real) 1.);
}

real MatCos2Phi(const Vector3r &w) {
    return MatCosPhi(w) * MatCosPhi(w);
}

real MatSin2Phi(const Vector3r &w) {
    return MatSinPhi(w) * MatSinPhi(w);
}

real MatCosDPhi(const Vector3r &wa, const Vector3r &wb) {
    return Clamp((wa.values.x * wb.values.x + wa.values.y * wb.values.y) /
                 std::sqrt((wa.values.x * wa.values.x + wa.values.y * wa.values.y) *
                           (wb.values.x * wb.values.x + wb.values.y * wb.values.y)), (real) -1., (real) 1.);
}

real FrDielectric(real cos_theta_i, real eta_i, real eta_t) {
    real clamped_cos_theta_i = Clamp(cos_theta_i, (real) -1., (real) 1.);

    bool entering = cos_theta_i > (real) 0.;
    if (!entering) {
        std::swap(eta_i, eta_t);
    }

    real sin_theta_i = std::sqrt(std::max((real) 0., 1 - clamped_cos_theta_i * clamped_cos_theta_i));
    real sin_theta_t = eta_i / eta_t * sin_theta_i;

    if (sin_theta_t >= (real) 1.) {
        return (real) 1.;
    }

    real cos_theta_t = std::sqrt(std::max((real) 0., (real) 1. - sin_theta_t * sin_theta_t));

    real r_par = ((eta_t * clamped_cos_theta_i) - (eta_i * cos_theta_t)) /
                 ((eta_t * clamped_cos_theta_i) + (eta_i * cos_theta_t));

    real r_per = ((eta_i * clamped_cos_theta_i) - (eta_t * cos_theta_t)) /
                 ((eta_i * clamped_cos_theta_i) + (eta_t * cos_theta_t));

    return (r_par * r_par + r_per * r_per) / (real) 2.;
}

bool MatSameHemisphere(const Vector3r &w, const Vector3r &wp) {
    return w.values.z * wp.values.z > 0;
}

Vector3r MatReflect(const Vector3r &wo, const Vector3r &n) {
    return wo * -1 + n * 2 * Dot(wo, n);
}

bool MatRefract(const Vector3r &wi, const Vector3r &n, real eta, Vector3r *wt) {
    real cos_theta_i = Dot(n, wi);
    real sin_2_theta_i = std::max((real) 0., (real) 1. - cos_theta_i * cos_theta_i);
    real sin_2_theta_t = eta * eta * sin_2_theta_i;

    if (sin_2_theta_t >= (real) 1.)
        return false;

    real cos_theta_t = std::sqrt((real) 1. - sin_2_theta_t);

    *wt = wi * -1 * eta + n * (eta * cos_theta_i - cos_theta_t);
    return true;
}

Spectrum Fresnel::Evaluate(real cos_i) {
    return Spectrum((real) 1.);
}

FresnelDielectric::FresnelDielectric(real eta_i, real eta_t)
        : eta_i(eta_i), eta_t(eta_t) {
}

Spectrum FresnelDielectric::Evaluate(real cos_theta_i) {
    return Spectrum(FrDielectric(cos_theta_i, eta_i, eta_t));
}


real MicrofacetDistribution::Pdf(const Vector3r &wo, const Vector3r &wh) {
    return D(wh) * MatAbsCosTheta(wh);
}

real MicrofacetDistribution::G1(const Vector3r &w) {
    return (real) 1. / ((real) 1. + Lambda(w));
}

real MicrofacetDistribution::G(const Vector3r &wo, const Vector3r &wi) {
    return 1 / (1 + Lambda(wo) + Lambda(wi));
}


BeckmannDistribution::BeckmannDistribution(real alpha) : alpha(alpha) {
}

real BeckmannDistribution::D(const Vector3r &wh) {
    real cos_2_theta = MatCos2Theta(wh);
    real cos_4_theta = cos_2_theta * cos_2_theta;
    real alpha_2 = alpha * alpha;

    real exp_term = std::exp((cos_2_theta - 1) / (alpha_2 * cos_2_theta));
    return exp_term / (3.141592f * alpha_2 * cos_4_theta);
}

real BeckmannDistribution::Lambda(const Vector3r &w) {
    real abs_tan_theta = std::abs(MatTanTheta(w));

    if (std::isinf(abs_tan_theta)) {
        return (real) 0.;
    }

    real a = 1 / (alpha * abs_tan_theta);

    if (a >= (real) 1.6) {
        return (real) 0.;
    } else {
        return (1 - 1.259f * a + 0.396f * a * a) /
               (3.535f * a + 2.181f * a * a);
    }
}

Vector3r BeckmannDistribution::SampleWh(const Vector3r &wo, const Point2r &u) {
    real tan_2_theta, phi;

    real log_sample = std::log(u.data[0]);
    if (std::isinf(log_sample)) log_sample = 0;

    tan_2_theta = -1 * alpha * alpha * log_sample;
    phi = 2 * u.data[1] * 3.1415926f;

    real cos_theta = 1 / std::sqrt(1 + tan_2_theta);
    real sin_theta = std::sqrt(std::max((real) 0., 1 - cos_theta * cos_theta));

    Vector3r wh = SphericalDirection(sin_theta, cos_theta, phi);
    if (!MatSameHemisphere(wo, wh)) wh = wh * -1;

    return wh;
}


// TrowbridgeReitzDistribution::TrowbridgeReitzDistribution(real alpha) : alpha(alpha) {
// }

// real TrowbridgeReitzDistribution::D(const Vector3r& wh) {
//     real tan_2_theta = MatTan2Theta(wh);

//     if (std::isinf(tan_2_theta)) {
//         return (real)0.;
//     }
//     real cos_4_theta = MatCos2Theta(wh) * MatCos2Theta(wh);
//     real e = tan_2_theta / (alpha * alpha);
//     return 1 / (3.1415926 * alpha * alpha * cos_4_theta * (1 + e) * (1 + e));
// }

// real TrowbridgeReitzDistribution::Lambda(const Vector3r& w) {
//     real abs_tan_theta = std::abs(MatTanTheta(w));

//     if (std::isinf(abs_tan_theta)) {
//         return (real)0.;
//     }

//     real sq = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
//     return (-1 + std::sqrt(1 + sq)) / 2;
// }

BxDF::BxDF(BxDFType type) : type(type) {
}

bool BxDF::MatchesFlags(BxDFType t) {
    return (type & t) == t;
}


LambertianReflection::LambertianReflection(const Spectrum &R) : BxDF(kBsdfReflection), R(R) {
}

Spectrum LambertianReflection::F(const Vector3r &wo, const Vector3r &wi) {
    return R * ((real) 1. / 3.141592f);
}

Point2r UnifromSampleDisk(const Point2r &u) {
    real r = std::sqrt(u.data[0]);
    real theta = 2 * 3.1415926f * u.data[1];

    return {r * std::cos(theta), r * std::sin(theta)};
}

Vector3r CosineSampleHemisphere(const Point2r &u) {
    Point2r d = UnifromSampleDisk(u);
    real z = std::sqrt(std::max((real) 0., 1 - d.values.x * d.values.y - d.values.y * d.values.y));
    return {d.values.x, d.values.y, z};
}

Vector3r UnifromSampleHemisphere(const Point2r &u) {
    real z = u.data[0];
    real r = std::sqrt(std::max((real)0, (real)1. - z * z));
    real phi = 2 * 3.141592 * u.data[1];
    return {r * std::cos(phi), r * std::sin(phi), z};
}

Spectrum LambertianReflection::SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) {
    //Implement cosine weighted sampling
    *wi = UnifromSampleHemisphere(u);
    if (wo.values.z < 0) wi->values.z *= -1;

    *pdf = Pdf(wo, *wi);

    return F(wo, *wi);
}

real LambertianReflection::Pdf(const Vector3r &wo, const Vector3r &wi) {
    return MatSameHemisphere(wo, wi) ? 1 / (2 * 3.1415926f) : (real) 0.;//MatSameHemisphere(wo, wi) ? MatAbsCosTheta(wi) / 3.1415926f : (real) 0.;
}

LambertianReflection::~LambertianReflection() {
}


MicrofacetReflection::MicrofacetReflection(const Spectrum &R, MicrofacetDistribution *distribution, Fresnel *fresnel)
        : BxDF(kBsdfReflection), R(R), distribution(distribution), fresnel(fresnel) {
}

Spectrum MicrofacetReflection::F(const Vector3r &wo, const Vector3r &wi) {
    real cos_theta_o = MatAbsCosTheta(wo), cos_theta_i = MatAbsCosTheta(wi);
    Vector3r wh = wi + wo;

    if (cos_theta_i == 0 || cos_theta_o == 0) return Spectrum((real) 0.);
    if (wh.values.x == 0 && wh.values.y == 0 && wh.values.z == 0) return Spectrum((real) 0.);

    wh = Normalize(wh);

    Spectrum fres = fresnel->Evaluate(Dot(wi, wh));
    Spectrum f = R * distribution->D(wh) * distribution->G(wo, wi) * fres *
                 Spectrum(4 * cos_theta_i * cos_theta_o).MulInverse();

    return f;
}

Spectrum
MicrofacetReflection::SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) {
    Vector3r wh = distribution->SampleWh(wo, u);
    *wi = MatReflect(wo, wh);

    if (!MatSameHemisphere(wo, *wi)) return Spectrum((real) 0.);

    *pdf = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));

    return F(wo, *wi);
}

real MicrofacetReflection::Pdf(const Vector3r &wo, const Vector3r &wi) {
    if (!MatSameHemisphere(wo, wi)) return (real) 0.;

    Vector3r wh = Normalize(wo + wi);
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
}

MicrofacetReflection::~MicrofacetReflection() {
    delete distribution;
    delete fresnel;
}


MicrofacetTransmission::MicrofacetTransmission(real eta_a, real eta_b, const Spectrum &R,
                                               MicrofacetDistribution *distribution, Fresnel *fresnel)
        : BxDF(kBsdfTransmission), eta_a(eta_a), eta_b(eta_b), R(R), distribution(distribution), fresnel(fresnel) {

}

Spectrum MicrofacetTransmission::F(const Vector3r &wo, const Vector3r &wi) {
    if (MatSameHemisphere(wo, wi)) return Spectrum(0);

    real cos_theta_o = MatCosTheta(wo), cos_theta_i = MatCosTheta(wi);
    if (cos_theta_i == 0 || cos_theta_o == 0) return Spectrum((real) 0.);

    real eta = (cos_theta_o > 0) ? (eta_b / eta_a) : (eta_a / eta_b);

    Vector3r wh = Normalize(wo + wi * eta);

    //if (wh.values.x == 0 && wh.values.y == 0 && wh.values.z == 0) return Spectrum((real)0.);
    if (!MatSameHemisphere(wo, wh)) wh = wh * -1;

    Spectrum fres = fresnel->Evaluate(Dot(wo, wh));

    real sqrt_denom = Dot(wo, wh) + eta * Dot(wi, wh);
    real factor = (1 / eta);

    Spectrum f =  (Spectrum(1) - fres) * R *
           std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
                    std::abs(Dot(wi, wh)) * std::abs(Dot(wo, wh)) * factor * factor /
                    (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom));
    return f;
}

Spectrum
MicrofacetTransmission::SampleF(const Vector3r &wo, Vector3r *wi, const Point2r &u, real *pdf) {
    if (wo.values.z == 0) return Spectrum(0);
    Vector3r wh = distribution->SampleWh(wo, u);

    real eta = MatCosTheta(wo) > 0 ? (eta_a / eta_b) : (eta_b / eta_a);

    if (!MatRefract(wo, wh, eta, wi)) {
        return Spectrum(0.);
    }


    // if (MatSameHemisphere(wo, *wi))
    //     return Spectrum(0.);

    *pdf = Pdf(wo, *wi);

    return F(wo, *wi);
}

real MicrofacetTransmission::Pdf(const Vector3r &wo, const Vector3r &wi) {
    if (MatSameHemisphere(wo, wi)) return 0;

    real eta = (MatCosTheta(wo) > 0) ? (eta_b / eta_a) : (eta_a / eta_b);

    Vector3r wh = Normalize(wo + wi * eta);
    if (!MatSameHemisphere(wo, wh)) wh = wh * -1;

    real sqrt_denominator = Dot(wo, wh) + Dot(wi, wh) * eta;
    real change_of_vars = std::abs((eta * eta * Dot(wi, wh)) / (sqrt_denominator * sqrt_denominator));
    return distribution->Pdf(wo, wh) * change_of_vars;
}

MicrofacetTransmission::~MicrofacetTransmission() {
    delete distribution;
    delete fresnel;
}


AshikhminShirleyReflection::AshikhminShirleyReflection(const Spectrum &r_d, const Spectrum &r_s,
                                                       MicrofacetDistribution *distribution)
        : BxDF(kBsdfReflection), r_d(r_d), r_s(r_s), distribution(distribution) {
}

Spectrum AshikhminShirleyReflection::SchlickFresnel(real cos_theta) {
    auto pow5 = [](real v) { return (v * v) * (v * v) * v; };
    return r_s + (Spectrum((real) 1.) - r_s) * pow5(1 - cos_theta);
}

Spectrum AshikhminShirleyReflection::F(const Vector3r &wo, const Vector3r &wi) {
    auto pow5 = [](real v) { return (v * v) * (v * v) * v; };

    Spectrum diffuse = (Spectrum((real) 1.) - r_s) *
                       (28.f / (23.f * 3.1415926)) * r_d *
                       (1 - pow5((real) 1. - (real) 0.5 * MatAbsCosTheta(wi))) *
                       (1 - pow5((real) 1. - (real) 0.5 * MatAbsCosTheta(wo)));

    Vector3r wh = wi + wo;
    if (wh.values.x == 0 && wh.values.y == 0 && wh.values.z == 0) return Spectrum((real) 0.);
    wh = Normalize(wh);

    Spectrum specular = SchlickFresnel(Dot(wh, wi)) * (distribution->D(wh) /
                                                       (4 * std::abs(Dot(wh, wi)) *
                                                        std::max(MatAbsCosTheta(wo), MatAbsCosTheta(wi))));

    return diffuse + specular;
}


BSDF::BSDF(const Vector3r &normal, real eta)
        : ns(normal), eta(eta) {
    ns.LocalCoordinateSystem(&ss, &ts);
}

Vector3r BSDF::WorldToLocal(const Vector3r &v) {
    return {Dot(v, ss), Dot(v, ts), Dot(v, ns)};
}

Vector3r BSDF::LocalToWorld(const Vector3r &v) {
    return {ss.values.x * v.values.x + ts.values.x * v.values.y + ns.values.x * v.values.z,
            ss.values.y * v.values.x + ts.values.y * v.values.y + ns.values.y * v.values.z,
            ss.values.z * v.values.x + ts.values.z * v.values.y + ns.values.z * v.values.z};
}

Spectrum BSDF::F(const Vector3r &wo_world, const Vector3r &wi_world) {
    Vector3r wo = WorldToLocal(wo_world);
    Vector3r wi = WorldToLocal(wi_world);

    Spectrum f(0);

    for (int i = 0; i < num_bxdfs; ++i) {
        f = f + bxdfs[i]->F(wo, wi);
    }

    return f;
}

Spectrum BSDF::SampleF(const Vector3r &wo_world, Vector3r *wi_world, const Point2r &u, real *pdf) {
    Vector3r wo = WorldToLocal(wo_world);
    Vector3r wi;

    Spectrum f;

    f = bxdfs[chosen_index]->SampleF(wo, &wi, u, pdf);
    auto chosen_type = bxdfs[chosen_index]->type;

    *wi_world = LocalToWorld(wi);

    int num_matching = 1;
    bool reflect = Dot(*wi_world, ns) * Dot(wo_world, ns) > 0;

    for (int i = 0; i < num_bxdfs; ++i) {
        if (i != chosen_index && ((reflect && bxdfs[i]->type == kBsdfReflection) || (!reflect && bxdfs[i]->type == kBsdfTransmission))) {
            *pdf += bxdfs[i]->Pdf(wo, wi);
            f = f + bxdfs[i]->F(wo, wi);
            num_matching++;
        }
    }
    *pdf /= num_matching;

    return f;
}

void BSDF::AddBxDF(BxDF *bxdf) {
    bxdfs[num_bxdfs++] = bxdf;
}

void BSDF::Choose(const Point2r &u) {
    if (num_bxdfs == 1) {
        chosen_index = 0;
    }
    else if (num_bxdfs == 2) {
        if (u.data[0] < 0.5) {
            chosen_index = 0;
        }
        else {
            chosen_index = 1;
        }
    }
    else {
        auto step = (real) 1. / num_bxdfs;
        for (int i = 0; i < num_bxdfs; ++i) {
            if ((i * step) <= u.data[0] && (i + 1) * step > u.data[0]) {
                chosen_index = i;
                break;
            }
        }
    }
}

BSDF::~BSDF() {
    for (int i = 0; i < num_bxdfs; ++i) {
        delete bxdfs[i];
    }
}


BSDF *Material::GetBSDF(const Interaction &interaction, real eta) {
    BSDF *bsdf = new BSDF(Normalize(Vector3r(interaction.normal)) * -1, eta);

    if (normal_map != nullptr) {
        Spectrum n = normal_map->Evaluate(interaction);
        Vector3r normal = Normalize((Vector3r(n.r, n.g, n.b) - Vector3r((real)0.5)) * (real)2.);

        Vector3r normal_tr =
                bsdf->ns + (bsdf->ts * -normal.values.x + bsdf->ss * normal.values.y +
                                    bsdf->ns * normal.values.z);

        bsdf = new (bsdf) BSDF(normal_tr, eta);
    }

    if (TestFlag(material_type, kMaterialDiffuse)) {
        LambertianReflection *bxdf = nullptr;

        if (color_map != nullptr) {
            bxdf = new LambertianReflection(diffuse_color * color_map->Evaluate(interaction));
        }
        else {
            bxdf = new LambertianReflection(diffuse_color);
        }

        bsdf->AddBxDF(bxdf);
    }

    if (TestFlag(material_type, kMaterialSpecular)) {
        MicrofacetReflection *bxdf = nullptr;
        BeckmannDistribution *distro = new BeckmannDistribution(alpha);
        Fresnel *fresnel = TestFlag(material_type, kMaterialTransparent) ? new FresnelDielectric(eta, ior) : new Fresnel;

        if (color_map != nullptr) {
            bxdf = new MicrofacetReflection (specular_color * color_map->Evaluate(interaction), distro, fresnel);
        }
        else {
            bxdf = new MicrofacetReflection (specular_color, distro, fresnel);
        }

        bsdf->AddBxDF(bxdf);
    }

    if (TestFlag(material_type, kMaterialTransparent)) {
        MicrofacetTransmission *bxdf = nullptr;
        BeckmannDistribution *distro = new BeckmannDistribution(alpha);
        Fresnel *fresnel = new FresnelDielectric(eta, ior);

        if (color_map != nullptr) {
            bxdf = new MicrofacetTransmission (eta, ior, transmission_color * color_map->Evaluate(interaction), distro, fresnel);
        }
        else {
            bxdf = new MicrofacetTransmission (eta, ior, transmission_color, distro, fresnel);
        }

        bsdf->AddBxDF(bxdf);
    }

    return bsdf;
}

bool TestFlag(unsigned int expression, unsigned int flag) {
    return (expression & flag) == flag;
}

/*
NormalMappedBSDF::NormalMappedBSDF(const Vector3r &normal, const Vector3r &map_normal, float eta) : BSDF(normal, eta) {
    Vector3r nnm = map_normal, snm, tnm;
    nnm.LocalCoordinateSystem(&snm, &tnm);
    ns = Normalize(nnm + ns);
    ss = Normalize(snm + ss);
    ts = Normalize(tnm + ts);
}
*/
