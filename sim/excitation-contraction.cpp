#include "all_bodies.h"
#include "all_body_relations.h"
#include "all_geometries.h"
#include "all_kernels.h"
#include "all_closures.h"
#include "all_particle_generators.h"
#include "all_particles.h"
#include "all_physical_dynamics.h"
#include "all_simbody.h"
#include "io_all.h"
#include "parameterization.h"
#include "sph_system.hpp"
using namespace SPH; 

// settings
// ------------------------------------------------------------
int geometry_flag = 1; // 1: ideal bi-ventricle, 2: atrium, 3: slab, 4: rabbit heart, 5: ideal bi-ventricle centered at the origin
Real end_time = 300; // simulation time, unit: ms
// ------------------------------------------------------------

std::string full_path_to_stl_file;
Real dp_0; // Initial particle spacing
Vec3d domain_lower_bound(3);
Vec3d domain_upper_bound(3);

// Material properties
// Unit:
// time t = ms
// length l = mm
// mass m = g
// density rho = g * (mm)^(-3)
// Pressure pa = g * (mm)^(-1) * (ms)^(-2)
// diffusion d = (mm)^(2) * (ms)^(-2)
Real rho0_s = 1.06e-3;
Real stress_scale = 1.0e-6; // Active stress factor
Real k_a = 100 * stress_scale;
std::array<Real, 4> a0 = {Real(496.0) * stress_scale, Real(15196.0) * stress_scale, Real(3283.0) * stress_scale, Real(662.0) * stress_scale};
std::array<Real, 4> b0 = {Real(7.209), Real(20.417), Real(11.176), Real(9.466)};
// reference stress to achieve weakly compressible condition 
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
// Electrophysiology parameters
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 0.8;
Real bias_coeff = 0.0;
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
Vec3d fiber_direction(1.0, 0.0, 0.0); 
Vec3d sheet_direction(0.0, 1.0, 0.0);

namespace SPH
{
// Define heart geometry
Vec3d translation(3);
class Heart : public ComplexShape
{
  public:
    explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
    {
        if (geometry_flag == 1) {
            // geometry: sample ventricle
            // Vec3d translation(-53.5, -70.0, -32.5); // original
            // NOTE: 
            // mean(geometry xyz) = [43.273098 34.759075 32.533287]
            // bounding box Vec3d domain_lower_bound(-55.0, -75.0, -35.0);
            // bounding box Vec3d domain_upper_bound(35.0, 5.0, 35.0);
            // mean(bounding box xyz) = (domain_lower_bound+domain_upper_bound)/2 = [-10.0 -35.0 0.0]
            // translation = - mean(geometry xyz) + mean(bounding box xyz), to move geometry to the center of the bounding box
            // (tried to set bounding box and geometry co-centered so no need to translate, but it did not work, thus keeping this original setup)
            translation << -43.273098-10.0, -34.759075-35.0, -32.533287+0.0;
            // translation << 0.0, 0.0, 0.0;
        } else if (geometry_flag == 2 || geometry_flag == 3 || geometry_flag == 4 || geometry_flag == 5) {
            // since the geometry is already set to centered in the bounding box, no translation is needed
            translation << 0.0, 0.0, 0.0;
        }

        Real length_scale = 1.0;
        add<TriangleMeshShapeSTL>(full_path_to_stl_file, translation, length_scale);
    }
};

// Imposing diffusion boundary condition 
class DiffusionBCs : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit DiffusionBCs(BodyPartByParticle &body_part, const std::string &species_name)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          pos_(particles_->getVariableDataByName<Vec3d>("Position")),
          phi_(particles_->registerStateVariable<Real>(species_name)) {};
    virtual ~DiffusionBCs() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Vec3d displacement = sph_body_.getInitialShape().findNormalDirection(pos_[index_i]);
        Vec3d face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vec3d center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);

        Real angle = face_norm.dot(center_norm);
        if (angle >= 0.0)
        {
            phi_[index_i] = 1.0;
        }
        else
        {
            if (pos_[index_i][1] < -sph_body_.getSPHAdaptation().ReferenceSpacing())
                phi_[index_i] = 0.0;
        }
    };

  protected:
    Vec3d *pos_;
    Real *phi_;
};

// Compute Fiber and Sheet direction after diffusion 
class ComputeFiberAndSheetDirections : public LocalDynamics
{
  protected:
    LocallyOrthotropicMuscle &muscle_material_;
    Vec3d *pos_;
    Real *phi_;
    Real beta_epi_, beta_endo_;
    Vec3d center_line_vector_; // parallel to the ventricular centerline and pointing  apex-to-base

  public:
    explicit ComputeFiberAndSheetDirections(SPHBody &sph_body, const std::string &species_name)
        : LocalDynamics(sph_body),
          muscle_material_(DynamicCast<LocallyOrthotropicMuscle>(this, sph_body_.getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vec3d>("Position")),
          phi_(particles_->registerStateVariable<Real>(species_name))
    {
        center_line_vector_ = Vec3d(0.0, 1.0, 0.0);
        beta_epi_ = -(70.0 / 180.0) * M_PI;
        beta_endo_ = (80.0 / 180.0) * M_PI;
    };
    virtual ~ComputeFiberAndSheetDirections() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        // Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
        // Present  doi.org/10.1016/j.cma.2016.05.031
        // Probe the face norm from level set field
        Vec3d displacement = sph_body_.getInitialShape().findNormalDirection(pos_[index_i]);
        Vec3d face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vec3d center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);
        if (face_norm.dot(center_norm) <= 0.0)
        {
            face_norm = -face_norm;
        }
        // Compute the centerline's projection on the plane orthogonal to face norm
        Vec3d circumferential_direction = getCrossProduct(center_line_vector_, face_norm);
        Vec3d cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
        // The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo
        Real beta = (beta_epi_ - beta_endo_) * phi_[index_i] + beta_endo_;
        // Compute the rotation matrix through Rodrigues rotation formulation
        Vec3d f_0 = cos(beta) * cd_norm + sin(beta) * getCrossProduct(face_norm, cd_norm) +
                   face_norm.dot(cd_norm) * (1.0 - cos(beta)) * face_norm;

        if (pos_[index_i][1] < -sph_body_.getSPHAdaptation().ReferenceSpacing())
        {
            muscle_material_.local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
            muscle_material_.local_s0_[index_i] = face_norm;
        }
        else
        {
            muscle_material_.local_f0_[index_i] = Vec3d::Zero();
            muscle_material_.local_s0_[index_i] = Vec3d::Zero();
        }
    };
};

// define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
  public:
    MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
    {
        Real l = domain_upper_bound[0] - domain_lower_bound[0];
        Real w = domain_upper_bound[2] - domain_lower_bound[2];
        halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
        resolution_ = 20;
        translation_ = Vec3d(-10.0, -1.0 * dp_0, 0.0);
    }
};

} // namespace SPH

int main(int ac, char *av[])
{
    Real space_buffer = 5.0; // space_buffer = 0.0 works well too, it seems no need for space buffer

    if (geometry_flag == 1) {
        full_path_to_stl_file = "./input/heart-new.stl"; // ideal bi-ventricle. x = 2.54:83.53, y = 0.0:70.0, z = 2.55:62.59
        domain_lower_bound << -55.0, -75.0, -35.0;
        domain_upper_bound << 35.0, 5.0, 35.0;
        // domain_lower_bound << 2.54-space_buffer,0.0-space_buffer, 2.55-space_buffer;
        // domain_upper_bound << 83.53+space_buffer, 70.0+space_buffer, 62.59+space_buffer;
        dp_0 = 2.0; 
    } else if (geometry_flag == 2) { // atrium
        full_path_to_stl_file = "./input/39_2-1-1-ReLA_212_thick_smooth.stl";
        domain_lower_bound << -107.2671-space_buffer,-34.7543-space_buffer, 68.1186-space_buffer;
        domain_upper_bound << 18.7329+space_buffer, 55.2457+space_buffer, 160.1186+space_buffer;
        dp_0 = 2.0;

        // full_path_to_stl_file = "./input/39_2-1-1-ReLA_212.stl";
        // domain_lower_bound << -104.0955-space_buffer,-31.5827-space_buffer, 71.2901-space_buffer;
        // domain_upper_bound << 14.9045+space_buffer, 52.4173+space_buffer, 156.2901+space_buffer;
        // dp_0 = 0.5; // for the 2*sqrt(2) thin atrium, dp_0 >= 0.5 will not work: will stuck before particle relaxation
    } else if (geometry_flag == 3) { // slab
        full_path_to_stl_file = "./input/slab.stl";
        domain_lower_bound << -35.0-space_buffer, -20.0-space_buffer, -25.0-space_buffer;
        domain_upper_bound << 35.0+space_buffer, 20.0+space_buffer, 25.0+space_buffer;
        dp_0 = 2.0;

        // NOTE: if all particles have positive coordinates, 
        // the resulted simulation will have propagating voltage and mechanical stress, 
        // but will not have mechanical movement
        // full_path_to_stl_file = "./input/slab_all_positive.stl";
        // domain_lower_bound << 10.0-space_buffer, 5.0-space_buffer, 7.0-space_buffer;
        // domain_upper_bound << 80.0+space_buffer, 45.0+space_buffer, 57.0+space_buffer;
        // dp_0 = 2.0;
    } else if (geometry_flag == 4) { // rabbit heart
        full_path_to_stl_file = "./input/2024-08-14-rec26-frame002.npy.stl";
        domain_lower_bound << -15.6976-space_buffer, -14.5304-space_buffer, -13.8727-space_buffer;
        domain_upper_bound << 16.6515+space_buffer, 14.5837+space_buffer, 8.7716+space_buffer;
        dp_0 = 0.5; // 0.5 works
                    // 2.0 is too large for the small geometry 
                    // 1.0, 1.5 will result in NaN in simulation 
                    // 0.8 some particle will be NaN
    } else if (geometry_flag == 5) { // ideal bi-ventricle centered at the origin
        full_path_to_stl_file = "./input/ideal_bi_ventricle.stl";
        domain_lower_bound << -40.714302-space_buffer, -34.77151-space_buffer, -29.979214-space_buffer;
        domain_upper_bound << 40.27483+space_buffer, 35.22849+space_buffer, 30.058975+space_buffer;
        dp_0 = 2.0;
    }
    BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); // domain bounds of the system

    //----------------------------------------------------------------------
    //	SPHSystem section
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);

    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment(); // handle command line arguments

    //----------------------------------------------------------------------
    //	SPH Particle relaxation for body fitted particle distribution.
    //----------------------------------------------------------------------
    int ite = 0;

    SolidBody herat_model(sph_system, makeShared<Heart>("HeartModel"));

    herat_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    herat_model.defineClosure<LocallyOrthotropicMuscle, IsotropicDiffusion>(
        ConstructArgs(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0),
        ConstructArgs(diffusion_species_name, diffusion_coeff));
    herat_model.generateParticles<BaseParticles, Lattice>();
    
    // Topology 
    InnerRelation herat_model_inner(herat_model);
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(herat_model);
    RelaxationStepInner relaxation_step_inner(herat_model_inner);

    // Relaxation output
    BodyStatesRecordingToVtp write_herat_model_state_to_vtp({herat_model});
    ReloadParticleIO write_particle_reload_files(herat_model);
    write_particle_reload_files.addToReload<Vec3d>(herat_model, "Fiber");
    write_particle_reload_files.addToReload<Vec3d>(herat_model, "Sheet");

    // Physics relaxation starts here.
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    write_herat_model_state_to_vtp.writeToFile(0.0);

    // From here the time stepping begins.        
    int relax_step = 1000; // original: 1000. NOTE: if the step is not enough, the simulation will result NaN after some time)
    while (ite < relax_step)
    {
        relaxation_step_inner.exec();
        ite++;
        if (ite % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
            write_herat_model_state_to_vtp.writeToFile(ite);
        }
    }

    // Diffusion process to initialize fiber direction
    GetDiffusionTimeStepSize get_time_step_size(herat_model);
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>> diffusion_relaxation_1(herat_model_inner);
    SimpleDynamics<ComputeFiberAndSheetDirections> compute_fiber_sheet(herat_model, diffusion_species_name);
    BodySurface surface_part(herat_model);
    SimpleDynamics<DiffusionBCs> impose_diffusion_bc(surface_part, diffusion_species_name);
    impose_diffusion_bc.exec();
    write_herat_model_state_to_vtp.addToWrite<Real>(herat_model, diffusion_species_name);
    write_herat_model_state_to_vtp.writeToFile(ite);

    int diffusion_step = 100;
    Real dt = get_time_step_size.exec();
    while (ite <= diffusion_step + relax_step)
    {
        diffusion_relaxation_1.exec(dt);
        impose_diffusion_bc.exec();
        if (ite % 10 == 0)
        {
            std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
            write_herat_model_state_to_vtp.writeToFile(ite);
        }
        ite++;
    }
    compute_fiber_sheet.exec();
    ite++;
    write_herat_model_state_to_vtp.writeToFile(ite);
    compute_fiber_sheet.exec();
    write_particle_reload_files.writeToFile(0);

    //----------------------------------------------------------------------
    // SPH simulation section
    //----------------------------------------------------------------------
    SolidBody mechanics_heart(sph_system, makeShared<Heart>("MechanicalHeart"));
    mechanics_heart.defineMaterial<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);

    mechanics_heart.generateParticles<BaseParticles, Reload>("HeartModel");

    SolidBody physiology_heart(sph_system, makeShared<Heart>("PhysiologyHeart"));
    AlievPanfilowModel aliev_panfilow_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
    physiology_heart.defineClosure<Solid, MonoFieldElectroPhysiology<LocalDirectionalDiffusion>>(
        Solid(), ConstructArgs(&aliev_panfilow_model, ConstructArgs(diffusion_coeff, bias_coeff, fiber_direction)));

    physiology_heart.generateParticles<BaseParticles, Reload>("HeartModel");

    //----------------------------------------------------------------------
    // SPHBody relation (topology) section
    //----------------------------------------------------------------------
    InnerRelation physiology_heart_inner(physiology_heart);
    InnerRelation mechanics_body_inner(mechanics_heart);
    ContactRelation physiology_heart_contact(physiology_heart, {&mechanics_heart});
    ContactRelation mechanics_body_contact(mechanics_heart, {&physiology_heart});

    //----------------------------------------------------------------------
    // SPH Method section
    //----------------------------------------------------------------------
    // Corrected configuration.
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_excitation(physiology_heart_inner);
    
    // Time step size calculation.
    GetDiffusionTimeStepSize get_physiology_time_step(physiology_heart);
    
    // Diffusion process for diffusion body.
    electro_physiology::ElectroPhysiologyDiffusionInnerRK2<LocalDirectionalDiffusion> diffusion_relaxation_2(physiology_heart_inner);

    // Solvers for ODE system.
    electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_heart, aliev_panfilow_model);
    electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_heart, aliev_panfilow_model);

    // Active mechanics.
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_contraction(mechanics_body_inner);
    InteractionDynamics<CorrectInterpolationKernelWeights> correct_kernel_weights_for_interpolation(mechanics_body_contact);
    
    // Interpolate the active contract stress from electrophysiology body. 
    InteractionDynamics<InterpolatingAQuantity<Real>> active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
    
    // Interpolate the particle position in physiology_heart  from mechanics_heart. 
    // TODO: this is a bug, we should interpolate displacement other than position.
    InteractionDynamics<InterpolatingAQuantity<Vec3d>> interpolation_particle_position(physiology_heart_contact, "Position", "Position");
    
    // active and passive stress relaxation. 
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(mechanics_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(mechanics_body_inner);

    // Time step size calculation. 
    ReduceDynamics<solid_dynamics::AcousticTimeStep> get_mechanics_time_step(mechanics_heart);
    
    // Constrain region of the inserted body. 
    MuscleBaseShapeParameters muscle_base_parameters;
    BodyRegionByParticle muscle_base(mechanics_heart, makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(muscle_base);

    //----------------------------------------------------------------------
    // SPH Output section
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(physiology_heart, "Voltage");
    write_states.addToWrite<Real>(physiology_heart, "GateVariable");
    write_states.addToWrite<Real>(physiology_heart, "ActiveContractionStress");
    // write_states.addToWrite<Real>(mechanics_heart, "ActiveContractionStress");

    //----------------------------------------------------------------------
    // Pre-simulation.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correct_configuration_excitation.exec();
    correct_configuration_contraction.exec();
    correct_kernel_weights_for_interpolation.exec();

    // Output initial states and observations
    write_states.writeToFile(0);
    
    //----------------------------------------------------------------------
    // Physical parameters for main loop.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int screen_output_interval = 10;
    ite = 0;
    int reaction_step = 2;
    Real Ouput_T = end_time / 200.0;
    Real Observer_time = 0.01 * Ouput_T;
    // Real dt = 0.0; // Default acoustic time step sizes for physiology
    dt = 0.0; // Default acoustic time step sizes for physiology
    Real dt_s = 0.0; // Default acoustic time step sizes for mechanics

    // Statistics for computing time. 
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    std::cout << "Main Loop Starts Here : " << "\n";
    
    std::vector<long unsigned int> s1_pacing_particle_id = {};
    s1_pacing_particle_id = {5040,5657,5659,5676,5677,5678};

    std::vector<int> s2_pacing_particle_id;
    s2_pacing_particle_id = {272,279,288,297,299,309,310,322,323,336,454,460,461,468,469,470,477,478,479,487,488,489,490,500,501,502,503,513,514,515,516,527,528,529,530,531,541,542,543,544,545,552,553,554,555,563,564,565,566,567,574,575,576,577,584,586,587,606,643,644,648,649,654,655,656,662,663,664,665,671,673,674,675,676,684,685,686,687,688,689,698,699,700,701,702,709,710,711,712,713,714,720,721,722,723,724,730,731,732,733,734,740,741,742,743,744,749,750,751,752,758,759,760,761,762,768,769,770,771,772,777,778,779,780,786,787,788,796,805,832,833,834,838,839,840,841,846,847,848,849,850,857,858,859,860,861,869,870,871,872,873,874,882,883,884,885,886,887,894,895,896,897,898,899,905,906,907,908,909,915,916,917,918,919,924,925,926,927,932,933,934,935,940,941,942,943,948,949,950,951,956,957,958,959,964,965,966,967,972,973,974,975,980,981,982,983,989,990,991,998,999,1005,1020,1021,1022,1026,1027,1028,1029,1034,1035,1036,1037,1038,1044,1045,1046,1047,1048,1049,1056,1057,1058,1059,1060,1061,1062,1069,1070,1071,1072,1073,1074,1080,1081,1082,1083,1089,1090,1091,1092,1093,1098,1099,1100,1101,1106,1107,1108,1109,1114,1115,1116,1117,1122,1123,1124,1125,1130,1131,1132,1133,1138,1139,1140,1141,1146,1147,1148,1149,1153,1154,1155,1159,1160,1161,1162,1167,1168,1169,1170,1175,1176,1177,1178,1183,1184,1185,1186,1194,1202,1207,1210,1211,1212,1216,1217,1218,1219,1220,1221,1227,1228,1229,1230,1231,1232,1238,1239,1240,1241,1242,1243,1249,1250,1251,1252,1253,1259,1260,1261,1262,1263,1268,1269,1270,1271,1276,1277,1278,1279,1284,1285,1286,1287,1292,1293,1294,1295,1300,1301,1302,1303,1308,1309,1310,1311,1316,1317,1318,1319,1323,1324,1325,1329,1330,1331,1336,1337,1338,1339,1344,1345,1346,1347,1351,1352,1353,1357,1358,1359,1363,1364,1365,1369,1370,1371,1372,1378,1379,1380,1385,1386,1387,1388,1392,1393,1396,1397,1398,1399,1400,1404,1405,1406,1407,1408,1409,1415,1416,1417,1418,1419,1420,1426,1427,1428,1429,1430,1431,1436,1437,1438,1439,1444,1445,1446,1447,1452,1453,1454,1455,1460,1461,1462,1463,1467,1468,1469,1474,1475,1476,1477,1481,1482,1483,1488,1489,1490,1491,1495,1496,1497,1501,1502,1503,1508,1509,1510,1511,1515,1516,1517,1521,1522,1523,1527,1528,1529,1534,1535,1536,1537,1542,1543,1544,1545,1549,1550,1551,1555,1556,1557,1561,1562,1563,1567,1568,1569,1573,1574,1575,1576,1578,1579,1580,1581,1585,1586,1587,1588,1589,1590,1596,1597,1598,1599,1600,1601,1607,1608,1609,1610,1611,1612,1617,1618,1619,1620,1625,1626,1627,1628,1633,1634,1635,1636,1641,1642,1643,1644,1649,1650,1651,1652,1656,1657,1658,1662,1663,1664,1665,1670,1671,1672,1673,1677,1678,1679,1684,1685,1686,1687,1691,1692,1693,1697,1698,1699,1703,1704,1705,1710,1711,1712,1713,1717,1718,1719,1723,1724,1725,1729,1730,1731,1735,1736,1737,1741,1742,1743,1748,1749,1750,1751,1756,1757,1758,1759,1762,1763,1765,1766,1767,1768,1771,1772,1773,1774,1775,1780,1781,1782,1783,1784,1790,1791,1792,1793,1794,1795,1800,1801,1802,1803,1808,1809,1810,1811,1815,1816,1817,1818,1823,1824,1825,1826,1831,1832,1833,1834,1838,1839,1840,1845,1846,1847,1848,1852,1853,1854,1858,1859,1860,1865,1866,1867,1868,1872,1873,1874,1881,1882,1883,1893,1894,1895,1896,1907,1908,1909,1920,1921,1922,1933,1934,1935,1948,1949,1950,1964,1965,1966,1967,1980,1981,1982,1995,1996,1997,2010,2011,2012,2025,2026,2027,2037,2038,2039,2040,2041,2044,2045,2046,2047,2048,2052,2053,2054,2055,2056,2057,2063,2064,2065,2066,2067,2068,2073,2074,2075,2076,2081,2082,2083,2084,2088,2089,2090,2095,2096,2097,2102,2103,2104,2105,2109,2110,2111,2116,2117,2118,2119,2127,2128,2129,2130,2141,2142,2143,2154,2155,2156,2168,2169,2170,2171,2183,2184,2185,2186,2187,2202,2203,2204,2219,2220,2221,2238,2239,2240,2257,2258,2259,2277,2278,2279,2296,2297,2298,2315,2316,2317,2335,2336,2337,2356,2357,2358,2377,2378,2379,2398,2399,2400,2416,2417,2418,2420,2421,2422,2423,2424,2428,2429,2430,2431,2432,2433,2439,2440,2441,2442,2443,2444,2450,2451,2452,2453,2454,2459,2460,2461,2462,2466,2467,2468,2473,2474,2475,2476,2481,2482,2483,2484,2485,2486,2494,2495,2496,2497,2498,2499,2509,2510,2511,2512,2513,2514,2515,2516,2525,2526,2527,2528,2529,2530,2543,2544,2545,2546,2547,2548,2560,2561,2562,2563,2564,2565,2578,2579,2580,2581,2582,2583,2584,2599,2600,2601,2602,2603,2604,2605,2621,2622,2623,2624,2625,2626,2642,2643,2644,2645,2646,2647,2664,2665,2666,2667,2668,2687,2688,2689,2690,2691,2711,2712,2713,2714,2735,2736,2737,2758,2759,2760,2781,2782,2783,2804,2805,2806,2827,2828,2829,2851,2852,2853,2854,2876,2877,2878,2879,2901,2902,2903,2904,2907,2908,2909,2910,2911,2915,2916,2917,2918,2919,2920,2926,2927,2928,2929,2930,2935,2936,2937,2938,2943,2944,2945,2946,2947,2948,2949,2956,2957,2958,2959,2960,2961,2962,2963,2964,2965,2966,2974,2975,2976,2977,2978,2979,2980,2981,2982,2983,2994,2995,2996,2997,2998,2999,3000,3001,3002,3013,3014,3015,3016,3017,3018,3019,3020,3021,3033,3034,3035,3036,3037,3038,3039,3040,3041,3055,3056,3057,3058,3059,3060,3061,3062,3077,3078,3079,3080,3081,3082,3083,3084,3101,3102,3103,3104,3105,3106,3107,3108,3125,3126,3127,3128,3129,3130,3131,3150,3151,3152,3153,3154,3155,3156,3157,3158,3176,3177,3178,3179,3180,3181,3182,3183,3200,3201,3202,3203,3204,3205,3206,3207,3227,3228,3229,3230,3231,3232,3233,3234,3255,3256,3257,3258,3259,3260,3261,3282,3283,3284,3285,3286,3287,3308,3309,3310,3311,3312,3335,3336,3337,3361,3362,3363,3364,3388,3389,3390,3391,3416,3417,3418,3419,3444,3445,3446,3447,3472,3496,3499,3500,3501,3502,3506,3507,3508,3509,3510,3511,3516,3517,3518,3519,3520,3521,3522,3527,3528,3529,3530,3531,3532,3533,3534,3535,3540,3541,3542,3543,3544,3545,3546,3547,3548,3549,3550,3556,3557,3558,3559,3560,3561,3562,3563,3564,3565,3573,3574,3575,3576,3577,3578,3579,3580,3581,3582,3593,3594,3595,3596,3597,3598,3599,3600,3601,3612,3613,3614,3615,3616,3617,3618,3619,3620,3633,3634,3635,3636,3637,3638,3639,3640,3654,3655,3656,3657,3658,3659,3660,3661,3676,3677,3678,3679,3680,3681,3682,3683,3700,3701,3702,3703,3704,3705,3706,3723,3724,3725,3726,3727,3728,3748,3749,3750,3751,3752,3753,3773,3774,3775,3776,3777,3778,3779,3780,3799,3800,3801,3802,3803,3804,3805,3806,3825,3826,3827,3828,3829,3830,3831,3832,3833,3854,3855,3856,3857,3858,3859,3860,3881,3882,3883,3884,3885,3886,3887,3909,3910,3911,3912,3913,3914,3936,3937,3938,3939,3940,3941,3965,3966,3967,3968,3969,3970,3994,3995,3996,3997,3998,3999,4024,4025,4026,4027,4028,4054,4055,4056,4084,4113,4140,4141,4145,4146,4147,4148,4152,4153,4154,4155,4156,4157,4162,4163,4164,4165,4166,4167,4168,4173,4174,4175,4176,4177,4178,4179,4180,4181,4187,4188,4189,4190,4191,4192,4193,4194,4195,4196,4203,4204,4205,4206,4207,4208,4209,4210,4211,4212,4213,4222,4223,4224,4225,4226,4227,4228,4229,4230,4231,4232,4242,4243,4244,4245,4246,4247,4248,4249,4250,4251,4264,4265,4266,4267,4268,4269,4270,4271,4272,4284,4285,4286,4287,4288,4289,4290,4291,4292,4307,4308,4309,4310,4311,4312,4313,4331,4332,4333,4334,4335,4354,4355,4356,4357,4358,4379,4380,4381,4382,4404,4405,4406,4407,4429,4430,4431,4432,4433,4457,4458,4459,4460,4461,4464,4483,4484,4485,4486,4487,4488,4489,4490,4491,4511,4512,4513,4514,4515,4516,4517,4518,4540,4541,4542,4543,4544,4545,4566,4567,4568,4569,4570,4571,4572,4573,4596,4597,4598,4599,4600,4601,4602,4626,4627,4628,4629,4630,4631,4655,4656,4657,4658,4659,4660,4685,4686,4687,4688,4716,4774,4807,4812,4813,4814,4815,4816,4820,4821,4822,4823,4824,4825,4826,4827,4832,4833,4834,4835,4836,4837,4838,4839,4840,4846,4847,4848,4849,4850,4851,4852,4853,4854,4855,4862,4863,4864,4865,4866,4867,4868,4869,4870,4871,4872,4879,4880,4881,4882,4883,4884,4885,4886,4887,4888,4889,4898,4899,4900,4901,4902,4903,4904,4905,4906,4907,4908,4920,4921,4922,4923,4924,4925,4926,4927,4940,4941,4942,4943,4944,4945,4946,4947,4948,4961,4962,4963,4964,4965,4966,4967,4968,4982,4983,4984,4985,4986,4987,4988,4989,5007,5008,5009,5010,5011,5031,5032,5033,5034,5057,5058,5059,5082,5083,5084,5085,5109,5110,5111,5112,5135,5137,5138,5139,5161,5162,5163,5164,5165,5166,5167,5182,5183,5184,5185,5186,5187,5188,5189,5205,5206,5207,5208,5209,5210,5211,5229,5230,5231,5232,5233,5234,5235,5252,5253,5254,5255,5256,5257,5275,5276,5277,5278,5279,5297,5298,5299,5300,5320,5321,5423,5429,5430,5431,5432,5433,5434,5440,5442,5443,5444,5445,5446,5447,5454,5455,5456,5457,5458,5459,5460,5461,5468,5469,5470,5471,5472,5473,5474,5475,5476,5484,5485,5486,5487,5488,5489,5490,5491,5492,5493,5502,5503,5504,5505,5506,5507,5508,5509,5510,5511,5512,5522,5523,5524,5525,5526,5527,5528,5529,5530,5531,5543,5544,5545,5546,5547,5548,5549,5550,5551,5552,5566,5567,5568,5569,5570,5571,5572,5573,5588,5589,5590,5591,5592,5593,5611,5612,5613,5614,5615,5616,5632,5633,5634,5651,5652,5653,5670,5673,5689,5691,5710,5730,5731,5732,5733,5749,5750,5751,5752,5753,5754,5755,5767,5768,5769,5770,5771,5772,5773,5787,5788,5789,5790,5791,5792,5793,5807,5808,5809,5810,5811,5812,5813,5828,5829,5830,5832,5846,5847,5849,5992,5993,5994,5995,6005,6006,6007,6008,6009,6010,6020,6021,6022,6023,6024,6025,6035,6036,6037,6038,6039,6040,6041,6042,6051,6053,6054,6055,6056,6057,6058,6059,6071,6072,6073,6074,6075,6076,6077,6078,6090,6092,6093,6094,6095,6096,6097,6110,6111,6112,6113,6114,6115,6116,6117,6118,6128,6129,6130,6131,6132,6133,6134,6135,6148,6149,6150,6151,6162,6163,6164,6165,6166,6179,6180,6181,6196,6229,6246,6264,6280,6281,6282,6283,6284,6285,6297,6298,6299,6300,6301,6302,6316,6317,6318,6319,6320,6334,6336,6543,6544,6545,6546,6563,6564,6566,6581,6583,6586,6588,6598,6600,6601,6602,6603,6604,6614,6615,6616,6617,6618,6619,6620,6621,6629,6630,6631,6632,6633,6634,6635,6644,6646,6648,6649,6650,6662,6663,6664,6675,6687,6774,6776,7084,7095,7096,7097,7124,7126};

    int rotor_flag = 1; // 1: apply s2 pacing, 0: do not apply s2 pacing

    // simulation computation
    double s1_t = 0.0; // ms
    double s2_t = s1_t + 50.0; // ms
    double ap_min = 0.001006643;
    double ap_max = 0.006296945;
    double h_min = 0.101714267;
    double h_max = 0.26044646;

    Real *voltage = physiology_heart.getBaseParticles().getVariableDataByName<Real>("Voltage");
    Real *voltage_prev = physiology_heart.getBaseParticles().registerStateVariable<Real>("Voltage_prev");
    Real *h = physiology_heart.getBaseParticles().getVariableDataByName<Real>("GateVariable");
    Real *h_prev = physiology_heart.getBaseParticles().registerStateVariable<Real>("h_prev");
    size_t n_particles = physiology_heart.SizeOfLoopRange();
    for (size_t i = 0; i < n_particles; ++i) voltage_prev[i] = voltage[i]; // initialize prev to current

    int debug_display = 1;
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Ouput_T)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observer_time)
            {
                for (size_t i = 0; i < n_particles; ++i) voltage_prev[i] = voltage[i];
                for (size_t i = 0; i < n_particles; ++i) h_prev[i] = h[i];

                if (ite % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
                              << physical_time
                              << "	dt = " << dt
                              << "	dt_s = " << dt_s << "\n";
                }
                
                // apply S1 pacing
                if (physical_time >= s1_t && physical_time <= s1_t + 0.5){
                    std::cout << "apply s1 pacing" << std::endl;

                    for (size_t k = 0; k < s1_pacing_particle_id.size(); ++k){
                        size_t pid = s1_pacing_particle_id[k];
                        voltage[pid] = 0.92;
                    }
                }

                // apply S2 pacing to induce spiral wave
                if (rotor_flag == 1 && physical_time >= s2_t &&  physical_time <= s2_t + 5){ 
                    std::cout << "apply s2 pacing" << std::endl;

                    // for (size_t k = 0; k < n_particles; ++k) {
                    //     if (position[k][0] >= 0.0 && position[k][0] <= 6.0){
                    //         if (position[k][1] >= -6.0){
                    //             if (position[k][2] >= 12.0){
                    //                 voltage[k] = 0.95;
                    //             }
                    //         }
                    //     }
                    // }

                    /*
                    std::vector<int> id1;
                    for (int i = 0; i < n_particles; ++i) {
                        if (voltage_prev[i] >= ap_min && voltage_prev[i] <= ap_max) {
                            id1.push_back(i);
                        }
                    }

                    std::vector<int> id2;
                    for (int i = 0; i < n_particles; ++i) {
                        if (h_prev[i] >= h_min && h_prev[i] <= h_max) {
                            id2.push_back(i);
                        }
                    }

                    // Find intersection of id1 and id2
                    std::vector<int> s2_pacing_particle_id;
                    std::sort(id1.begin(), id1.end());
                    std::sort(id2.begin(), id2.end());
                    std::set_intersection(
                        id1.begin(), id1.end(),
                        id2.begin(), id2.end(),
                        std::back_inserter(s2_pacing_particle_id)
                    );

                    debug_display = 0;
                    if (debug_display == 1){
                        std::cout << "Intersection: ";
                        for (int x : s2_pacing_particle_id) std::cout << x << " ";
                        std::cout << std::endl;
                    }
                    */
                    
                    for (size_t k = 0; k < s2_pacing_particle_id.size(); ++k){
                        size_t pid = s2_pacing_particle_id[k];
                        voltage[pid] = 0.92;
                    }
                }

                // Strong splitting method. 
                // forward reaction
                int ite_forward = 0;
                while (ite_forward < reaction_step){
                    reaction_relaxation_forward.exec(0.5 * dt / Real(reaction_step));
                    ite_forward++;
                }
                
                // 2nd Runge-Kutta scheme for diffusion. 
                diffusion_relaxation_2.exec(dt);

                // backward reaction
                int ite_backward = 0;
                while (ite_backward < reaction_step){
                    reaction_relaxation_backward.exec(0.5 * dt / Real(reaction_step));
                    ite_backward++;
                }

                active_stress_interpolation.exec();

                Real dt_s_sum = 0.0;
                while (dt_s_sum < dt){
                    dt_s = get_mechanics_time_step.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    stress_relaxation_first_half.exec(dt_s);
                    constraint_holder.exec(dt_s);
                    stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                }

                ite++;
                dt = get_physiology_time_step.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
        }

        TickCount t2 = TickCount::now();
        interpolation_particle_position.exec();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
