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
int geometry_flag = 2; // 1: ideal bi-ventricle, 2: atrium, 3: slab, 4: rabbit heart, 5: ideal bi-ventricle centered at the origin
Real end_time = 500; // simulation time, unit: ms. 
// The simulation saving function has bug: if simulation time is too long (about > 2000 ms), it will create file name with negative number (PhysiologyHeart_-2147483648.vtp) and stop saving.
// And simulation results of 1500 ms long has problem: the activation waves does not look correct. Maybe because of the naming of the .vtp causing confusion, and some files are being overwritten? 
// This should be a bug of giving .vtp names.
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

std::vector<std::vector<size_t>> find_nearest_neighbors(
    const Vec3d* nodes,
    size_t n_nodes,
    size_t n_neighbors)
{
    std::vector<std::vector<size_t>> neighbor_id_2d;
    neighbor_id_2d.reserve(n_nodes);

    for (size_t i = 0; i < n_nodes; ++i)
    {
        std::vector<std::pair<Real, size_t>> dists;
        dists.reserve(n_nodes);

        // compute distance from node i to all other nodes
        for (size_t j = 0; j < n_nodes; ++j)
        {
            Vec3d diff = nodes[i] - nodes[j];
            dists.push_back({diff.squaredNorm(), j});
        }

        // sort by distance
        std::sort(dists.begin(), dists.end());

        // get nearest ids, skipping the first one (itself)
        std::vector<size_t> nearest_ids;
        nearest_ids.reserve(n_neighbors);
        for (size_t k = 1; k <= n_neighbors && k < dists.size(); ++k)
        {
            nearest_ids.push_back(dists[k].second);
        }
        neighbor_id_2d.push_back(nearest_ids);
    }

    return neighbor_id_2d;
}

int main(int ac, char *av[])
{
    int debug_display = 0;

    Real space_buffer = 5.0; // space_buffer = 0.0 works well too, it seems no need for space buffer

    if (geometry_flag == 1) {
        full_path_to_stl_file = "./input/heart-new.stl"; // sample ventricle. x = 2.54:83.53, y = 0.0:70.0, z = 2.55:62.59
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

    Vec3d *node = physiology_heart.getBaseParticles().getVariableDataByName<Vec3d>("Position");
    size_t n_particles = physiology_heart.SizeOfLoopRange();
    size_t n_neighbors = 6;
    std::vector<std::vector<size_t>> neighbor_id_2d = find_nearest_neighbors(node, n_particles, n_neighbors);

    // print for verification
    debug_display = 0;
    if (debug_display == 1) {
        std::cout << "Neighbor check (first 5 particles):\n";
        for (size_t i = 0; i < 5 && i < n_particles; ++i) {
            std::cout << "  Particle " << i << " neighbors: ";
            for (size_t neighbor_idx : neighbor_id_2d[i]) {
                std::cout << neighbor_idx << " ";
            }
            std::cout << "\n";
        }
    }

    std::vector<long unsigned int> s1_pacing_particle_id = {};
    if (geometry_flag == 1) { // ventricle
        s1_pacing_particle_id = {8262}; // assign a particle id, it will be the s1 pacing site
    } else if (geometry_flag == 2) { // atrium
        s1_pacing_particle_id = {23403};
    } else if (geometry_flag == 3) { // slab
        s1_pacing_particle_id = {1, 2, 3, 4, 5};
    } else if (geometry_flag == 4) { // rabbit heart
        s1_pacing_particle_id = {24710, 24720, 24721, 24732, 25743, 25744, 25755};
    }

    // simulation computation
    double s1_t = 0.0; // ms
    double s2_t = s1_t + 50.0; // ms
    double ap_min = 0.016783728;
    double ap_max = 0.645235435;
    double h_min = 0.431867775;
    double h_max = 2.380168937;
    double s2_region_size_factor = 0.3; // Controls how large the final region is

    Real *voltage = physiology_heart.getBaseParticles().getVariableDataByName<Real>("Voltage");
    Real *voltage_prev = physiology_heart.getBaseParticles().registerStateVariable<Real>("Voltage_prev");
    Real *h = physiology_heart.getBaseParticles().getVariableDataByName<Real>("GateVariable");
    Real *h_prev = physiology_heart.getBaseParticles().registerStateVariable<Real>("h_prev");

    int rotor_flag = 1;
    while (physical_time < end_time)
    {
        for (size_t i = 0; i < n_particles; ++i) voltage_prev[i] = voltage[i]; // assign the previous value
        for (size_t i = 0; i < n_particles; ++i) h_prev[i] = h[i]; // assign the previous value

        Real integration_time = 0.0;
        while (integration_time < Ouput_T)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observer_time)
            {
                if (ite % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = " << physical_time
                        << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
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

                    std::vector<int> s2_pacing_particle_id_auto;
                    std::sort(id1.begin(), id1.end());
                    std::sort(id2.begin(), id2.end());
                    std::set_intersection( // intersection of id1 and id2
                        id1.begin(), id1.end(),
                        id2.begin(), id2.end(),
                        std::back_inserter(s2_pacing_particle_id_auto)
                    );

                    debug_display = 0;
                    if (debug_display == 1){
                        std::cout << "s2_pacing_particle_id_auto: ";
                        for (int x : s2_pacing_particle_id_auto) std::cout << x << " ";
                        std::cout << std::endl;
                    }
                    
                    // reduce the s2 pacing region
                    std::vector<int> s2_pacing_particle_id;
                    if (!s2_pacing_particle_id_auto.empty())
                    {
                        std::set<int> id_set;
                        id_set.insert(s2_pacing_particle_id_auto[0]); // find one voxel to start

                        // For efficient lookups, create a set of the original auto-detected IDs
                        std::set<int> s2_auto_set(s2_pacing_particle_id_auto.begin(), s2_pacing_particle_id_auto.end());

                        size_t previous_size = 0;
                        while (id_set.size() < s2_pacing_particle_id_auto.size() * s2_region_size_factor)
                        {
                            previous_size = id_set.size();
                            std::vector<int> neighbors_to_add;
                            for (int particle_id : id_set)
                            {
                                if (particle_id < neighbor_id_2d.size()) {
                                    for (size_t neighbor : neighbor_id_2d[particle_id])
                                    {
                                        neighbors_to_add.push_back(neighbor);
                                    }
                                }
                            }

                            for (int neighbor_id : neighbors_to_add)
                            {
                                // Only add the neighbor if it was in the original auto-detected region
                                if (s2_auto_set.count(neighbor_id))
                                {
                                    id_set.insert(neighbor_id);
                                }
                            }
                        }
                        // Convert the final set of IDs back to a vector
                        s2_pacing_particle_id.assign(id_set.begin(), id_set.end());
                    }

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
