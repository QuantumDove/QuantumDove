
#include <cassert>

// Before doing anything else, grab in the definitions of the test
// space and trial spaces that make up the finite element system.
#include <qdove/base/test_space.h>
#include <qdove/base/trial_space.h>
#include <qdove/materials/constants.h>

// Use the *solution* to a fick equation as a material function
#include <qdove/models/fick.h>

// Start by solving Schroedinger's problem
#include <qdove/models/schroedinger.h>

// and next by solving Poissons's problem
#include <qdove/models/poisson.h>

// Let's use a function from the qdove library for a psudo potential.
#include <qdove/psuedopotentials/function_library.h>

// Also make use of the predefined generalised eigenspectrum system -
// we need this for the Schroedinger problem.
#include <qdove/generic_linear_algebra/eigenspectrum_system.h>

// Also make use of the predefined linear algebra system - we need
// this for the Poisson problem.
#include <qdove/generic_linear_algebra/linear_algebra_system.h>

// Next up are some deal.II objects that have not been generalised
// away yet...
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/data_out.h>

// The purpose of this example is to have a working (dimensionless)
// Schroedinger-Poisson solver that can be generalised piecewise to
// the qdove library.
//
template<int dim>
class SelfConsistentProblem
{
public:
  SelfConsistentProblem (dealii::Triangulation<dim> &triangulation);
  ~SelfConsistentProblem ();

  void write_gnuplot (const dealii::PETScWrappers::Vector &vector,
          const std::string                   &name,
          const unsigned int                   cycle);
  void run ();

private:
  // The description of the finite element basis
  qdove::TrialSpace<dim> trial_space;

  // Geometry description and boundary constraints
  qdove::TestSpace<dim> test_space;

  // Fick's problem, which is used to obtain a material profile.
  qdove::Fick::Solution<1> fick_solution;

  // The eigenspectrum system that will be solved by SLEPc
  qdove::Schroedinger::Problem<dim> schroedinger_problem;

  // The system of equations that will be solved by PETSc
  qdove::Poisson::Problem<dim> poisson_problem;

  // The eigenpairs from schroedinger's problem
  std::vector<dealii::PETScWrappers::Vector> eigenvectors;
  std::vector<double>                        eigenvalues;

  // The solution from poisson's problem
  dealii::PETScWrappers::Vector solution;

};

template<int dim>
SelfConsistentProblem<dim>::SelfConsistentProblem (dealii::Triangulation<dim> &triangulation)
  :
  trial_space (triangulation),
  test_space (trial_space),
  fick_solution (trial_space, test_space),
  schroedinger_problem (trial_space, test_space, 10),
  poisson_problem (trial_space, test_space)
{}

template<int dim>
SelfConsistentProblem<dim>::~SelfConsistentProblem ()
{}

template<int dim>
void
SelfConsistentProblem<dim>::write_gnuplot (const dealii::PETScWrappers::Vector &vector,
             const std::string                   &name,
             const unsigned int                   cycle)
{
  // Output a vector to gnuplot style file.
  std::ostringstream filename;
  filename << "solution-" << name << "-" << cycle << ".gpl";
  std::ofstream output (filename.str ().c_str ());

  dealii::DataOut<dim> data_out;
  data_out.attach_dof_handler (test_space.dofs ());
  data_out.add_data_vector (vector, name);

  // generate default patches and output.
  data_out.build_patches ();
  data_out.write_gnuplot (output);
}


template<int dim>
void
SelfConsistentProblem<dim>::run ()
{
  std::cout << "Test space:" << std::endl
      << "   Finite element type:          "
      << test_space.fe ().get_name ()
      << std::endl
      << "   Number of degrees of freedom: "
      << test_space.n_dofs ()
      << std::endl;

  const double energy_height = 0.5 * qdove::E0;
  const double doping = 1e25;

  //
  // Setup a material function
  //
  fick_solution.reinit ();
  fick_solution.set_initial_height (1.);
  fick_solution.set_initial_rate(2e-10);
  fick_solution.set_symmetric_profile(false);

  dealii::PETScWrappers::Vector transition_1 (test_space.n_dofs ());
  fick_solution.set_initial_length (50e-10);
  fick_solution.interpolate_analytic_solution (transition_1);

  dealii::PETScWrappers::Vector transition_2 (test_space.n_dofs ());
  fick_solution.set_initial_length (100e-10);
  fick_solution.interpolate_analytic_solution (transition_2);

  dealii::PETScWrappers::Vector transition_3 (test_space.n_dofs ());
  fick_solution.set_initial_length (150e-10);
  fick_solution.interpolate_analytic_solution (transition_3);


  dealii::PETScWrappers::Vector material_function (test_space.n_dofs ());
  for (std::size_t i=0; i<material_function.size(); ++i)
    material_function[i] =  (1.0 - transition_1[i])
                           + transition_2[i] * 0.20/0.35  // for Al_{0.20}Ga_{0.80}As
                           + transition_3[i] * 0.15/0.35;      // for Al_{0.35}Ga_{0.65}As
  write_gnuplot (material_function, "material_function", 0);


  // set up the effective mass
  dealii::PETScWrappers::Vector kinetic_energy_prefactor (test_space.n_dofs ());
  for (std::size_t i=0; i<kinetic_energy_prefactor.size(); ++i)
    kinetic_energy_prefactor[i] = (qdove::HBAR*qdove::HBAR) / (2.*(qdove::mstar_GaAs + 0.35*material_function[i]*0.083)*qdove::M0);
  write_gnuplot (kinetic_energy_prefactor, "kinetic_energy_prefactor", 0);

  // set up the initial guess for the electron potential (i.e. electrostatic potential multiplied by elementary charge)
  dealii::PETScWrappers::Vector potential (test_space.n_dofs ());
  potential = material_function;
  potential *= 0.9 * energy_height; //start with no space charge region by selecting the potential right below the the Fer

  double E_bonding = 0.005 * qdove::E0; //5meV bonding energy of the ions
  double E_fermi = energy_height - E_bonding;

  //
  // Main iteration Schroedinger <-> Poisson
  //
  for (std::size_t cycle = 0; cycle < 100; ++cycle)
  {
    write_gnuplot (potential, "potential", cycle);

    // Get started on Schroedinger's problem
    std::cout << "Schroedinger's problem:" << std::endl;
    schroedinger_problem.reinit ();

    schroedinger_problem.assemble (kinetic_energy_prefactor, potential);
    schroedinger_problem.solve ();
    schroedinger_problem.get_solution_eigenpairs (eigenvalues, eigenvectors);
    write_gnuplot (eigenvectors[0], "electron_function", cycle);

    // output
    std::cout << "   Eigenvalues:                  ";
    for (unsigned int i=0; i<eigenvalues.size (); ++i)
      std::cout << eigenvalues[i] << " ";
    std::cout << std::endl;

    std::cout << "   Scaled values (meV):           ";
    for (unsigned int i=0; i<eigenvalues.size (); ++i)
      std::cout << eigenvalues[i]/(1e-03*qdove::E0) << " ";
    std::cout << std::endl;

    // Compute density from the eigenfunctions and eigenvalues:
    double kBT = 300.0 * qdove::KB;
    dealii::PETScWrappers::Vector density(test_space.n_dofs());
    for (unsigned int i=0; i<eigenvectors.size (); ++i)
    {
      double E_wavefunction = eigenvalues[i];
      double occupancy      =  kBT * std::log(1.0+std::exp((E_fermi - E_wavefunction)/kBT)); //integrated Fermi-Dirac statistic
      for (unsigned int j=0; j<eigenvectors[0].size (); ++j)
      {
        double density_of_states = (qdove::mstar_GaAs + 0.35*material_function[i]*0.083) * qdove::M0 / (qdove::HBAR*qdove::HBAR*qdove::PI);
        density[j] += density_of_states * eigenvectors[i][j] * eigenvectors[i][j] * occupancy;
      }
    }
    write_gnuplot (density, "density", cycle);

    // Set up RHS for Poisson's equation
    dealii::PETScWrappers::Vector rho(density.size());
    for (std::size_t i=0; i<solution.size(); ++i)
      if (potential[i] > E_fermi)
        rho[i] = -doping;
    rho += density;
    rho *= qdove::E0 * qdove::E0 / (qdove::permittivity_GaAs);
    write_gnuplot (rho, "rho", cycle);

    // Solve Poisson:
    std::cout << "Poisson's problem:" << std::endl;
    poisson_problem.reinit ();
    poisson_problem.assemble (rho);
    poisson_problem.solve ();
    poisson_problem.get_solution_vector (solution);
    write_gnuplot (solution, "poisson_solution", cycle);

    // Update potential:
    double max_update = 0;
    for (std::size_t i=0; i<solution.size(); ++i)
    {
      double old_value = potential[i] - 0.9 * energy_height * material_function[i];
      max_update = std::max(max_update, std::fabs(solution[i] - old_value));
    }

    double cutoff_update_value = 0.02 * qdove::E0;
    double alpha = (max_update > cutoff_update_value) ? cutoff_update_value / max_update : 1.0; //update by at most 10 mV for stability reasons

    alpha = 0.2;
    std::cout << "alpha: " << alpha << std::endl;

    for (std::size_t i=0; i<potential.size(); ++i)
    {
      double old_value = potential[i] - 0.9 * energy_height * material_function[i];
      double update = solution[i] - old_value;
      potential[i] += alpha * update;
    }
  }
}

int main (int argc, char **argv)
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
      {
  // Create a grid
  dealii::Triangulation<1> triangulation;
  dealii::GridGenerator::hyper_cube (triangulation, 0, 200e-10);
  triangulation.refine_global (8);

  // Run Schroedinger's problem on that grid
  SelfConsistentProblem<1> self_consistent_problem (triangulation);
  self_consistent_problem.run ();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
