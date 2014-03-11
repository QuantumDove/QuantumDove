
// Before doing anything else, grab in the definitions of the test
// space and trial spaces that make up the finite element system.
#include <qdove/base/test_space.h>
#include <qdove/base/trial_space.h>
#include <qdove/materials/constants.h>

// Start by solving Schroedinger's problem 
#include <qdove/models/schroedinger.h>

// Let's use a function from the qdove library for a psudo potential.
#include <qdove/psuedopotentials/function_library.h>

// Also make use of the predefined generalised eigenspectrum system -
// we need this for the Schroedinger problem.
#include <qdove/generic_linear_algebra/eigenspectrum_system.h>

// Next up are some deal.II objects that have not been generalised
// away yet...
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/data_out.h>

// C++
#include <cassert>
#include <iostream>
#include <fstream>

// The purpose of this example is to solve a problem with a known
// analytical solution.
template<int dim>
class SchroedingerProblem
{
public:
  SchroedingerProblem (dealii::Triangulation<dim> &triangulation);
  ~SchroedingerProblem ();

  void write_gnuplot (const dealii::PETScWrappers::Vector &vector,
		      const std::string                   &name);
  void run ();

private:
  // The description of the finite element basis
  qdove::TrialSpace<dim> trial_space;

  // Geometry description and boundary constraints
  qdove::TestSpace<dim> test_space;

  // The eigenspectrum system that will be solved by SLEPc
  qdove::Schroedinger::Problem<dim> schroedinger_problem;

  // The eigenpairs from schroedinger's problem
  std::vector<dealii::PETScWrappers::Vector> eigenvectors;
  std::vector<double>                        eigenvalues;
};

template<int dim>
SchroedingerProblem<dim>::SchroedingerProblem (dealii::Triangulation<dim> &triangulation)
  :
  trial_space (triangulation),
  test_space (trial_space),
  schroedinger_problem (trial_space, test_space)
{}

template<int dim>
SchroedingerProblem<dim>::~SchroedingerProblem ()
{}

template<int dim>
void
SchroedingerProblem<dim>::write_gnuplot (const dealii::PETScWrappers::Vector &vector,
					 const std::string                   &name)
{
  // Output a vector to gnuplot style file. 
  std::ostringstream filename;
  filename << "solution-" << name << ".gpl";
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
SchroedingerProblem<dim>::run ()
{
  std::cout << "Test space:" << std::endl
	    << "   Finite element type:          "
	    << test_space.fe ().get_name ()
	    << std::endl
	    << "   Number of degrees of freedom: " 
	    << test_space.n_dofs ()
	    << std::endl;

  const double radius = 50e-10;

  dealii::PETScWrappers::Vector material_function (test_space.n_dofs ());
  for (unsigned int i=0; i<material_function.size (); ++i)
    material_function[i] = 1.;
  write_gnuplot (material_function, "material_function");

  // Get started on Schroedinger's problem  
  std::cout << "Schroedinger's problem:" 
	    << std::endl;
  schroedinger_problem.reinit ();

  // set up the effective mass
  dealii::PETScWrappers::Vector eff_mass (test_space.n_dofs ());
  const double ke_term = (qdove::HBAR*qdove::HBAR) / (2.*qdove::mstar_p_CdTe*qdove::M0);
  std::cout << "   Kinetic energy term:          " 
	    << ke_term << std::endl;
  eff_mass.add (ke_term, material_function);
  write_gnuplot (eff_mass, "effective_mass");

  // set up the potential
  dealii::PETScWrappers::Vector potential (test_space.n_dofs ());

  schroedinger_problem.assemble (eff_mass, potential);
  schroedinger_problem.solve ();
  schroedinger_problem.get_solution_eigenpairs (eigenvalues, eigenvectors);
  write_gnuplot (eigenvectors[0], "electron_function");

  // output
  std::cout << "   Eigenvalues:                  ";
  for (unsigned int i=0; i<eigenvalues.size (); ++i)
    std::cout << eigenvalues[i] << " ";
  std::cout << std::endl;

  std::cout << "   Scaled values (eV):           ";
  for (unsigned int i=0; i<eigenvalues.size (); ++i)
    std::cout << eigenvalues[i]/(1e-03*qdove::E0) << " ";
  std::cout << std::endl;

  std::cout << "   Scaled analytic values:       ";
  for (unsigned int i=0; i<eigenvalues.size (); ++i)
    {
      // This equation is equation 3.13 in Harrisson's book. It reads:
      // E_n = (hbar*pi*n)^2 / (2mL^2)
      const unsigned int n          = i+1;
      const double L                = 2.*radius;
      const double analytical_value = (qdove::HBAR*qdove::HBAR*qdove::PI*qdove::PI*n*n) / (2*qdove::mstar_p_CdTe*qdove::M0*L*L);
      std::cout << analytical_value/(1e-03*qdove::E0) << " ";
    }
  std::cout << std::endl;
}

int main (int argc, char **argv)
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
      {
	// Create a grid
	dealii::Triangulation<1> triangulation;
	dealii::GridGenerator::hyper_cube (triangulation, -50e-10, 50e-10);
	triangulation.refine_global (9);

	// Run Schroedinger's problem on that grid
	SchroedingerProblem<1> schroedinger_problem (triangulation);
	schroedinger_problem.run ();
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
