## Generic linear algebra clases.
set (src   
    eigenspectrum_system
    generic_eigenspectrum_solver
    generic_linear_algebra_solver
    linear_algebra_system
  )

add_library (generic_linear_algebra OBJECT ${src})
