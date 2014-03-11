## Base clases.
set (src
  ## Main models
  fick
  poisson
  schroedinger
  ## Auxillary models
  statistics
  )

add_library (models OBJECT ${src})
