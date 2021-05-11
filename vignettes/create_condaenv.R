
# Create conda environment with pyscenic
reticulate::conda_create(
  envname = "pyscenic",
  packages = c("numpy"),
  python_version = "3.7"
)

# indicate that we want to use a specific condaenv
reticulate::use_condaenv("pyscenic")

# Install additional packages
reticulate::conda_install(
  envname = "pyscenic",
  packages = "cytoolz",
  channel = "anaconda",
)

# Install pyscenic  with pip
reticulate::conda_install(
  envname = "pyscenic",
  packages = "pyscenic",
  pip = TRUE
)

