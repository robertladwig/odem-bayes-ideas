Robert: Some ODEM news fyi: updated code is now on Github: https://github.com/LynetteGao/Limno_DataScience/tree/dev (dev branch), and I attached the model equations plus experiments (calibration, mass-balance, NTL-LTER lake results) to this mail.

# rstan quickstart

Getting started (once per R session)
```r
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
```

# One-time configuration steps

## Once per rstan installation

Optional makevars that can result in compiled Stan programs that execute much faster than they otherwise would:
```
Makevars_file <- file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
file.edit(Makevars_file)
```
New contents of ~/.R/Makevars (was empty before):
```
CXX14FLAGS=-O3 -march=native -mtune=native
CXX14FLAGS += -arch x86_64 -ftemplate-depth-256
```

## Once per RStudio project

In Tools | Project Options | Code Editing, check the boxes for `Ensure that files end with newline` and `Strip trailing whitespace`.
