# .init.R
# Functions to initialize this Exercise session
# Boris Steipe
# ====================================================================

if (! file.exists("myScript.R")) {
    file.copy("tmp.R", "myScript.R")
}

file.edit("R_EDA-HypothesisTesting.R")

# [End]
