# Package
version     = "1.3.1"
author      = "Edoardo Giacopuzzi"
description = "GREEN-VARAN:Annotate variants in regulatory regions from GREEN-DB"
license     = "MIT"

# Deps
requires "nim >= 0.10.0", "hts >= 0.3.4"
requires "argparse 0.10.1"

srcDir = "src"
bin = @["greenvaran"]

skipDirs = @["test"]
