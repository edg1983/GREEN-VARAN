# GREEN-VARAN
# Author: Edoardo Giacopuzzi
# Annotate VCF files with regulatory annotations from GREEN-DB

import argparse
import os
import tables
import strformat
import ./greenvaran/smallvars
import ./greenvaran/version_logo
import ./greenvaran/sv
import ./greenvaran/make_querytab

proc main*() =
  type pair = object
    f: proc(dropfirst:bool)
    description: string

  var dispatcher = {
    "smallvars": pair(f:smallvars.main, description:"Annotate small variants VCF"),
    "sv": pair(f:sv.main, description:"Annotate structural variants VCF"),
    "querytab": pair(f:make_querytab.main, description: "Generate a table for greendb_query from VCF")
    }.toOrderedTable

  var args = commandLineParams()
  
  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""
  
  if args[0] == "version":
      quit fmt"GREEN-VARAN version {VERSION}", QuitSuccess

  echo LOGO
  dispatcher[args[0]].f(false)

when isMainModule:
    main()