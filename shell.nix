{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  name = "cshto-dev-shell";

  buildInputs = [
    pkgs.python3
    pkgs.python3Packages.pip
    pkgs.python3Packages.pyyaml
    pkgs.python3Packages.pytest
    pkgs.python3Packages.setuptools
  ];

  shellHook = ''
    # Add the src directory to PYTHONPATH so the modules can be found.
    export PYTHONPATH=./:$PYTHONPATH
    echo "Development shell for cshto loaded. PYTHONPATH is set to: $PYTHONPATH"
  '';
}