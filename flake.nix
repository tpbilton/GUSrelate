{
  description = "GUSrelate";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    GUSbase = {
      url = github:tpbilton/GUSbase;
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs = { self, nixpkgs, flake-utils, GUSbase }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
          };

          flakePkgs = {
            GUSbase = GUSbase.packages.${system}.default;
          };

          GUSrelate = with pkgs;
            rPackages.buildRPackage {
              name = "GUSrelate";
              src = ./.;
              propagatedBuildInputs = with rPackages;
                [ data_table, R6, Rdpack, ggplot2, plotly, flakePkgs.GUSbase ];
            };

          R-with-GUSrelate = with pkgs;
            rWrapper.override {
              packages = with rPackages;
                [ GUSrelate ];
            };
        in
          with pkgs;
          {
            devShells.default = mkShell {
              buildInputs = [ R-with-GUSrelate ];
              shellHook = ''
            mkdir -p "$(pwd)/_libs"
            export R_LIBS_USER="$(pwd)/_libs"
          '';
            };

            packages.default = GUSrelate;
          });
}
