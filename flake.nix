{
  description = "Flake for the SparseVectors library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-darwin";
      pkgs = import nixpkgs { inherit system; };
    in {
      devShells.default = pkgs.mkShell {
        buildInputs = with pkgs; [
          gcc      # g++, gcc, etc.
          mpfr     # Multiple Precision Floating-Point library
        ];
      };
    };
}
