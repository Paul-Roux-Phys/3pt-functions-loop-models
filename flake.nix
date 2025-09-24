{
  description = "Flake for the SparseVectors library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }: {
    devShells.x86_64-darwin.default =
      let
        system = "x86_64-darwin";
        pkgs = import nixpkgs { inherit system; };
      in
        pkgs.mkShell {
          name = "SparseVectors";
          nativeBuildInputs = with pkgs; [
            gcc
            gnumake
            pkg-config
            mpfr     # Multiple Precision Floating-Point library
            libmpc   # Multiple precision complex floats
            arpack   # eigensolver
            cln
          ];
          buildInputs = with pkgs; [
            doxygen
          ];
        };
  };
}
