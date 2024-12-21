{
  description = "Nix environment as a flake";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs = { self, nixpkgs }: {
    devShells.default = 
      let
        pkgs = import nixpkgs { system = "x86_64-darwin"; };
      in
      pkgs.mkShell {
        name = "TransferMatrices";

        # Tools available in the shell environment
        nativeBuildInputs = with pkgs; [
          gcc
          pkg-config
        ];

        # Libraries or other runtime dependencies
        buildInputs = with pkgs; [
          arpack
          mpfr
          mpc
          boost
          blas
        ];
      };
  };
}

