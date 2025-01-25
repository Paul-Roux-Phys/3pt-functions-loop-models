{
  description = "Transfer Matrices Development Environment";

  # Inputs for the flake
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  # Outputs for the flake
  outputs = { self, nixpkgs }: {
    # Define the devShell for the system architecture you are working with
    devShells.x86_64-darwin.default =
      let
        pkgs = import nixpkgs { system = "x86_64-darwin"; };
      in
        pkgs.mkShell {
	        name = "TransferMatrices";

          # Development tools available in the shell environment
	        nativeBuildInputs = with pkgs; [
	          gcc
	          arpack
	          mpfr
	          libmpc
	          boost
	          blas
	        ];
        };
  };
}
