{ pkgs ? import <nixpkgs> {} }:

pkgs.stdenv.mkDerivation {
  name = "test nix environment";
  nativeBuildInputs = with pkgs; [
    gcc
    pkg-config
    llvmPackages_19.libcxx
  ];
  buildInputs = with pkgs; [
    arpack
    mpfr
    mpc
    boost
    blas
  ];
}
