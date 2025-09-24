#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

#include "BigFloats.hpp"
#include "On3pt.hpp"

using namespace OnLoops;

unsigned int L;
int k1 = 0, k2 = 0, k3 = 0;
const unsigned int prec = 1024;
using Weight = mp::mpfloat<prec>;
Weight lambda;
Parameters<Weight> params;

OnState<Weight> v(1e3), buf(1e3);

void set_bottom_state(LinkPattern<Weight> &bot) {
  for (int i = 0; i < params.k1; i++) {
    bot.set_site(i, bot_defect_labels[i]);
  }
}

void compute_diagram(int M) {
  LinkPattern<Weight> bottom_state(0);
  set_bottom_state(bottom_state);

  v.add(bottom_state);

  for (int i = 0; i < M; i++) {
    v.transfer(buf, params);
  }

  assert(k2 == 2 || k2 == 1 || k2 == 0);
  if (params.k2 == 2 || params.k2 == 1)
    v.insert_mid_op(buf, params);

  for (int i = 0; i < M; i++) {
    v.transfer(buf, params);
  }
}

LinkPattern<Weight> top_state(int sigma3) {
  LinkPattern<Weight> top(0);
  int contracted = (k1 + k2 - k3) / 2;
  int nb_mid_def = k2 - contracted;
  int nb_bot_def = k1 - contracted;
  params.sigma1 *= -1;
  for (int i = 0; i < nb_mid_def; i++) {
    top.set_site(i, PERMUTE_MID(mid_defect_labels[contracted + i], params));
  }
  for (int i = 0; i < nb_bot_def; i++) {
    top.set_site(i + nb_mid_def,
                 PERMUTE_BOT(bot_defect_labels[contracted + i], params));
  }
  for (int i = 0; i < sigma3; i++) {
    int last = top[k3 - 1];
    for (int j = k3 - 1; j > 0; j--) {
      top.set_site(j, top[j - 1]);
    }
    top.set_site(0, last);
  }
  params.sigma1 *= -1;
  return top;
}

std::vector<Weight> compute_3pt_diagrams(int M, std::ostream &out) {
  std::vector<Weight> res(k1 * k2 * k3);
  LinkPattern<Weight> top;

  for (int sigma1 = 0; sigma1 < k1; sigma1++) {
    for (int sigma2 = 0; sigma2 < k2; sigma2++) {
      params = Parameters<Weight>(lambda, k1, k2, k3, sigma1, sigma2);
      compute_diagram(M);
      if (sigma1 == 0 && sigma2 == 0)
        std::cout << "Space dimension = " << v.size() << std::endl;
      for (int sigma3 = 0; sigma3 < k3; sigma3++) {
        top = top_state(sigma3);
        res[sigma1 * k2 * k3 + sigma2 * k3 + sigma3] = v[top];
        out << res[sigma1 * k2 * k3 + sigma2 * k3 + sigma3] << ",";
      }
      v.clear();
    }
  }
  return res;
}

std::vector<Weight> compute_2pt_diagrams(int M, std::ostream &out) {
  params = Parameters<Weight>(lambda, k1, k2, k3, 0, 0);
  std::vector<Weight> res(k1);
  LinkPattern<Weight> top;

  compute_diagram(M);
  std::cout << "Dimension = " << v.size() << std::endl;

  for (int sigma = 0; sigma < k1; sigma++) {
    top = top_state(sigma);
    res[sigma] = v[top];
    out << res[sigma] << ",";
  }

  v.clear();
  return res;
}

std::vector<Weight> compute_2pt_diagrams_mid(int M, std::ostream &out) {
  std::vector<Weight> res(k1);
  LinkPattern<Weight> top = top_state(0);;

  for (int sigma = 0; sigma < k1; sigma++) {
    params = Parameters<Weight>(lambda, k1, k2, k3, 0, sigma);
    compute_diagram(M);
    if (sigma == 0)
      std::cout << "Dimension = " << v.size() << std::endl;
    res[sigma] = v[top];
    out << res[sigma] << ",";
    v.clear();
  }

  return res;
}

std::vector<Weight> compute_partition_function(int M, std::ostream &out) {
  params = Parameters(lambda, 0, 0, 0, 0, 0);
  compute_diagram(M);
  LinkPattern<Weight> top = top_state(0);
  std::vector<Weight> res(1);
  res[0] = v[top];
  out << res[0];
  v.clear();
  return res;
}

std::vector<Weight> compute_all_diagrams(int M, std::ostream &out) {
  out << std::fixed << std::setprecision(8);
  if (k1 == 0 && k2 == 0)
    return compute_partition_function(M, out);
  else if (k2 == 0)
    return compute_2pt_diagrams(M, out);
  else if (k1 == 0 || k3 == 0)
    return compute_2pt_diagrams_mid(M, out);
  else
    return compute_3pt_diagrams(M, out);
  out << std::endl;
}

mpc::mpcomplex<prec> compute_3pt_function(std::vector<Weight> vals, int rs1,
                                          int rs2, int rs3) {
  mpc::mpcomplex<prec> res = 0;
  mpc::mpcomplex<prec> im(0, 1);
  static const Weight pi = mp::pi<prec>();
  for (int sigma1 = 0; sigma1 < k1; sigma1++)
    for (int sigma2 = 0; sigma2 < k2; sigma2++)
      for (int sigma3 = 0; sigma3 < k3; sigma3++) {
        res += exp(2 * im * pi *
                   (rs1 * sigma1 / (Weight)k1 + rs2 * sigma2 / (Weight)k2 +
                    rs3 * sigma3 / (Weight)k3)) *
	  vals[sigma1 * k2 * k3 + sigma2 * k3 + sigma3];
      }
  return res;
}

std::vector<mpc::mpcomplex<prec>>
compute_3pt_functions(std::vector<Weight> vals, std::ostream &out) {
  std::vector<mpc::mpcomplex<prec>> res(k1 * k2 * k3);
  for (int rs1 = 0; rs1 < k1; rs1++)
    for (int rs2 = 0; rs2 < k2; rs2++)
      for (int rs3 = 0; rs3 < k3; rs3++) {
	res[rs1 * k2 * k3 + rs2 * k3 + rs3] = compute_3pt_function(vals, rs1, rs2, rs3);
	out << res[rs1 * k2 * k3 + rs2 * k3 + rs3] << ",";
      }
  return res;
}

mpc::mpcomplex<prec> compute_2pt_function(std::vector<Weight> vals, int rs) {
  mpc::mpcomplex<prec> res = 0;
  mpc::mpcomplex<prec> im(0, 1);
  static const Weight pi = mp::pi<prec>();
  for (int sigma = 0; sigma < k1; sigma++) {
    res += exp(2 * im * pi * (rs * sigma / (Weight)k1)) * vals[sigma];
  }
  return res;
}

std::vector<mpc::mpcomplex<prec>>
compute_2pt_functions(std::vector<Weight> vals, std::ostream &out) {
  std::vector<mpc::mpcomplex<prec>> res(k1);
  for (int rs = 0; rs < k1; rs++) {
    res[rs] = compute_2pt_function(vals, rs);
    out << res[rs] << ",";
  }
  return res;
}

std::vector<mpc::mpcomplex<prec>>
compute_2_3pt_functions(std::vector<Weight> vals, std::ostream &out) {
  if (k1 == 0 || k2 == 0 || k3 == 0)
    return compute_2pt_functions(vals, out);
  return compute_3pt_functions(vals, out);
}

void write_csv_header(std::ostream &out) {
  out << "L,lambda,";
  if (k1 == 0 && k2 == 0) {
    out << "Z";
  } else if (k1 == 0 || k2 == 0 || k3 == 0) {
    for (int sigma = 0; sigma < k1; sigma++)
      out << "d[" << sigma << "],";
    for (int sigma = 0; sigma < k1; sigma++)
      out << "Z[" << sigma << "],";
  } else {
    for (int sigma1 = 0; sigma1 < k1; sigma1++)
      for (int sigma2 = 0; sigma2 < k2; sigma2++)
        for (int sigma3 = 0; sigma3 < k3; sigma3++) {
          out << "d[" << sigma1 << sigma2 << sigma3 << "],";
    }
    for (int sigma1 = 0; sigma1 < k1; sigma1++)
      for (int sigma2 = 0; sigma2 < k2; sigma2++)
        for (int sigma3 = 0; sigma3 < k3; sigma3++) {
          out << "Z[" << sigma1 << sigma2 << sigma3 << "],";
        }
  }
  out << std::endl;
}

int process_command_line_args(int argc, char *argv[]);
bool open_results_file(std::ofstream &out);

int main(int argc, char *argv[]) {
  if (!process_command_line_args(argc, argv))
    return 1;

  std::ofstream out;
  if (!open_results_file(out))
    return 1;

  write_csv_header(out);
  out << L << "," << lambda << ",";
  auto diags = compute_all_diagrams(10 * L, out);
  auto parts = compute_2_3pt_functions(diags, out);
  out << std::endl;

  return 0;
}

int process_command_line_args(int argc, char *argv[]) {
  const char usage[] = "Usage: %s L lambda k1 k2 k3\n";
  int opt;
  // parse command-line options
  while ((opt = getopt(argc, argv, "h")) != -1) {
    switch (opt) {
    case 'h':
      printf(usage, argv[0]);
      return 1;
    case '?':
      printf(usage, argv[0]);
      break;
    }
  }

  // parse other arguments
  if (argc - optind != 5) {
    printf(usage, argv[0]);
    return 0;
  }
  L = std::atoi(argv[optind++]);
  lambda = Weight(std::atof(argv[optind++]));
  k1 = std::atoi(argv[optind++]);
  k2 = std::atoi(argv[optind++]);
  k3 = std::atoi(argv[optind++]);

  std::cout << "L = " << L << ", lambda = " << lambda << ", k1 = " << k1
            << ", k2 = " << k2 << ", k3 = " << k3 << std::endl;

  return 1;
}

bool open_results_file(std::ofstream &out) {
  std::string dirname = std::to_string(k1) + std::to_string(k2) + std::to_string(k3);
  fs::path src_dir = fs::path(__FILE__).parent_path();
  fs::path dir = src_dir / "results" / dirname;
  if (!fs::exists(dir)) {
    if (!fs::create_directories(dir)) {
      std::cerr << "Error: could not create directory " << dir << "\n";
      return 0;
    }
  }

  // Format lambda nicely (avoid too many decimals in filename)
  std::ostringstream lambda_str;
  lambda_str << std::fixed << std::setprecision(3) << lambda;

  // Build filename
  std::string filename =
      "L=" + std::to_string(L) + "_lambda=" + lambda_str.str() + ".csv";

  fs::path filepath = dir / filename;

  out.open(filepath);

  if (!out.is_open()) {
    std::cerr << "Error: could not open " << filepath << " for writing\n";
    return 0;
  } else {
    std::cout << "Writing in " << filepath << std::endl;
  }
  return 1;
}
