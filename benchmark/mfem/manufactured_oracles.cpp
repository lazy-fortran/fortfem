// Independent MFEM implementations of the shared Poisson and Ampere cases.
#include "mfem.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

using mfem::Array;
using mfem::BilinearForm;
using mfem::Coefficient;
using mfem::ConstantCoefficient;
using mfem::CurlCurlIntegrator;
using mfem::DiffusionIntegrator;
using mfem::DomainLFIntegrator;
using mfem::Element;
using mfem::FiniteElementSpace;
using mfem::FunctionCoefficient;
using mfem::Geometry;
using mfem::GridFunction;
using mfem::GSSmoother;
using mfem::H1_FECollection;
using mfem::IntegrationRule;
using mfem::IntRules;
using mfem::LinearForm;
using mfem::Mesh;
using mfem::ND_FECollection;
using mfem::OperatorPtr;
using mfem::PCG;
using mfem::SparseMatrix;
using mfem::Vector;
using mfem::VectorCoefficient;
using mfem::VectorFEDomainLFIntegrator;
using mfem::VectorFEMassIntegrator;
using mfem::VectorFunctionCoefficient;

namespace {
using Clock = std::chrono::steady_clock;

struct Record {
  int n;
  int cells;
  int dofs;
  double l2_error;
  double graph_error;
  double assembly_seconds;
  double solve_seconds;
  double total_seconds;
};

double elapsed(const Clock::time_point &start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

double poisson_exact(const Vector &point) {
  return std::sin(M_PI * point(0)) * std::sin(M_PI * point(1));
}

double poisson_source(const Vector &point) {
  return 2.0 * M_PI * M_PI * poisson_exact(point);
}

void poisson_gradient(const Vector &point, Vector &gradient) {
  gradient(0) = M_PI * std::cos(M_PI * point(0)) * std::sin(M_PI * point(1));
  gradient(1) = M_PI * std::sin(M_PI * point(0)) * std::cos(M_PI * point(1));
}

void ampere_exact(const Vector &point, Vector &field) {
  field(0) = std::sin(M_PI * point(1));
  field(1) = std::sin(M_PI * point(0));
}

void ampere_source(const Vector &point, Vector &source) {
  ampere_exact(point, source);
  source *= M_PI * M_PI + 1.0;
}

void ampere_curl(const Vector &point, Vector &curl) {
  curl(0) = M_PI * (std::cos(M_PI * point(0)) - std::cos(M_PI * point(1)));
}

Array<int> essential_dofs(Mesh &mesh, FiniteElementSpace &space) {
  Array<int> essential;
  Array<int> marker(mesh.bdr_attributes.Max());
  marker = 1;
  space.GetEssentialTrueDofs(marker, essential);
  return essential;
}

void solve_system(BilinearForm &form, LinearForm &load,
                  const Array<int> &essential, GridFunction &solution,
                  double &assembly_seconds, double &solve_seconds) {
  const auto assembly_start = Clock::now();
  load.Assemble();
  form.Assemble();
  OperatorPtr system;
  Vector right_hand_side;
  Vector unknown;
  form.FormLinearSystem(essential, solution, load, system, unknown,
                        right_hand_side);
  assembly_seconds = elapsed(assembly_start);

  const auto solve_start = Clock::now();
  GSSmoother smoother(static_cast<SparseMatrix &>(*system));
  PCG(*system, smoother, right_hand_side, unknown, 0, 2000, 1.0e-12, 0.0);
  solve_seconds = elapsed(solve_start);

  Vector residual(system->Height());
  system->Mult(unknown, residual);
  residual -= right_hand_side;
  const double relative_residual =
      residual.Norml2() / std::max(1.0, right_hand_side.Norml2());
  if (relative_residual > 1.0e-5) {
    throw std::runtime_error("MFEM linear solve residual is " +
                             std::to_string(relative_residual));
  }
  form.RecoverFEMSolution(unknown, load, solution);
}

Record solve_poisson(int n) {
  const auto total_start = Clock::now();
  Mesh mesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1.0, 1.0);
  H1_FECollection collection(1, 2);
  FiniteElementSpace space(&mesh, &collection);
  const Array<int> essential = essential_dofs(mesh, space);
  FunctionCoefficient exact(poisson_exact);
  FunctionCoefficient source(poisson_source);
  VectorFunctionCoefficient gradient(2, poisson_gradient);
  LinearForm load(&space);
  auto *load_integrator = new DomainLFIntegrator(source);
  load_integrator->SetIntRule(&IntRules.Get(Geometry::TRIANGLE, 8));
  load.AddDomainIntegrator(load_integrator);
  BilinearForm form(&space);
  form.AddDomainIntegrator(new DiffusionIntegrator);
  GridFunction solution(&space);
  solution = 0.0;
  double assembly_seconds;
  double solve_seconds;
  solve_system(form, load, essential, solution, assembly_seconds,
               solve_seconds);
  const IntegrationRule *rules[Geometry::NumGeom] = {};
  rules[Geometry::TRIANGLE] = &IntRules.Get(Geometry::TRIANGLE, 10);
  return {n,
          mesh.GetNE(),
          space.GetTrueVSize(),
          solution.ComputeL2Error(exact, rules),
          solution.ComputeH1Error(&exact, &gradient, rules),
          assembly_seconds,
          solve_seconds,
          elapsed(total_start)};
}

Record solve_ampere(int n) {
  const auto total_start = Clock::now();
  Mesh mesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 1.0, 1.0);
  ND_FECollection collection(1, 2);
  FiniteElementSpace space(&mesh, &collection);
  const Array<int> essential = essential_dofs(mesh, space);
  VectorFunctionCoefficient exact(2, ampere_exact);
  VectorFunctionCoefficient source(2, ampere_source);
  VectorFunctionCoefficient curl(1, ampere_curl);
  LinearForm load(&space);
  auto *load_integrator = new VectorFEDomainLFIntegrator(source);
  load_integrator->SetIntRule(&IntRules.Get(Geometry::TRIANGLE, 8));
  load.AddDomainIntegrator(load_integrator);
  ConstantCoefficient one(1.0);
  BilinearForm form(&space);
  form.AddDomainIntegrator(new CurlCurlIntegrator(one));
  form.AddDomainIntegrator(new VectorFEMassIntegrator(one));
  GridFunction solution(&space);
  solution = 0.0;
  double assembly_seconds;
  double solve_seconds;
  solve_system(form, load, essential, solution, assembly_seconds,
               solve_seconds);
  const IntegrationRule *rules[Geometry::NumGeom] = {};
  rules[Geometry::TRIANGLE] = &IntRules.Get(Geometry::TRIANGLE, 10);
  return {n,
          mesh.GetNE(),
          space.GetTrueVSize(),
          solution.ComputeL2Error(exact, rules),
          solution.ComputeHCurlError(&exact, &curl, rules),
          assembly_seconds,
          solve_seconds,
          elapsed(total_start)};
}

void require_decreasing(const std::vector<Record> &records,
                        const std::string &case_name) {
  for (std::size_t index = 1; index < records.size(); ++index) {
    if (records[index].l2_error >= records[index - 1].l2_error ||
        records[index].graph_error >= records[index - 1].graph_error) {
      throw std::runtime_error(case_name + " errors did not decrease");
    }
  }
}

void write_records(const std::filesystem::path &path,
                   const std::vector<Record> &records,
                   const std::string &graph_name) {
  std::ofstream stream(path);
  stream << "mesh_id,cells,dofs,h,l2_error," << graph_name
         << ",assembly_seconds,solve_seconds,total_seconds\n";
  stream << std::setprecision(17);
  for (const Record &record : records) {
    stream << "unit-square-n" << record.n << ',' << record.cells << ','
           << record.dofs << ',' << 1.0 / record.n << ',' << record.l2_error
           << ',' << record.graph_error << ',' << record.assembly_seconds << ','
           << record.solve_seconds << ',' << record.total_seconds << '\n';
  }
}
} // namespace

int main(int argc, char **argv) {
  const std::filesystem::path output =
      argc > 1 ? argv[1] : "benchmark-results/mfem";
  std::filesystem::create_directories(output);
  const std::vector<int> levels = {4, 8, 16, 32, 64};
  std::vector<Record> poisson;
  std::vector<Record> ampere;
  for (const int n : levels) {
    poisson.push_back(solve_poisson(n));
    ampere.push_back(solve_ampere(n));
  }
  require_decreasing(poisson, "Poisson");
  require_decreasing(ampere, "Ampere");
  write_records(output / "poisson.csv", poisson, "h1_error");
  write_records(output / "ampere.csv", ampere, "hcurl_error");
}
