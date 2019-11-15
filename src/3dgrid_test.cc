#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
using namespace dealii;

int main()
{
  Triangulation<3> triangulation;
  const int nx = 190;
  const int ny = 190;
  const int nz = 150;

  // Build the points
  std::vector<Point<3>> vertices((nx+1)*(ny+1)*(nz+1));
  for (int i = 0; i < nx + 1; ++i) 
    for (int j = 0; j < ny + 1; ++j)
      for (int k = 0; k < nz + 1; ++k)
        vertices[k*(nx+1)*(ny+1) + j*(nx+1) + i] = 
        {(double)i, (double)j, (double)k};

  // Build the cells
  std::vector<CellData<3>> cells(nx*ny*nz, CellData<3>());
  int n = 0;
  for (int i = 0; i < nx; ++i) 
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
        cells[n].vertices[0] = k*(nx+1)*(ny+1) + j*(nx+1) + i;
        cells[n].vertices[1] = k*(nx+1)*(ny+1) + j*(nx+1) + i + 1;
        cells[n].vertices[2] = k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i;
        cells[n].vertices[3] = k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i + 1;
        cells[n].vertices[4] = (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i;
        cells[n].vertices[5] = (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i + 1;
        cells[n].vertices[6] = (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i;
        cells[n].vertices[7] = (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i + 1;
        cells[n].material_id = 0;
        ++n;
      }
  auto t1 = Clock::now();
  triangulation.create_triangulation(vertices, cells, SubCellData());
  auto t2 = Clock::now();
  std::cout << "Grid created in "
    << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
    << " s" << std::endl;
  
  std::ofstream out("3dgrid.eps");
  GridOut       grid_out;
  grid_out.write_eps(triangulation, out);
  
  // Output grid
  /*
  t1 = Clock::now();
  DataOut<3> data_out;
  DoFHandler<3> dof_handler(triangulation);
  FE_Q<3> fe(1);
  dof_handler.distribute_dofs(fe);
  data_out.attach_dof_handler(dof_handler);
  Vector<double> solution(dof_handler.n_dofs()); 
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::ofstream output("3dgrid.vtk");
  data_out.write_vtk(output);
  t2 = Clock::now();
  std::cout << "Grid written in "
    << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
    << " s" << std::endl;
  
  data_out.clear();
  dof_handler.clear();
  */
}
