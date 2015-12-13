/**
 Copyright (c) 2014, Carlos Osuna (cosunae@gmail.com)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "netcdf.h"
#include <assert.h>
#include <vector>
#include <fstream> 
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

#define MAX_NDIMS 40


void gmsh_render(int ncells_per_diamond, std::vector<double>& cartx_vert, std::vector<double>& carty_vert, std::vector<double>& cartz_vert, std::vector<double>& vertex_cell, bool doLabel, int nmpi=0)
{
    assert( cartx_vert.size() == carty_vert.size() && cartx_vert.size() == cartz_vert.size() );
    assert(vertex_cell.size()%3==0);

    std::ofstream pp;
    pp.open("outgrid.msh");
    pp << "$MeshFormat" << std::endl;
    pp << "2.2 0 8" << std::endl;
    pp << "$EndMeshFormat" << std::endl;
    pp << "$Nodes" << std::endl;
    pp << cartx_vert.size() << std::endl;

    int factor=10000000;
    for(int i=0; i < cartx_vert.size(); ++i){
        pp << i+1 << " " << (int)(cartx_vert[i]*factor) << " " << (int)(carty_vert[i]*factor) << " " << (int)(cartz_vert[i]*factor) << std::endl;
    }

    pp << "$Elements" << std::endl;
    int ncells = vertex_cell.size()/3;

    int ncell_partition = ncells/nmpi;
    if(ncell_partition%4)
        ncell_partition += (4-ncell_partition%4);

    pp << ncells << std::endl;

    for(int i=0; i < ncells; ++i)
    {
        int partition = i/ncell_partition;
        pp << i+1 << " " << 2 << " " << 4 << " " << 1 << " " << 1 << " " << 1  << " " << partition << " " << vertex_cell[0*ncells + i] << " " << 
            vertex_cell[1*ncells + i] << " " <<
            vertex_cell[2*ncells + i] << " " << std::endl;
    }

    pp << "$EndElements" << std::endl;
}

int main(int argc, char** argv)
{  
   namespace po = boost::program_options;

   int ndiamonds_label=0;
   int nmpi=0;
   bool doLabels= false;

   // Declare the supported options.
   po::options_description desc("Allowed options");
   desc.add_options()
      ("help", "produce help message")
      ("label", "label all cells and vertices")
      ("ndiamonds-label", po::value<int>(&ndiamonds_label)->default_value(1),  "number of diamonds to label (-1 for all)")
      ("input-file", po::value< std::vector<std::string> >(), "input file")
      ("nmpi", po::value<int>(&nmpi)->default_value(1),  "number of mpi")
   ;
   
   po::positional_options_description p;
   p.add("input-file", -1);

   po::variables_map vm;
   po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
   po::notify(vm);    

   if (vm.count("help")) {
       std::cout << desc << "\n";
       return 1;
   }

   if(!vm.count("input-file"))
   { 
       std::cout << "Error: need to specify input file" << std::endl;
       std::cout << desc << std::endl;
       return 1;
   }

   const std::string ncfilename = vm["input-file"].as< std::vector<std::string> >()[0];
   
   if ( !boost::filesystem::exists( ncfilename ) )
   {
       std::cout << "Can't find netcdf file!" << std::endl;
       return 1;
   }

   if(vm.count("label")) doLabels=true;


   int ncid;
   int ierror = nc_open(ncfilename.c_str(), NC_NOWRITE, & ncid);
   assert(ierror==NC_NOERR);

   int cell_dimid;
   ierror = nc_inq_dimid(ncid, "cell", &cell_dimid);
   assert(ierror==NC_NOERR);

   size_t cell_dimlen;
   ierror = nc_inq_dimlen(ncid, cell_dimid, &cell_dimlen);
   assert(ierror==NC_NOERR);

   assert(cell_dimlen%20==0);

   int vertex_dimid;
   ierror = nc_inq_dimid(ncid, "vertex", &vertex_dimid);
   assert(ierror==NC_NOERR);
 
   size_t vertex_dimlen;
   ierror = nc_inq_dimlen(ncid, vertex_dimid, &vertex_dimlen);
   assert(ierror==NC_NOERR);

   int cartx_vert_id, carty_vert_id, cartz_vert_id;
   ierror= nc_inq_varid(ncid, "cartesian_x_vertices", &cartx_vert_id);
   assert(ierror==NC_NOERR);
   ierror= nc_inq_varid(ncid, "cartesian_y_vertices", &carty_vert_id);
   assert(ierror==NC_NOERR);
   ierror= nc_inq_varid(ncid, "cartesian_z_vertices", &cartz_vert_id);
   assert(ierror==NC_NOERR);


   int cart_vert_numdim;
   ierror = nc_inq_varndims(ncid, cartx_vert_id, &cart_vert_numdim);
   assert(ierror==NC_NOERR && cart_vert_numdim==1);
   ierror = nc_inq_varndims(ncid, carty_vert_id, &cart_vert_numdim);
   assert(ierror==NC_NOERR && cart_vert_numdim==1);
   ierror = nc_inq_varndims(ncid, cartz_vert_id, &cart_vert_numdim);
   assert(ierror==NC_NOERR && cart_vert_numdim==1);


   int cart_vert_dimids[1];
   ierror = nc_inq_vardimid(ncid, cartx_vert_id, cart_vert_dimids);
   assert(ierror==NC_NOERR && cart_vert_dimids[0] == vertex_dimid);
   ierror = nc_inq_vardimid(ncid, carty_vert_id, cart_vert_dimids);
   assert(ierror==NC_NOERR && cart_vert_dimids[0] == vertex_dimid);
   ierror = nc_inq_vardimid(ncid, cartz_vert_id, cart_vert_dimids);
   assert(ierror==NC_NOERR && cart_vert_dimids[0] == vertex_dimid);

   std::vector<double> cartx_vert, carty_vert, cartz_vert;

   cartx_vert.resize( vertex_dimlen );
   carty_vert.resize( vertex_dimlen );
   cartz_vert.resize( vertex_dimlen );

   ierror = nc_get_var_double(ncid, cartx_vert_id, &cartx_vert[0]);
   assert(ierror==NC_NOERR);
   ierror= nc_get_var_double(ncid, carty_vert_id, &carty_vert[0]);
   assert(ierror==NC_NOERR);
   ierror= nc_get_var_double(ncid, cartz_vert_id, &cartz_vert[0]);
   assert(ierror==NC_NOERR);

   int vertex_cell_id;
   ierror= nc_inq_varid(ncid, "vertex_of_cell", &vertex_cell_id);
   assert(ierror==NC_NOERR);
   int vertex_cell_numdim;
   ierror = nc_inq_varndims(ncid, vertex_cell_id, &vertex_cell_numdim);
   assert(ierror==NC_NOERR && vertex_cell_numdim==2);

   int vertex_cell_dimids[2];
   size_t vertex_cell_dimlens[2];
   ierror = nc_inq_vardimid(ncid, vertex_cell_id, vertex_cell_dimids);
   assert(ierror==NC_NOERR);

   ierror = nc_inq_dimlen(ncid, vertex_cell_dimids[0], &vertex_cell_dimlens[0]);
   assert(ierror==NC_NOERR && vertex_cell_dimlens[0] == 3);
   ierror = nc_inq_dimlen(ncid, vertex_cell_dimids[1], &vertex_cell_dimlens[1]);
   assert(ierror==NC_NOERR);

   std::vector<double> vertex_cell;
   vertex_cell.resize( vertex_cell_dimlens[0] * vertex_cell_dimlens[1] );
   ierror = nc_get_var_double(ncid, vertex_cell_id, &vertex_cell[0]);
   assert(ierror==NC_NOERR);

   gmsh_render(cell_dimlen/20, cartx_vert, carty_vert, cartz_vert, vertex_cell, doLabels, nmpi);

   nc_close(ncid);
};
