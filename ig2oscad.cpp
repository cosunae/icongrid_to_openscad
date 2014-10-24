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

#define MAX_NDIMS 40

void dump_header(std::ofstream& ofs)
{
 ofs << "use <polyhedra_functions.scad>;" << std::endl;
 ofs<< " use <Write.scad> " << std::endl;

}

void dump_scad_grid(std::vector<double>& cartx_vert, std::vector<double>& carty_vert, std::vector<double>& cartz_vert, std::vector<double>& vector_cell)
{
    std::vector<std::string> colors;
    colors.push_back("NavajoWhite");
    colors.push_back("YellowGreen");
    colors.push_back("Gold");
    colors.push_back("Crimson");
    colors.push_back("Purple");
    colors.push_back("DarkSlateGray");
    colors.push_back("OrangeRed");
    colors.push_back("Lime");
    colors.push_back("Olive");
    colors.push_back("Plum");
    colors.push_back("Maroon");

    assert( cartx_vert.size() == carty_vert.size() && cartx_vert.size() == cartz_vert.size() );
    assert(vector_cell.size()%3==0);
    std::ofstream file("output.scad");

    dump_header(file);

//    file << "color(\"Blue\",0.5) { \n  polyhedron(" << std::endl;
    file << "  allpoints = [" << std::endl;

    for(int i=0; i < cartx_vert.size()-1; ++i){
          file << "    [" << cartx_vert[i] << ", " << carty_vert[i] << ", " << cartz_vert[i] << "],"<< std::endl;
    }
    file << "    [" << cartx_vert[cartx_vert.size()-1] << ", " << carty_vert[cartx_vert.size()-1] << ", " << cartz_vert[cartx_vert.size()-1] << "],"<< std::endl;

    file << "  ];" << std::endl;
    //,\n  triangles = [" << std::endl;
    int triangles_cnt = 0;
    
//    file << std::string("color(\"") + colors[color_cnt]+"\",0.7) {\n polyhedron( \n  points = allpoints,\n faces = [" << std::endl;
    file << "faces" << boost::lexical_cast<std::string>(triangles_cnt) << " = [ " << std::endl;

    int ncells = vector_cell.size()/3;
    for(int i=0; i < ncells; ++i)
    {
        if((i+1)%256==0) {
            file<<"    [" << vector_cell[0*ncells + i]-1 << "," << vector_cell[1*ncells+i]-1 << "," << vector_cell[2*ncells+i]-1 << "]" << std::endl;
            file << "  ];\n" << std::endl;
            
            file << std::string("color(\"") + colors[triangles_cnt % colors.size()]+"\",0.7) {\n difference() { \n polyhedron( " << 
            " points = allpoints, faces = faces" << boost::lexical_cast<std::string>(triangles_cnt) << ");\n" << 
            " engrave_face_word(faces"<< boost::lexical_cast<std::string>(triangles_cnt) << ",allpoints,\n" << 
            "font=\"orbitron\",\nword=\"123456789\",\nratio=1,thickness=0.3);\n}\n}" << std::endl;

//            file<<"    [" << vector_cell[0*ncells + i]-1 << "," << vector_cell[1*ncells+i]-1 << "," << vector_cell[2*ncells+i]-1 << "]" << std::endl;
//            file << "  ]\n); }" << std::endl;
            ++triangles_cnt;
            if(i != ncells-1) {
                file << "faces" << boost::lexical_cast<std::string>(triangles_cnt) << " = [ " << std::endl;
//                file << std::string("color(\"") + colors[color_cnt]+"\",0.7) {\n polyhedron( \n  points = allpoints,\n faces = [" << std::endl;
            }
        }
        else {
            file<<"    [" << vector_cell[0*ncells + i]-1 << "," << vector_cell[1*ncells+i]-1 << "," << vector_cell[2*ncells+i]-1 << "]," << std::endl;
        }      
    }

//        file<<"    [" << vector_cell[i*3+0]-1 << "," << vector_cell[i*3+1]-1 << "," << vector_cell[i*3+2]-1 << "]," << std::endl;
//              file<<"    [" << vector_cell[0*ncells + i]-1 << "," << vector_cell[1*ncells+i]-1 << "," << vector_cell[2*ncells+i]-1 << "]," << std::endl;
  //  file<<"    [" << vector_cell[(ncells-3)*3+0]-1 << "," << vector_cell[(ncells-3)*3+1]-1 << "," << vector_cell[(ncells-3)*3+2]-1 << "]" << std::endl;
//    file<< "  ]\n);" << std::endl;
//    file << "}" << std::endl;

    file.close();
}

int main(int argc, char** argv)
{  
   int ncid;
   int ierror = nc_open("iconR2B03-grid.nc", NC_NOWRITE, & ncid);
   assert(ierror==NC_NOERR);

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

   dump_scad_grid(cartx_vert, carty_vert, cartz_vert, vertex_cell);

   nc_close(ncid);
};
