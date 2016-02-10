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
#include <limits>

#define MAX_NDIMS 40
#define PI 3.1415

struct neighbours_of
{
    neighbours_of(unsigned int neigh_dim, unsigned int n_nodes) :
        m_neigh_dim(neigh_dim), m_data(n_nodes*neigh_dim) {}

    double* data_ptr() { return &(m_data[0]);}

    double& at(unsigned int nnode_id, unsigned int neigh_id) {
        assert(nnode_id > 0);
        assert(nnode_id -1 < n_nodes());
        assert(neigh_id < m_neigh_dim);
        return m_data[index(nnode_id, neigh_id)];
    }

    unsigned int index(unsigned int node_id, unsigned int neigh_id){
        return (node_id-1)*neigh_dim() + neigh_id;
    }

    std::vector<double> neighbours_of_node(unsigned int node_id)
    {
        assert(node_id>0);
        std::vector<double> res;
        std::copy(m_data.begin() + index(node_id,0), m_data.begin() + index(node_id, neigh_dim()), std::back_inserter(res));
        return res;
    }

    void insert_node(std::vector<double> neighbours)
    {
        assert(neighbours.size() == m_neigh_dim);

        std::copy(neighbours.begin(), neighbours.end(), std::back_inserter(m_data));
    }

    bool node_has_neighbour_in(unsigned int node_id, std::vector<unsigned int> list_id)
    {
        assert(node_id>0);
        for(size_t i=0; i < neigh_dim(); ++i) {
            if(std::find(list_id.begin(), list_id.end(), at(node_id, i)) != list_id.end()) return true;
        }
        return false;
    }

    unsigned int n_nodes() const {
        assert((m_data.size() % m_neigh_dim)==0);
        return m_data.size() / m_neigh_dim;}
    unsigned int neigh_dim() const { return m_neigh_dim;}
private: 
    std::vector<double> m_data;
    unsigned int m_neigh_dim;
};

class mesh
{
public:
    mesh(unsigned int ncells, unsigned int nvertices, unsigned int n_vertices_of_cell, unsigned int n_vertices_of_vertex) :
        m_n_vertices_of_cell(n_vertices_of_cell), m_n_vertices_of_vertex(n_vertices_of_vertex),
        m_cartx_vert(nvertices), m_carty_vert(nvertices), m_cartz_vert(nvertices),
        m_vertex_of_cell(n_vertices_of_cell, ncells), m_vertices_of_vertex(n_vertices_of_vertex, nvertices)
    {}

    std::vector<double >& cart_vert(unsigned int dim) {
        assert(dim < 3);
        if(dim == 0) return m_cartx_vert;
        else if(dim == 1) return m_cartx_vert;
        else if(dim == 2) return m_cartx_vert;
    }
    neighbours_of& vertex_of_cell() {return m_vertex_of_cell;}
    neighbours_of& vertices_of_vertex() { return m_vertices_of_vertex;}

    void insert_cart_vert(std::vector<double>& cartx, std::vector<double>& carty, std::vector<double>& cartz) {
        assert(cartx.size() == carty.size() && carty.size() == cartz.size());
        for(size_t i=0; i != cartx.size(); ++i)
        {
            (m_cartx_vert[i])= cartx[i];
            (m_carty_vert[i]) = carty[i];
            (m_cartz_vert[i]) = cartz[i];
        }
    }

    void gmsh_render(bool doLabel, int nmpi=0)
    {
        std::ofstream pp;
        pp.open("outgrid.msh");
        pp << "$MeshFormat" << std::endl;
        pp << "2.2 0 8" << std::endl;
        pp << "$EndMeshFormat" << std::endl;
        pp << "$Nodes" << std::endl;
        pp << m_cartx_vert.size() << std::endl;

        int factor=10000000;
        for(int i=0; i < m_cartx_vert.size(); ++i){
            pp << i+1 << " " << (int)((m_cartx_vert[i])*factor) << " " << (int)((m_carty_vert[i])*factor) << " " <<
                                (int)((m_cartz_vert[i])*factor) << std::endl;
        }

        pp << "$Elements" << std::endl;

        int ncell_partition = ncells()/nmpi;
        if(ncell_partition%4)
            ncell_partition += (4-ncell_partition%4);

        pp << ncells() << std::endl;

        for(int i=0; i < ncells(); ++i)
        {
            int partition = i/ncell_partition;
            pp << i+1 << " " << 2 << " " << 4 << " " << 1 << " " << 1 << " " << 1  << " " << partition << " " <<
                  m_vertex_of_cell.at(i+1, 0) << " " <<
                m_vertex_of_cell.at(i+1,1) << " " <<
                m_vertex_of_cell.at(i+1,2) << " " << std::endl;
        }

        pp << "$EndElements" << std::endl;
    }

    void add_ghost_nodes()
    {
        auto north_pole = std::max_element(m_carty_vert.begin(), m_carty_vert.end());
        auto south_pole = std::min_element(m_carty_vert.begin(), m_carty_vert.end());

        int north_pole_idx = std::distance(m_carty_vert.begin(), north_pole);
        int south_pole_idx = std::distance(m_carty_vert.begin(), south_pole);

        std::vector<unsigned int> nodes_id_east_with_ghost, nodes_id_west_with_ghost;
        recur_find_east(north_pole_idx, south_pole_idx, nodes_id_east_with_ghost);
        recur_find_west(north_pole_idx, south_pole_idx, nodes_id_west_with_ghost);

        insert_ghost_nodes(nodes_id_east_with_ghost, nodes_id_west_with_ghost);
    }

    void fill_from_netcdf(const int ncid, const int vertex_cell_id, const int vertices_of_vertex_id)
    {

        const unsigned int buffsize = std::max(m_vertex_of_cell.n_nodes(), m_vertices_of_vertex.n_nodes()) *
                std::max(m_vertex_of_cell.neigh_dim(), m_vertices_of_vertex.neigh_dim());
        std::vector<double> tmp(buffsize);

        int ierror = nc_get_var_double(ncid, vertex_cell_id, &tmp[0]);
        assert(ierror==NC_NOERR);
        unsigned int n_nodes = m_vertex_of_cell.n_nodes();
        unsigned int neigh_dim = m_vertex_of_cell.neigh_dim();

        for(unsigned int node_it = 0; node_it != n_nodes; ++node_it)
        {
            for(unsigned int neigh_it = 0; neigh_it != neigh_dim; ++neigh_it)
            {
                assert(node_it + neigh_it*n_nodes < buffsize);
                m_vertex_of_cell.at(node_it+1, neigh_it) = tmp[node_it + neigh_it*n_nodes];
            }
        }

        ierror = nc_get_var_double(ncid, vertices_of_vertex_id, &tmp[0]);
        assert(ierror==NC_NOERR);
        n_nodes = m_vertices_of_vertex.n_nodes();
        neigh_dim = m_vertices_of_vertex.neigh_dim();

        for(unsigned int node_it = 0; node_it != n_nodes; ++node_it)
        {
            for(unsigned int neigh_it = 0; neigh_it != neigh_dim; ++neigh_it)
            {
                assert(node_it + neigh_it*n_nodes < buffsize);
                m_vertices_of_vertex.at(node_it+1, neigh_it) = tmp[node_it + neigh_it*n_nodes];
            }
        }
    }

    unsigned int ncells() const { m_vertex_of_cell.n_nodes();}
private:

    void recur_find_west(unsigned int root_idx, unsigned int south_pole_idx, std::vector<unsigned int>& nodes_id_with_ghost)
    {
        double lat = m_carty_vert[root_idx];
        const unsigned int root_id = root_idx+1;
        const unsigned int south_pole_id = south_pole_idx+1;

        double min_lon=std::numeric_limits<double>::max();
        int next_id =0;
        for(int vertex_cnt = 0; vertex_cnt != m_vertices_of_vertex.neigh_dim(); ++vertex_cnt)
        {
            int vertex_id = m_vertices_of_vertex.at(root_id, vertex_cnt);
            if(vertex_id == 0) continue;

            unsigned int vertex_idx = vertex_id-1;

            if(m_carty_vert[vertex_idx] > lat) continue;

            if(vertex_id == south_pole_id) return;

            if(min_lon > m_cartx_vert[vertex_idx]) {
                min_lon = m_cartx_vert[vertex_idx];
                next_id = vertex_id;
            }
        }

        nodes_id_with_ghost.push_back(next_id);
        recur_find_west(next_id-1, south_pole_idx, nodes_id_with_ghost);
    }

    void recur_find_east(unsigned int root_idx, unsigned int south_pole_idx, std::vector<unsigned int>& nodes_id_with_ghost)
    {
        double lat = m_carty_vert[root_idx];
        const unsigned int root_id = root_idx+1;
        const unsigned int south_pole_id = south_pole_idx+1;

        double max_lon=-std::numeric_limits<double>::max();
        int next_id =0;
        for(int vertex_cnt = 0; vertex_cnt != m_vertices_of_vertex.neigh_dim(); ++vertex_cnt)
        {
            int vertex_id = m_vertices_of_vertex.at(root_id, vertex_cnt);
            if(vertex_id == 0) continue;

            unsigned int vertex_idx = vertex_id-1;

            if(m_carty_vert[vertex_idx] > lat) continue;

            if(vertex_id == south_pole_id) return;

            if(max_lon < m_cartx_vert[vertex_idx]) {
                max_lon = m_cartx_vert[vertex_idx];
                next_id = vertex_id;
            }
        }

        nodes_id_with_ghost.push_back(next_id);
        recur_find_east(next_id-1, south_pole_idx, nodes_id_with_ghost);
    }

    void replace_and_insert_ghost_elements(std::vector<unsigned int>& east_nodes_id_with_ghost, std::vector<unsigned int>& west_nodes_id_with_ghost,
                                           std::vector<unsigned int>& east_ghost_nodes_id, std::vector<unsigned int>& west_ghost_nodes_id)
    {
        std::vector< std::shared_ptr<std::vector<double> > > new_elements;
        for(size_t element_cnt = 0; element_cnt != m_vertex_of_cell.n_nodes(); ++element_cnt)
        {
            if(m_vertex_of_cell.node_has_neighbour_in(element_cnt+1, east_nodes_id_with_ghost) &&
                    m_vertex_of_cell.node_has_neighbour_in(element_cnt+1, west_nodes_id_with_ghost))
            {

                std::shared_ptr<std::vector<double> > orig_neigh = std::make_shared<std::vector<double> >(m_vertex_of_cell.neighbours_of_node(element_cnt+1));
                for(int c=0; c != m_vertex_of_cell.neigh_dim(); ++c){
                    std::vector<unsigned int>::iterator iter;
                    while( (iter = std::find(west_nodes_id_with_ghost.begin(), west_nodes_id_with_ghost.end(), m_vertex_of_cell.at(element_cnt+1, c))) !=
                            west_nodes_id_with_ghost.end())
                    {
                        unsigned int dist = std::distance(west_nodes_id_with_ghost.begin(), iter);
                        m_vertex_of_cell.at(element_cnt+1, c) = east_ghost_nodes_id[dist];
                    }

                    while( (iter = std::find(east_nodes_id_with_ghost.begin(), east_nodes_id_with_ghost.end(), (*orig_neigh)[c])) !=
                            east_nodes_id_with_ghost.end())
                    {
                        unsigned int dist = std::distance(east_nodes_id_with_ghost.begin(), iter);
                        (*orig_neigh)[c] = west_ghost_nodes_id[dist];
                    }
                }
                new_elements.push_back(orig_neigh);
            }
        }
        for(int c=0; c != new_elements.size(); ++c)
        {
            m_vertex_of_cell.insert_node(*(new_elements[c]));
        }

    }

    void insert_ghost_nodes(std::vector<unsigned int>& east_nodes_id_with_ghost,
                            std::vector<unsigned int>& west_nodes_id_with_ghost)
    {
        std::vector<unsigned int> east_ghost_nodes_id;
        std::vector<unsigned int> west_ghost_nodes_id;

        for(size_t i= 0; i != east_nodes_id_with_ghost.size(); ++i)
        {
            unsigned int node_id = east_nodes_id_with_ghost[i];
            double lat = m_carty_vert[node_id-1];

            double newlon = m_cartx_vert[node_id-1]-2*PI;
            unsigned int newnode_id = m_cartx_vert.size()+1;

            m_cartx_vert.push_back(newlon);
            m_carty_vert.push_back(lat);
            m_cartz_vert.push_back(0);
            west_ghost_nodes_id.push_back(newnode_id);
        }

        for(size_t i= 0; i != west_nodes_id_with_ghost.size(); ++i)
        {
            unsigned int node_id = west_nodes_id_with_ghost[i];
            double lat = m_carty_vert[node_id-1];

            double newlon = m_cartx_vert[node_id-1] + 2*PI;
            unsigned int newnode_id = m_cartx_vert.size()+1;

            m_cartx_vert.push_back(newlon);
            m_carty_vert.push_back(lat);
            m_cartz_vert.push_back(0);
            east_ghost_nodes_id.push_back(newnode_id);
        }

        for(size_t i= 0; i != east_nodes_id_with_ghost.size(); ++i)
        {
            unsigned int node_id = east_nodes_id_with_ghost[i];

            for(size_t c=0; c != m_vertices_of_vertex.neigh_dim(); ++c)
            {
                unsigned int neigh_node_id = m_vertices_of_vertex.at(node_id, c);

                std::vector<unsigned int>::iterator iter;
                if((iter = std::find(west_nodes_id_with_ghost.begin(), west_nodes_id_with_ghost.end(), neigh_node_id))
                        != west_nodes_id_with_ghost.end() ) {
                    m_vertices_of_vertex.at(node_id, c) = east_ghost_nodes_id[std::distance(west_nodes_id_with_ghost.begin(), iter)];
                }
            }
        }

        for(size_t i= 0; i != west_nodes_id_with_ghost.size(); ++i)
        {
            unsigned int node_id = west_nodes_id_with_ghost[i];

            for(size_t c=0; c != m_vertices_of_vertex.neigh_dim(); ++c)
            {
                unsigned int neigh_node_id = m_vertices_of_vertex.at(node_id, c);

                std::vector<unsigned int>::iterator iter;
                if((iter = std::find(east_nodes_id_with_ghost.begin(), east_nodes_id_with_ghost.end(), neigh_node_id))
                        != east_nodes_id_with_ghost.end() ) {
                    m_vertices_of_vertex.at(node_id, c) = west_ghost_nodes_id[std::distance(east_nodes_id_with_ghost.begin(), iter)];
                }
            }
        }

        replace_and_insert_ghost_elements(east_nodes_id_with_ghost, west_nodes_id_with_ghost, east_ghost_nodes_id, west_ghost_nodes_id);
    }

    std::vector< double > m_cartx_vert;
    std::vector< double > m_carty_vert;
    std::vector< double > m_cartz_vert;

    neighbours_of m_vertex_of_cell;
    neighbours_of m_vertices_of_vertex;
    unsigned int m_n_vertices_of_cell, m_n_vertices_of_vertex;
};

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
      ("ndims", po::value<int>(&nmpi)->default_value(3),  "number of dimension for the gmesh")
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

   unsigned int ndims = 3;
   if(vm.count("ndims"))
   {
       ndims = vm["ndims"].as<int>();
       if(ndims != 3 && ndims != 2) {
           std::cout << "Error: wrong number of dimensions" << std::endl;
           return 1;
       }
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

   int lon_vert_id, lat_vert_id;
   ierror= nc_inq_varid(ncid, "longitude_vertices", &lon_vert_id);
   assert(ierror==NC_NOERR);
   ierror= nc_inq_varid(ncid, "latitude_vertices", &lat_vert_id);
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

   assert(vertex_cell_dimlens[1] == cell_dimlen);

   int vertices_of_vertex_id;
   ierror= nc_inq_varid(ncid, "vertices_of_vertex", &vertices_of_vertex_id);
   assert(ierror==NC_NOERR);
   int vertices_of_vertex_numdim;
   ierror = nc_inq_varndims(ncid, vertices_of_vertex_id, &vertices_of_vertex_numdim);
   assert(ierror==NC_NOERR && vertices_of_vertex_numdim==2);

   int vertices_of_vertex_dimids[2];
   size_t vertices_of_vertex_dimlens[2];
   ierror = nc_inq_vardimid(ncid, vertices_of_vertex_id, vertices_of_vertex_dimids);
   assert(ierror==NC_NOERR);

   ierror = nc_inq_dimlen(ncid, vertices_of_vertex_dimids[0], &vertices_of_vertex_dimlens[0]);
   assert(ierror==NC_NOERR);
   ierror = nc_inq_dimlen(ncid, vertices_of_vertex_dimids[1], &vertices_of_vertex_dimlens[1]);

   assert(vertex_dimlen == vertices_of_vertex_dimlens[1]);

   assert(ierror==NC_NOERR);

   // extract data from netcdf

   mesh icon_mesh(vertex_cell_dimlens[1], vertices_of_vertex_dimlens[1], vertex_cell_dimlens[0], vertices_of_vertex_dimlens[0]);

   std::vector<double> cartx_vert, carty_vert, cartz_vert;

   cartx_vert.resize( vertex_dimlen );
   carty_vert.resize( vertex_dimlen );
   cartz_vert.resize( vertex_dimlen );

   if(ndims == 3) {

       ierror = nc_get_var_double(ncid, cartx_vert_id, &cartx_vert[0]);
       assert(ierror==NC_NOERR);
       ierror= nc_get_var_double(ncid, carty_vert_id, &carty_vert[0]);
       assert(ierror==NC_NOERR);
       ierror= nc_get_var_double(ncid, cartz_vert_id, &cartz_vert[0]);
       assert(ierror==NC_NOERR);
   }
   else {
       ierror = nc_get_var_double(ncid, lon_vert_id, &cartx_vert[0]);
       assert(ierror==NC_NOERR);

       ierror = nc_get_var_double(ncid, lat_vert_id, &carty_vert[0]);
       assert(ierror==NC_NOERR);

       for(size_t i=0; i < vertex_dimlen; ++i)
           cartz_vert[i] = 0;
   }
   icon_mesh.insert_cart_vert(cartx_vert, carty_vert, cartz_vert);

   icon_mesh.fill_from_netcdf(ncid, vertex_cell_id, vertices_of_vertex_id);

   if(ndims==2) {
       icon_mesh.add_ghost_nodes();
   }

   icon_mesh.gmsh_render(doLabels, nmpi);

   nc_close(ncid);
};
