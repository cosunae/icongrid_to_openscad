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

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCubeSource.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkActor2D.h>
#include <vtkIdFilter.h>
#include <vtkLabeledDataMapper.h>
#include <vtkCellCenters.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkTextProperty.h>

#define MAX_NDIMS 40


void vtk_render(int ncells_per_diamond, std::vector<double>& cartx_vert, std::vector<double>& carty_vert, std::vector<double>& cartz_vert, std::vector<double>& vertex_cell, int ndiamonds_to_label=0)
{
    assert( cartx_vert.size() == carty_vert.size() && cartx_vert.size() == cartz_vert.size() );
    assert(vertex_cell.size()%3==0);


    // Create a cube for coloring cells
    vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();

    // Setup vertex points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


    for(int i=0; i < cartx_vert.size(); ++i){
        points->InsertNextPoint(cartx_vert[i], carty_vert[i], cartz_vert[i] );
    }

    // Setup vertex points for labels
    vtkSmartPointer<vtkPoints> pointsLabels = vtkSmartPointer<vtkPoints>::New();

    for(int i=0; i < ncells_per_diamond*ndiamonds_to_label; ++i){
        pointsLabels->InsertNextPoint(cartx_vert[i], carty_vert[i], cartz_vert[i] );
    }

    int ncells = vertex_cell.size()/3;

    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    // Add the polygon to a list of polygons
    vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();

    std::vector<vtkIdType> lineId(ncells*3*2);
    for(int i=0; i < ncells; ++i)
    {
        lineId[i*3*2+0] = vertex_cell[0*ncells + i]-1;
        lineId[i*3*2+1] = vertex_cell[1*ncells + i]-1;
        lineId[i*3*2+2] = vertex_cell[1*ncells + i]-1;
        lineId[i*3*2+3] = vertex_cell[2*ncells + i]-1;
        lineId[i*3*2+4] = vertex_cell[2*ncells + i]-1;
        lineId[i*3*2+5] = vertex_cell[0*ncells + i]-1;

    }
    for(int i=0; i < ncells; ++i)
    {
        lines->InsertNextCell(2, &(lineId[i*6]));
        lines->InsertNextCell(2, &(lineId[i*6+2]));
        lines->InsertNextCell(2, &(lineId[i*6+4]));
    }

    std::vector<vtkIdType> faces(ncells*3);
    for(int i=0; i < ncells; ++i)
    {
        faces[i*3+0]=vertex_cell[0*ncells + i]-1;
        faces[i*3+1]=vertex_cell[1*ncells + i]-1;
        faces[i*3+2]=vertex_cell[2*ncells + i]-1;
    }

    for(int i=0; i < ncells; ++i)
    {
      polygons->InsertNextCell(3, &(faces[i*3]));
    }

    vtkSmartPointer<vtkCellArray> polygonsLabels = vtkSmartPointer<vtkCellArray>::New();

    for(int i=0; i < ncells/20; ++i)
    {
          polygonsLabels->InsertNextCell(3, &(faces[i*3]));
    }

    // Create a PolyData to draw faces
    vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(polygons);

    // Create a PolyData to draw edges
    vtkSmartPointer<vtkPolyData> polyDataLines =  vtkSmartPointer<vtkPolyData>::New();
    polyDataLines->SetPoints(points);
    polyDataLines->SetLines(lines);
   
   
    // Create a PolyData to draw labels
    vtkSmartPointer<vtkPolyData> polyDataLabels =vtkSmartPointer<vtkPolyData>::New();
    polyDataLabels->SetPoints(points);
    polyDataLabels->SetPolys(polygonsLabels);
    
    vtkUnsignedCharArray *faceColors = vtkUnsignedCharArray::New(); 
    faceColors->SetNumberOfComponents(3); 
    // Setup colors
    unsigned char colors[10][3] ={
        {220,220,220}, // gainsboro
        {227, 207, 87}, // Banana
        {255, 99, 71}, // Tomato
        {255, 99, 71}, // Wheat
        {230, 230, 250}, // Lavender
        {255, 125, 64}, // Flesh
        {135, 38, 89}, // Raspberry
        {250, 128, 114}, // Salmon
        {189, 255, 201}, // Mint
        {51, 161, 201} // Peacock
    };

    unsigned char gainsboro[3] = {220,220,220};
    unsigned char white[3] = {255, 255, 255};
    unsigned char black[3] = {0, 0, 0};
    unsigned char banana[3] = {227, 207, 87}; // Banana
    unsigned char tomato[3] = {255, 99, 71}; // Tomato
    unsigned char wheat[3] = {245, 222, 179}; // Wheat
    unsigned char lavender[3] = {230, 230, 250}; // Lavender
    unsigned char flesh[3]= {255, 125, 64}; // Flesh
    unsigned char raspberry[3] = {135, 38, 89}; // Raspberry
    unsigned char salmon[3] = {250, 128, 114}; // Salmon
    unsigned char mint[3] = {189, 255, 201}; // Mint
    unsigned char peacock[3] = {51, 161, 201}; // Peacock
                    
    int color_cnt=0;

    for(int j=0; j < ncells/ncells_per_diamond; ++j) {
        for(int i=0 ; i < ncells_per_diamond ; ++i) {
            faceColors->InsertNextTupleValue(colors[color_cnt%10]);
        }
        color_cnt++;
    }
   
    cubeSource->GetOutput()->DeepCopy(polygonPolyData);                                   
    polygonPolyData->GetCellData()->SetScalars(faceColors);
    
// Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(polygonPolyData);
#else
    mapper->SetInputData(polygonPolyData);
#endif

    vtkSmartPointer<vtkPolyDataMapper> mapperLines =    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapperLines->SetInput(polyDataLines);
#else
    mapperLines->SetInputData(polyDataLines);
#endif

  // Generate data arrays containing point and cell ids
    vtkSmartPointer<vtkIdFilter> ids =
    vtkSmartPointer<vtkIdFilter>::New();
    ids->SetInputData( polyDataLabels );
    ids->PointIdsOn();
    ids->CellIdsOn();
    ids->FieldDataOn();

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
 // Create labels for points
//  visPts = vtkSmartPointer<vtkSelectVisiblePoints>::New();
//  visPts->SetInputConnection( ids->GetOutputPort() );
//  visPts->SetRenderer( ren1 );
//  visPts->SelectionWindowOn();
//  visPts->SetSelection( xmin, xmin + xLength, ymin, ymin + yLength );


  // Create the mapper to display the point ids.  Specify the
    // format to use for the labels.  Also create the associated actor.
      vtkSmartPointer<vtkLabeledDataMapper> ldm = vtkSmartPointer<vtkLabeledDataMapper>::New();
    ldm->SetInputConnection( ids->GetOutputPort() );
    ldm->SetLabelModeToLabelFieldData();


  // Create labels for cells
    vtkSmartPointer<vtkCellCenters> cc =   vtkSmartPointer<vtkCellCenters>::New();
    cc->SetInputConnection( ids->GetOutputPort() );
           
    vtkSmartPointer<vtkSelectVisiblePoints> visCells = vtkSmartPointer<vtkSelectVisiblePoints>::New();
    visCells->SetInputConnection( cc->GetOutputPort() );
    visCells->SetRenderer( renderer );
   // Create the mapper to display the cell ids.  Specify the
     // format to use for the labels.  Also create the associated actor.
    vtkSmartPointer<vtkLabeledDataMapper> cellMapper =  vtkSmartPointer<vtkLabeledDataMapper>::New();
    cellMapper->SetInputConnection( visCells->GetOutputPort() );
    cellMapper->SetLabelModeToLabelFieldData();
    cellMapper->GetLabelTextProperty()->SetColor( 0, 0.4, 0 );


    vtkSmartPointer<vtkActor2D> pointLabels = vtkSmartPointer<vtkActor2D>::New();
    pointLabels->SetMapper( ldm );

    vtkSmartPointer<vtkActor2D> cellLabels = vtkSmartPointer<vtkActor2D>::New();
    cellLabels->SetMapper( cellMapper );

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
 

    vtkSmartPointer<vtkActor> actorLines = vtkSmartPointer<vtkActor>::New();
    actorLines->SetMapper(mapperLines);


    // Visualize
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
 
    renderer->AddActor(actor);
    renderer->AddActor(actorLines);
//  renderer->AddActor2D(pointLabels);
    renderer->AddActor2D(cellLabels);
//    renderer->SetBackground(.5,.3,.31); // Background color salmon
    renderer->SetBackground(0.36, 0.36, 0.36);
 
    renderWindow->Render();
    renderWindowInteractor->Start();

}

int main(int argc, char** argv)
{  

   if(argc < 2) {
     printf("Error, need to specify input file");
     exit(1);
   }
   std::string ncfilename = argv[1];

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

   vtk_render(cell_dimlen/20, cartx_vert, carty_vert, cartz_vert, vertex_cell);


   nc_close(ncid);
};
