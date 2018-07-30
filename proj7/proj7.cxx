/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include <bitset>
#include <cmath>

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "TriangleList.h"
#include "tricase.cxx"

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);
}

// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));
}

int GetCellPoint(int cellX, int cellY, int cellZ,
                 int edge, int* end1, int* end2)
{
  switch(edge)
  {
    case 0:
      end1[0] = cellX; end1[1] = cellY; end1[2] = cellZ;
      end2[0] = cellX + 1; end2[1] = cellY; end2[2] = cellZ;
    break;
    case 1:
      end1[0] = cellX + 1; end1[1] = cellY; end1[2] = cellZ;
      end2[0] = cellX + 1; end2[1] = cellY; end2[2] = cellZ + 1;
    break;
    case 2:
      end1[0] = cellX; end1[1] = cellY; end1[2] = cellZ + 1;
      end2[0] = cellX + 1; end2[1] = cellY; end2[2] = cellZ + 1;
    break;
    case 3:
      end1[0] = cellX; end1[1] = cellY; end1[2] = cellZ;
      end2[0] = cellX; end2[1] = cellY; end2[2] = cellZ + 1;
    break;
    case 4:
      end1[0] = cellX; end1[1] = cellY + 1; end1[2] = cellZ;
      end2[0] = cellX + 1; end2[1] = cellY + 1; end2[2] = cellZ;
    break;
    case 5:
      end1[0] = cellX + 1; end1[1] = cellY + 1; end1[2] = cellZ;
      end2[0] = cellX + 1; end2[1] = cellY + 1; end2[2] = cellZ + 1;
    break;
    case 6:
      end1[0] = cellX; end1[1] = cellY + 1; end1[2] = cellZ + 1;
      end2[0] = cellX + 1; end2[1] = cellY + 1; end2[2] = cellZ + 1;
    break;
    case 7:
      end1[0] = cellX; end1[1] = cellY + 1; end1[2] = cellZ;
      end2[0] = cellX; end2[1] = cellY + 1; end2[2] = cellZ + 1;
    break;
    case 8:
      end1[0] = cellX; end1[1] = cellY; end1[2] = cellZ;
      end2[0] = cellX; end2[1] = cellY + 1; end2[2] = cellZ;
    break;
    case 9:
      end1[0] = cellX + 1; end1[1] = cellY; end1[2] = cellZ;
      end2[0] = cellX + 1; end2[1] = cellY + 1; end2[2] = cellZ;
    break;
    case 10:
      end1[0] = cellX; end1[1] = cellY; end1[2] = cellZ + 1;
      end2[0] = cellX; end2[1] = cellY + 1; end2[2] = cellZ + 1;
    break;
    case 11:
      end1[0] = cellX + 1; end1[1] = cellY; end1[2] = cellZ + 1;
      end2[0] = cellX + 1; end2[1] = cellY + 1; end2[2] = cellZ + 1;
    break;
  }
  return 0;
}

int InterpolateForPoint(int* end1, int* end2,
                        float *xCoords, float *yCoords, float *zCoords,
                        int* dims, float* field, float isoValue, float* point)
{
  int index1 = GetPointIndex(end1, dims);
  int index2 = GetPointIndex(end2, dims);
  float proportion = (isoValue - field[index1]) / (field[index2] - field[index1]);
  point[0] = xCoords[end1[0]] + proportion * (xCoords[end2[0]] - xCoords[end1[0]]);
  point[1] = yCoords[end1[1]] + proportion * (yCoords[end2[1]] - yCoords[end1[1]]);
  point[2] = zCoords[end1[2]] + proportion * (zCoords[end2[2]] - zCoords[end1[2]]);
  return 0;
}

bool IsOutOfBounds(float* xCoords, float* yCoords, float* zCoords,
                   int* dims, float* point)
{
  if(point[0] < xCoords[0] || point[0] > xCoords[dims[0]-1]
     || point[1] < yCoords[0] || point[1] > yCoords[dims[1]-1]
     || point[2] < zCoords[0] || point[2] > zCoords[dims[2]-1])
    return true;
  return false;
}

int AddTrianglesForCase(TriangleList& triangleList,
                       int cellX, int cellY, int cellZ,
                       float* xCoords, float* yCoords, float* zCoords,
                       int* dims, float* field, int caseId, float isoValue)
{
  int edge1, edge2, edge3;
  int end1[3], end2[3];
  float point1[3], point2[3], point3[3];
  int triangleOffset = 0;
  //edge 1
  while(!(triCase[caseId][triangleOffset] == -1))
  {
    //Loop for all triangles
    edge1 = triCase[caseId][triangleOffset++];
    edge2 = triCase[caseId][triangleOffset++];
    edge3 = triCase[caseId][triangleOffset++];

    GetCellPoint(cellX, cellY, cellZ, edge1, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, zCoords, dims, field, isoValue, point1);

    GetCellPoint(cellX, cellY, cellZ, edge2, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, zCoords, dims, field, isoValue, point2);

    GetCellPoint(cellX, cellY, cellZ, edge3, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, zCoords, dims, field, isoValue, point3);

    bool _1check = IsOutOfBounds(xCoords, yCoords, zCoords, dims, point1);
    bool _2check = IsOutOfBounds(xCoords, yCoords, zCoords, dims, point2);
    bool _3check = IsOutOfBounds(xCoords, yCoords, zCoords, dims, point3);

    if(_1check || _2check || _3check)
      continue;

    triangleList.AddTriangle(point1[0], point1[1], point1[2],
                             point2[0], point2[1], point2[2],
                             point3[0], point3[1], point3[2]);
  }
  return 0;
}

int GenerateIsoLines(TriangleList& triangleList,
                     float* xCoords, float* yCoords, float* zCoords,
                     float* field, int* dims, float isoValue=3.2)
{
  // Iterate over all cells.
  for(int xind = 0; xind < dims[0] - 1; xind++)
    for(int yind = 0; yind < dims[1] - 1; yind++)
      for(int zind = 0; zind < dims[2] - 1; zind++)
      {
        int logicalPointIndex[3] = {xind, yind, zind};
        int vert0 = GetPointIndex(logicalPointIndex, dims);
        ++logicalPointIndex[0];
        int vert1 = GetPointIndex(logicalPointIndex, dims);
        ++logicalPointIndex[1];
        int vert5 = GetPointIndex(logicalPointIndex, dims);
        --logicalPointIndex[0];
        int vert4 = GetPointIndex(logicalPointIndex, dims);
        ++logicalPointIndex[2];
        int vert6 = GetPointIndex(logicalPointIndex, dims);
        --logicalPointIndex[1];
        int vert2 = GetPointIndex(logicalPointIndex, dims);
        ++logicalPointIndex[0];
        int vert3 = GetPointIndex(logicalPointIndex, dims);
        ++logicalPointIndex[1];
        int vert7 = GetPointIndex(logicalPointIndex, dims);
        // Gather the caseID.
        std::bitset<8> casebits;
        casebits[0] = field[vert0] > isoValue ? 1 : 0;
        casebits[1] = field[vert1] > isoValue ? 1 : 0;
        casebits[2] = field[vert2] > isoValue ? 1 : 0;
        casebits[3] = field[vert3] > isoValue ? 1 : 0;
        casebits[4] = field[vert4] > isoValue ? 1 : 0;
        casebits[5] = field[vert5] > isoValue ? 1 : 0;
        casebits[6] = field[vert6] > isoValue ? 1 : 0;
        casebits[7] = field[vert7] > isoValue ? 1 : 0;
        int caseId = (int)casebits.to_ulong();

        // Get the Segments that are affected.
        AddTrianglesForCase(triangleList, xind, yind, zind,
                            xCoords, yCoords, zCoords,
                            dims, field, caseId, isoValue);
      }
  return 0;
}

int main(int argc, char** argv)
{
    int  i, j;
    float isoValue = 3.2f;
    if(argc == 3)
      isoValue = atof(argv[2]);

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName(argc>=2?argv[1]:"proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    TriangleList triangleList;

    // YOUR CODE TO GENERATE ISOSURFACE SHOULD GO HERE!

    GenerateIsoLines(triangleList, X, Y, Z, F, dims, isoValue);

    vtkPolyData *pd = triangleList.MakePolyData();

    //This can be useful for debugging
    /*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
    */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
