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
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
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
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
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
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
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
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
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
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    //Two end points of the segment (x1,y1) and (x2,y2).
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

const int EdgeLookUp[16][4] =
{ {-1, -1, -1, -1}, //0
  {0, 3, -1, -1},   //1
  {0, 1, -1, -1},   //2
  {3, 1, -1, -1},   //3
  {2, 3, -1, -1},   //4
  {0, 2, -1, -1},   //5
  {0, 1, 2, 3},     //6
  {1, 2 , -1, -1},  //7
  {1, 2, -1, -1},   //8
  {0, 3, 1, 2},     //9
  {0, 2, -1, -1},   //10
  {2, 3, -1, -1},   //11
  {1, 3, -1, -1},   //12
  {0, 1, -1, -1},   //13
  {0, 3, -1, -1},   //14
  { -1, -1, -1, -1} //15
};

int GetCellPoint(int cellX, int cellY, int edge, int* end1, int* end2)
{
  switch(edge)
  {
    case 0:
      end1[0] = cellX; end1[1] = cellY;
      end2[0] = cellX + 1; end2[1] = cellY;
    break;
    case 1:
      end1[0] = cellX + 1; end1[1] = cellY;
      end2[0] = cellX + 1; end2[1] = cellY + 1;
    break;
    case 2:
      end1[0] = cellX; end1[1] = cellY + 1;
      end2[0] = cellX + 1; end2[1] = cellY + 1;
    break;
    case 3:
      end1[0] = cellX; end1[1] = cellY;
      end2[0] = cellX; end2[1] = cellY + 1;
    break;
  }
  return 0;
}

int InterpolateForPoint(int* end1, int* end2,
                        float *xCoords, float *yCoords,
                        int* dims, float* field, float isoValue, float* point)
{
  int index1 = GetPointIndex(end1, dims);
  int index2 = GetPointIndex(end2, dims);
  float proportion = (isoValue - field[index1]) / (field[index2] - field[index1]);
  point[0] = xCoords[end1[0]] + proportion * (xCoords[end2[0]] - xCoords[end1[0]]);
  point[1] = yCoords[end1[1]] + proportion * (yCoords[end2[1]] - yCoords[end1[1]]);
  return 0;
}

int AddSegmentsForCase(SegmentList& segmentList,
                       int cellX, int cellY,
                       float* xCoords, float* yCoords, int* dims,
                       float* field, int caseId, float isoValue)
{
  int edge1, edge2;
  int end1[2], end2[2];
  float point1[2], point2[2];
  //edge 1
  if(!(EdgeLookUp[caseId][0] == -1))
  {
    edge1 = EdgeLookUp[caseId][0];
    edge2 = EdgeLookUp[caseId][1];
    GetCellPoint(cellX, cellY, edge1, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, dims, field, isoValue, point1);
    GetCellPoint(cellX, cellY, edge2, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, dims, field, isoValue, point2);
    segmentList.AddSegment(point1[0], point1[1], point2[0], point2[1]);
  }
  //edge 2
  if(!(EdgeLookUp[caseId][2] == -1))
  {
    edge1 = EdgeLookUp[caseId][2];
    edge2 = EdgeLookUp[caseId][3];
    GetCellPoint(cellX, cellY, edge1, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, dims, field, isoValue, point1);
    GetCellPoint(cellX, cellY, edge2, end1, end2);
    InterpolateForPoint(end1, end2, xCoords, yCoords, dims, field, isoValue, point2);
    segmentList.AddSegment(point1[0], point1[1], point2[0], point2[1]);
  }
  return 0;
}

int GenerateIsoLines(SegmentList& segmentList,
                     float* xCoords, float* yCoords,
                     float* field, int* dims, float isoValue=3.2)
{
  // Iterate over all cells.
  for(int xind = 0; xind < dims[0] - 1; xind++)
    for(int yind = 0; yind < dims[1] - 1; yind++)
    {
      // Get the vertices that we need.

      /*     E2
         2 ------ 3
         |        |
      E3 |        | E1   We have point 1, point 2 is +1 in X.
         |   E0   |      Point 3 is + 1 in Y and Point 4 +1 in X and Y.
         0 ------ 1
      */
      int logicalPointIndex[2] = {xind, yind};
      int vert0 = GetPointIndex(logicalPointIndex, dims);
      ++logicalPointIndex[0];
      int vert1 = GetPointIndex(logicalPointIndex, dims);
      ++logicalPointIndex[1];
      int vert3 = GetPointIndex(logicalPointIndex, dims);
      --logicalPointIndex[0];
      int vert2 = GetPointIndex(logicalPointIndex, dims);

      // Gather the caseID.
      std::bitset<4> casebits;
      casebits[0] = field[vert0] > isoValue ? 1 : 0;
      casebits[1] = field[vert1] > isoValue ? 1 : 0;
      casebits[2] = field[vert2] > isoValue ? 1 : 0;
      casebits[3] = field[vert3] > isoValue ? 1 : 0;
      int caseId = (int)casebits.to_ulong();
      // Get the Segments that are affected.
      AddSegmentsForCase(segmentList, xind, yind,
                         xCoords, yCoords, dims, field, caseId, isoValue);
    }
  return 0;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!

    GenerateIsoLines(sl, X, Y, F, dims);

    vtkPolyData *pd = sl.MakePolyData();

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
