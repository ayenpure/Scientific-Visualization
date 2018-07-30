#include <string>

#include "vtkArrowSource.h"
#include "vtkLineSource.h"
#include "vtkStreamTracer.h"
#include "vtkExtractRectilinearGrid.h"
#include "vtkGlyph3D.h"
#include "vtkActor.h"
#include "vtkAlgorithmOutput.h"
#include "vtkAppendPolyData.h"
#include "vtkContourFilter.h"
#include "vtkCutter.h"
#include "vtkDataSetReader.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRectilinearGrid.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

int main(int argc, char **argv) {
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(800, 800);
  vtkSmartPointer<vtkRenderWindowInteractor> renderInteract =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderInteract->SetRenderWindow(renderWindow);

  // Define viewport ranges
  double xmins[4] = {0, .5, 0, .5};
  double xmaxs[4] = {0.5, 1, 0.5, 1};
  double ymins[4] = {0, 0, .5, .5};
  double ymaxs[4] = {0.5, 0.5, 1, 1};

  // If need be, parameterize this using command line args.
  std::string filename("proj8.vtk");
  vtkDataSetReader *reader = vtkDataSetReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkRectilinearGrid *mesh = (vtkRectilinearGrid *)reader->GetOutput();
  double *range = mesh->GetPointData()->GetScalars()->GetRange();
  mesh->GetPointData()->SetActiveAttribute("grad", vtkDataSetAttributes::VECTORS);
  mesh->GetPointData()->SetActiveAttribute("hardyglobal", vtkDataSetAttributes::SCALARS);


  // Read input file here.
  for (unsigned i = 0; i < 4; i++) {
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow->AddRenderer(renderer);
    renderer->SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i]);
    // Create a mapper and actor.
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkPolyDataAlgorithm> filter = nullptr;

    if (i == 0) {

      vtkSmartPointer<vtkContourFilter> contour =
          vtkSmartPointer<vtkContourFilter>::New();
      contour->SetInputConnection(reader->GetOutputPort());
      contour->SetValue(0, 2.5);
      contour->SetValue(1, 5.0);

      filter = contour;
    } else if (i == 1) {

      vtkSmartPointer<vtkExtractRectilinearGrid> extractFilter = vtkSmartPointer<vtkExtractRectilinearGrid>::New();
      extractFilter->SetInputData(mesh);
      extractFilter->SetSampleRate(5, 5, 5);
      //extractFilter->IncludeBoundaryOn();
      extractFilter->Update();
      vtkRectilinearGrid* subMesh = extractFilter->GetOutput();

      vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
      /*arrow->SetGlyphTypeToArrow();
      arrow->SetScale(0.3);*/

      vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
      //glyph->SetInputConnection(reader->GetOutputPort());
      glyph->SetInputData(subMesh);
      glyph->SetSourceConnection(arrow->GetOutputPort());
      glyph->ScalingOn();
      //glyph->OrientOn();
      glyph->SetVectorModeToUseVector();
      /*vtkSmartPointer<vtkHedgeHog> hedgehog = vtkSmartPointer<vtkHedgeHog>::New();
      hedgehog->SetVectorModeToUseVector();
      hedgehog->SetInputData(subMesh);
      hedgehog->SetScaleFactor(5.0);*/
      filter = glyph;

    } else if (i == 2) {

      vtkSmartPointer<vtkAppendPolyData> appendCuts =
          vtkSmartPointer<vtkAppendPolyData>::New();

      vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
      plane1->SetOrigin(0, 0, 0);
      plane1->SetNormal(1, 0, 0);
      vtkSmartPointer<vtkCutter> cutter1 = vtkSmartPointer<vtkCutter>::New();
      cutter1->SetInputConnection(reader->GetOutputPort());
      cutter1->SetCutFunction(plane1);
      appendCuts->AddInputConnection(cutter1->GetOutputPort());

      vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
      plane2->SetOrigin(0, 0, 0);
      plane2->SetNormal(0, 1, 0);
      vtkSmartPointer<vtkCutter> cutter2 = vtkSmartPointer<vtkCutter>::New();
      cutter2->SetInputConnection(reader->GetOutputPort());
      cutter2->SetCutFunction(plane2);
      appendCuts->AddInputConnection(cutter2->GetOutputPort());

      vtkSmartPointer<vtkPlane> plane3 = vtkSmartPointer<vtkPlane>::New();
      plane3->SetOrigin(0, 0, 0);
      plane3->SetNormal(0, 0, 1);
      vtkSmartPointer<vtkCutter> cutter3 = vtkSmartPointer<vtkCutter>::New();
      cutter3->SetInputConnection(reader->GetOutputPort());
      cutter3->SetCutFunction(plane3);
      appendCuts->AddInputConnection(cutter3->GetOutputPort());

      filter = appendCuts;
    } else {
      vtkSmartPointer<vtkLineSource> seeds = vtkSmartPointer<vtkLineSource>::New();
      seeds->SetPoint1(-9, 0, 0);
      seeds->SetPoint2(9, 0, 0);
      seeds->SetResolution(18);
      vtkSmartPointer<vtkStreamTracer> streamTracer = vtkSmartPointer<vtkStreamTracer>::New();
      streamTracer->SetIntegratorTypeToRungeKutta4();
      streamTracer->SetInputData(mesh);
      streamTracer->SetSourceConnection(seeds->GetOutputPort());
      streamTracer->SetMaximumPropagation(100);
      streamTracer->SetInitialIntegrationStep(.2);

      filter = streamTracer;
    }

    if (filter == nullptr) {
      std::cout << "Error in retrieving data to render" << std::endl;
      return EXIT_FAILURE;
    }

    mapper->SetInputConnection(filter->GetOutputPort());
    mapper->SetScalarRange(range[0], range[1]);
    actor->SetMapper(mapper);
    renderer->AddActor(actor);
    renderer->ResetCamera();

  }

  renderWindow->Render();
  renderWindow->SetWindowName("Multiple ViewPorts");

  renderInteract->Start();
  return EXIT_SUCCESS;
}
