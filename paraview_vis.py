# trace generated using paraview version 5.9.0

import glob, os
vtus = glob.glob("*.vtu")

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

meshes = []
for vtu in vtus:
	# create a new 'XML Unstructured Grid Reader'
	mesh = XMLUnstructuredGridReader(registrationName=vtu, FileName=[vtu])
	mesh.PointArrayStatus = ['q']

	mesh.TimeArray = 'None'
	meshes.append(mesh)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=meshes)

# show data in view
groupDatasets1Display = Show(groupDatasets1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
groupDatasets1Display.Representation = 'Surface'
groupDatasets1Display.ColorArrayName = [None, '']
groupDatasets1Display.SelectTCoordArray = 'None'
groupDatasets1Display.SelectNormalArray = 'None'
groupDatasets1Display.SelectTangentArray = 'None'
groupDatasets1Display.OSPRayScaleArray = 'q'
groupDatasets1Display.OSPRayScaleFunction = 'PiecewiseFunction'
groupDatasets1Display.SelectOrientationVectors = 'None'
groupDatasets1Display.SelectScaleArray = 'None'
groupDatasets1Display.GlyphType = 'Arrow'
groupDatasets1Display.GlyphTableIndexArray = 'None'
groupDatasets1Display.GaussianRadius = 0.05
groupDatasets1Display.SetScaleArray = ['POINTS', 'q']
groupDatasets1Display.ScaleTransferFunction = 'PiecewiseFunction'
groupDatasets1Display.OpacityArray = ['POINTS', 'q']
groupDatasets1Display.OpacityTransferFunction = 'PiecewiseFunction'
groupDatasets1Display.DataAxesGrid = 'GridAxesRepresentation'
groupDatasets1Display.PolarAxes = 'PolarAxesRepresentation'
groupDatasets1Display.ScalarOpacityUnitDistance = 0.1919383103666485
groupDatasets1Display.OpacityArrayName = ['POINTS', 'q']
groupDatasets1Display.ExtractedBlockIndex = 1

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(groupDatasets1Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
groupDatasets1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# set scalar coloring
ColorBy(groupDatasets1Display, ('POINTS', 'q'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
groupDatasets1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
groupDatasets1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'q'
qLUT = GetColorTransferFunction('q')

# get opacity transfer function/opacity map for 'q'
qPWF = GetOpacityTransferFunction('q')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(896, 1171)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.8799999952316284, 0.0, 10000.0]
renderView1.CameraFocalPoint = [0.8799999952316284, 0.0, 0.0]
renderView1.CameraParallelScale = 8.46722565490744

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# input()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
