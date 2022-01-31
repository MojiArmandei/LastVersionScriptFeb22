#importing main packages
from abaqus import *
from abaqusConstants import *
import regionToolset
import numpy as np 
import mesh


##############################
## Checking errors due to Compatability of Commands
backwardCompatibility.setValues(includeDeprecated=True, reportDeprecated=False)
#### This command forces Abaqus to generate Reply file with coordinates and not only index of geometry entity
cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)""")

#######################################
# JUST TO RESET ALL PREVIOUS MODELS
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
###############################################
cae_no_parts_input_file=ON 
##############################################################
####### Most important variabls => Independent Variables 
#############################################################
#############################################
#### Geometric Nonlinearity
 
##### Large & Small Deformation

Non_Lin_Geom = 'ON'
# Non_Lin_Geom = 'OFF'

###  0) Radius of Pipe (for D/thick)
##### EXTERNAL RADIUS OF PIPE
radius = 103
# radius = 103 * 2
# radius = 103 * 3
# radius = 103 * 4
# radius = 103 * 5

###############################
########################
LoadingType = 'Bending'
# LoadingType = 'Int_Pressure'
# LoadingType = 'Tensile' ### Tensile applied by Pressure/
# LoadingType = 'Torsion'
# LoadingType = 'ALL'

Applied_Bending = 0.1
Applied_Int_Pressure = 10
Applied_Tensile_Pressure = 5
Applied_Torsion = 0.1

#### 1 )  NUMBER OF CPU USED FOR ANALYSIS (Depends on your computer)
NCPUs = 6

#### 2)  SELECT TYPE OF THE MATERIAL DEFINITION USED IN THE MODEL
MaterialType = 'Plastic-Deformation'  # DEFAULT
# MaterialType = 'Elastic-Plastic'




# a= 2.06 ### MAXI 2.07 for clad 3

radius = 103.0
thick = 20.6
clad = 3.0
# clad = 2.0

a= clad* 0.69
A1 = a/45
A2 = 2* A1
A3 = 3 * A1
A4 = 4 * A1
A5 = 5 * A1
A6 = 6 * A1
A7 = 7 * A1
A8 = 8 * A1
A9 = 9 * A1
A10 = 10 * A1
lbox = 15 * A1
h3 = (thick-lbox)/2 + lbox
Wedge = a/450

Epsilon = A1/2

theta_pi = 0.04
# theta_pi = 0.12
# theta_pi = 0.20


D_ext = radius * 2
L = 3 * D_ext 
C = theta_pi * 2 * D_ext
LxPositive=((2*pi*radius)/2.0)-(thick+C)
LDiv1 = thick

LDiv2 = thick + L * 0.2
LDiv3 = thick + L * 0.6 

Weld_Angle = 0
# Weld_Angle = 30


NseedSpider = 3
Nseedcurve = 3 
NseedC = 12 
Nseeda = 3
Nseedlbox = 4
Nseedthick = 5
Nseeedclad = 3
NseeedLDiv1 = 7
NseeedLDiv2 = 10
NseeedLDiv3 = 9
NseeedL = 7
NseedthickZ = 6
NseedLxPos = 8

###############################################################
####### Materials  
################################################################

## 1) BaseSteel

BaseSteelDens = 7800.0
BaseSteelE= 2.06E5
BaseSteelPratio= 0.3
N_Base = 5 
Yield_Stress_Base= 275


Yield_Strain_Base = Yield_Stress_Base / BaseSteelE
Alpha= 1.0

SS_Curve_Base = [[Yield_Stress_Base,0.0],
                [Yield_Stress_Base*1.1,Yield_Stress_Base]
                ]


##2)  WeldSteel
WeldSteelDens= 7800.0
WeldSteelE=2.06E5
WeldSteelPratio= 0.3
N_Weld = 5 
Yield_Stress_Weld= 275 

Yield_Strain_Weld= Yield_Stress_Weld / WeldSteelE
Alpha=1.0

SS_Curve_Weld = [[Yield_Strain_Weld,0.0],
                 [Yield_Strain_Weld,Yield_Stress_Weld]
                 ]

##3)  CladSteelWeld
CladSteelWeldDens= 7800.0
CladSteelWeldE=2.06E5
CladSteelWeldPratio= 0.3
N_CladWeld = 5 
Yield_Stress_CladSteelWeld= 275 


Yield_Strain_CladSteelWeld = Yield_Stress_CladSteelWeld / CladSteelWeldE
Alpha=1.0

SS_Curve_CladWeld = [[Yield_Strain_CladSteelWeld,0.0],
                 [Yield_Strain_Weld,Yield_Stress_CladSteelWeld]
                ]


##4)  CladSteelBase
CladSteelBaseDens= 7800.0
CladSteelBaseE=2.06E5
CladSteelBasePratio= 0.3
N_clad = 5 
Yield_Stress_CladSteelBase= 275 

Yield_Strain_CladSteelBase = Yield_Stress_CladSteelBase / CladSteelBaseE
N_CladSteelBase= N_Base
Alpha=1.0

SS_Curve_Clad = [[Yield_Strain_CladSteelBase,0.0],
                 [Yield_Strain_CladSteelBase,Yield_Stress_CladSteelBase]
                ]


#---------
#create the model
mdb.models.changeKey(fromName='Model-1',toName='Pipe_Model_Internal')
PipeModel=mdb.models['Pipe_Model_Internal']

#---------
#Create the part
import sketch
import part

# a) sketch the pipe plate before rolling using rectangle tool
Orig=[0.0,0.0,0.0]
PipeSketch=PipeModel.ConstrainedSketch(name='Pipe Plate Sketch', sheetSize=100.0)
PipeSketch.setPrimaryObject(option=STANDALONE)
PipeSketch.rectangle(point1=(Orig[0], Orig[1]), point2=(thick, thick)) 

#b) creating the Main part
PipePart = PipeModel.Part(name='PIPE', dimensionality=THREE_D,type=DEFORMABLE_BODY)
PipePart = PipeModel.parts['PIPE']
PipePart.BaseSolidExtrude(sketch=PipeSketch, depth= thick )
PipeSketch.unsetPrimaryObject()

Faces=PipePart.faces
Edges=PipePart.edges
Cells=PipePart.cells

# a-1) make a sketch tranformation plane in YZ plane

Face_Sketch_plane=Faces.findAt(coordinates=(thick/2, thick/2, 0.0))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, thick/2, 0.0))
Origin_Sketch_plane=Orig

##sketch of crack face
PipePart = PipeModel.parts['PIPE']
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

sketch_1 = PipeModel.ConstrainedSketch(name='PIPESketch',   sheetSize=47.7, gridSpacing=1.19, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.ArcByCenterEnds(center=(0.0,thick), point1=(0.0,thick-a), point2=(-a,thick), direction=CLOCKWISE)

pickedFaces=Faces.findAt(coordinates=(thick/2, thick/2, 0.0))
Edge_For_PartitionBySketch=Edges.findAt(coordinates=(0.0, clad, 0.0))

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_PartitionBySketch, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()


 #############################################################
# ##### Sketching Spider Mesh at crack face
# ############################################################
Face_Sketch_plane=Faces.findAt(coordinates=(0.0, thick/2, thick/2))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, Epsilon, 0.0))
Origin_Sketch_plane= (0.0, 0.0, 0.0)


#sketch of crack face
PipePart = PipeModel.parts['PIPE']
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

sketch_1 = PipeModel.ConstrainedSketch(name='PIPESketch',   sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.ArcByCenterEnds(center=(0.0, thick- a), point1=(0.0, thick-a-A1), point2=(0.0,thick-a+A1), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A2), point2=(0.0,thick-a+A2), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A3), point2=(0.0,thick-a+A3), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A4), point2=(0.0,thick-a+A4), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A5), point2=(0.0,thick-a+A5), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A6), point2=(0.0,thick-a+A6), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A7), point2=(0.0,thick-a+A7), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A8), point2=(0.0,thick-a+A8), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A9), point2=(0.0,thick-a+A9), direction=COUNTERCLOCKWISE)
sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-A10), point2=(0.0,thick-a+A10), direction=COUNTERCLOCKWISE)


PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=Face_Sketch_plane, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['PIPESketch']

################# Sketching clad

PipePart = PipeModel.parts['PIPE']
Faces = PipePart.faces
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick*0.1, thick*0.01))
Edge_For_Sketch =  Edges.findAt(coordinates=(0.0, thick*0.3, thick))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, origin=(0.0, 0.0, thick))
sketch = PipeModel.ConstrainedSketch(name='__profile__',sheetSize=58.26, gridSpacing=1.45, transform=Transform)
sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch, filter=COPLANAR_EDGES)

sketch.Line(point1=(-thick, thick-clad), point2=(0.0, thick-clad))

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch,  faces=Face_Sketch_plane, sketch=sketch)
sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

######## Extruding Clad partion

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick * 0.1, thick * 0.1, thick * 0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick - clad, thick*0.3)), 
    )
ExtrudeLine = Edges.findAt(coordinates=(thick*0.3, thick, thick))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)



##### Sketching vertical clad
PipePart = PipeModel.parts['PIPE']
Faces = PipePart.faces
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(thick*0.1, thick*0.1, 0.0))
Edge_For_Sketch =  Edges.findAt(coordinates=(0.0, thick*0.3, 0.0))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
sketch = PipeModel.ConstrainedSketch(name='__profile__',sheetSize=71.36, gridSpacing=1.78, transform=Transform)
sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch, filter=COPLANAR_EDGES)

sketch.Line(point1=(-(clad), thick), point2=(-(clad), 0.0))

pickedFaces = Faces.findAt(
    ((thick*0.1, thick*0.1, 0.0), ), 
    ((thick*0.1, thick*0.99, 0.0), )
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch,  faces=pickedFaces, sketch=sketch)
sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

##### Extruding vertical clad line

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.2), ),
    ((thick*0.01, thick*0.1, thick*0.2), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(clad, thick*0.99, 0.0)), 
    Edges.findAt(coordinates=(clad, thick*0.09, 0.0))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, 0.0, thick*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


###### Extruding spider pattern curves at crack face
## A1
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A1)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A2
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A2)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A3
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A3)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)


## A4
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A4)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A5
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A5)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

### A6
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A6)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

# A7
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A7)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A8
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A8)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A9
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A9)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)

## A10
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((thick*0.01, thick*0.99, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a,  A10)), 
    )
SweepLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
PipePart.PartitionCellBySweepEdge(sweepPath=SweepLine, cells=pickedCells, edges=pickedEdges)


########## extrude crack curve 

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(
    ((a*cos(pi/4), thick-a*cos(pi/4), A1*0.6), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A1*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A2*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A3*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A4*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A5*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A6*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A7*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A8*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A9*1.1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A10*1.1), ),
    )

pickedEdges =(
    Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0)), 
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, thick, thick*0.3))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

############### Sketching thick half partition

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(thick*0.9, thick*0.1, 0.0))
Edge_For_Sketch = Edges.findAt(coordinates=(0.0, thick*0.3, 0.0))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,  sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))

Sketch = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=55.5, gridSpacing=1.38, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(-thick, (thick-clad)/2), point2=(0.0, (thick-clad)/2))

pickedFaces = Faces.findAt(
    ((thick*0.9, thick*0.1, 0.0), ),
     ((thick*0.01, thick*0.1, 0.0), )
     )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch,  faces=pickedFaces, sketch=Sketch)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

#### Extrusion of thick half partition
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.findAt(
    ((thick*0.9, thick*0.3, thick*0.1), ),
    ((thick*0.01, thick*0.3, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(thick*0.9, (thick-clad)/2, 0.0)),
    Edges.findAt(coordinates=(thick*0.01, (thick-clad)/2, 0.0))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, 0.0, thick*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)



##### sketching lbox
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick*0.1,  thick*0.2))
Edge_For_Sketch = Edges.findAt(coordinates=(0.0, thick*0.3, 0.0))

Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))

Sketch = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=69.12, gridSpacing=1.72, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart = PipeModel.parts['PIPE']
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(lbox, thick), point2=(lbox, 0.0))

pickedFaces = Faces.findAt(
    ((0.0, thick*0.1,  thick*0.2), ),
    ((0.0, thick*0.65,  thick*0.2), ),
    ((0.0, thick-a*1.1,  thick*0.2), ),
    ((0.0, thick-a*0.91,  thick*0.2), ),
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=Sketch)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


######## Extrude lbox

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces
pickedCells = Cells.findAt(
    ((thick*0.01, thick-a*0.1, lbox*0.99), ),
    ((thick*0.01, thick-a*1.01, lbox*0.99), ),
    ((thick*0.9, thick-a*1.01, lbox*0.99), ),
    ((thick*0.9, thick*0.7, thick*0.99), ),
    ((thick*0.01, thick*0.7, thick*0.99), ),
    ((thick*0.9, thick*0.1, thick*0.99), ),
    ((thick*0.01, thick*0.1, thick*0.99), ),
  
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a*0.1, lbox)),
    Edges.findAt(coordinates=(0.0,thick-a*1.01, lbox)),
    Edges.findAt(coordinates=(0.0, thick*0.7,lbox)), 
    Edges.findAt(coordinates=(0.0, thick*0.1, lbox))

    )
ExtrudeLine = Edges.findAt(coordinates=(a*1.1, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=REVERSE)

#######################################
####### sketching V-groove 
if Weld_Angle > 0:   
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges
    Faces = PipePart.faces
    THETA = Weld_Angle/2 * pi/180
    RR = thick/cos(THETA)
    rr = RR * sin(THETA)
    Face_Sketch_plane = Faces.findAt(coordinates=(thick, thick*0.1, thick*0.3))
    Edge_For_Sketch = Edges.findAt(coordinates=(thick, thick*0.1, lbox))
    Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,sketchPlaneSide=SIDE1, origin=(thick, 0.0, thick))
    sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=60.86, gridSpacing=1.52, transform=Transform)
    sketch_1.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)
    sketch_1.Line(point1=(thick-lbox, thick), point2=(thick-lbox-rr, 0.0))
    pickedFaces = Faces.findAt(
        ((thick, thick*0.1, thick*0.3), ), 
        ((thick, 1.1*(thick-clad)/2, thick*0.3), ),
        ((thick, thick-clad*0.4, thick*0.3), )
        )

    PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=sketch_1)
    sketch_1.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']

##### Extruding V-groove line 
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells

    pickedCells = Cells.findAt(
        ((thick*0.01, thick-a*0.1, lbox*1.1), ),
        ((thick*0.01, thick-a*1.01, lbox*1.1), ),
        ((thick*0.9, thick-a*1.01, lbox*1.1), ),
        ((thick*0.9, thick*0.7, thick*0.9), ),
        ((thick*0.01, thick*0.7, thick*0.9), ),
        ((thick*0.9, thick*0.1, thick*0.9), ),
        ((thick*0.01, thick*0.1, thick*0.9), ),
        )
    pickedEdges =(
        Edges.findAt(coordinates=(thick, thick-a*0.2, lbox + a*0.2 * tan(THETA))), 
        Edges.findAt(coordinates=(thick, thick - clad*1.1, lbox + clad*1.1*tan(THETA))),
        Edges.findAt(coordinates=(thick, thick*0.1, lbox+thick*0.9*tan(THETA)))
        )
    ExtrudeLine = Edges.findAt(coordinates=(thick*0.3, thick, thick))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

######## Solid Extrusion for C (LxNegative)

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick*0.1,thick*0.3))
Edge_For_Sketch = Edges.findAt(coordinates=(0.0, thick*0.1, 0.0))

Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))

Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=45.2, gridSpacing=1.13, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(0.0, 0.0), point2=(0.0, thick))
Sketch.Line(point1=(0.0, thick), point2=(thick, thick))
Sketch.Line(point1=(thick, thick), point2=(thick, 0.0))
Sketch.Line(point1=(thick, 0.0), point2=(0.0, 0.0))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,  sketchPlaneSide=SIDE1, 
    sketchOrientation=RIGHT, sketch=Sketch, depth=C, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

############## Extruding A1 in LxNegatvie (C)

##### A1
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A1*cos(pi/4), A1*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(pi/4), A1*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ##### A2

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A2*cos(pi/4), A2*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(pi/4), A2*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A3

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A3*cos(pi/4), A3*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(pi/4), A3*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A4

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A4*cos(pi/4), A4*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(pi/4), A4*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### A5

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A5*cos(pi/4), A5*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(pi/4), A5*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A6

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A6*cos(pi/4), A6*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(pi/4), A6*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A7

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A7*cos(pi/4), A7*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(pi/4), A7*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A8

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A8*cos(pi/4), A8*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(pi/4), A8*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A9

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A9*cos(pi/4), A9*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(pi/4), A9*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### A10
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(((-C*0.3, thick*0.3, thick*0.9), ))
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A10*cos(pi/4), A10*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A10*cos(pi/4), A10*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, 0.))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


####### Extrude Center of Crack in LxNegative (C)

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(
    ((-C*0.76, thick-a, A1*0.99), ),
    ((-C*0.76, thick-a, A2*0.99), ),
    ((-C*0.76, thick-a, A3*0.99), ),
    ((-C*0.76, thick-a, A4*0.99), ),
    ((-C*0.76, thick-a, A5*0.99), ),
    ((-C*0.76, thick-a, A6*0.99), ),
    ((-C*0.76, thick-a, A7*0.99), ),
    ((-C*0.76, thick-a, A8*0.99), ),
    ((-C*0.76, thick-a, A9*0.99), ),
    ((-C*0.76, thick-a, A10*0.99), ),
    ((-C*0.76, thick-a, A10*1.99), ),
    )
if Weld_Angle >0 :
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-a, A1*0.9)),
        Edges.findAt(coordinates=(0.0, thick-a, A2*0.9)),
        Edges.findAt(coordinates=(0.0, thick-a, A3*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A4*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A5*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A6*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A7*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A8*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A9*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A9*1.01)), 
        Edges.findAt(coordinates=(0.0, thick-a, A10*1.01)),
        Edges.findAt(coordinates=(0.0, thick-a, lbox*1.1)),
        Edges.findAt(coordinates=(0.0, thick-a, thick*0.8)),
        )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-a, A1*0.9)),
        Edges.findAt(coordinates=(0.0, thick-a, A2*0.9)),
        Edges.findAt(coordinates=(0.0, thick-a, A3*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A4*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A5*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A6*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A7*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A8*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A9*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-a, A9*1.01)), 
        Edges.findAt(coordinates=(0.0, thick-a, A10*1.01)),
        Edges.findAt(coordinates=(0.0, thick-a, thick*0.8)),
        )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.4, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

########### extruding Weld Fixed Width line in Lxnegative (C)

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.75, thick-a*0.99, thick*0.2), ),
    ((-C*0.75, thick-a*1.01, thick*0.2), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a*0.99, lbox)),
    Edges.findAt(coordinates=(0.0, thick-a*1.01, lbox)), 
    Edges.findAt(coordinates=(0.0, thick-clad*1.1, lbox)), 
    Edges.findAt(coordinates=(0.0, thick*0.01, lbox))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.4, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#########################################
###### extrude V-groove in LxNegative (C)
if Weld_Angle > 0:
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    pickedCells = Cells.findAt(
        ((-C*0.3, thick-a*0.2, thick*0.6), ),
        ((-C*0.3, thick-a*1.1,thick*0.6), ),

        )
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-a*0.2, lbox+a*0.2*tan(THETA))), 
        Edges.findAt(coordinates=(0.0, thick-a*1.01, lbox+a*1.01*tan(THETA))), 
        Edges.findAt(coordinates=(0.0, thick-clad*1.01, lbox+clad*1.01*tan(THETA))), 
        Edges.findAt(coordinates=(0.0, thick*0.1, lbox+thick*0.9*tan(THETA))), 

        )
    ExtrudeLine = Edges.findAt(coordinates=(-C*0.3, 0.0, thick))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### extrude thick-clad lin in Lxnegative (C)
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
            ((-C*0.4, thick*0.2, lbox*0.99), ), 
            ((-C*0.4, thick*0.2, lbox*1.01), ), 
            ((-C*0.4, thick*0.2, thick*0.9), )
            )
if Weld_Angle > 0 : 
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-clad, lbox*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-clad, lbox*1.1)),
        Edges.findAt(coordinates=(0.0, thick-clad, thick*0.9))
        )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-clad, lbox*0.9)), 
        Edges.findAt(coordinates=(0.0, thick-clad, thick*0.9))
        )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.4, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### extrude (thick-clad)/2 lin in Lxnegative (C)

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
            ((-C*0.4, thick*0.2, lbox*0.99), ), 
            ((-C*0.4, thick*0.2, lbox*1.1), ), 
            ((-C*0.4, thick*0.2, thick*0.8), )
            )
if Weld_Angle > 0 : 
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, (thick-clad)/2, lbox*0.9)), 
        Edges.findAt(coordinates=(0.0, (thick-clad)/2, lbox*1.1)),
        Edges.findAt(coordinates=(0.0, (thick-clad)/2, thick*0.9)),
    )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, (thick-clad)/2, lbox*0.9)), 
        Edges.findAt(coordinates=(0.0, (thick-clad)/2, thick*0.9)),
    )

ExtrudeLine = Edges.findAt(coordinates=(-C*0.4, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##########################################
#### Solid extrude of LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges


Face_Sketch_plane = Faces.findAt(coordinates=(thick, thick*0.15,thick*0.4))
Edge_For_Sketch = Edges.findAt(coordinates=(thick, thick*0.33, thick))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(thick, 0.0, thick))
Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=214.08, gridSpacing=5.35, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(-thick, 0.0), point2=(0.0, 0.0))
Sketch.Line(point1=(0.0, 0.0), point2=(0.0, -thick))
Sketch.Line(point1=(0.0, -thick), point2=(-thick, -thick))
Sketch.Line(point1=(-thick, -thick), point2=(-thick, 0.0))

PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=Sketch, depth=LxPositive,   flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

######## Extrude lbox in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((LxPositive*0.8, thick, thick*0.3), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.99, lbox)), 
    Edges.findAt(coordinates=(thick, thick*0.7, lbox)), 
    Edges.findAt(coordinates=(thick, thick*0.1, lbox))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

######### Extrude Clad in LxPositive
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((LxPositive*0.8, thick*0.7, lbox*0.3), ), 
    ((LxPositive*0.8, thick*0.7, lbox*1.01), ), 
    ((LxPositive*0.8, thick*0.7, thick*0.8), )
    )
if Weld_Angle > 0 :
    pickedEdges =(
        Edges.findAt(coordinates=(thick, thick-clad, lbox *0.3)), 
        Edges.findAt(coordinates=(thick, thick-clad, lbox *1.01)), 
        Edges.findAt(coordinates=(thick, thick-clad, thick*0.8))
    )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(thick, thick-clad, lbox *0.3)), 
        Edges.findAt(coordinates=(thick, thick-clad, thick*0.8))
    )

ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

######### Extrude (thick-Clad)/2 in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((LxPositive*0.34, thick*0.3, lbox*0.3), ), 
    ((LxPositive*0.34, thick*0.3, lbox*1.01), ), 
    ((LxPositive*0.34, thick*0.3, thick*0.8), )
    )
if Weld_Angle > 0 :    
    pickedEdges =(
        Edges.findAt(coordinates=(thick, (thick-clad)/2, lbox*0.3)), 
        Edges.findAt(coordinates=(thick, (thick-clad)/2, lbox*1.01)), 
        Edges.findAt(coordinates=(thick, (thick-clad)/2, thick*0.8))
        )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(thick, (thick-clad)/2, lbox*0.3)), 
        Edges.findAt(coordinates=(thick, (thick-clad)/2, thick*0.8))
    )

ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### Extrude V-groove in Lx Positive:

if Weld_Angle > 0 :
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    pickedCells = Cells.findAt(
        ((LxPositive*0.3, thick-a*0.3, lbox*1.1), ),
        ((LxPositive*0.3, thick-clad*1.1, lbox*1.1), ),
        ((LxPositive*0.3, thick*0.1, lbox*1.1), )
        )

    pickedEdges =(
        Edges.findAt(coordinates=(thick, thick-a*0.3, lbox+a*0.3*tan(THETA))),
        Edges.findAt(coordinates=(thick, thick-clad *1.1, lbox+clad*1.1*tan(THETA))),
        Edges.findAt(coordinates=(thick, thick*0.1,lbox+thick*0.9*tan(THETA)))
        )
    ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.3, 0.0, thick))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,cells=pickedCells, edges=pickedEdges, sense=FORWARD)


############# Sketching h3 

if Weld_Angle == 0:
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    Face_Sketch_plane = Faces.findAt(coordinates=(LxPositive+thick, thick*0.1, thick*0.4))
    Edge_For_Sketch = Edges.findAt(coordinates=(LxPositive+thick, thick*0.1, thick))

    Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,sketchPlaneSide=SIDE1, origin=(LxPositive+thick, 0.0, thick))
    Sketch = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=443.14, gridSpacing=11.07, transform=Transform)
    Sketch.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

    Sketch.Line(point1=(h3, thick), point2=(h3, 0.0))

    pickedFaces = Faces.findAt(
        ((LxPositive+thick, thick*0.1, thick*0.4), ),
        ((LxPositive+thick, thick*0.7, thick*0.4), ), 
        ((LxPositive+thick, thick*0.99, thick*0.4), )
        )
    PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=Sketch)

    ######### Extruding h3 line in all length

    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=thick)

    pickedEdges =(
        Edges.findAt(coordinates=(thick+LxPositive, thick*0.99, (thick-lbox)/2 )), 
        Edges.findAt(coordinates=(thick+LxPositive, thick*0.7, (thick-lbox)/2)), 
        Edges.findAt(coordinates=(thick+LxPositive, thick*0.1, (thick-lbox)/2))
        )
    ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, thick))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)


################ Solid Extrude for L in Z direction

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(LxPositive*0.7, thick*0.3, thick))
Edge_For_Sketch = Edges.findAt(coordinates=(thick+LxPositive, thick*0.2, thick))

Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge= Edge_For_Sketch, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(thick+LxPositive, 0.0,thick))
Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=443.14, gridSpacing=11.07, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(0.0, 0.0), point2=(0.0, thick))
Sketch.Line(point1=(0.0, thick), point2=(-C-thick-LxPositive, thick))
Sketch.Line(point1=(-C-thick-LxPositive, thick), point2=(-C-thick-LxPositive, 0.0))
Sketch.Line(point1=(-C-thick-LxPositive, 0.0), point2=(0.0, 0.0))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane,  sketchUpEdge=Edge_For_Sketch,sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=Sketch, depth=L-thick, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

#### extrude thick line throguh L
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.2, thick*0.3, L*0.3), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.99, thick)), 
    Edges.findAt(coordinates=(thick, thick-clad*1.1, thick)), 
    Edges.findAt(coordinates=(thick, thick*0.01, thick)))

ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, L*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,cells=pickedCells, edges=pickedEdges, sense=FORWARD)

################################ Weld_Fixed_Width
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.2, thick*0.3, L*0.3), )
    )
## revised for mesh purposes, fixed width changed to clad to be closer to the last curve
pickedEdges =(
    Edges.findAt(coordinates=(clad, thick*0.99, thick)), 
    Edges.findAt(coordinates=(clad, thick-clad*1.1, thick)), 
    Edges.findAt(coordinates=(clad, thick*0.01, thick)))

ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, L*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,cells=pickedCells, edges=pickedEdges, sense=FORWARD)

############################ Zero
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.2, thick*0.3, L*0.3), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a*0.99, thick)), 
    Edges.findAt(coordinates=(0.0, thick-a*1.01, thick)), 
    Edges.findAt(coordinates=(0.0, thick-clad*1.01, thick)), 
   Edges.findAt(coordinates=(0.0, thick*0.3, thick))
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, L*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,cells=pickedCells, edges=pickedEdges, sense=FORWARD)

########### Extrude of (thick-clad)/2 in L
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((LxPositive*0.7, thick*0.3, L*0.3), ),
    ((thick*0.7, thick*0.3, L*0.3), ),
    ((clad*0.7, thick*0.3, L*0.3), ),
    ((-C*0.7, thick*0.3, L*0.3), ),


    )
pickedEdges =(
    Edges.findAt(coordinates=(LxPositive*0.56, (thick-clad)/2, thick)), 
    Edges.findAt(coordinates=(thick*0.3, (thick-clad)/2, thick)), 
    Edges.findAt(coordinates=(thick*0.01, (thick-clad)/2, thick)),
    Edges.findAt(coordinates=(-C*0.2, (thick-clad)/2, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, L*0.3))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#######################

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((a*0.3, thick*0.9, thick*1.5), ), 
    ((-C*0.3, thick*0.9, thick*1.5), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(-C*0.3, thick-a, thick)), 
    Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*2))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

pickedCells = Cells.findAt(
    ((-C*0.3, thick-clad*0.99, thick*1.99), ),
    ((a*0.3, thick-clad*0.99, thick*1.99), ),
    ((thick*0.8, thick*0.99, thick*1.99), ),
    ((LxPositive*0.7,thick*0.99, thick*1.99), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(-C*0.3, thick-clad, thick)), 
    Edges.findAt(coordinates=(a*0.3, thick-clad, thick)), 
    Edges.findAt(coordinates=(thick*0.8, thick-clad, thick)),
    Edges.findAt(coordinates=(LxPositive*0.7, thick-clad, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, (thick-clad)/2, thick*1.2))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

###################### Sketching LDiv Lines

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(thick+LxPositive, thick-a*0.99, thick*1.5))
Edge_For_Sketch = Edges.findAt(coordinates=(thick+LxPositive, thick*0.1, thick))

Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,sketchPlaneSide=SIDE1, origin=(thick+LxPositive, 0.0, thick))

Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=1195.51, gridSpacing=29.88, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)

PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(-LDiv1, thick), point2=(-LDiv1, 0.0))
Sketch.Line(point1=(-LDiv2, thick), point2=(-LDiv2, 0.0))
Sketch.Line(point1=(-LDiv3, thick), point2=(-LDiv3, 0.0))

pickedFaces = Faces.findAt(
    ((thick+LxPositive, thick-a*0.99, thick*1.5), ), 
    ((thick+LxPositive, thick-clad*1.01, thick*1.5), ), 
    ((thick+LxPositive, thick*0.1, thick*1.5), )
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=Sketch)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


############################
##### Extruding LDiv1 
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=L)

pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*0.99, thick+LDiv1)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*1.1, thick+LDiv1)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick*0.1, thick+LDiv1))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, thick))

PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)


#### Extruding LDiv2
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=L)

pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*0.99, thick+LDiv2)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*1.1, thick+LDiv2)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick*0.1, thick+LDiv2))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, thick))

PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)


# ##### Extruding LDiv3
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=L)

pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*0.99, thick+LDiv3)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick-clad*1.1, thick+LDiv3)), 
    Edges.findAt(coordinates=(thick+LxPositive, thick*0.1, thick+LDiv3))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive*0.34, 0.0, thick))

PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

### Sketching LxPositive*0.35 Line

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(thick+LxPositive/2, thick-a*0.99,L))
Edge_For_Sketch = Edges.findAt(coordinates=(thick, thick*0.1,  L))

Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch,sketchPlaneSide=SIDE1, origin=(thick, 0.0, L))

Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=1274.31, gridSpacing=31.85, transform=Transform)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(LxPositive*0.35, thick), point2=(LxPositive*0.35, 0.0))

pickedFaces = Faces.findAt(
    ((thick+LxPositive/2, thick-a*0.99,L), ), 
    ((thick+LxPositive/2, thick-clad*1.01,L), ), 
    ((thick+LxPositive/2, thick*0.1,L), )
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=Sketch)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

########## Extruding LxPositive *0.35

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=L)

pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive*0.35, thick-clad*0.99, L)), 
    Edges.findAt(coordinates=(thick+LxPositive*0.35, thick-clad*1.1, L)),
    Edges.findAt(coordinates=(thick+LxPositive*0.35, thick*0.1, L))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick, 0.0, L*0.75))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

################### Material
################ Create Material



if MaterialType == 'Elastic-Plastic':
    PipeModel.Material('BaseSteel')
    PipeModel.materials['BaseSteel'].Density(table=((BaseSteelDens, ), ))
    PipeModel.materials['BaseSteel'].Elastic(table=((BaseSteelE,BaseSteelPratio), ))
    PipeModel.materials['BaseSteel'].Plastic(table=((Yield_Stress_Base, 0.0),  (Yield_Stress_Base, Yield_Strain_Base)))
    PipeModel.HomogeneousSolidSection(name='BaseSteel_Section', material='BaseSteel', thickness=None)
    ##


if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='BaseSteel')
    PipeModel.materials['BaseSteel'].Density(table=((BaseSteelDens, ), ))
    PipeModel.materials['BaseSteel'].Elastic(table=((BaseSteelE,  BaseSteelPratio), ))
    PipeModel.materials['BaseSteel'].Plastic(table=(((Yield_Stress_Base, 0.0),  (Yield_Stress_Base, Yield_Strain_Base)),  
                                                    (SS_Curve_Base[1][0], SS_Curve_Base[1][1]),
                                                    ))
    PipeModel.HomogeneousSolidSection(name='BaseSteelSection', material='BaseSteel', thickness=None)

if MaterialType == 'Plastic-Deformation':
    BaseSteelMaterial=PipeModel.Material('BaseSteel')
    BaseSteelMaterial.DeformationPlasticity(table=((BaseSteelE, BaseSteelPratio, Yield_Stress_Base, N_Base, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='BaseSteelSection', material='BaseSteel', thickness=None)

cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0,yMax=thick-clad,zMin=lbox,zMax=L)

region = PipePart.Set(cells=cells1, name='SteelBase_Set')
PipePart.SectionAssignment(region=region, sectionName='BaseSteelSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='WeldSteel')
    PipeModel.materials['WeldSteel'].Density(table=((WeldSteelDens, ), ))
    PipeModel.materials['WeldSteel'].Elastic(table=((WeldSteelE,  WeldSteelPratio), ))
    PipeModel.materials['WeldSteel'].Plastic(table=((SS_Curve_Weld[0][0], SS_Curve_Weld[0][1]),  (SS_Curve_Weld[1][0], SS_Curve_Weld[1][1])))
    PipeModel.HomogeneousSolidSection(name='WeldSteelSection', material='WeldSteel', thickness=None)

if MaterialType == 'Plastic-Deformation':
    WeldSteelMaterial=PipeModel.Material('WeldSteel')
    WeldSteelMaterial.DeformationPlasticity(table=((WeldSteelE, WeldSteelPratio, Yield_Stress_Weld, N_Weld, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='WeldSteelSection', material='WeldSteel', thickness=None) 

cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick-clad,zMin=0.0,zMax=lbox)
region = PipePart.Set(cells=cells, name='WeldSteel_Set')
PipePart.SectionAssignment(region=region, sectionName='WeldSteelSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='CladSteelWeld')
    PipeModel.materials['CladSteelWeld'].Density(table=((CladSteelWeldDens, ), ))
    PipeModel.materials['CladSteelWeld'].Elastic(table=((CladSteelWeldE,  CladSteelWeldPratio), ))
    PipeModel.materials['CladSteelWeld'].Plastic(table=((SS_Curve_CladWeld[0][0], SS_Curve_CladWeld[0][1]),  (SS_Curve_CladWeld[1][0], SS_Curve_CladWeld[1][1])))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeldSection', material='CladSteelWeld', thickness=None)

if MaterialType == 'Plastic-Deformation':
    CladSteelWeldMaterial=PipeModel.Material('CladSteelWeld')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelWeldE, CladSteelWeldPratio, Yield_Stress_CladSteelWeld, N_CladWeld, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeldSection', material='CladSteelWeld', thickness=None) 

cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=thick-clad,yMax=thick,zMin=0.0,zMax=lbox)
region = PipePart.Set(cells=cells, name='CladSteelWeld_Set')
PipePart.SectionAssignment(region=region, sectionName='CladSteelWeldSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

if MaterialType == 'Elastic-Plastic':
    PipeModel.Material('CladSteelBase')
    PipeModel.materials['CladSteelBase'].Density(table=((CladSteelBaseDens, ), ))
    PipeModel.materials['CladSteelBase'].Elastic(table=((CladSteelBaseE,CladSteelBasePratio), ))
    PipeModel.materials['CladSteelBase'].Plastic(table=((SS_Curve_Clad[0][0], SS_Curve_Clad[0][1]),  (SS_Curve_Clad[1][0], SS_Curve_Clad[1][1])))
    PipeModel.HomogeneousSolidSection(name='CladSteelBaseSection', material='CladSteelBase', thickness=None)

if MaterialType == 'Plastic-Deformation':
    CladSteelWeldMaterial=PipeModel.Material('CladSteelBase')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelBaseE, CladSteelBasePratio, Yield_Stress_CladSteelBase, N_clad, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelBaseSection', material='CladSteelBase', thickness=None)   


cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=thick-clad,yMax=thick,zMin=lbox,zMax=L)

region = PipePart.Set(cells=cells1, name='CladSteelBase_Set')
PipePart.SectionAssignment(region=region, sectionName='CladSteelBaseSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

###### Sketching A11 (a1)  to improve mesh

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

a1 = 12 * A1

Face_Sketch_plane=Faces.findAt(coordinates=(-C*0.3, thick*0.98, L))
Edge_For_Sketch = Edges.findAt(coordinates=(0.0, thick*0.98, L))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_For_Sketch, sketchPlaneSide=SIDE1, origin=(0.0, thick, L))
sketch_1 = mdb.models['Pipe_Model_Internal'].ConstrainedSketch(name='__profile__',  sheetSize=84.83, gridSpacing=2.12, transform=Transform)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.ArcByCenterEnds(center=(0.0, 0.0), point1=((a-a1), 0.0), point2=(0.0, -(a-a1)), direction=CLOCKWISE)

sketch_1.Line(point1=(0.0, -(a-a1)), point2=(-C, -(a-a1)))

pickedFaces = Faces.findAt(
    ((-C*0.3, thick*0.98, L), ), 
    ((a*0.2,    thick*0.98, L), ))
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_For_Sketch, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

##### Extruding A11 (a1)
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick,yMin=0.,yMax=thick,zMin=0.0,zMax= L )

pickedEdges =(
    Edges.findAt(coordinates=((a-a1)*cos(pi/3), thick-(a-a1)*sin(pi/3), L)), 
    Edges.findAt(coordinates=(-C*0.3, thick-(a-a1), L))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, thick, L*0.75))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

##############################################################
############ SEEDs MESH
################################################################
##### partial global meshing
PipePart.seedPart(size=11.0, deviationFactor=0.1, minSizeFactor=0.1)
PipePart.generateMesh()



pickedEdges = Edges.findAt(
    #### at -C
    ((-C, thick-a+A2*cos(pi/4), A2*sin(pi/4)), ),
    ((-C, thick-a+A3*cos(pi/4), A3*sin(pi/4)), ),
    ((-C, thick-a+A4*cos(pi/4), A4*sin(pi/4)), ),
    ((-C, thick-a+A5*cos(pi/4), A5*sin(pi/4)), ),
    ((-C, thick-a+A6*cos(pi/4), A6*sin(pi/4)), ),
    ((-C, thick-a+A7*cos(pi/4), A7*sin(pi/4)), ),
    ((-C, thick-a+A8*cos(pi/4), A8*sin(pi/4)), ),
    ((-C, thick-a+A9*cos(pi/4), A9*sin(pi/4)), ),
    ((-C, thick-a+A10*cos(pi/4), A10*sin(pi/4)), ),

    ((-C, thick-a-A2*cos(pi/4), A2*sin(pi/4)), ),
    ((-C, thick-a-A3*cos(pi/4), A3*sin(pi/4)), ),
    ((-C, thick-a-A4*cos(pi/4), A4*sin(pi/4)), ),
    ((-C, thick-a-A5*cos(pi/4), A5*sin(pi/4)), ),
    ((-C, thick-a-A6*cos(pi/4), A6*sin(pi/4)), ),
    ((-C, thick-a-A7*cos(pi/4), A7*sin(pi/4)), ),
    ((-C, thick-a-A8*cos(pi/4), A8*sin(pi/4)), ),
    ((-C, thick-a-A9*cos(pi/4), A9*sin(pi/4)), ),
    ((-C, thick-a-A10*cos(pi/4), A10*sin(pi/4)), ),

     #### at X=0
    ((0.0, thick-a+A2*cos(pi/4), A2*sin(pi/4)), ),
    ((0.0, thick-a+A3*cos(pi/4), A3*sin(pi/4)), ),
    ((0.0, thick-a+A4*cos(pi/4), A4*sin(pi/4)), ),
    ((0.0, thick-a+A5*cos(pi/4), A5*sin(pi/4)), ),
    ((0.0, thick-a+A6*cos(pi/4), A6*sin(pi/4)), ),
    ((0.0, thick-a+A7*cos(pi/4), A7*sin(pi/4)), ),
    ((0.0, thick-a+A8*cos(pi/4), A8*sin(pi/4)), ),
    ((0.0, thick-a+A9*cos(pi/4), A9*sin(pi/4)), ),
    ((0.0, thick-a+A10*cos(pi/4), A10*sin(pi/4)), ),

    ((0.0, thick-a-A2*cos(pi/4), A2*sin(pi/4)), ),
    ((0.0, thick-a-A3*cos(pi/4), A3*sin(pi/4)), ),
    ((0.0, thick-a-A4*cos(pi/4), A4*sin(pi/4)), ),
    ((0.0, thick-a-A5*cos(pi/4), A5*sin(pi/4)), ),
    ((0.0, thick-a-A6*cos(pi/4), A6*sin(pi/4)), ),
    ((0.0, thick-a-A7*cos(pi/4), A7*sin(pi/4)), ),
    ((0.0, thick-a-A8*cos(pi/4), A8*sin(pi/4)), ),
    ((0.0, thick-a-A9*cos(pi/4), A9*sin(pi/4)), ),
    ((0.0, thick-a-A10*cos(pi/4), A10*sin(pi/4)), ),

     #### at X=a
    ((a+A2*cos(pi/4), thick, A2*sin(pi/4)), ),
    ((a+A3*cos(pi/4), thick, A3*sin(pi/4)), ),
    ((a+A4*cos(pi/4), thick, A4*sin(pi/4)), ),
    ((a+A5*cos(pi/4), thick, A5*sin(pi/4)), ),
    ((a+A6*cos(pi/4), thick, A6*sin(pi/4)), ),
    ((a+A7*cos(pi/4), thick, A7*sin(pi/4)), ),
    ((a+A8*cos(pi/4), thick, A8*sin(pi/4)), ),
    ((a+A9*cos(pi/4), thick, A9*sin(pi/4)), ),
    ((a+A10*cos(pi/4), thick, A10*sin(pi/4)), ),

    ((a-A2*cos(pi/4), thick, A2*sin(pi/4)), ),
    ((a-A3*cos(pi/4), thick, A3*sin(pi/4)), ),
    ((a-A4*cos(pi/4), thick, A4*sin(pi/4)), ),
    ((a-A5*cos(pi/4), thick, A5*sin(pi/4)), ),
    ((a-A7*cos(pi/4), thick, A7*sin(pi/4)), ),
    ((a-A8*cos(pi/4), thick, A8*sin(pi/4)), ),
    ((a-A9*cos(pi/4), thick, A9*sin(pi/4)), ),
    ((a-A10*cos(pi/4), thick, A10*sin(pi/4)), ),

    ((a-A1*cos(pi/4), thick, A1*sin(pi/4)), ),
    ((a+A1*cos(pi/4), thick, A1*sin(pi/4)), ),

    ### ((a*cos(pi/4), thick-a*cos(pi/4), Wedge), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A1), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A2), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A3), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A4), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A5), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A6), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A7), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A8), ),
    ((a*cos(pi/4), thick-a*cos(pi/4), A9), ),


    ((0.0, thick-(a-A1*cos(pi/4)), A1*sin(pi/4)), ),
    ((0.0, thick-(a+A1*cos(pi/4)), A1*sin(pi/4)), ),

    ((-C, thick-(a-A1*cos(pi/4)), A1*sin(pi/4)), ),
    ((-C, thick-(a+A1*cos(pi/4)), A1*sin(pi/4)), ),

    # ((a, thick, Wedge), ),
    # ((0.0, thick-a, Wedge), ),


    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedSpider, constraint=FINER)
# # ##################

##########  nSeeda
if Weld_Angle == 0:
    pickedEdges = Edges.findAt(

        ((clad, thick-Epsilon, 0.0), ),
        ((clad, thick-Epsilon, lbox), ),
        ((clad, thick-Epsilon, thick-h3), ),
        ((clad, thick-Epsilon, thick), ),
        ((clad, thick-Epsilon, thick+LDiv1), ),
        ((clad, thick-Epsilon, thick+LDiv2), ),
        ((clad, thick-Epsilon, thick+LDiv3), ),
        ((clad, thick-Epsilon, L), ),

        ((thick, thick-Epsilon, 0.0), ),
        ((thick, thick-Epsilon, lbox), ),
        ((thick, thick-Epsilon, thick-h3), ),
        ((thick, thick-Epsilon, thick), ),
        ((thick, thick-Epsilon, thick+LDiv1), ),
        ((thick, thick-Epsilon, thick+LDiv2), ),
        ((thick, thick-Epsilon, thick+LDiv3), ),
        ((thick, thick-Epsilon, L), ),

        ((thick+LxPositive*0.35, thick-Epsilon, 0.0), ),
        ((thick+LxPositive*0.35, thick-Epsilon, lbox), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick-h3), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv1), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv2), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv3), ),
        ((thick+LxPositive*0.35, thick-Epsilon, L), ),

        ((thick+LxPositive, thick-Epsilon, 0.0), ),
        ((thick+LxPositive, thick-Epsilon, lbox), ),
        ((thick+LxPositive, thick-Epsilon, thick-h3), ),
        ((thick+LxPositive, thick-Epsilon, thick), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv1), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv2), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv3), ),
        ((thick+LxPositive, thick-Epsilon, L), ),

    )
else:
        pickedEdges = Edges.findAt(
 
        ((clad, thick-Epsilon, 0.0), ),
        ((clad, thick-Epsilon, lbox), ),
        ((clad, thick-Epsilon, thick), ),
        ((clad, thick-Epsilon, thick+LDiv1), ),
        ((clad, thick-Epsilon, thick+LDiv2), ),
        ((clad, thick-Epsilon, thick+LDiv3), ),
        ((clad, thick-Epsilon, L), ),

        ((thick, thick-Epsilon, 0.0), ),
        ((thick, thick-Epsilon, lbox), ),
        ((thick, thick-Epsilon, thick), ),
        ((thick, thick-Epsilon, thick+LDiv1), ),
        ((thick, thick-Epsilon, thick+LDiv2), ),
        ((thick, thick-Epsilon, thick+LDiv3), ),
        ((thick, thick-Epsilon, L), ),

        ((thick+LxPositive*0.35, thick-Epsilon, 0.0), ),
        ((thick+LxPositive*0.35, thick-Epsilon, lbox), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv1), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv2), ),
        ((thick+LxPositive*0.35, thick-Epsilon, thick+LDiv3), ),
        ((thick+LxPositive*0.35, thick-Epsilon, L), ),

        ((thick+LxPositive, thick-Epsilon, 0.0), ),
        ((thick+LxPositive, thick-Epsilon, lbox), ),
        ((thick+LxPositive, thick-Epsilon, thick), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv1), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv2), ),
        ((thick+LxPositive, thick-Epsilon, thick+LDiv3), ),
        ((thick+LxPositive, thick-Epsilon, L), ),
    )

PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseeda, constraint=FINER)

#######################
#######
RR5 = ((a*sin(pi/4))) /cos(THETA)
rr5 = RR5 * sin(THETA)

pickedEdges = Edges.findAt(

    (((a-sqrt(A10**2)) * cos(pi/4) , thick-((a-sqrt(A10**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A9**2)) * cos(pi/4) , thick-((a-sqrt(A9**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A8**2)) * cos(pi/4) , thick-((a-sqrt(A8**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A7**2)) * cos(pi/4) , thick-((a-sqrt(A7**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A6**2)) * cos(pi/4) , thick-((a-sqrt(A6**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A5**2)) * cos(pi/4) , thick-((a-sqrt(A5**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A4**2)) * cos(pi/4) , thick-((a-sqrt(A4**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A3**2)) * cos(pi/4) , thick-((a-sqrt(A3**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A2**2)) * cos(pi/4) , thick-((a-sqrt(A2**2))*sin(pi/4)), 0.0), ), 
    (((a-sqrt(A1**2)) * cos(pi/4) , thick-((a-sqrt(A1**2))*sin(pi/4)), 0.0), ), 

    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A1), ), 

    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A2), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A3), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A4), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A5), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A6), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A7), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A8), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A8), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), A10), ), 

    (((a+A1)*cos(pi/4), thick-(a+A1)*sin(pi/4), 0.0), ),
    (((a+A2)*cos(pi/4), thick-(a+A2)*sin(pi/4), 0.0), ),
    (((a+A3)*cos(pi/4), thick-(a+A3)*sin(pi/4), 0.0), ),
    (((a+A4)*cos(pi/4), thick-(a+A4)*sin(pi/4), 0.0), ),
    (((a+A5)*cos(pi/4), thick-(a+A5)*sin(pi/4), 0.0), ),
    (((a+A6)*cos(pi/4), thick-(a+A6)*sin(pi/4), 0.0), ),
    (((a+A7)*cos(pi/4), thick-(a+A7)*sin(pi/4), 0.0), ),
    (((a+A8)*cos(pi/4), thick-(a+A8)*sin(pi/4), 0.0), ),
    (((a+A9)*cos(pi/4), thick-(a+A9)*sin(pi/4), 0.0), ),
    (((a+A10)*cos(pi/4), thick-(a+A10)*sin(pi/4), 0.0), ),

    ((a * cos(pi/4) , thick-(a*sin(pi/4)), lbox), ), 
    # ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick-h3), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv1), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv2), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv3), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), L), ), 


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseedcurve, constraint=FINER)



Sintheta = 0.0 
XA10 = sqrt(A10**2)
pickedEdges = Edges.findAt(
    (( -C*0.2, thick-a+sqrt(A10**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A9**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A8**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A7**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A6**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A5**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A4**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A3**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A2**2), 0.0), ),
    (( -C*0.2, thick-a+sqrt(A1**2), 0.0), ),


    (( -C*0.2, thick-a-A1, 0.0), ),
    (( -C*0.2, thick-a-A2, 0.0), ),
    (( -C*0.2, thick-a-A3, 0.0), ),
    (( -C*0.2, thick-a-A4, 0.0), ),
    (( -C*0.2, thick-a-A5, 0.0), ),
    (( -C*0.2, thick-a-A6, 0.0), ),
    (( -C*0.2, thick-a-A7, 0.0), ),
    (( -C*0.2, thick-a-A8, 0.0), ),
    (( -C*0.2, thick-a-A9, 0.0), ),
    (( -C*0.2, thick-a-A10, 0.0), ),

     (( -C*0.2, thick, 0.0), ),
    (( -C*0.2, thick-clad, 0.0), ),
    (( -C*0.2, (thick-clad)/2, 0.0), ),
    (( -C*0.2, 0.0, 0.0), ),

    (( -C*0.2, thick, lbox), ),
    (( -C*0.2, thick-a, lbox), ),
    (( -C*0.2, thick-clad, lbox), ),
    (( -C*0.2, (thick-clad)/2, lbox), ),
    (( -C*0.2, 0.0, lbox), ),

    # (( -C*0.2, thick, thick-h3), ),
    # (( -C*0.2, thick-a, thick-h3), ),
    # (( -C*0.2, thick-clad, thick-h3), ),
    # (( -C*0.2, (thick-clad)/2, thick-h3), ),
    # (( -C*0.2, 0.0, thick-h3), ),
   
    (( -C*0.2, thick, thick), ),
    (( -C*0.2, thick-a, thick), ),
    (( -C*0.2, thick-clad, thick), ),
    (( -C*0.2, (thick-clad)/2, thick), ),
    (( -C*0.2, 0.0, thick), ),

    (( -C*0.2, thick, thick+LDiv1), ),
    (( -C*0.2, thick-clad, thick+LDiv1), ),
    (( -C*0.2, thick-a, thick+LDiv1), ),
    (( -C*0.2, (thick-clad)/2, thick+LDiv1), ),
    (( -C*0.2, 0.0, thick+LDiv1), ),

    (( -C*0.2, thick, thick+LDiv2), ),
    (( -C*0.2, thick-clad, thick+LDiv2), ),
    (( -C*0.2, thick-a, thick+LDiv2), ),
    (( -C*0.2, (thick-clad)/2, thick+LDiv2), ),
    (( -C*0.2, 0.0, thick+LDiv2), ),

    (( -C*0.2, thick, thick+LDiv3), ),
    (( -C*0.2, thick-a, thick+LDiv3), ),
    (( -C*0.2, thick-clad, thick+LDiv3), ),
    (( -C*0.2, (thick-clad)/2, thick+LDiv3), ),
    (( -C*0.2, 0.0, thick+LDiv3), ),

    (( -C*0.2, thick, L), ),
    (( -C*0.2, thick-a, L), ),
    (( -C*0.2, thick-clad, L), ),
    (( -C*0.2, (thick-clad)/2, L), ),
    (( -C*0.2, 0.0, L), ),
    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedC, constraint=FINER)

# # ########################


pickedEdges = Edges.findAt(
        (( thick*1.1, thick, 0), ),
        (( thick*1.1, thick, lbox), ),
        (( thick*1.1, thick, thick), ),
        # (( thick*1.1, thick, thick-h3), ),
        (( thick*1.1, thick, thick+LDiv1), ),
        (( thick*1.1, thick, thick+LDiv2), ),
        (( thick*1.1, thick, thick+LDiv3), ),
        (( thick*1.1, thick, L), ),

        
        (( LxPositive*0.9, thick, 0), ),
        (( LxPositive*0.9, thick, lbox), ),
        (( LxPositive*0.9, thick, thick), ),
        # (( LxPositive*0.9, thick, thick-h3), ),
        (( LxPositive*0.9, thick, thick+LDiv1), ),
        (( LxPositive*0.9, thick, thick+LDiv2), ),
        (( LxPositive*0.9, thick, thick+LDiv3), ),
        (( LxPositive*0.9, thick, L), ),


        (( thick*1.1, thick-clad, 0), ),
        (( thick*1.1, thick-clad, lbox), ),
        (( thick*1.1, thick-clad, thick), ),
        # (( thick*1.1, thick-clad, thick-h3), ),
        (( thick*1.1, thick-clad, thick+LDiv1), ),
        (( thick*1.1, thick-clad, thick+LDiv2), ),
        (( thick*1.1, thick-clad, thick+LDiv3), ),
        (( thick*1.1, thick-clad, L), ),

        
        (( LxPositive*0.9, thick-clad, 0), ),
        (( LxPositive*0.9, thick-clad, lbox), ),
        (( LxPositive*0.9, thick-clad, thick), ),
        # (( LxPositive*0.9, thick-clad, thick-h3), ),
        (( LxPositive*0.9, thick-clad, thick+LDiv1), ),
        (( LxPositive*0.9, thick-clad, thick+LDiv2), ),
        (( LxPositive*0.9, thick-clad, thick+LDiv3), ),
        (( LxPositive*0.9, thick-clad, L), ),

        (( thick*1.1, (thick-clad)/2, 0), ),
        (( thick*1.1, (thick-clad)/2, lbox), ),
        (( thick*1.1, (thick-clad)/2, thick), ),
        # (( thick*1.1, (thick-clad)/2, thick-h3), ),
        (( thick*1.1, (thick-clad)/2, thick+LDiv1), ),
        (( thick*1.1, (thick-clad)/2, thick+LDiv2), ),
        (( thick*1.1, (thick-clad)/2, thick+LDiv3), ),
        (( thick*1.1, (thick-clad)/2, L), ),

        
        (( LxPositive*0.9, (thick-clad)/2, 0), ),
        (( LxPositive*0.9, (thick-clad)/2, lbox), ),
        (( LxPositive*0.9, (thick-clad)/2, thick), ),
        # (( LxPositive*0.9, (thick-clad)/2, thick-h3), ),
        (( LxPositive*0.9, (thick-clad)/2, thick+LDiv1), ),
        (( LxPositive*0.9, (thick-clad)/2, thick+LDiv2), ),
        (( LxPositive*0.9, (thick-clad)/2, thick+LDiv3), ),
        (( LxPositive*0.9, (thick-clad)/2, L), ),

        (( thick*1.1, 0.0, 0), ),
        (( thick*1.1, 0.0, lbox), ),
        (( thick*1.1, 0.0, thick), ),
        # (( thick*1.1, 0.0, thick-h3), ),
        (( thick*1.1, 0.0, thick+LDiv1), ),
        (( thick*1.1, 0.0, thick+LDiv2), ),
        (( thick*1.1, 0.0, thick+LDiv3), ),
        (( thick*1.1, 0.0, L), ),

        
        (( LxPositive*0.9, 0.0, 0), ),
        (( LxPositive*0.9, 0.0, lbox), ),
        (( LxPositive*0.9, 0.0, thick), ),
        # (( LxPositive*0.9, 0.0, thick-h3), ),
        (( LxPositive*0.9, 0.0, thick+LDiv1), ),
        (( LxPositive*0.9, 0.0, thick+LDiv2), ),
        (( LxPositive*0.9, 0.0, thick+LDiv3), ),
        (( LxPositive*0.9, 0.0, L), ),

        )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedLxPos, constraint=FINER)


pickedEdges = Edges.findAt(
        (( -C, thick, lbox*0.9), ),
        (( -C, thick-a, lbox*0.9), ),
        (( -C, thick-clad, lbox*0.9), ),
        (( -C, (thick-clad)/2, lbox*0.9), ),
        (( -C, 0.0, lbox*0.9), ),

        (( 0, thick, lbox*0.9), ),
        (( 0, thick-a, lbox*0.9), ),
        (( 0, thick-clad, lbox*0.9), ),
        (( 0, (thick-clad)/2, lbox*0.9), ),
        (( 0, 0.0, lbox*0.9), ),

        (( clad, thick, lbox*0.9), ),
        (( clad, thick-clad, lbox*0.9), ),
        (( clad, (thick-clad)/2, lbox*0.9), ),
        (( clad, 0.0, lbox*0.9), ),

        (( thick, thick, lbox*0.9), ),
        (( thick, thick-clad, lbox*0.9), ),
        (( thick, (thick-clad)/2, lbox*0.9), ),
        (( thick, 0.0, lbox*0.9), ),

        (( thick+LxPositive*0.35, thick, lbox*0.9), ),
        (( thick+LxPositive*0.35, thick-clad, lbox*0.9), ),
        (( thick+LxPositive*0.35, (thick-clad)/2, lbox*0.9), ),
        (( thick+LxPositive*0.35, 0.0, lbox*0.9), ),

    
        (( thick+LxPositive, thick, lbox*0.9), ),
        (( thick+LxPositive, thick-clad, lbox*0.9), ),
        (( thick+LxPositive, (thick-clad)/2, lbox*0.9), ),
        (( thick+LxPositive, 0.0, lbox*0.9), ),

        )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseedlbox, constraint=FINER)

# # ##################


pickedEdges = Edges.findAt(

        (( -C, thick, lbox*1.1), ),
        (( -C, thick-a, lbox*1.1), ),
        (( -C, thick-clad, lbox*1.1), ),
        (( -C, (thick-clad)/2, lbox*1.1), ),
        (( -C, 0.0, lbox*1.1), ),

        (( 0.0, thick, lbox*1.1), ),
        (( 0.0, thick-a, lbox*1.1), ),
        (( 0.0, thick-clad, lbox*1.1), ),
        (( 0.0, (thick-clad)/2, lbox*1.1), ),
        (( 0.0, 0.0, lbox*1.1), ),

        (( clad, thick, lbox*1.1), ),
        (( clad, thick-clad, lbox*1.1), ),
        (( clad, (thick-clad)/2, lbox*1.1), ),
        (( clad, 0.0, lbox*1.1), ),

        (( thick, thick, lbox*1.1), ),
        (( thick, thick-clad, lbox*1.1), ),
        (( thick, (thick-clad)/2, lbox*1.1), ),
        (( thick, 0.0, lbox*1.1), ),

        (( thick+LxPositive*0.35, thick, lbox*1.1), ),
        (( thick+LxPositive*0.35, thick-clad, lbox*1.1), ),
        (( thick+LxPositive*0.35, (thick-clad)/2, lbox*1.1), ),
        (( thick+LxPositive*0.35, 0.0, lbox*1.1), ),

        (( thick+LxPositive, thick, lbox*1.1), ),
        (( thick+LxPositive, thick-clad, lbox*1.1), ),
        (( thick+LxPositive, (thick-clad)/2, lbox*1.1), ),
        (( thick+LxPositive, 0.0, lbox*1.1), ),

        (( -C, thick, thick*0.9), ),
        (( -C, thick-a, thick*0.9), ),
        (( -C, thick-clad, thick*0.9), ),
        (( -C, (thick-clad)/2, thick*0.9), ),
        (( -C, 0.0, thick*0.9), ),

        (( 0.0, thick, thick*0.9), ),
        (( 0.0, thick-a, thick*0.9), ),
        (( 0.0, thick-clad, thick*0.9), ),
        (( 0, (thick-clad)/2, thick*0.9), ),
        (( 0, 0.0, thick*0.9), ),

        (( clad, thick, thick*0.9), ),
        (( clad, thick-clad, thick*0.9), ),
        (( clad, (thick-clad)/2, thick*0.9), ),
        (( clad, 0.0, thick*0.9), ),
    
        (( thick, thick, thick*0.9), ),
        (( thick, thick-clad, thick*0.9), ),
        (( thick, (thick-clad)/2, thick*0.9), ),
        (( thick, 0.0, thick*0.9), ),

        (( thick+LxPositive*0.35, thick, thick*0.9), ),
        (( thick+LxPositive*0.35, thick-clad, thick*0.9), ),
        (( thick+LxPositive*0.35, (thick-clad)/2, thick*0.9), ),
        (( thick+LxPositive*0.35, 0.0, thick*0.9), ),

        (( thick+LxPositive, thick, thick*0.9), ),
        (( thick+LxPositive, thick-clad, thick*0.9), ),
        (( thick+LxPositive, (thick-clad)/2, thick*0.9), ),
        (( thick+LxPositive, 0.0, thick*0.9), ),      

    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedthickZ, constraint=FINER)

#####################

pickedEdges = Edges.findAt(
    (( -C, thick, thick*1.01), ), 
    (( -C, thick-a, thick*1.01), ),      
    (( -C, thick-clad, thick*1.01), ),      
    (( -C, (thick-clad)/2, thick*1.01), ),      
    (( -C, 0.0, thick*1.01), ),   

    (( 0.0, thick, thick*1.01), ), 
    (( 0.0, thick-a, thick*1.01), ),      
    (( 0.0, thick-clad, thick*1.01), ),      
    (( 0.0, (thick-clad)/2, thick*1.01), ),      
    (( 0.0, 0.0, thick*1.01), ), 

    (( clad, thick, thick*1.01), ), 
    (( clad, thick-clad, thick*1.01), ),      
    (( clad, (thick-clad)/2, thick*1.01), ),      
    (( clad, 0.0, thick*1.01), ), 

    (( thick, thick, thick*1.01), ), 
    (( thick, thick-clad, thick*1.01), ),      
    (( thick, (thick-clad)/2, thick*1.01), ),      
    (( thick, 0.0, thick*1.01), ),

    (( thick+LxPositive*0.35, thick, thick*1.01), ), 
    (( thick+LxPositive*0.35, thick-clad, thick*1.01), ),      
    (( thick+LxPositive*0.35, (thick-clad)/2, thick*1.01), ),      
    (( thick+LxPositive*0.35, 0.0, thick*1.01), ),

    (( thick+LxPositive, thick, thick*1.01), ), 
    (( thick+LxPositive, thick-clad, thick*1.01), ),      
    (( thick+LxPositive, (thick-clad)/2, thick*1.01), ),      
    (( thick+LxPositive, 0.0, thick*1.01), ),
    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv1, constraint=FINER)

pickedEdges = Edges.findAt(

    (( -C, thick, thick+LDiv1*1.01), ), 
    (( -C, thick-a, thick+LDiv1*1.01), ),      
    (( -C, thick-clad, thick+LDiv1*1.01), ),      
    (( -C, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( -C, 0.0, thick+LDiv1*1.01), ),   

    (( 0.0, thick, thick+LDiv1*1.01), ), 
    (( 0.0, thick-a, thick+LDiv1*1.01), ),      
    (( 0.0, thick-clad, thick+LDiv1*1.01), ),      
    (( 0.0, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( 0.0, 0.0, thick+LDiv1*1.01), ), 

    (( clad, thick, thick+LDiv1*1.01), ), 
    (( clad, thick-clad, thick+LDiv1*1.01), ),      
    (( clad, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( clad, 0.0, thick+LDiv1*1.01), ), 

    (( thick, thick, thick+LDiv1*1.01), ), 
    (( thick, thick-clad, thick+LDiv1*1.01), ),      
    (( thick, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( thick, 0.0, thick+LDiv1*1.01), ),

    (( thick+LxPositive*0.35, thick, thick+LDiv1*1.01), ), 
    (( thick+LxPositive*0.35, thick-clad, thick+LDiv1*1.01), ),      
    (( thick+LxPositive*0.35, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( thick+LxPositive*0.35, 0.0, thick+LDiv1*1.01), ),

    (( thick+LxPositive, thick, thick+LDiv1*1.01), ), 
    (( thick+LxPositive, thick-clad, thick+LDiv1*1.01), ),      
    (( thick+LxPositive, (thick-clad)/2, thick+LDiv1*1.01), ),      
    (( thick+LxPositive, 0.0, thick+LDiv1*1.01), ),
    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv2, constraint=FINER)


pickedEdges = Edges.findAt(

    (( -C, thick, thick+LDiv2*1.01), ), 
    (( -C, thick-a, thick+LDiv2*1.01), ),      
    (( -C, thick-clad, thick+LDiv2*1.01), ),      
    (( -C, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( -C, 0.0, thick+LDiv2*1.01), ),   

    (( 0.0, thick, thick+LDiv2*1.01), ), 
    (( 0.0, thick-a, thick+LDiv2*1.01), ),      
    (( 0.0, thick-clad, thick+LDiv2*1.01), ),      
    (( 0.0, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( 0.0, 0.0, thick+LDiv2*1.01), ), 

    (( clad, thick, thick+LDiv2*1.01), ), 
    (( clad, thick-clad, thick+LDiv2*1.01), ),      
    (( clad, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( clad, 0.0, thick+LDiv2*1.01), ), 

    (( thick, thick, thick+LDiv2*1.01), ), 
    (( thick, thick-clad, thick+LDiv2*1.01), ),      
    (( thick, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( thick, 0.0, thick+LDiv2*1.01), ),

    (( thick+LxPositive*0.35, thick, thick+LDiv2*1.01), ), 
    (( thick+LxPositive*0.35, thick-clad, thick+LDiv2*1.01), ),      
    (( thick+LxPositive*0.35, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( thick+LxPositive*0.35, 0.0, thick+LDiv2*1.01), ),

    (( thick+LxPositive, thick, thick+LDiv2*1.01), ), 
    (( thick+LxPositive, thick-clad, thick+LDiv2*1.01), ),      
    (( thick+LxPositive, (thick-clad)/2, thick+LDiv2*1.01), ),      
    (( thick+LxPositive, 0.0, thick+LDiv2*1.01), ),
    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv3, constraint=FINER)

pickedEdges = Edges.findAt(

    (( -C, thick, thick+LDiv3*1.01), ), 
    (( -C, thick-a, thick+LDiv3*1.01), ),      
    (( -C, thick-clad, thick+LDiv3*1.01), ),      
    (( -C, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( -C, 0.0, thick+LDiv3*1.01), ),   

    (( 0.0, thick, thick+LDiv3*1.01), ), 
    (( 0.0, thick-a, thick+LDiv3*1.01), ),      
    (( 0.0, thick-clad, thick+LDiv3*1.01), ),      
    (( 0.0, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( 0.0, 0.0, thick+LDiv3*1.01), ), 

    (( clad, thick, thick+LDiv3*1.01), ), 
    (( clad, thick-clad, thick+LDiv3*1.01), ),      
    (( clad, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( clad, 0.0, thick+LDiv3*1.01), ), 

    (( thick, thick, thick+LDiv3*1.01), ), 
    (( thick, thick-clad, thick+LDiv3*1.01), ),      
    (( thick, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( thick, 0.0, thick+LDiv3*1.01), ),

    (( thick+LxPositive*0.35, thick, thick+LDiv3*1.01), ), 
    (( thick+LxPositive*0.35, thick-clad, thick+LDiv3*1.01), ),      
    (( thick+LxPositive*0.35, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( thick+LxPositive*0.35, 0.0, thick+LDiv3*1.01), ),

    (( thick+LxPositive, thick, thick+LDiv3*1.01), ), 
    (( thick+LxPositive, thick-clad, thick+LDiv3*1.01), ),      
    (( thick+LxPositive, (thick-clad)/2, thick+LDiv3*1.01), ),      
    (( thick+LxPositive, 0.0, thick+LDiv3*1.01), ),
    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedL, constraint=FINER)

pickedEdges = Edges.findAt(
    (( -C, thick*0.1, 0.0), ), 
    (( -C, thick-clad*1.1, 0.0), ), 
    (( -C, thick*0.1, lbox), ), 
    (( -C, thick-clad*1.1, lbox), ), 
    (( -C, thick*0.1, thick), ), 
    (( -C, thick-clad*1.1, thick), ), 
    (( -C, thick*0.1, thick+LDiv1), ), 
    (( -C, thick-clad*1.1, thick+LDiv1), ), 
    (( -C, thick*0.1, thick+LDiv2), ), 
    (( -C, thick-clad*1.1, thick+LDiv2), ), 
    (( -C, thick*0.1, thick+LDiv3), ), 
    (( -C, thick-clad*1.1, thick+LDiv3), ), 
    (( -C, thick*0.1, L), ), 
    (( -C, thick-clad*1.1, L), ), 

    (( 0.0, thick*0.1, 0.0), ), 
    (( 0.0, thick-clad*1.1, 0.0), ), 
    (( 0.0, thick*0.1, lbox), ), 
    (( 0.0, thick-clad*1.1, lbox), ), 
    (( 0.0, thick*0.1, thick), ), 
    (( 0.0, thick-clad*1.1, thick), ), 
    (( 0.0, thick*0.1, thick+LDiv1), ), 
    (( 0.0, thick-clad*1.1, thick+LDiv1), ), 
    (( 0.0, thick*0.1, thick+LDiv2), ), 
    (( 0.0, thick-clad*1.1, thick+LDiv2), ), 
    (( 0.0, thick*0.1, thick+LDiv3), ), 
    (( 0.0, thick-clad*1.1, thick+LDiv3), ), 
    (( 0.0, thick*0.1, L), ), 
    (( 0.0, thick-clad*1.1, L), ), 

    (( clad, thick*0.1, 0.0), ), 
    (( clad, thick-clad*1.1, 0.0), ), 
    (( clad, thick*0.1, lbox), ), 
    (( clad, thick-clad*1.1, lbox), ), 
    (( clad, thick*0.1, thick), ), 
    (( clad, thick-clad*1.1, thick), ), 
    (( clad, thick*0.1, thick+LDiv1), ), 
    (( clad, thick-clad*1.1, thick+LDiv1), ), 
    (( clad, thick*0.1, thick+LDiv2), ), 
    (( clad, thick-clad*1.1, thick+LDiv2), ), 
    (( clad, thick*0.1, thick+LDiv3), ), 
    (( clad, thick-clad*1.1, thick+LDiv3), ), 
    (( clad, thick*0.1, L), ), 
    (( clad, thick-clad*1.1, L), ), 

     (( thick, thick*0.1, 0.0), ), 
    (( thick, thick-clad*1.1, 0.0), ), 
    (( thick, thick*0.1, lbox), ), 
    (( thick, thick-clad*1.1, lbox), ), 
    (( thick, thick*0.1, thick), ), 
    (( thick, thick-clad*1.1, thick), ), 
    (( thick, thick*0.1, thick+LDiv1), ), 
    (( thick, thick-clad*1.1, thick+LDiv1), ), 
    (( thick, thick*0.1, thick+LDiv2), ), 
    (( thick, thick-clad*1.1, thick+LDiv2), ), 
    (( thick, thick*0.1, thick+LDiv3), ), 
    (( thick, thick-clad*1.1, thick+LDiv3), ), 
    (( thick, thick*0.1, L), ), 
    (( thick, thick-clad*1.1, L), ), 

    (( thick+LxPositive*0.35, thick*0.1, 0.0), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, 0.0), ), 
    (( thick+LxPositive*0.35, thick*0.1, lbox), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, lbox), ), 
    (( thick+LxPositive*0.35, thick*0.1, thick), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, thick), ), 
    (( thick+LxPositive*0.35, thick*0.1, thick+LDiv1), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, thick+LDiv1), ), 
    (( thick+LxPositive*0.35, thick*0.1, thick+LDiv2), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, thick+LDiv2), ), 
    (( thick+LxPositive*0.35, thick*0.1, thick+LDiv3), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, thick+LDiv3), ), 
    (( thick+LxPositive*0.35, thick*0.1, L), ), 
    (( thick+LxPositive*0.35, thick-clad*1.1, L), ), 

     (( thick+LxPositive, thick*0.1, 0.0), ), 
    (( thick+LxPositive, thick-clad*1.1, 0.0), ), 
    (( thick+LxPositive, thick*0.1, lbox), ), 
    (( thick+LxPositive, thick-clad*1.1, lbox), ), 
    (( thick+LxPositive, thick*0.1, thick), ), 
    (( thick+LxPositive, thick-clad*1.1, thick), ), 
    (( thick+LxPositive, thick*0.1, thick+LDiv1), ), 
    (( thick+LxPositive, thick-clad*1.1, thick+LDiv1), ), 
    (( thick+LxPositive, thick*0.1, thick+LDiv2), ), 
    (( thick+LxPositive, thick-clad*1.1, thick+LDiv2), ), 
    (( thick+LxPositive, thick*0.1, thick+LDiv3), ), 
    (( thick+LxPositive, thick-clad*1.1, thick+LDiv3), ), 
    (( thick+LxPositive, thick*0.1, L), ), 
    (( thick+LxPositive, thick-clad*1.1, L), ), 
)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseedthick, constraint=FINER)

# ##################################
# #################################
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,  secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)

cells = Cells.getByBoundingBox(xMin=-20000, xMax=20000,yMin=-20000,yMax=20000,zMin=-20000,zMax=20000)
pickedRegions =(cells,)
PipePart.setElementType(regions=pickedRegions, elemTypes=(elemType1,elemType2))
PipePart.generateMesh()

# ###################

 ### STEP

if LoadingType == 'Bending':    
    if Non_Lin_Geom == 'ON':
        PipeModel.StaticStep(name='Bending_Step', previous='Initial', nlgeom=ON,  description='Remote Bending', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    else:
        PipeModel.StaticStep(name='Bending_Step', previous='Initial',   description='Remote Bending', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)

elif LoadingType == 'Int_Pressure':
    if Non_Lin_Geom == 'ON':
        PipeModel.StaticStep(name='Int_Pressure_Step', previous='Initial', nlgeom=ON,  description='Internal Pressure', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    else:
        PipeModel.StaticStep(name='Int_Pressure_Step', previous='Initial',   description='Internal Pressure', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    
elif LoadingType == 'Tensile':
    if Non_Lin_Geom == 'ON':
        PipeModel.StaticStep(name='Tensile_Step', previous='Initial', nlgeom=ON ,  description='Axial Tensile', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    else:
        PipeModel.StaticStep(name='Tensile_Step', previous='Initial' ,  description='Axial Tensile', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)

elif LoadingType == 'Torsion':
    if Non_Lin_Geom == 'ON':
        PipeModel.StaticStep(name='Torsion_Step', previous='Initial', nlgeom=ON ,   description='Torsion', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    else:
        PipeModel.StaticStep(name='Torsion_Step', previous='Initial',   description='Torsion', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
elif LoadingType == 'ALL':
    if Non_Lin_Geom == 'ON':
        PipeModel.StaticStep(name='ALL_Step', previous='Initial', nlgeom=ON ,  description='ALL together', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
    else:
        PipeModel.StaticStep(name='ALL_Step', previous='Initial',   description='ALL together', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)



CR1 = a 

PipeAssambly = PipeModel.rootAssembly
PipeAssambly.DatumCsysByDefault(CARTESIAN)
PipePart = PipeModel.parts['PIPE']
PipeAssambly.Instance(name='Pipe_Instance', part=PipePart, dependent=ON)

n1 = PipeAssambly.instances['Pipe_Instance'].nodes
nodes1 = n1[0:10000000000]
PipeAssambly.Set(nodes=nodes1, name='ALLN')

# Movement of Instance- Assembly

PipeAssambly.translate(instanceList=('Pipe_Instance', ), vector=(C, -thick, 0))
PipeAssambly.rotate(instanceList=('Pipe_Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.2), angle=-90.0)
PipeAssambly.rotate(instanceList=('Pipe_Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(-1.2, 0.0, 0.0), angle=180.0)

# # #### Reference Point
PipeAssambly = PipeModel.rootAssembly
PipeAssambly.ReferencePoint(point=(0, 0.0, -L))

PipeModel.setValues(noPartsInputFile=ON)

########### Surface Set for Internal Pressure

BALA = C + thick + LxPositive

PipeAssambly = PipeModel.rootAssembly
Faces = PipeAssambly.instances['Pipe_Instance'].faces
side1Faces1 = Faces.findAt(

    ((0.0, BALA - Epsilon, -Epsilon), ), 
    ((0.0, BALA - Epsilon, -lbox-Epsilon), ), 
    ((0.0, BALA - Epsilon, -thick-Epsilon), ), 
    ((0.0, BALA - Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0, BALA - Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0, BALA - Epsilon, -thick-LDiv3-Epsilon), ), 

    ((0.0, C + thick + Epsilon, -Epsilon), ), 
    ((0.0, C + thick + Epsilon, -lbox-Epsilon), ), 
    ((0.0, C + thick + Epsilon, -thick-Epsilon), ), 
    ((0.0, C + thick + Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0, C + thick + Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0, C + thick + Epsilon, -thick-LDiv3-Epsilon), ), 

    ((0.0, C + thick - Epsilon, -Epsilon), ), 
    ((0.0, C + thick - Epsilon, -lbox-Epsilon), ), 
    ((0.0, C + thick - Epsilon, -thick-Epsilon), ), 
    ((0.0, C + thick - Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0, C + thick - Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0, C + thick - Epsilon, -thick-LDiv3-Epsilon), ), 

    ((0.0, C +a+ Epsilon, -lbox+Epsilon), ), 
    ((0.0, C +a+ Epsilon, -lbox-Epsilon), ), 
    ((0.0, C +a+ Epsilon, -thick-Epsilon), ), 
    ((0.0, C +a+ Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0, C +a+ Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0, C +a+ Epsilon, -thick-LDiv3-Epsilon), ), 


     ((0.0, C + Epsilon, -Epsilon), ), 
    ((0.0, C + Epsilon, -lbox-Epsilon), ), 
    ((0.0, C + Epsilon, -thick-Epsilon), ), 
    ((0.0, C + Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0, C + Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0, C + Epsilon, -thick-LDiv3-Epsilon), ), 

     ((0.0,  Epsilon, -Epsilon), ), 
    ((0.0,  Epsilon, -lbox-Epsilon), ), 
    ((0.0,  Epsilon, -thick-Epsilon), ), 
    ((0.0,  Epsilon, -thick-LDiv1-Epsilon), ), 
    ((0.0,  Epsilon, -thick-LDiv2-Epsilon), ), 
    ((0.0,  Epsilon, -thick-LDiv3-Epsilon), ), 
    )
PipeAssambly.Surface(side1Faces=side1Faces1, name='SurfacePressureSet')


side1Faces2 = Faces.findAt(
        ((-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-thick+Epsilon, C*0.1, -L), ), 

        ((-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-thick+Epsilon, C*1.01, -L), ), 

        ((-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick-Epsilon, -L), ),

        ((-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+Epsilon, -L), ), 

        ((-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        
)
PipeAssambly.Surface(side1Faces=side1Faces2, name='Tensile_Faces_Set')

########## MPC

PipeAssembly = PipeModel.rootAssembly
RefPoints = PipeAssembly.referencePoints
refPoints1=(RefPoints[5], )
region1=PipeAssembly.Set(referencePoints=refPoints1, name='Ref_Point_Set')
Faces = PipeAssembly.instances['Pipe_Instance'].faces
Edges = PipeAssembly.instances['Pipe_Instance'].edges

faces1 = Faces.findAt(
        ((-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-thick+Epsilon, C*0.1, -L), ), 

        ((-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-thick+Epsilon, C*1.01, -L), ), 

        ((-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick-Epsilon, -L), ),

        ((-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+Epsilon, -L), ), 

        ((-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
       
)

if LoadingType == 'Bending':    
    region2=PipeAssembly.Set( faces=faces1, name='MPC_Faces_Set')
    PipeModel.MultipointConstraint(name='MPC_Ref_Point', controlPoint=region1,surface=region2, mpcType=BEAM_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)

if LoadingType == 'Torsion':    
    region2=PipeAssembly.Set( faces=faces1, name='MPC_Faces_Set')
    PipeModel.MultipointConstraint(name='MPC_Ref_Point', controlPoint=region1,surface=region2, mpcType=BEAM_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)
if LoadingType == 'ALL':    
    region2=PipeAssembly.Set( faces=faces1, name='MPC_Faces_Set')
    PipeModel.MultipointConstraint(name='MPC_Ref_Point', controlPoint=region1,surface=region2, mpcType=BEAM_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)
  
###############################
################################
PipeAssembly = PipeModel.rootAssembly
Faces = PipeAssembly.instances['Pipe_Instance'].faces

BALA = C+thick+LxPositive
faces1 = Faces.findAt(

        ((-a*0.01, 0.0, -lbox*0.01), ), 
        ((-a*1.01, 0.0, -lbox*0.99), ), 
        ((-clad*1.01, 0.0, -lbox*0.99), ), 
        ((-thick*0.99, 0.0, -lbox*0.99), ), 
    
        ((-a+A1*0.01, 0.0, -A1*0.9), ),        
        ((-a-A1*0.9, 0.0, -Epsilon), ), 

        ((-a+A2*0.01, 0.0, -A2*0.9), ), 
        ((-a-A2*0.9, 0.0, -Epsilon), ), 

        ((-a+A3*0.01, 0.0, -A3*0.9), ),
         ((-a-A3*0.9, 0.0, -Epsilon), ), 

        ((-a+A4*0.01, 0.0, -A4*0.9), ),
        ((-a-A4*0.9, 0.0, -Epsilon), ), 

        ((-a+A5*0.01, 0.0, -A5*0.9), ), 
        ((-a-A5*0.9, 0.0, -Epsilon), ), 

        ((-a+A6*0.01, 0.0, -A6*0.9), ), 
        ((-a-A6*0.9, 0.0, -Epsilon), ), 

        ((-a+A7*0.01, 0.0, -A7*0.9), ), 
        ((-a-A7*0.9, 0.0, -Epsilon), ), 

        ((-a+A8*0.01, 0.0, -A8*0.9), ),
         ((-a-A8*0.9, 0.0, -Epsilon), ), 

        ((-a+A9*0.01, 0.0, -A9*0.9), ), 
         ((-a-A9*0.9, 0.0, -Epsilon), ), 

        
        # ######
          ((-a*0.99, 0.0, -lbox*1.01), ), 
        ((-a*1.01, 0.0, -lbox*1.01), ), 
        ((-clad*1.01, 0.0, -lbox*1.01), ), 
        ((-thick*0.99, 0.0, -lbox*1.01), ), 

        ((-a*0.99, 0.0, -thick*0.99), ), 
        ((-a*1.01, 0.0, -thick*0.99), ), 
        ((-clad*1.01, 0.0, -thick*0.99), ), 
        ((-thick*0.99, 0.0, -thick*0.99), ), 

        ((-a*0.99, 0.0, -thick-Epsilon), ), 
        ((-a*1.01, 0.0, -thick-Epsilon), ), 
        ((-clad*1.01, 0.0, -thick-Epsilon), ), 
        ((-thick*0.99, 0.0, -thick-Epsilon), ), 

        ((-a*0.99, 0.0, -thick-LDiv1-Epsilon), ), 
        ((-a*1.01, 0.0, -thick-LDiv1-Epsilon), ), 
        ((-clad*1.01, 0.0, -thick-LDiv1-Epsilon), ), 
        ((-thick*0.99, 0.0, -thick-LDiv1-Epsilon), ), 

        ((-a*0.99, 0.0, -thick-LDiv2-Epsilon), ), 
        ((-a*1.01, 0.0, -thick-LDiv2-Epsilon), ), 
        ((-clad*1.01, 0.0, -thick-LDiv2-Epsilon), ), 
        ((-thick*0.99, 0.0, -thick-LDiv2-Epsilon), ), 

        ((-a*0.99, 0.0, -thick-LDiv3-Epsilon), ), 
        ((-a*1.01, 0.0, -thick-LDiv3-Epsilon), ), 
        ((-clad*1.01, 0.0, -thick-LDiv3-Epsilon), ), 
        ((-thick*0.99, 0.0, -thick-LDiv3-Epsilon), ), 
     
#         ############################################

        ((-a*0.01, BALA, -lbox*0.01), ), 
        ((-a*1.01, BALA, -lbox*0.99), ), 
        ((-clad*1.01, BALA, -lbox*0.99), ), 
        ((-thick*0.99, BALA, -lbox*0.99), ), 

        ((-a*0.99, BALA, -lbox*1.01), ), 
        ((-a*1.01, BALA, -lbox*1.01), ), 
        ((-clad*1.01, BALA, -lbox*1.01), ), 
        ((-thick*0.99, BALA, -lbox*1.01), ), 

        ((-a*0.99, BALA, -thick*0.99), ), 
        ((-a*1.01, BALA, -thick*0.99), ), 
        ((-clad*1.01, BALA, -thick*0.99), ), 
        ((-thick*0.99, BALA, -thick*0.99), ), 

        ((-a*0.99, BALA, -thick-Epsilon), ), 
        ((-a*1.01, BALA, -thick-Epsilon), ), 
        ((-clad*1.01, BALA, -thick-Epsilon), ), 
        ((-thick*0.99, BALA, -thick-Epsilon), ), 

        ((-a*0.99, BALA, -thick-LDiv1-Epsilon), ), 
        ((-a*1.01, BALA, -thick-LDiv1-Epsilon), ), 
        ((-clad*1.01, BALA, -thick-LDiv1-Epsilon), ), 
        ((-thick*0.99, BALA, -thick-LDiv1-Epsilon), ), 

        ((-a*0.99, BALA, -thick-LDiv2-Epsilon), ), 
        ((-a*1.01, BALA, -thick-LDiv2-Epsilon), ), 
        ((-clad*1.01, BALA, -thick-LDiv2-Epsilon), ), 
        ((-thick*0.99, BALA, -thick-LDiv2-Epsilon), ), 

        ((-a*0.99, BALA, -thick-LDiv3-Epsilon), ), 
        ((-a*1.01, BALA, -thick-LDiv3-Epsilon), ), 
        ((-clad*1.01, BALA, -thick-LDiv3-Epsilon), ), 
        ((-thick*0.99, BALA, -thick-LDiv3-Epsilon), ), 
        
      
)
faces2 = Faces.findAt(
        ((-a-A1*0.01, 0.0, -A1*0.9), ), 
        ((-a+A10*0.01, 0.0, -A10*0.9), ), 
        ((-a-A10*0.01, 0.0, -A10*0.9), ), 

)
Faces1 = faces1 + faces2

region = PipeAssembly.Set(faces=Faces1, name='Set_Faces_XSYM')

if LoadingType == 'Bending':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Bending_Step', region=region, localCsys=None)
elif LoadingType == 'Int_Pressure':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Int_Pressure_Step', region=region, localCsys=None)
elif LoadingType == 'Tensile':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Tensile_Step', region=region, localCsys=None)
elif LoadingType == 'Torsion':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Torsion_Step', region=region, localCsys=None)
elif LoadingType == 'ALL':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='ALL_Step', region=region, localCsys=None)




PipeAssembly = PipeModel.rootAssembly
Faces = PipeAssembly.instances['Pipe_Instance'].faces
KK = 0

faces1 = Faces.findAt(
    ((-thick*0.99, C*0.01, 0.0), ), 
    ((-thick*0.99, C+lbox*0.01, 0.0), ), 
    ((-thick*0.99, C+lbox*1.01, 0.0), ), 
    ((-thick*0.99, C+thick*0.99, 0.0), ), 
    ((-thick*0.99, C+thick*1.01, 0.0), ), 
    ((-thick*0.99, C+thick+LxPositive*0.99, 0.0), ), 

    ((-clad*1.01, C*0.01, 0.0), ), 
    ((-clad*1.01, C+lbox*0.01, 0.0), ), 
    ((-clad*1.01, C+lbox*1.01, 0.0), ), 
    ((-clad*1.01, C+thick*0.99, 0.0), ), 
    ((-clad*1.01, C+thick*1.01, 0.0), ), 
    ((-clad*1.01, C+thick+LxPositive*0.99, 0.0), ), 

    ((-clad*0.99, C*0.01, 0.0), ), 
    ((-clad*0.99, C+lbox*0.01, 0.0), ), 
    ((-clad*0.99, C+lbox*1.01, 0.0), ), 
    ((-clad*0.99, C+thick*0.99, 0.0), ), 
    ((-clad*0.99, C+thick*1.01, 0.0), ), 
    ((-clad*0.99, C+thick+LxPositive*0.99, 0.0), ), 

    ((-a-A10*0.99, C*0.01, 0.0), ), 
    ((-a-A10*0.99, C+Epsilon, 0.0), ), 

    
    ((-a-A9*0.99, C*0.01, 0.0), ), 
    ((-a-A9*0.99, C+Epsilon, 0.0), ), 

    ((-a-A8*0.99, C*0.01, 0.0), ), 
    ((-a-A8*0.99, C+Epsilon, 0.0), ), 

    ((-a-A7*0.99, C*0.01, 0.0), ), 
    ((-a-A7*0.99, C+Epsilon, 0.0), ), 

    ((-a-A6*0.99, C*0.01, 0.0), ), 
    ((-a-A6*0.99, C+Epsilon, 0.0), ), 

    ((-a-A5*0.99, C*0.01, 0.0), ), 
    ((-a-A5*0.99, C+Epsilon, 0.0), ), 

    ((-a-A4*0.99, C*0.01, 0.0), ), 
    ((-a-A4*0.99, C+Epsilon, 0.0), ), 

    ((-a-A3*0.99, C*0.01, 0.0), ), 
    ((-a-A3*0.99, C+Epsilon, 0.0), ), 

    ((-a-A2*0.99, C*0.01, 0.0), ), 
    ((-a-A2*0.99, C+Epsilon, 0.0), ),     
  
    )
faces2 = Faces.findAt(
    ((-a-A1*0.99, C*0.01, 0.0), ), 
    ((-a-A1*0.99, C+Epsilon, 0.0), ),     
)
Faces1 = faces1+faces2


region = PipeAssembly.Set(faces=Faces1, name='Set_Faces_ZSYM')

if LoadingType == 'Bending':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Bending_Step', region=region, localCsys=None)
elif LoadingType == 'Int_Pressure':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Int_Pressure_Step', region=region, localCsys=None)
elif LoadingType == 'Tensile':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Tensile_Step', region=region, localCsys=None)
elif LoadingType == 'Torsion':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Torsion_Step', region=region, localCsys=None)
elif LoadingType == 'ALL':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='ALL_Step', region=region, localCsys=None)



# ####################### Applied LOADs
if LoadingType == 'Bending':
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='Y_Fixed_RefP', createStepName='Bending_Step', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)

    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='Bending', createStepName='Bending_Step', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=-Applied_Bending, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        fieldName='', localCsys=None)

if LoadingType == 'Int_Pressure':
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.surfaces['SurfacePressureSet']
    PipeModel.Pressure(name='PressureInternal', createStepName='Int_Pressure_Step', region=region, 
    distributionType=UNIFORM, field='', magnitude=Applied_Int_Pressure, amplitude=UNSET)

if LoadingType == 'Tensile':
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.surfaces['Tensile_Faces_Set']
    PipeModel.Pressure(name='Tensile_Pressure', createStepName='Tensile_Step', region=region, distributionType=UNIFORM, 
    field='', magnitude= -Applied_Tensile_Pressure, amplitude=UNSET)

if LoadingType == 'Torsion':
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='Y_Fixed_RefP', createStepName='Torsion_Step', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)

    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='TorsionDisp', createStepName='Torsion_Step', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3= -Applied_Torsion, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        fieldName='', localCsys=None)

if LoadingType == 'ALL':
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='Y_Fixed_RefP', createStepName='ALL_Step', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    ### Bending
    
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='Bending', createStepName='ALL_Step', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=-Applied_Bending, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        fieldName='', localCsys=None)
    ### Torsion
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.sets['Ref_Point_Set']
    PipeModel.DisplacementBC(name='TorsionDisp', createStepName='ALL_Step', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3= -Applied_Torsion, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        fieldName='', localCsys=None)

    ### Pressure Tensile
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.surfaces['Tensile_Faces_Set']
    PipeModel.Pressure(name='Tensile_Pressure', createStepName='ALL_Step', region=region, distributionType=UNIFORM, 
    field='', magnitude= -Applied_Tensile_Pressure, amplitude=UNSET)
    ### Internal Pressure
    PipeAssembly = PipeModel.rootAssembly
    region = PipeAssembly.surfaces['SurfacePressureSet']
    PipeModel.Pressure(name='PressureInternal', createStepName='ALL_Step', region=region, 
    distributionType=UNIFORM, field='', magnitude=Applied_Int_Pressure, amplitude=UNSET)


# ######## Crack definition

PipeAssembly = PipeModel.rootAssembly
PipeAssembly.makeIndependent(instances=(PipeAssembly.instances['Pipe_Instance'], ))

PipeAssembly = PipeModel.rootAssembly
Edges = PipeAssembly.instances['Pipe_Instance'].edges
edges1 = Edges.findAt(
	((-CR1, C*0.1, 0.0), ), 
	((-CR1*cos(pi/8), C+CR1*sin(pi/8), 0.0), ),
    ((-CR1*cos(3*pi/8), C+CR1*sin(3*pi/8), 0.0), )

	)
crackFront = regionToolset.Region(edges=edges1)
crackTip = regionToolset.Region(edges=edges1)
PipeAssembly.engineeringFeatures.ContourIntegral(name='Crack_Interface', symmetric=ON,  crackFront=crackFront, crackTip=crackTip, 
    extensionDirectionMethod=CRACK_NORMAL, crackNormal=((0.0, 0.0, 0.0), (0.0,  0.0, -1.0)), midNodePosition=0.25, collapsedElementAtTip=DUPLICATE_NODES) 
####################################
# #### Creating a set of all nodes for NMAP function
n1 = PipeAssambly.instances['Pipe_Instance'].nodes
nodes1 = n1[0:10000000000]
PipeAssambly.Set(nodes=nodes1, name='ALLN1')



############ OutPut J-Integral

if LoadingType == 'Bending':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Bending_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=40)  
    ###History OUTPUT 
    ###OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='Bending_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)  
    #### History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Bending_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Bending_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'Int_Pressure':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Int_Pressure_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)  
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Int_Pressure_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)
    # ###History OUTPUT 
    ##OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='Int_Pressure_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ##History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Int_Pressure_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Int_Pressure_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'Tensile':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Tensile_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30) 
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Tensile_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)
    # ###History OUTPUT 
    ###OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='Tensile_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ##History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Tensile_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Tensile_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'Torsion':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Torsion_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30) 
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Torsion_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)
    # ###History OUTPUT 
    ###OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='Torsion_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ##History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Torsion_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Torsion_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'ALL':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='ALL_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30) 
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='ALL_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)
    # ###History OUTPUT 
    ###OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='ALL_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ##History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='ALL_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='ALL_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)


del PipeModel.fieldOutputRequests['F-Output-1']
del PipeModel.historyOutputRequests['H-Output-1']


# # # # ###############################

PipeAssembly = PipeModel.rootAssembly
PipeAssembly.rotate(instanceList=('Pipe_Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)
PipeAssembly.translate(instanceList=('Pipe_Instance', ), vector=(-thick, BALA, 0.0))

radvar = str(radius)
vari = str(180/(radius*pi))
PipeModel.keywordBlock.synchVersions(storeNodesAndElements=False)
PipeModel.keywordBlock.insert(35, 
  '\n*NMAP, NSET=ALLN1, TYPE=RECTANGULAR\n%s,0.,0.,999999.,0.,0.\n%s,1.,0.\n1.,1.,1.\n*NMAP, NSET=ALLN1, TYPE=CYLINDRICAL\n0.,0.,0.,0.,0.,1.\n0.,1.,0.\n1.,%s,1.\n*NMAP, NSET=ALLN1, TYPE=RECTANGULAR\n0.,0.,0., 1.,0.,0.\n1.,1.,0.' %(radvar,radvar,vari) )

mdb.Job(name='Int_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick))), model='Pipe_Model_Internal', description='', type=ANALYSIS, 
     atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
     memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
     explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
     modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
     scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=NCPUs, 
     activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=NCPUs)
PipeModel.setValues(noPartsInputFile=ON)
mdb.jobs['Int_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].writeInput(consistencyChecking=OFF)
#mdb.jobs['Int_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].submit(consistencyChecking=OFF)
#mdb.jobs['Int_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].waitForCompletion()

