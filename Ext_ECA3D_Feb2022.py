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
###  0) Radius of Pipe
##### EXTERNAL RADIUS OF PIPE
# radius = 103
radius = 103 * 2
# radius = 103 * 3
# radius = 103 * 4
# radius = 103 * 5

###############################
########################
LoadingType = 'Bending'
# LoadingType = 'Int_Pressure'
# LoadingType = 'Tensile'
# LoadingType = 'Torsion'


Applied_Bending = 0.24
Applied_Int_Pressure = 10
Applied_Tensile_Pressure = 5
#### 1 )  NUMBER OF CPU USED FOR ANALYSIS (Depends on your computer)
NCPUs = 6

#### 2)  SELECT TYPE OF THE MATERIAL DEFINITION USED IN THE MODEL
MaterialType = 'Plastic-Deformation'  
# MaterialType = 'Elastic-Plastic'

##### 8 ) Select weld angle
N_Mat = 5 
# N_Mat = 10
# N_Mat = 20 


###################

# a= 2.06
# a= 4.12
# a=6.18
# a = 8.24
# a = 10.3
a = 12.36

radius = 103.0
thick = 20.6
clad = 2


# theta_pi = 0.04
# theta_pi = 0.12
theta_pi = 0.20


D_ext = radius * 2
L = 3 * D_ext 
C = theta_pi * 2 * D_ext
LxPositive=((2*pi*radius)/2.0)-(thick+C)
LDiv1 = thick

LDiv2 = thick + L * 0.2
LDiv3 = thick + L * 0.6 


Weld_Angle = 0
# Weld_Angle = 20
Weld_Angle = 45
# Weld_Angle = 60

# KeyHole = 0.05
KeyHole = 0.0

A1 = a/75
A2 = 2* A1
A3 = 3 * A1
A4 = 4 * A1
A5 = 5 * A1
A6 = 6 * A1
A7 = 7 * A1
A8 = 8 * A1
A9 = 9 * A1
A10 = 10 * A1

Weld_Fixed_Width = 12 * A1
a2 = Weld_Fixed_Width
Wedge = a/250

########### Just a dummy variable to be used in coordinates
Epsilon = A1/2
############################
NseedSpider = 4
Nseedcurve = 3
NseedC = 10 
Nseeda = 4
NseedWeldFixed = 4
Nseedthick = 4
Nseeedclad = 3
NseeedLDiv1 = 5
NseeedLDiv2 = 5
NseeedLDiv3 = 5
NseeedL = 5
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
                [Yield_Stress_Base,Yield_Stress_Base]
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


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################


#---------
#create the model
mdb.models.changeKey(fromName='Model-1',toName='Pipe_Model_External')
PipeModel=mdb.models['Pipe_Model_External']

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


#######
##################################################################################################################
########################### Cut Extrude to make a KeyHole ########################################################
##################################################################################################################

if KeyHole > 0.0:
    Transform= (0.0, -0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, thick/2, thick/2)
    sketch_2 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=112.92, gridSpacing=2.82, transform=Transform)
    sketch_2.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=sketch_2, filter=COPLANAR_EDGES)

    KK= thick/2-KeyHole
    sketch_2.Line(point1=(-thick/2, thick/2), point2=(-thick/2+KeyHole, thick/2))
    sketch_2.Line(point1=(-thick/2+KeyHole, thick/2), point2=(-thick/2+KeyHole, thick/2-a))
    sketch_2.ArcByCenterEnds(center=(-thick/2, thick/2-a), point1=(-thick/2+KeyHole, thick/2-a), point2=(-thick/2,  thick/2-a-KeyHole), direction=CLOCKWISE)
    sketch_2.Line(point1=(-thick/2, thick/2-a-KeyHole), point2=(-thick/2, thick/2))
    sketch_2.unsetPrimaryObject()

    PipePart = PipeModel.parts['PIPE']
    Edges = PipePart.edges

    edges = Edges.findAt((( a * np.cos(pi/4), thick - a * np.cos(pi/4), 0.0), ))
    pathEdges=edges
    SketchPlane=Faces.findAt(coordinates=(0.0 , thick  - a , thick - a ))
    SketchEdge=Edges.findAt(coordinates=( 0.0 , thick / 4 , thick))

    PipePart.CutSweep(path=pathEdges, sketchPlane=SketchPlane, sketchUpEdge=SketchEdge,  sketchOrientation=RIGHT, profile=sketch_2, flipSweepDirection=ON)
    del PipeModel.sketches['__profile__']

# #############################################################
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

# ####################################
##### Extruding a crack curve

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
pickedCells = Cells.findAt(((thick/2, thick/2, thick/3), ))

pickedEdges =(
    Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), KeyHole)), 
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, thick, thick/3))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,   cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### Extruding Spider mesh curves

### A1
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A1*cos(pi/4), A1*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(pi/4), A1*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

## A2
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A2*cos(pi/4), A2*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(pi/4), A2*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A3
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A3*cos(pi/4), A3*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(pi/4), A3*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A4
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A4*cos(pi/4), A4*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(pi/4), A4*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A5

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A5*cos(pi/4), A5*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(pi/4), A5*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A6
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A6*cos(pi/4), A6*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(pi/4), A6*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A7
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A7*cos(pi/4), A7*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(pi/4), A7*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)
## A8

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A8*cos(pi/4), A8*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(pi/4), A8*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A9

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A9*cos(pi/4), A9*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(pi/4), A9*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# ## A10

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, Epsilon, Epsilon), ), 
    ((Epsilon, thick-Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A10*cos(pi/4), A10*cos(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A10*cos(pi/4), A10*cos(pi/4)))
    )

ExtrudeEdge = Edges.findAt(coordinates=(a*cos(pi/4),thick-a*cos(pi/4), KeyHole))
PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeEdge, cells=pickedCells, edges=pickedEdges)

# #######################################################################
# ##### Sketching Clad line $ Weld_Fixed_Width


################### External Crack
PipePart = PipeModel.parts['PIPE']
Faces = PipePart.faces
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick-Epsilon, thick-Epsilon))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, Epsilon, 0.0))
Origin_Sketch_plane = (0.0, 0.0, 0.0)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)

PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(0.0, clad), point2=(thick, clad))
sketch_1.Line(point1=(Weld_Fixed_Width, thick), point2=(Weld_Fixed_Width, 0.0))


pickedFaces = Faces.findAt(
    ((0.0, thick-Epsilon, thick-Epsilon), ), 
    ((0.0, Epsilon, thick-Epsilon), )
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


#############################################

### Extruding Clad line

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((0.0, Epsilon, thick-Epsilon), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(0.0, clad, thick-Epsilon)), 
    Edges.findAt(coordinates=(0.0, clad, KeyHole+ Epsilon)))
ExtrudeLine = Edges.findAt(coordinates=(thick*0.75, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

# # ## Extruding Fixed_Weld_Width line

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((Epsilon, thick-Epsilon, thick-Epsilon), ), 
    ((Epsilon, clad+Epsilon, thick-Epsilon), ), 
    ((Epsilon, clad-Epsilon, thick-Epsilon), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-Epsilon, Weld_Fixed_Width)), 
    Edges.findAt(coordinates=(0.0, clad-Epsilon, Weld_Fixed_Width)), 
    Edges.findAt(coordinates=(0.0, clad+Epsilon, Weld_Fixed_Width))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick*0.75, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)


##########################################################
##### Skethcin a2 to improve mesh

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(thick-Epsilon,clad+Epsilon, 0.0))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, Epsilon, 0.0))
Origin_Sketch_plane = (0.0, 0.0, 0.0) 
transformXY = PipePart.MakeSketchTransform(sketchPlane= Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(-(a+a2), thick), point2=(-((a+a2)), 0.0))


pickedFaces = Faces.findAt(
    ((thick-Epsilon,clad+Epsilon, 0.0), ), 
    ((thick-Epsilon,clad-Epsilon, 0.0), )
    )

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


# # ## Extruding clad a2

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(
    ((thick-Epsilon, thick-Epsilon, thick-Epsilon), ), 
    ((thick-Epsilon, thick-Epsilon, Weld_Fixed_Width-Epsilon), ), 
    ((thick-Epsilon, Epsilon, thick-Epsilon), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(a+a2, clad+Epsilon, 0.0)), 
    Edges.findAt(coordinates=(a+a2, clad-Epsilon, 0.0))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, 0.0, Epsilon))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

pickedCells = Cells.findAt(
    ((thick*0.95, clad*0.75, 0.0), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(a+a2, clad, Weld_Fixed_Width*0.3)), 
    )
ExtrudeLine = Edges.findAt(coordinates=(a+a2, clad*0.75, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#######################################################################################
##################### V-groove angle sketch

THETA = Weld_Angle/2 * pi/180
if Weld_Angle > 0 :
    PipePart = PipeModel.parts['PIPE']
    Faces = PipePart.faces
    Edges = PipePart.edges
    RR = thick/cos(THETA)
    rr = RR * sin(THETA)

    Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick-Epsilon, Weld_Fixed_Width+Epsilon ))
    Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, Epsilon, 0.0))
    Origin_Sketch_plane = (0.0, 0.0, 0.0)
    transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,  sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
    sketch_1 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=55.42, gridSpacing=1.38, transform= transformXY)
    sketch_1.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

    sketch_1.Line(point1=(Weld_Fixed_Width , 0.0), point2=(Weld_Fixed_Width+rr, thick))

    pickedFaces = Faces.findAt(
        ((0.0, thick-Epsilon, Weld_Fixed_Width+Epsilon), ), 
        ((0.0, clad-Epsilon, Weld_Fixed_Width+Epsilon), ), 
        ((0.0, clad+Epsilon, Weld_Fixed_Width+Epsilon), )
        )
    PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane,  faces=pickedFaces, sketch=sketch_1)
    sketch_1.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']

    ########################################################################
    ##### Extruding V-groove line
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    pickedCells = Cells.getByBoundingBox(xMin=0,xMax=thick,yMin=0.,yMax=thick,zMin=0.0,zMax=thick)

    RR1 = (clad*0.3) /cos(THETA)
    rr1 = RR1 * sin(THETA)
    RR2 = (clad*1.3) /cos(THETA)
    rr2 = RR2 * sin(THETA)
    RR3 = (thick-Epsilon)/cos(THETA)
    rr3 = RR3 * sin(THETA)

    pickedEdges =(
        Edges.findAt(coordinates=(0.0, clad*0.3, Weld_Fixed_Width+rr1)), 
        Edges.findAt(coordinates=(0.0, clad*1.3, Weld_Fixed_Width+rr2)), 
        Edges.findAt(coordinates=(0.0, thick-Epsilon, Weld_Fixed_Width+rr3))
        )
    ExtrudeLine = Edges.findAt(coordinates=(Epsilon, 0.0, 0.0))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#####################################################################
######### Sketching wedge

if KeyHole == 0:    
    PipePart = PipeModel.parts['PIPE']
    Faces = PipePart.faces
    Edges= PipePart.edges

    Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick-a-Wedge/3, Wedge*0.5))
    Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, thick-a-Wedge/3, 0.0))
    Origin_Sketch_plane = (0.0, thick-a, 0.0)
    transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
    sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=16.5, gridSpacing=0.41, transform=transformXY)
    sketch_1.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)
    sketch_1.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, Wedge), point2=(0.0, -Wedge), direction=CLOCKWISE)
    pickedFaces = Faces.findAt(
        ((0.0, thick-a+Wedge/3, Wedge*0.5), ),
        ((0.0, thick-a-Wedge/3, Wedge*0.5), )
        )
    ExtrudeLine = Edges.findAt(coordinates=(0.0, thick-a+Wedge/3, 0.0))
    PipePart.PartitionFaceBySketch(sketchUpEdge=ExtrudeLine, faces=pickedFaces, sketch=sketch_1)
    sketch_1.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']
    ### Extrude Wedge
    PipePart = PipeModel.parts['PIPE']
    Faces = PipePart.faces
    Edges= PipePart.edges
    pickedCells = Cells.findAt(
        (((a-A1*0.2)*cos(pi/4), thick-(a-A1*0.2)*cos(pi/4), Wedge*0.2), ), 
        (((a+A1*0.2)*cos(pi/4), thick-(a+A1*0.2)*cos(pi/4), Wedge*0.2 ), )
        )
    
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-a-Wedge*cos(pi/4), Wedge*cos(pi/4))), 
        Edges.findAt(coordinates=(0.0, thick-a+Wedge*cos(pi/4), Wedge*cos(pi/4)))
        )
    ExtrudeLine = Edges.findAt(coordinates=(a*cos(pi/4), thick-a*cos(pi/4), 0.0))
    PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)


##################################################################################
### Solid Extrude C 

PipePart = PipeModel.parts['PIPE']
Faces = PipePart.faces
Edges = PipePart.edges

Face_Sketch_plane = Faces.findAt(coordinates=(0.0, thick-Epsilon, thick-Epsilon))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, Epsilon, thick))
Origin_Sketch_plane = (0.0, 0.0, 0.0)
transformXY = PipePart.MakeSketchTransform(sketchPlane= Face_Sketch_plane , sketchUpEdge= Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=Origin_Sketch_plane)
sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart = PipeModel.parts['PIPE']
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(0.0, 0.0), point2=(0.0, thick))
sketch_1.Line(point1=(0.0, thick), point2=(thick, thick))
sketch_1.Line(point1=(thick, thick), point2=(thick, 0.0))
sketch_1.Line(point1=(thick, 0.0), point2=(0.0, 0.0))

PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=sketch_1, depth= C, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)

sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

######################################################################
#### Extrude Cut for KeyHole

if KeyHole > 0:
    PipePart = PipeModel.parts['PIPE']
    Faces = PipePart.faces
    Edges = PipePart.edges

    Face_Sketch_plane = Faces.findAt(coordinates=(-C, thick*0.3,thick*0.3))
    Edge_Sketch_plane = Edges.findAt(coordinates=(-C, thick*0.3, 0.0))
    Origin_Sketch_plane = (-C, 0.0, 0.0)
    transformXY = PipePart.MakeSketchTransform(sketchPlane= Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=Origin_Sketch_plane)

    sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=174.79, gridSpacing=4.36, transform=transformXY)
    sketch_1.setPrimaryObject(option=SUPERIMPOSE)

    PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

    sketch_1.ArcByCenterEnds(center=(0.0, thick-a), point1=(0.0, thick-a-KeyHole), point2=(KeyHole, thick-a),direction=COUNTERCLOCKWISE)

    sketch_1.Line(point1=(0.0, thick-a-KeyHole), point2=(0.0, thick))
    sketch_1.Line(point1=(0.0, thick), point2=(KeyHole, thick))
    sketch_1.Line(point1=(KeyHole, thick), point2=(KeyHole, thick-a))

    PipePart.CutExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=sketch_1, depth=C, flipExtrudeDirection=OFF)
    sketch_1.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']

#########################################################################
######### Extrude of crack face curves Spider pattern

#### Wedge
if KeyHole == 0:
    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells

    pickedCells = Cells.findAt(
        ((-Epsilon, Epsilon, KeyHole), )
        )

    Edges = PipePart.edges
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, thick-a-Wedge*cos(pi/4), Wedge*sin(pi/4))),
        Edges.findAt(coordinates=(0.0, thick-a+Wedge*cos(pi/4), Wedge*sin(pi/4))),
        )
    ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


######### A1
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A1*cos(pi/4), A1*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(pi/4), A1*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

######### A2
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A2*cos(pi/4), A2*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(pi/4), A2*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A3
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A3*cos(pi/4), A3*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(pi/4), A3*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A4
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A4*cos(pi/4), A4*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(pi/4), A4*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A5
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A5*cos(pi/4), A5*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(pi/4), A5*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A6
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A6*cos(pi/4), A6*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(pi/4), A6*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


# ######### A7
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A7*cos(pi/4), A7*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(pi/4), A7*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A8
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A8*cos(pi/4), A8*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(pi/4), A8*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A9
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )

Edges = PipePart.edges

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A9*cos(pi/4), A9*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(pi/4), A9*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ######### A10
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a-A10*cos(pi/4), A10*sin(pi/4))),
    Edges.findAt(coordinates=(0.0, thick-a+A10*cos(pi/4), A10*sin(pi/4))),
    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

################ extruding Central Line of crack

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=0,yMin=0.,yMax=thick,zMin=0.0,zMax=thick)
if Weld_Angle > 0 :
    if KeyHole ==0:
        pickedEdges =(
            Edges.findAt(coordinates=(0.0, thick-a, Wedge*0.9)),
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
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width-Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width+Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, thick-Epsilon)),
            )
        ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
        PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)
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
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width-Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width+Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, thick-Epsilon)),
            )
        ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
        PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

else:
    if KeyHole ==0:
        pickedEdges =(
            Edges.findAt(coordinates=(0.0, thick-a, Wedge*0.9)),
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
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width-Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width+Epsilon)),
            )
        ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
        PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)
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
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width-Epsilon)),
            Edges.findAt(coordinates=(0.0, thick-a, Weld_Fixed_Width+Epsilon)),
            )
        ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
        PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)



################# extrude clad line in C
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), )
    )
Edges = PipePart.edges

if Weld_Angle > 0:    
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, clad,Weld_Fixed_Width-Epsilon )),
        Edges.findAt(coordinates=(0.0, clad,Weld_Fixed_Width+Epsilon )),
        Edges.findAt(coordinates=(0.0, clad,thick-Epsilon )),
        )
    ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
else:
    pickedEdges =(
        Edges.findAt(coordinates=(0.0, clad,Weld_Fixed_Width-Epsilon )),
        Edges.findAt(coordinates=(0.0, clad,thick-Epsilon )),
        )
    ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))

PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ################# Weld Fixed Width clad line in C
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells

pickedCells = Cells.findAt(
    ((-Epsilon, Epsilon, KeyHole), ),
    ((-Epsilon, clad+Epsilon, KeyHole), ),
    ((-Epsilon, thick-Epsilon, KeyHole), ),
    )
Edges = PipePart.edges
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick - Epsilon,Weld_Fixed_Width )),
    Edges.findAt(coordinates=(0.0, clad+Epsilon,Weld_Fixed_Width )),
    Edges.findAt(coordinates=(0.0, clad-Epsilon,Weld_Fixed_Width )),

    )
ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ##########################################################################
# ######## Extruding V-groove in C 

if Weld_Angle > 0 :

    PipePart = PipeModel.parts['PIPE']

    Cells = PipePart.cells
    pickedCells = Cells.findAt(
        ((-Epsilon, thick-Epsilon, Weld_Fixed_Width+Epsilon), ), 
        ((-Epsilon, clad+Epsilon, Weld_Fixed_Width+Epsilon), ),
        ((-Epsilon, clad-Epsilon, Weld_Fixed_Width+Epsilon), )
        )
    RR1 = (clad*0.3) /cos(THETA)
    rr1 = RR1 * sin(THETA)
    RR2 = (clad*1.3) /cos(THETA)
    rr2 = RR2 * sin(THETA)
    RR3 = (thick-Epsilon)/cos(THETA)
    rr3 = RR3 * sin(THETA)

    pickedEdges =(
        Edges.findAt(coordinates=(0.0, clad*0.3, Weld_Fixed_Width+rr1)), 
        Edges.findAt(coordinates=(0.0, clad*1.3, Weld_Fixed_Width+rr2)), 
        Edges.findAt(coordinates=(0.0, thick-Epsilon, Weld_Fixed_Width+rr3))
    )
    ExtrudeLine = Edges.findAt(coordinates=(-Epsilon, 0.0, 0.0))
    PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

###########################################################################
########################## Solid Extruding LxPositive
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(thick, clad*0.2, Weld_Fixed_Width*0.3))
Edge_Sketch_plane = Edges.findAt(coordinates=(thick, 0.0, Weld_Fixed_Width*0.2))
Origin_Sketch_plane = (thick, 0.0, 0.0)
transformXY = PipePart.MakeSketchTransform(sketchPlane= Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=Origin_Sketch_plane)
sketch_1 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=206.38, gridSpacing=5.15, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)


sketch_1.Line(point1=(0.0, thick), point2=(thick, thick))
sketch_1.Line(point1=(thick, thick), point2=(thick, 0.0))
sketch_1.Line(point1=(thick, 0.0), point2=(0.0, 0.0))
sketch_1.Line(point1=(0.0, 0.0), point2=(0.0, thick))

PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=sketch_1, depth=LxPositive, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

#######################################################################
############# Extrude Weld_Fixed_Width in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(((thick+Epsilon, Epsilon, 0.0), ))

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick-Epsilon, Weld_Fixed_Width)), 
    Edges.findAt(coordinates=(thick, Epsilon, Weld_Fixed_Width))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick+Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

############ Extrude clad in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells = Cells.findAt(
    ((thick+Epsilon, 0.0, thick-Epsilon), ), 
    ((thick+Epsilon, 0.0, Epsilon), )
    )
if Weld_Angle > 0:
    pickedEdges =(
        Edges.findAt(coordinates=(thick, clad, thick-Epsilon)), 
        Edges.findAt(coordinates=(thick, clad, Epsilon)), 
        Edges.findAt(coordinates=(thick, clad, Weld_Fixed_Width + Epsilon))
        )
else:
    pickedEdges =(
        Edges.findAt(coordinates=(thick, clad, thick-Epsilon)), 
        Edges.findAt(coordinates=(thick, clad, Epsilon)), 
    )

ExtrudeLine = Edges.findAt(coordinates=(thick+Epsilon, 0.0, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

############### Extrude V-groove in LxPositive

if Weld_Angle > 0 :

    PipePart = PipeModel.parts['PIPE']
    Cells = PipePart.cells
    Edges = PipePart.edges

    pickedCells = Cells.findAt(
        ((thick+Epsilon, clad * 1.1, thick-Epsilon), ), 
        ((thick+Epsilon, clad*0.2, thick-Epsilon), )
        )

    RR1 = (clad*0.3) /cos(THETA)
    rr1 = RR1 * sin(THETA)
    RR2 = (clad*1.3) /cos(THETA)
    rr2 = RR2 * sin(THETA) 

    pickedEdges =(
        Edges.findAt(coordinates=(thick, clad*0.3, Weld_Fixed_Width+rr1)), 
        Edges.findAt(coordinates=(thick, clad*1.3, Weld_Fixed_Width+rr2))
        )
    ExtrudeLine = Edges.findAt(coordinates=(thick+Epsilon, 0.0, 0.0))
    PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ##########################################################################################
# ############### Solid Extrude of L

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(thick*0.45,clad*1.3, thick))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, clad*1.3, thick))
Origin_Sketch_plane = (0.0 ,0.0, thick)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin= Origin_Sketch_plane)
sketch_1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)

sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)


sketch_1.Line(point1=(-C, thick), point2=(-C, 0.0))
sketch_1.Line(point1=(-C, 0.0), point2=(thick+LxPositive, 0.0))
sketch_1.Line(point1=(thick+LxPositive, 0.0), point2=(thick+LxPositive, thick))
sketch_1.Line(point1=(thick+LxPositive, thick), point2=(-C, thick))

PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane,sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=sketch_1, depth= L - thick, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

# ##################################################################
# ########### Extruding clad in L

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

PipePart = PipeModel.parts['PIPE']
pickedCells = Cells.findAt(
    ((0, Epsilon, thick+Epsilon), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(LxPositive-Epsilon, clad, thick)), 
    Edges.findAt(coordinates=(thick-Epsilon, clad, thick)), 
    Edges.findAt(coordinates=(Epsilon, clad, thick)), 
    Edges.findAt(coordinates=(-C + Epsilon, clad, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, thick+Epsilon))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ##################################################################
# ########### Extruding clad in X = 0

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.findAt(
    ((0, clad + Epsilon, thick + Epsilon), ), 
    ((0, clad - Epsilon, thick + Epsilon), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, clad- Epsilon, thick)), 
    Edges.findAt(coordinates=(0.0, clad + Epsilon, thick)), 
    Edges.findAt(coordinates=(0.0, thick - Epsilon, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, 0.0, thick+Epsilon))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ########### Extruding clad in X = a2 along L

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.findAt(
    ((thick * 0.2, clad+Epsilon, L * 0.7), ), 
    ((thick * 0.2, clad - Epsilon, L * 0.7), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(a + a2, clad * 1.2, thick)), 
    Edges.findAt(coordinates=(a + a2, clad * 0.2, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, 0.0, L * 0.7))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

########### Extruding clad in X = thick along L

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.findAt(
    ((thick , clad+Epsilon, L * 0.7), ), 
    ((thick , clad - Epsilon, L * 0.7), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(thick, clad * 1.2, thick)), 
    Edges.findAt(coordinates=(thick, clad * 0.2, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, 0.0, L * 0.7))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

################################ Extruding Curve line of crack

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.findAt(
    ((a * 0.1 , thick - Epsilon, L * 0.7), ), 
    ((-C * 0.3, thick - Epsilon, L * 0.7), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(-C * 0.2, thick - a, thick)), 
    Edges.findAt(coordinates=(a * cos(pi/4), thick - a * cos(pi/4), thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, thick, L * 0.2))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#####################################################################
########### Sketching LXPositive/2 partition


PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(thick * 1.2, thick * 0.6, L))
Edge_Sketch_plane = Edges.findAt(coordinates=(thick, clad * 0.25, L))
Origin_Sketch_plane = (thick, 0.0, L)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

sketch_1 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=1289.15, gridSpacing=32.22, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(LxPositive/2, thick), point2=(LxPositive/2, 0.0))

pickedFaces = Faces.findAt(
    ((thick * 1.2, thick * 0.99, L), ), 
    ((thick * 1.2, clad * 0.99,L), )
    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

########### Extruding LxPositive/2 along L 

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

pickedCells = Cells.getByBoundingBox(xMin=0,xMax=thick + LxPositive,yMin=0.,yMax=thick,zMin=0.0,zMax= L )

pickedEdges =(
    Edges.findAt(coordinates=(thick + LxPositive/2, clad*0.6, L)), 
    Edges.findAt(coordinates=(thick + LxPositive/2, thick * 0.99, L))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick, thick, L * 0.7))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

###########################################################################################
################ Sketching LDivs

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(thick + LxPositive,clad*0.6,L*0.7))
Edge_Sketch_plane = Edges.findAt(coordinates=(thick + LxPositive, 0.5, thick))
Origin_Sketch_plane = (thick + LxPositive, 0.0, thick)

transformXY = PipePart.MakeSketchTransform(sketchPlane= Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane , sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
sketch_1 = mdb.models['Pipe_Model_External'].ConstrainedSketch(name='__profile__', sheetSize=1274.31, gridSpacing=31.85, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(-LDiv1, 0.0), point2=(-LDiv1, thick))
sketch_1.Line(point1=(-LDiv2, 0.0), point2=(-LDiv2, thick))
sketch_1.Line(point1=(-LDiv3, 0.0), point2=(-LDiv3, thick))

pickedFaces = Faces.findAt(
    ((thick + LxPositive, clad*0.6, L*0.7), ), 
    ((thick + LxPositive, thick*0.99, L*0.7), )
    )
ExtrudeLine = Edges.findAt(coordinates=(thick + LxPositive, clad * 0.25,thick))
PipePart.PartitionFaceBySketch(sketchUpEdge=ExtrudeLine, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

################# Extruding LDiv1 in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick + LxPositive,yMin=0.,yMax=thick,zMin=0.0,zMax= L )
pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick * 0.9,thick+LDiv1)), 
    Edges.findAt(coordinates=(thick+LxPositive, clad * 0.25, thick +LDiv1))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive * 0.7, thick, thick))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

################# Extruding LDiv2 in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick + LxPositive,yMin=0.,yMax=thick,zMin=0.0,zMax= L )
pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick * 0.9,thick+LDiv2)), 
    Edges.findAt(coordinates=(thick+LxPositive, clad * 0.25, thick +LDiv2))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive * 0.7, thick, thick))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

################# Extruding LDiv3 in LxPositive

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick + LxPositive,yMin=0.,yMax=thick,zMin=0.0,zMax= L )
pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive, thick * 0.9,thick+LDiv3)), 
    Edges.findAt(coordinates=(thick+LxPositive, clad * 0.25, thick +LDiv3))
    )
ExtrudeLine = Edges.findAt(coordinates=(LxPositive * 0.7, thick, thick))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

##########################################################################################################
############ Sketching Lthick under the crack to improve mesh

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(-C * 0.9, clad*1.01, 0.0))
Edge_Sketch_plane = Edges.findAt(coordinates=(-C, clad*1.01, 0.0))
Origin_Sketch_plane = (-C, clad, 0.0)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
sketch_1 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=649.78, gridSpacing=16.24, transform=transformXY)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(0.0, (thick-clad)-a-a2), point2=(-(C+thick+LxPositive),  (thick-clad)-a-a2))
pickedFaces = Faces.findAt(
    ((-C * 0.9, clad*1.01, 0.0),),
    ((a*0.2, clad*1.01, 0.0),),  
    ((thick*0.9, clad*1.01, 0.0), ),
    ((thick*1.1, clad*1.01, 0.0), ),  
    ((thick+LxPositive*.8, clad*1.01, 0.0),)   

    )
PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

####################### Extruding Lthick along L
PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick + LxPositive,yMin=0.,yMax=thick,zMin=0.0,zMax= L )

pickedEdges =(
    Edges.findAt(coordinates=(-C * 0.7, thick-a-a2, 0.0)), 
    Edges.findAt(coordinates=(a * 0.3, thick-a-a2, 0.0)),
    Edges.findAt(coordinates=(thick*0.9, thick-a-a2, 0.0)),
    Edges.findAt(coordinates=(thick * 1.3, thick-a-a2, 0.0)),
    Edges.findAt(coordinates=(thick + LxPositive*0.9, thick-a-a2, 0.0))
    )

ExtrudeLine = Edges.findAt(coordinates=(thick+LxPositive, thick, Weld_Fixed_Width*0.6))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

#####################
### Creating diagonal line in crack 

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces

Face_Sketch_plane = Faces.findAt(coordinates=(a*0.1, thick*0.9, L))
Edge_Sketch_plane = Edges.findAt(coordinates=(a+a2, thick*0.89, L))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,  sketchPlaneSide=SIDE1, origin=(0.0, thick, L))
sketch_1 = PipeModel.ConstrainedSketch(name='__profile__',  sheetSize=1250.38, gridSpacing=31.25, transform=Transform)
sketch_1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch_1, filter=COPLANAR_EDGES)

sketch_1.Line(point1=(0.0, 0.0), point2=(a+a2, -(a+a2)))
pickedFaces = Faces.findAt(
    ((a*0.01, thick-a*0.1, L), ), 
    ((a*1.01, thick*0.9, L), )
    )

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch_1)
sketch_1.unsetPrimaryObject()
del mdb.models['Pipe_Model_External'].sketches['__profile__']

### Extruding Diagnoal line in crack 

pickedCells = Cells.getByBoundingBox(xMin=0,xMax=thick,yMin=0.,yMax=thick,zMin=0.0,zMax= L )

pickedEdges =(
    Edges.findAt(coordinates=(a*0.2*cos(pi/4), thick-a*0.2*cos(pi/4), L)), 
    Edges.findAt(coordinates=(a*1.01*cos(pi/4), thick-a*1.01*cos(pi/4), L)), 

    )
ExtrudeLine = Edges.findAt(coordinates=(0.0, thick, L*0.7))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=REVERSE)

#### Sketcing a1

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
Edges = PipePart.edges
Faces = PipePart.faces


if a/thick > 0.1:
    a1 = a2
else:
    a1 = A10+A1

Face_Sketch_plane = Faces.findAt(coordinates=(a*0.3, thick-Epsilon, L))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, thick-a*0.4, L))
Transform = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, origin=(0.0, thick, L))
sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=177.22, gridSpacing=4.43, transform=Transform)
sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=sketch1, filter=COPLANAR_EDGES)

sketch1.ArcByCenterEnds(center=(0.0, 0.0), point1=(-(a-a1), 0.0), point2=(0.0, (a-a1)),  direction=CLOCKWISE)
sketch1.Line(point1=(0.0, (a-a1)), point2=(C, (a-a1)))

pickedFaces = Faces.findAt(
    ((a*0.3, thick-Epsilon, L), ),
    ((-C*0.3,thick-Epsilon, L), ), 
    ((Epsilon, thick-a*0.4, L), )
    )

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=sketch1)
sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


#### Extruding L10 through all Z direction length to improve mesh

PipePart = PipeModel.parts['PIPE']
Cells = PipePart.cells
pickedCells = Cells.getByBoundingBox(xMin=-C,xMax=thick,yMin=0.,yMax=thick,zMin=0.0,zMax= L )

pickedEdges =(
    Edges.findAt(coordinates=(-C*0.3, thick-a+a1, L)), 
    Edges.findAt(coordinates=((a-a1)*cos(pi/3), thick-(a-a1)*sin(pi/3), L)),
    Edges.findAt(coordinates=((a-a1)*cos(pi/6), thick-(a-a1)*sin(pi/6), L)),

    )
ExtrudeLine = Edges.findAt(coordinates=(-C, thick, L*0.75))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=REVERSE)

######################################################################################################
######################################################################################################
############### MESH & SEED ###########################################################################
########################################################################################################
#######################################################################################################

# ##### partial global meshing
PipePart.seedPart(size=11.0, deviationFactor=0.1, minSizeFactor=0.1)
PipePart.generateMesh()

# # ###### Seeding Cruved lines in spider pattern

PipePart = PipeModel.parts['PIPE']
Edges = PipePart.edges

if KeyHole == 0:
    KK = Wedge
    KK1 = KK * sin(pi/4)
else:
    KK = KeyHole
    KK1 = KK

pickedEdges = Edges.findAt(

     ((a+KK*cos(pi/4), thick, KK*sin(pi/4)), ), 

     ((a+A1*cos(pi/4), thick, A1*sin(pi/4)), ), 
    ((a+A2*cos(pi/4), thick, A2*sin(pi/4)), ), 
    ((a+A3*cos(pi/4), thick, A3*sin(pi/4)), ), 
    ((a+A4*cos(pi/4), thick, A4*sin(pi/4)), ), 
    ((a+A5*cos(pi/4), thick, A5*sin(pi/4)), ), 
    ((a+A6*cos(pi/4), thick, A6*sin(pi/4)), ), 
    ((a+A7*cos(pi/4), thick, A7*sin(pi/4)), ), 
    ((a+A8*cos(pi/4), thick, A8*sin(pi/4)), ), 
    ((a+A9*cos(pi/4), thick, A9*sin(pi/4)), ), 
    ((a+A10*cos(pi/4), thick, A10*sin(pi/4)), ), 

    ((a-KK*cos(pi/4), thick, KK1), ), 
    ((a-A1*cos(pi/4), thick, A1*sin(pi/4)), ), 
    ((a-A2*cos(pi/4), thick, A2*sin(pi/4)), ), 
    ((a-A3*cos(pi/4), thick, A3*sin(pi/4)), ), 
    ((a-A4*cos(pi/4), thick, A4*sin(pi/4)), ), 
    ((a-A5*cos(pi/4), thick, A5*sin(pi/4)), ), 
    ((a-A6*cos(pi/4), thick, A6*sin(pi/4)), ), 
    ((a-A7*cos(pi/4), thick, A7*sin(pi/4)), ), 
    ((a-A8*cos(pi/4), thick, A8*sin(pi/4)), ), 
    ((a-A9*cos(pi/4), thick, A9*sin(pi/4)), ), 
    ((a-A10*cos(pi/4), thick, A10*sin(pi/4)), ), 

#     ##############################################
    ((0.0, thick-(a+KK*cos(pi/4)), KK*sin(pi/4)), ), 
    ((0.0, thick-(a+A1*cos(pi/4)), A1*sin(pi/4)), ), 
    ((0.0, thick-(a+A2*cos(pi/4)), A2*sin(pi/4)), ), 
    ((0.0, thick-(a+A3*cos(pi/4)), A3*sin(pi/4)), ), 
    ((0.0, thick-(a+A4*cos(pi/4)), A4*sin(pi/4)), ), 
    ((0.0, thick-(a+A5*cos(pi/4)), A5*sin(pi/4)), ), 
    ((0.0, thick-(a+A6*cos(pi/4)), A6*sin(pi/4)), ), 
    ((0.0, thick-(a+A7*cos(pi/4)), A7*sin(pi/4)), ), 
    ((0.0, thick-(a+A8*cos(pi/4)), A8*sin(pi/4)), ), 
    ((0.0, thick-(a+A9*cos(pi/4)), A9*sin(pi/4)), ), 
    ((0.0, thick-(a+A10*cos(pi/4)), A10*sin(pi/4)), ), 

    ((0.0, thick-(a-KK*cos(pi/4)), KK1), ), 
    ((0.0, thick-(a-A1*cos(pi/4)), A1*sin(pi/4)), ), 
    ((0.0, thick-(a-A2*cos(pi/4)), A2*sin(pi/4)), ), 
    ((0.0, thick-(a-A3*cos(pi/4)), A3*sin(pi/4)), ), 
    ((0.0, thick-(a-A4*cos(pi/4)), A4*sin(pi/4)), ), 
    ((0.0, thick-(a-A5*cos(pi/4)), A5*sin(pi/4)), ), 
    ((0.0, thick-(a-A6*cos(pi/4)), A6*sin(pi/4)), ), 
    ((0.0, thick-(a-A7*cos(pi/4)), A7*sin(pi/4)), ), 
    ((0.0, thick-(a-A8*cos(pi/4)), A8*sin(pi/4)), ), 
    ((0.0, thick-(a-A9*cos(pi/4)), A9*sin(pi/4)), ), 
    ((0.0, thick-(a-A10*cos(pi/4)), A10*sin(pi/4)), ), 

#     ##############################################
    ((-C, thick-(a+KK*cos(pi/4)), KK*sin(pi/4)), ), 
    ((-C, thick-(a+A1*cos(pi/4)), A1*sin(pi/4)), ), 
    ((-C, thick-(a+A2*cos(pi/4)), A2*sin(pi/4)), ), 
    ((-C, thick-(a+A3*cos(pi/4)), A3*sin(pi/4)), ), 
    ((-C, thick-(a+A4*cos(pi/4)), A4*sin(pi/4)), ), 
    ((-C, thick-(a+A5*cos(pi/4)), A5*sin(pi/4)), ), 
    ((-C, thick-(a+A6*cos(pi/4)), A6*sin(pi/4)), ), 
    ((-C, thick-(a+A7*cos(pi/4)), A7*sin(pi/4)), ), 
    ((-C, thick-(a+A8*cos(pi/4)), A8*sin(pi/4)), ), 
    ((-C, thick-(a+A9*cos(pi/4)), A9*sin(pi/4)), ), 
    ((-C, thick-(a+A10*cos(pi/4)), A10*sin(pi/4)), ), 

    ((-C, thick-a+KK*cos(pi/4), KK1), ), 
    ((-C, thick-(a-A1*cos(pi/4)), A1*sin(pi/4)), ), 
    ((-C, thick-(a-A2*cos(pi/4)), A2*sin(pi/4)), ), 
    ((-C, thick-(a-A3*cos(pi/4)), A3*sin(pi/4)), ), 
    ((-C, thick-(a-A4*cos(pi/4)), A4*sin(pi/4)), ), 
    ((-C, thick-(a-A5*cos(pi/4)), A5*sin(pi/4)), ), 
    ((-C, thick-(a-A6*cos(pi/4)), A6*sin(pi/4)), ), 
    ((-C, thick-(a-A7*cos(pi/4)), A7*sin(pi/4)), ), 
    ((-C, thick-(a-A8*cos(pi/4)), A8*sin(pi/4)), ), 
    ((-C, thick-(a-A9*cos(pi/4)), A9*sin(pi/4)), ), 
    ((-C, thick-(a-A10*cos(pi/4)), A10*sin(pi/4)), ), 

    )
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedSpider, constraint=FINER)

####################################################
#### Curves of crack face in thickness


RR5 = (thick-(a*sin(pi/4))) /cos(THETA)
rr5 = RR5 * sin(THETA)

pickedEdges = Edges.findAt(

    (((a-sqrt(A10**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A10**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A9**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A9**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A8**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A8**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A7**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A7**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A6**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A6**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A5**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A5**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A4**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A4**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A3**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A3**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A2**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A2**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 
    (((a-sqrt(A1**2-KeyHole**2)) * cos(pi/4) , thick-((a-sqrt(A1**2-KeyHole**2))*sin(pi/4)), KeyHole), ), 

    ((a * cos(pi/4) , thick-(a*sin(pi/4)), KeyHole), ), 

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
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), Weld_Fixed_Width), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), Weld_Fixed_Width+rr5), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv1), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv2), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), thick+LDiv3), ), 
    ((a * cos(pi/4) , thick-(a*sin(pi/4)), L), ), 

    (((a+A1) * cos(pi/4) , thick-((a+A1)*sin(pi/4)), 0.0), ), 
    (((a+A2) * cos(pi/4) , thick-((a+A2)*sin(pi/4)), 0.0), ), 
    (((a+A3) * cos(pi/4) , thick-((a+A3)*sin(pi/4)), 0.0), ), 
    (((a+A4) * cos(pi/4) , thick-((a+A4)*sin(pi/4)), 0.0), ), 
    (((a+A5) * cos(pi/4) , thick-((a+A5)*sin(pi/4)), 0.0), ), 
    (((a+A6) * cos(pi/4) , thick-((a+A6)*sin(pi/4)), 0.0), ), 
    (((a+A7) * cos(pi/4) , thick-((a+A7)*sin(pi/4)), 0.0), ), 
    (((a+A8) * cos(pi/4) , thick-((a+A8)*sin(pi/4)), 0.0), ), 
    (((a+A9) * cos(pi/4) , thick-((a+A9)*sin(pi/4)), 0.0), ), 
    (((a+A10) * cos(pi/4) , thick-((a+A10)*sin(pi/4)), 0.0), ), 

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseedcurve, constraint=FINER)

# ##########################################################
# ########### Seeding C

RR6 = (clad) /cos(THETA)
rr6 = RR6 * sin(THETA)
RR7 = (thick-a-Weld_Fixed_Width) /cos(THETA)
rr7 = RR7 * sin(THETA)
RR8 = (thick-a) /cos(THETA)
rr8 = RR8 * sin(THETA)
RR9 = (thick) /cos(THETA)
rr9 = RR9 * sin(THETA)
RR10 = (thick-a*0.2) /cos(THETA)
rr10 = RR10 * sin(THETA)
pickedEdges = Edges.findAt(

    (( -C * 0.2 , 0.0, 0.0), ), 
    (( -C * 0.2 , clad, 0.0), ), 
    (( -C * 0.2 , thick - a- Weld_Fixed_Width, 0.0), ), 
    (( -C * 0.2 , thick - a- A10, 0.0), ), 
    (( -C * 0.2 , thick -a- A9, 0.0), ), 
    (( -C * 0.2 , thick -a- A8, 0.0), ), 
    (( -C * 0.2 , thick -a- A7, 0.0), ), 
    (( -C * 0.2 , thick -a- A6, 0.0), ), 
    (( -C * 0.2 , thick -a- A5, 0.0), ),
    (( -C * 0.2 , thick -a- A4, 0.0), ), 
    (( -C * 0.2 , thick -a- A3, 0.0), ), 
    (( -C * 0.2 , thick -a- A2, 0.0), ), 
    (( -C * 0.2 , thick -a- A1, 0.0), ), 
    (( -C * 0.2 , thick -a- KeyHole, 0.0), ), 

    (( -C * 0.2 , thick -a, KeyHole), ), 
    (( -C * 0.2 , thick -a, A1), ), 
    (( -C * 0.2 , thick -a, A2), ), 
    (( -C * 0.2 , thick -a, A3), ), 
    (( -C * 0.2 , thick -a, A4), ), 
    (( -C * 0.2 , thick -a, A5), ), 
    (( -C * 0.2 , thick -a, A6), ), 
    (( -C * 0.2 , thick -a, A7), ), 
    (( -C * 0.2 , thick -a, A8), ), 
    (( -C * 0.2 , thick -a, A9), ), 
    (( -C * 0.2 , thick -a, A10), ), 
    (( -C * 0.2 , thick -a, Weld_Fixed_Width), ), 

    (( -C * 0.2 , thick - a+ sqrt(A10**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A9**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A8**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A7**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A6**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A5**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A4**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A3**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A2**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A1**2-KeyHole**2), KeyHole), ), 
    (( -C * 0.2 , thick - a+ sqrt(A1**2-KeyHole**2), KeyHole), ), 
     (( -C * 0.2 , thick , KeyHole), ), 

    (( -C * 0.2 , thick , Weld_Fixed_Width), ), 
    (( -C * 0.2 , thick-a , Weld_Fixed_Width), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , Weld_Fixed_Width), ), 
    (( -C * 0.2 , clad , Weld_Fixed_Width), ), 
    (( -C * 0.2 , 0.0 , Weld_Fixed_Width), ), 

    (( -C * 0.2 , clad , Weld_Fixed_Width+rr6), ), 
    (( -C * 0.2 , thick -a - Weld_Fixed_Width , Weld_Fixed_Width+rr7), ), 
    (( -C * 0.2 , thick -a  , Weld_Fixed_Width+rr8), ), 
    (( -C * 0.2 , thick , Weld_Fixed_Width+rr9), ), 

    (( -C * 0.2 , thick , thick), ), 
    (( -C * 0.2 , thick-a , thick), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , thick), ), 
    (( -C * 0.2 , clad , thick), ), 
    (( -C * 0.2 , 0 , thick), ), 

    (( -C * 0.2 , thick , thick+LDiv1), ), 
    (( -C * 0.2 , thick-a , thick+LDiv1), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , thick+LDiv1), ), 
    (( -C * 0.2 , clad , thick+LDiv1), ), 
    (( -C * 0.2 , 0 , thick+LDiv1), ),

    (( -C * 0.2 , thick , thick+LDiv2), ), 
    (( -C * 0.2 , thick-a , thick+LDiv2), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , thick+LDiv2), ), 
    (( -C * 0.2 , clad , thick+LDiv2), ), 
    (( -C * 0.2 , 0 , thick+LDiv2), ), 
 

    (( -C * 0.2 , thick , thick+LDiv2), ), 
    (( -C * 0.2 , thick-a , thick+LDiv2), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , thick+LDiv2), ), 
    (( -C * 0.2 , clad , thick+LDiv2), ), 
    (( -C * 0.2 , 0 , thick+LDiv2), ), 

    (( -C * 0.2 , thick , thick+LDiv3), ), 
    (( -C * 0.2 , thick-a , thick+LDiv3), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , thick+LDiv3), ), 
    (( -C * 0.2 , clad , thick+LDiv3), ), 
    (( -C * 0.2 , 0 , thick+LDiv3), ), 

    (( -C * 0.2 , thick , L), ), 
    (( -C * 0.2 , thick-a , L), ), 
    (( -C * 0.2 , thick-a-Weld_Fixed_Width , L), ), 
    (( -C * 0.2 , clad , L), ), 
    (( -C * 0.2 , 0 , L), ), 


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedC, constraint=FINER)

pickedEdges = Edges.findAt(

    (( a * 0.2 , thick , KeyHole), ), 
    (( a * 0.2 , thick , Weld_Fixed_Width), ), 
    (( a * 0.2 , thick , Weld_Fixed_Width+rr9), ), 
    (( a * 0.2 , thick , thick), ), 
    (( a * 0.2 , thick , thick+LDiv1), ), 
    (( a * 0.2 , thick , thick+LDiv2), ), 
    (( a * 0.2 , thick , thick+LDiv3), ), 
    (( a * 0.2 , thick , L), ), 

    (( a * 1.2 , thick , Weld_Fixed_Width), ), 
    (( a * 1.2 , thick , Weld_Fixed_Width+rr9), ), 
    (( a * 1.2 , thick , thick), ), 
    (( a * 1.2 , thick , thick+LDiv1), ), 
    (( a * 1.2 , thick , thick+LDiv2), ), 
    (( a * 1.2 , thick , thick+LDiv3), ), 
    (( a * 1.2 , thick , L), ), 


    (( 0.0 , thick-a*0.2 , KeyHole), ), 
    (( 0.0 , thick-a*0.2 , Weld_Fixed_Width), ), 
    (( 0.0 , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 
    (( 0.0 , thick-a*0.2 , thick), ), 
    (( 0.0 , thick-a*0.2 , thick+LDiv1), ), 
    (( 0.0 , thick-a*0.2 , thick+LDiv2), ), 
    (( 0.0 , thick-a*0.2 , thick+LDiv3), ), 
    (( 0.0 , thick-a*0.2 , L), ), 

    (( -C , thick-a*0.2 , KeyHole), ), 
    (( -C , thick-a*0.2 , Weld_Fixed_Width), ), 
    (( -C , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 
    (( -C , thick-a*0.2 , thick), ), 
    (( -C , thick-a*0.2 , thick+LDiv1), ), 
    (( -C , thick-a*0.2 , thick+LDiv2), ), 
    (( -C , thick-a*0.2 , thick+LDiv3), ), 
    (( -C , thick-a*0.2 , L), ), 

    ((a+ a2 , thick-a*0.2, 0.0), ),  
    ((thick , thick-a*0.2, 0.0), ),  
    ((thick+LxPositive/2 , thick-a*0.2, 0.0), ),  
    ((thick+LxPositive , thick-a*0.2, 0.0), ),  

    ((a+ a2 , thick-a*0.2, Weld_Fixed_Width), ),  
    ((thick , thick-a*0.2, Weld_Fixed_Width), ),  
    ((thick+LxPositive/2 , thick-a*0.2, Weld_Fixed_Width), ),  
    ((thick+LxPositive , thick-a*0.2, Weld_Fixed_Width), ),

    ((a+a2 , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 
    ((thick , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 
    ((thick+LxPositive/2 , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 
    ((thick+LxPositive , thick-a*0.2 , Weld_Fixed_Width+rr10), ), 

    ((a+ a2 , thick-a*0.2, thick), ),  
    ((thick , thick-a*0.2, thick), ),  
    ((thick+LxPositive/2 , thick-a*0.2, thick), ),  
    ((thick+LxPositive , thick-a*0.2, thick), ), 

    ((a+ a2 , thick-a*0.2, thick+LDiv1), ),  
    ((thick , thick-a*0.2, thick+LDiv1), ),  
    ((thick+LxPositive/2 , thick-a*0.2, thick+LDiv1), ),  
    ((thick+LxPositive , thick-a*0.2, thick+LDiv1), ),

    ((a+ a2 , thick-a*0.2, thick+LDiv2), ),  
    ((thick , thick-a*0.2, thick+LDiv2), ),  
    ((thick+LxPositive/2 , thick-a*0.2, thick+LDiv2), ),  
    ((thick+LxPositive , thick-a*0.2, thick+LDiv2), ),

    ((a+ a2 , thick-a*0.2, thick+LDiv3), ),  
    ((thick , thick-a*0.2, thick+LDiv3), ),  
    ((thick+LxPositive/2 , thick-a*0.2, thick+LDiv3), ),  
    ((thick+LxPositive , thick-a*0.2, thick+LDiv3), ),   

    ((a+ a2 , thick-a*0.2, L), ),  
    ((thick , thick-a*0.2, L), ),  
    ((thick+LxPositive/2 , thick-a*0.2, L), ),  
    ((thick+LxPositive , thick-a*0.2, L), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseeda, constraint=FINER)

pickedEdges = Edges.findAt(

    (( -C , 0, Weld_Fixed_Width *0.3), ), 
    (( -C , clad, Weld_Fixed_Width *0.3), ), 
    (( -C , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    (( -C , thick, Weld_Fixed_Width *0.3), ), 

    (( 0.0 , 0, Weld_Fixed_Width *0.3), ), 
    (( 0.0 , clad, Weld_Fixed_Width *0.3), ), 
    (( 0.0 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    (( 0.0 , thick, Weld_Fixed_Width *0.3), ), 

    ((a+Weld_Fixed_Width , 0, Weld_Fixed_Width *0.3), ), 
    ((a+Weld_Fixed_Width , clad, Weld_Fixed_Width *0.3), ), 
    ((a+Weld_Fixed_Width , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    ((a+Weld_Fixed_Width , thick, Weld_Fixed_Width *0.3), ), 

    ((thick , 0, Weld_Fixed_Width *0.3), ), 
    ((thick , clad, Weld_Fixed_Width *0.3), ), 
    ((thick , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    ((thick , thick, Weld_Fixed_Width *0.3), ),

    ((thick+LxPositive/2 , 0, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive/2 , clad, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive/2 , thick, Weld_Fixed_Width *0.3), ),  

    ((thick+LxPositive , 0, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive , clad, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive , thick-a-Weld_Fixed_Width, Weld_Fixed_Width *0.3), ), 
    ((thick+LxPositive , thick, Weld_Fixed_Width *0.3), ),  

    ((a*0.2 , thick-a-Weld_Fixed_Width, 0.0), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, thick), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, thick+LDiv1), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, thick+LDiv2), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, thick+LDiv3), ),  
    ((a*0.2 , thick-a-Weld_Fixed_Width, L), ),  

    ((a*0.2 , clad, 0.0), ),  
    ((a*0.2 , clad, Weld_Fixed_Width), ),  
    ((a*0.2 , clad, thick), ),  
    ((a*0.2 , clad, thick+LDiv1), ),  
    ((a*0.2 , clad, thick+LDiv2), ),  
    ((a*0.2 , clad, thick+LDiv3), ),  
    ((a*0.2 , clad, L), ),  

    ((a*0.2 , 0.0, 0.0), ),  
    ((a*0.2 , 0.0, Weld_Fixed_Width), ),  
    ((a*0.2 , 0.0, thick), ),  
    ((a*0.2 , 0.0, thick+LDiv1), ),  
    ((a*0.2 , 0.0, thick+LDiv2), ),  
    ((a*0.2 , 0.0, thick+LDiv3), ),  
    ((a*0.2 , 0.0, L), ),  

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedWeldFixed, constraint=FINER)

if a/thick >=0.5:
    Nseedthick = 3

RR11 = (thick-a-Weld_Fixed_Width - Epsilon) /cos(THETA)
rr11 = RR11 * sin(THETA)
pickedEdges = Edges.findAt(

    ((-C , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((-C , thick-a-Weld_Fixed_Width - Epsilon, L), ),  

    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((0.0 , thick-a-Weld_Fixed_Width - Epsilon, L), ),

#     ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((a+a2 , thick-a-Weld_Fixed_Width - Epsilon, L), ),  

    ((thick , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((thick , thick-a-Weld_Fixed_Width - Epsilon, L), ),  

    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width - Epsilon, L), ),

    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, 0.0), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, Weld_Fixed_Width+rr11), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, thick), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv1), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv2), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, thick+LDiv3), ),  
    ((thick+LxPositive , thick-a-Weld_Fixed_Width - Epsilon, L), ),
)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseedthick, constraint=FINER)

RR12 = (clad*0.2) /cos(THETA)
rr12 = RR12 * sin(THETA)

pickedEdges = Edges.findAt(

    ((-C , clad*0.2, 0.0), ),
    ((-C , clad*0.2, Weld_Fixed_Width), ),
    ((-C , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((-C , clad*0.2, thick), ),
    ((-C , clad*0.2, thick+LDiv1), ),
    ((-C , clad*0.2, thick+LDiv2), ),
    ((-C , clad*0.2, thick+LDiv3), ),
    ((-C , clad*0.2, L), ),

    ((0.0 , clad*0.2, 0.0), ),
    ((0.0 , clad*0.2, Weld_Fixed_Width), ),
    ((0.0 , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((0.0 , clad*0.2, thick), ),
    ((0.0 , clad*0.2, thick+LDiv1), ),
    ((0.0 , clad*0.2, thick+LDiv2), ),
    ((0.0 , clad*0.2, thick+LDiv3), ),
    ((0.0 , clad*0.2, L), ),

    ((a+a2 , clad*0.2, 0.0), ),
    ((a+a2 , clad*0.2, Weld_Fixed_Width), ),
    ((a+a2 , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((a+a2 , clad*0.2, thick), ),
    ((a+a2 , clad*0.2, thick+LDiv1), ),
    ((a+a2 , clad*0.2, thick+LDiv2), ),
    ((a+a2 , clad*0.2, thick+LDiv3), ),
    ((a+a2 , clad*0.2, L), ),

    ((thick , clad*0.2, 0.0), ),
    ((thick , clad*0.2, Weld_Fixed_Width), ),
    ((thick , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((thick , clad*0.2, thick), ),
    ((thick , clad*0.2, thick+LDiv1), ),
    ((thick , clad*0.2, thick+LDiv2), ),
    ((thick , clad*0.2, thick+LDiv3), ),
    ((thick , clad*0.2, L), ),

    ((thick+LxPositive/2 , clad*0.2, 0.0), ),
    ((thick+LxPositive/2 , clad*0.2, Weld_Fixed_Width), ),
    ((thick+LxPositive/2 , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((thick+LxPositive/2 , clad*0.2, thick), ),
    ((thick+LxPositive/2 , clad*0.2, thick+LDiv1), ),
    ((thick+LxPositive/2 , clad*0.2, thick+LDiv2), ),
    ((thick+LxPositive/2 , clad*0.2, thick+LDiv3), ),
    ((thick+LxPositive/2 , clad*0.2, L), ),

    ((thick+LxPositive , clad*0.2, 0.0), ),
    ((thick+LxPositive , clad*0.2, Weld_Fixed_Width), ),
    ((thick+LxPositive , clad*0.2, Weld_Fixed_Width+rr12), ),
    ((thick+LxPositive , clad*0.2, thick), ),
    ((thick+LxPositive , clad*0.2, thick+LDiv1), ),
    ((thick+LxPositive , clad*0.2, thick+LDiv2), ),
    ((thick+LxPositive , clad*0.2, thick+LDiv3), ),
    ((thick+LxPositive , clad*0.2, L), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=Nseeedclad, constraint=FINER)

if Weld_Angle ==0:
    pickedEdges = Edges.findAt(
        ((-C , thick, thick-Epsilon), ),
        ((-C , thick-a, thick-Epsilon), ),
        ((-C , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((-C , clad, thick-Epsilon), ),
        ((-C , 0.0, thick-Epsilon), ),

        ((0.0 , thick, thick-Epsilon), ),
        ((0.0 , thick-a, thick-Epsilon), ),
        ((0.0 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((0.0 , clad, thick-Epsilon), ),
        ((0.0 , 0.0, thick-Epsilon), ),

        ((a+a2 , thick, thick-Epsilon), ),
        ((a+a2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((a+a2 , clad, thick-Epsilon), ),
        ((a+a2 , 0.0, thick-Epsilon), ),

        ((thick , thick, thick-Epsilon), ),
        ((thick , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick , clad, thick-Epsilon), ),
        ((thick , 0.0, thick-Epsilon), ),

        ((thick+LxPositive/2 , thick, thick-Epsilon), ),
        ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive/2 , clad, thick-Epsilon), ),
        ((thick+LxPositive/2 , 0.0, thick-Epsilon), ),

        ((thick+LxPositive , thick, thick-Epsilon), ),
        ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive , clad, thick-Epsilon), ),
        ((thick+LxPositive , 0.0, thick-Epsilon), ),

    )
    PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedthickZ, constraint=FINER)

elif Weld_Angle == 20:
    NseedthickZ1 = NseedthickZ - 1
    NseedthickZ2 = 3
    pickedEdges = Edges.findAt(
        ((-C , thick, thick-Epsilon), ),
        ((-C , thick-a, thick-Epsilon), ),
        ((-C , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((-C , clad, thick-Epsilon), ),
        ((-C , 0.0, thick-Epsilon), ),

        ((0.0 , thick, thick-Epsilon), ),
        ((0.0 , thick-a, thick-Epsilon), ),
        ((0.0 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((0.0 , clad, thick-Epsilon), ),
        ((0.0 , 0.0, thick-Epsilon), ),

        ((a+a2 , thick, thick-Epsilon), ),
        ((a+a2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((a+a2 , clad, thick-Epsilon), ),
        ((a+a2 , 0.0, thick-Epsilon), ),

        ((thick , thick, thick-Epsilon), ),
        ((thick , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick , clad, thick-Epsilon), ),
        ((thick , 0.0, thick-Epsilon), ),

        ((thick+LxPositive/2 , thick, thick-Epsilon), ),
        ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive/2 , clad, thick-Epsilon), ),
        ((thick+LxPositive/2 , 0.0, thick-Epsilon), ),

        ((thick+LxPositive , thick, thick-Epsilon), ),
        ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive , clad, thick-Epsilon), ),
        ((thick+LxPositive , 0.0, thick-Epsilon), ),

    )
    PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedthickZ1, constraint=FINER)
else:
    NseedthickZ1 = NseedthickZ - 2
    NseedthickZ2 = 4
    pickedEdges = Edges.findAt(
        ((-C , thick, thick-Epsilon), ),
        ((-C , thick-a, thick-Epsilon), ),
        ((-C , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((-C , clad, thick-Epsilon), ),
        ((-C , 0.0, thick-Epsilon), ),

        ((0.0 , thick, thick-Epsilon), ),
        ((0.0 , thick-a, thick-Epsilon), ),
        ((0.0 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((0.0 , clad, thick-Epsilon), ),
        ((0.0 , 0.0, thick-Epsilon), ),

        ((a , thick, thick-Epsilon), ),

        ((a+a2 , thick, thick-Epsilon), ),
        ((a+a2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((a+a2 , clad, thick-Epsilon), ),
        ((a+a2 , 0.0, thick-Epsilon), ),

        ((thick , thick, thick-Epsilon), ),
        ((thick , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick , clad, thick-Epsilon), ),
        ((thick , 0.0, thick-Epsilon), ),

        ((thick+LxPositive/2 , thick, thick-Epsilon), ),
        ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive/2 , clad, thick-Epsilon), ),
        ((thick+LxPositive/2 , 0.0, thick-Epsilon), ),

        ((thick+LxPositive , thick, thick-Epsilon), ),
        ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick-Epsilon), ),
        ((thick+LxPositive , clad, thick-Epsilon), ),
        ((thick+LxPositive , 0.0, thick-Epsilon), ),

    )
    PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedthickZ1, constraint=FINER)

    pickedEdges = Edges.findAt(
        ((-C , thick, Weld_Fixed_Width+Epsilon), ),
        ((-C , thick-a, Weld_Fixed_Width+Epsilon), ),
        ((-C , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((-C , clad, Weld_Fixed_Width+Epsilon), ),
        ((-C , 0.0, Weld_Fixed_Width+Epsilon), ),

        ((0.0 , thick, Weld_Fixed_Width+Epsilon), ),
        ((0.0 , thick-a, Weld_Fixed_Width+Epsilon), ),
        ((0.0 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((0.0 , clad, Weld_Fixed_Width+Epsilon), ),
        ((0.0 , 0.0, Weld_Fixed_Width+Epsilon), ),

        ((a+a2 , thick, Weld_Fixed_Width+Epsilon), ),
        ((a+a2 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((a+a2 , clad, Weld_Fixed_Width+Epsilon), ),
        ((a+a2 , 0.0, Weld_Fixed_Width+Epsilon), ),

        ((thick , thick, Weld_Fixed_Width+Epsilon), ),
        ((thick , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((thick , clad, Weld_Fixed_Width+Epsilon), ),
        ((thick , 0.0, Weld_Fixed_Width+Epsilon), ),

        ((thick+LxPositive/2 , thick, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive/2 , clad, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive/2 , 0.0, Weld_Fixed_Width+Epsilon), ),

        ((thick+LxPositive , thick, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive , clad, Weld_Fixed_Width+Epsilon), ),
        ((thick+LxPositive , 0.0, Weld_Fixed_Width+Epsilon), ),

    )
    PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedthickZ2, constraint=FINER)



pickedEdges = Edges.findAt(

    ((-C , thick, thick + Epsilon), ),
    ((-C , thick-a, thick + Epsilon), ),
    ((-C , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((-C , clad, thick + Epsilon), ),
    ((-C , 0.0, thick + Epsilon), ),

    ((0.0 , thick, thick + Epsilon), ),
    ((0.0 , thick-a, thick + Epsilon), ),
    ((0.0 , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((0.0 , clad, thick + Epsilon), ),
    ((0.0 , 0.0, thick + Epsilon), ),

    ((a , thick, thick + Epsilon), ),


    ((a+a2 , thick, thick + Epsilon), ),
    ((a+a2 , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((a+a2 , clad, thick + Epsilon), ),
    ((a+a2 , 0.0, thick + Epsilon), ),

    ((thick , thick, thick + Epsilon), ),
    ((thick , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((thick , clad, thick + Epsilon), ),
    ((thick , 0.0, thick + Epsilon), ),

    ((thick+LxPositive/2 , thick, thick + Epsilon), ),
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((thick+LxPositive/2 , clad, thick + Epsilon), ),
    ((thick+LxPositive/2 , 0.0, thick + Epsilon), ),

    ((thick+LxPositive , thick, thick + Epsilon), ),
    ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick + Epsilon), ),
    ((thick+LxPositive , clad, thick + Epsilon), ),
    ((thick+LxPositive , 0.0, thick + Epsilon), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv1, constraint=FINER)

pickedEdges = Edges.findAt(

    ((-C , thick, thick + LDiv1+ Epsilon), ),
    ((-C , thick-a, thick + LDiv1+ Epsilon), ),
    ((-C , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((-C , clad, thick + LDiv1+ Epsilon), ),
    ((-C , 0.0, thick + LDiv1+ Epsilon), ),

    ((0.0 , thick, thick + LDiv1+ Epsilon), ),
    ((0.0 , thick-a, thick + LDiv1+ Epsilon), ),
    ((0.0 , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((0.0 , clad, thick + LDiv1+ Epsilon), ),
    ((0.0 , 0.0, thick + LDiv1+ Epsilon), ),

    ((a , thick, thick + LDiv1+ Epsilon), ),


    ((a+a2 , thick, thick + LDiv1+ Epsilon), ),
    ((a+a2 , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((a+a2 , clad, thick + LDiv1+ Epsilon), ),
    ((a+a2 , 0.0, thick + LDiv1+ Epsilon), ),

    ((thick , thick, thick + LDiv1+ Epsilon), ),
    ((thick , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((thick , clad, thick + LDiv1+ Epsilon), ),
    ((thick , 0.0, thick + LDiv1+ Epsilon), ),

    ((thick+LxPositive/2 , thick, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive/2 , clad, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive/2 , 0.0, thick + LDiv1+ Epsilon), ),

    ((thick+LxPositive , thick, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive , clad, thick + LDiv1+ Epsilon), ),
    ((thick+LxPositive , 0.0, thick + LDiv1+ Epsilon), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv2, constraint=FINER)

pickedEdges = Edges.findAt(

    ((-C , thick, thick + LDiv2+ Epsilon), ),
    ((-C , thick-a, thick + LDiv2+ Epsilon), ),
    ((-C , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((-C , clad, thick + LDiv2+ Epsilon), ),
    ((-C , 0.0, thick + LDiv2+ Epsilon), ),

    ((0.0 , thick, thick + LDiv2+ Epsilon), ),
    ((0.0 , thick-a, thick + LDiv2+ Epsilon), ),
    ((0.0 , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((0.0 , clad, thick + LDiv2+ Epsilon), ),
    ((0.0 , 0.0, thick + LDiv2+ Epsilon), ),

    ((a , thick, thick + LDiv2+ Epsilon), ),


    ((a+a2 , thick, thick + LDiv2+ Epsilon), ),
    ((a+a2 , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((a+a2 , clad, thick + LDiv2+ Epsilon), ),
    ((a+a2 , 0.0, thick + LDiv2+ Epsilon), ),

    ((thick , thick, thick + LDiv2+ Epsilon), ),
    ((thick , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((thick , clad, thick + LDiv2+ Epsilon), ),
    ((thick , 0.0, thick + LDiv2+ Epsilon), ),

    ((thick+LxPositive/2 , thick, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive/2 , clad, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive/2 , 0.0, thick + LDiv2+ Epsilon), ),

    ((thick+LxPositive , thick, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive , clad, thick + LDiv2+ Epsilon), ),
    ((thick+LxPositive , 0.0, thick + LDiv2+ Epsilon), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedLDiv3, constraint=FINER)

pickedEdges = Edges.findAt(

    ((-C , thick, thick + LDiv3+ Epsilon), ),
    ((-C , thick-a, thick + LDiv3+ Epsilon), ),
    ((-C , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((-C , clad, thick + LDiv3+ Epsilon), ),
    ((-C , 0.0, thick + LDiv3+ Epsilon), ),

    ((0.0 , thick, thick + LDiv3+ Epsilon), ),
    ((0.0 , thick-a, thick + LDiv3+ Epsilon), ),
    ((0.0 , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((0.0 , clad, thick + LDiv3+ Epsilon), ),
    ((0.0 , 0.0, thick + LDiv3+ Epsilon), ),

    ((a , thick, thick + LDiv3+ Epsilon), ),

    ((a+a2 , thick, thick + LDiv3+ Epsilon), ),
    ((a+a2 , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((a+a2 , clad, thick + LDiv3+ Epsilon), ),
    ((a+a2 , 0.0, thick + LDiv3+ Epsilon), ),

    ((thick , thick, thick + LDiv3+ Epsilon), ),
    ((thick , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((thick , clad, thick + LDiv3+ Epsilon), ),
    ((thick , 0.0, thick + LDiv3+ Epsilon), ),

    ((thick+LxPositive/2 , thick, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive/2 , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive/2 , clad, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive/2 , 0.0, thick + LDiv3+ Epsilon), ),

    ((thick+LxPositive , thick, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive , thick-a-Weld_Fixed_Width, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive , clad, thick + LDiv3+ Epsilon), ),
    ((thick+LxPositive , 0.0, thick + LDiv3+ Epsilon), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseeedL, constraint=FINER)

RR13 = (thick-a-Weld_Fixed_Width) /cos(THETA)
rr13 = RR13 * sin(THETA)

pickedEdges = Edges.findAt(

    ((thick+Epsilon , thick, 0.0), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, 0.0), ),
    ((thick+Epsilon , clad, 0.0), ),
    ((thick+Epsilon , 0.0, 0.0), ),

    ((thick+LxPositive-Epsilon , thick, 0.0), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, 0.0), ),
    ((thick+LxPositive-Epsilon , clad, 0.0), ),
    ((thick+LxPositive-Epsilon , 0.0, 0.0), ),

    #######

    ((thick+Epsilon , thick, Weld_Fixed_Width), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, Weld_Fixed_Width), ),
    ((thick+Epsilon , clad, Weld_Fixed_Width), ),
    ((thick+Epsilon , 0.0, Weld_Fixed_Width), ),

    ((thick+LxPositive-Epsilon , thick, Weld_Fixed_Width), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, Weld_Fixed_Width), ),
    ((thick+LxPositive-Epsilon , clad, Weld_Fixed_Width), ),
    ((thick+LxPositive-Epsilon , 0.0, Weld_Fixed_Width), ),

    ##########

    ((thick+Epsilon , thick, Weld_Fixed_Width+rr9), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+rr13), ),
    ((thick+Epsilon , clad, Weld_Fixed_Width+rr6), ),
    ((thick+Epsilon , 0.0, Weld_Fixed_Width), ),

    ((thick+LxPositive-Epsilon , thick, Weld_Fixed_Width+rr9), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, Weld_Fixed_Width+rr13), ),
    ((thick+LxPositive-Epsilon , clad, Weld_Fixed_Width+rr6), ),
    ((thick+LxPositive-Epsilon , 0.0, Weld_Fixed_Width), ),

    ##########

    ((thick+Epsilon , thick, thick), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, thick), ),
    ((thick+Epsilon , clad, thick), ),
    ((thick+Epsilon , 0.0, thick), ),

    ((thick+LxPositive-Epsilon , thick, thick), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, thick), ),
    ((thick+LxPositive-Epsilon , clad, thick), ),
    ((thick+LxPositive-Epsilon , 0.0, thick), ),

    ##########

    ((thick+Epsilon , thick, thick+LDiv1), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv1), ),
    ((thick+Epsilon , clad, thick+LDiv1), ),
    ((thick+Epsilon , 0.0, thick+LDiv1), ),

    ((thick+LxPositive-Epsilon , thick, thick+LDiv1), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv1), ),
    ((thick+LxPositive-Epsilon , clad, thick+LDiv1), ),
    ((thick+LxPositive-Epsilon , 0.0, thick+LDiv1), ),

#     ##########

    ((thick+Epsilon , thick, thick+LDiv2), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv2), ),
    ((thick+Epsilon , clad, thick+LDiv2), ),
    ((thick+Epsilon , 0.0, thick+LDiv2), ),

    ((thick+LxPositive-Epsilon , thick, thick+LDiv2), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv2), ),
    ((thick+LxPositive-Epsilon , clad, thick+LDiv2), ),
    ((thick+LxPositive-Epsilon , 0.0, thick+LDiv2), ),

    ##########

    ((thick+Epsilon , thick, thick+LDiv3), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv3), ),
    ((thick+Epsilon , clad, thick+LDiv3), ),
    ((thick+Epsilon , 0.0, thick+LDiv3), ),

    ((thick+LxPositive-Epsilon , thick, thick+LDiv3), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, thick+LDiv3), ),
    ((thick+LxPositive-Epsilon , clad, thick+LDiv3), ),
    ((thick+LxPositive-Epsilon , 0.0, thick+LDiv3), ),

    ##########

    ((thick+Epsilon , thick, L), ),
    ((thick+Epsilon , thick-a-Weld_Fixed_Width, L), ),
    ((thick+Epsilon , clad, L), ),
    ((thick+Epsilon , 0.0, L), ),

    ((thick+LxPositive-Epsilon , thick, L), ),
    ((thick+LxPositive-Epsilon , thick-a-Weld_Fixed_Width, L), ),
    ((thick+LxPositive-Epsilon , clad, L), ),
    ((thick+LxPositive-Epsilon , 0.0, L), ),
)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=NseedLxPos, constraint=FINER)

elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,  secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6 , elemLibrary=STANDARD)

cells = Cells.getByBoundingBox(xMin=-20000, xMax=20000,yMin=-20000,yMax=20000,zMin=-20000,zMax=20000)
pickedRegions =(cells,)
PipePart.setElementType(regions=pickedRegions, elemTypes=(elemType1,elemType2,))

if KeyHole == 0:
    PipePart.setSweepPath(region=Cells.findAt(coordinates=(-C, thick-a+Wedge/2, 0.0)), 
        edge=Edges.findAt(coordinates=(-C+0.1, thick-a, 0.00)), sense=REVERSE)
        
PipePart.generateMesh()


# # #####################################################################################
# # ############## MATERIALs  ##########################################################
# # #######################################################################################


import material
if MaterialType==None:
    MaterialType = 'Plastic Deformation'

if MaterialType != 'Elastic-Plastic':
    if MaterialType != 'Plastic Deformation':
        MaterialType = 'Plastic Deformation'

RR14 = (thick) /cos(THETA)
rr14 = RR14 * sin(THETA)

if MaterialType == 'Elastic-Plastic':
    PipeModel.Material('BaseSteel')
    PipeModel.materials['BaseSteel'].Density(table=((BaseSteelDens, ), ))
    PipeModel.materials['BaseSteel'].Elastic(table=((BaseSteelE,BaseSteelPratio), ))
    PipeModel.materials['BaseSteel'].Plastic(table=((Yield_Stress_Base, 0.0),  (Yield_Stress_Base, Yield_Strain_Base)))
    PipeModel.HomogeneousSolidSection(name='BaseSteel_Section', material='BaseSteel', thickness=None)
    ##
    if Weld_Angle == 0:
        CELLS = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=Weld_Fixed_Width,zMax=thick+L)
    else:
        cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=Weld_Fixed_Width+rr14,zMax=thick+L) 
        cells2 = Cells.findAt(
        ((-C + Epsilon, clad+Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-a+Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-Epsilon, thick-Epsilon), ),


        (( a-Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( a+Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( a*0.4, thick-Epsilon, thick-Epsilon), ),


        (( Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( a+a2+Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( a+a2+Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( a+a2+Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( thick+Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( thick+Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( thick+Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( thick+LxPositive-Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( thick+LxPositive-Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( thick+LxPositive-Epsilon, thick-a+Epsilon, thick-Epsilon), ),
        )
        CELLS = cells+cells2    
    region = PipePart.Set(cells=CELLS, name='BaseSteel_Set')
    PipePart.SectionAssignment(region=region, sectionName='BaseSteel_Section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

    PipeModel.Material('WeldSteel')
    PipeModel.materials['WeldSteel'].Density(table=((WeldSteelDens, ), ))
    PipeModel.materials['WeldSteel'].Elastic(table=((WeldSteelE,WeldSteelPratio), ))
    PipeModel.materials['WeldSteel'].Plastic(table=((Yield_Stress_Weld, 0.0),  (Yield_Stress_Weld , Yield_Strain_Weld)))
    PipeModel.HomogeneousSolidSection(name='WeldSteel_Section', material='WeldSteel', thickness=None)
   ##
    cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=0.0,zMax=thick*0.99)
    region = PipePart.Set(cells=cells, name='WeldSteel_Set')
    PipePart.SectionAssignment(region=region, sectionName='WeldSteel_Section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)



    PipeModel.Material('CladSteelWeld')
    PipeModel.materials['CladSteelWeld'].Density(table=((CladSteelWeldDens, ), ))
    PipeModel.materials['CladSteelWeld'].Elastic(table=((CladSteelWeldE,CladSteelWeldPratio), ))
    PipeModel.materials['CladSteelWeld'].Plastic(table=((Yield_Stress_CladSteelWeld, 0.0),  (Yield_Stress_CladSteelWeld , Yield_Strain_CladSteelWeld)))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeld_section', material='CladSteelWeld', thickness=None)
    cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=0.0,zMax=thick*0.99)
    region = PipePart.Set(cells=cells, name='CladSteelWeld_Set')
    PipePart.SectionAssignment(region=region, sectionName='CladSteelWeld_section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)


    PipeModel.Material('CladSteelBase')
    PipeModel.materials['CladSteelBase'].Density(table=((CladSteelBaseDens, ), ))
    PipeModel.materials['CladSteelBase'].Elastic(table=((CladSteelBaseE,CladSteelBasePratio), ))
    PipeModel.materials['CladSteelBase'].Plastic(table=((Yield_Stress_CladSteelBase, 0.0),  (Yield_Stress_CladSteelBase , Yield_Strain_CladSteelBase)))
    PipeModel.HomogeneousSolidSection(name='CladSteelBase_section', material='CladSteelBase', thickness=None)
    if Weld_Angle == 0:
        cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=Weld_Fixed_Width,zMax=thick+L)
    else:
        cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=thick,zMax=thick+L)
        cells2 = Cells.findAt(
        ((-C + Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick - Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick + Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick + LxPositive - Epsilon, clad-Epsilon, thick-Epsilon), ),
        )
        cells = cells1+ cells2
    region = PipePart.Set(cells=cells, name='CladSteelBase_Set')
    PipePart.SectionAssignment(region=region, sectionName='CladSteelBase_section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)



elif MaterialType == 'Plastic Deformation':

    BaseSteelMaterial=PipeModel.Material('BaseSteel')
    BaseSteelMaterial.DeformationPlasticity(table=((BaseSteelE, BaseSteelPratio, Yield_Stress_Base, N_Base, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='BaseSteel_Section', material='BaseSteel', thickness=None)
    ##PART 1 - 
    if Weld_Angle == 0:
        CELLS = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=Weld_Fixed_Width,zMax=thick+L)
    else:
        cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=Weld_Fixed_Width+rr14,zMax=thick+L) 
        cells2 = Cells.findAt(
        ((-C + Epsilon, clad+Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-a+Epsilon, thick-Epsilon), ),
        ((-C + Epsilon, thick-Epsilon, thick-Epsilon), ),


        (( a-Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( a+Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-Epsilon, thick-Epsilon), ),
        (( a*0.4, thick-Epsilon, thick-Epsilon), ),



        (( Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( a+a2+Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( a+a2+Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( a+a2+Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( thick+Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( thick+Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( thick+Epsilon, thick-a+Epsilon, thick-Epsilon), ),

        (( thick+LxPositive-Epsilon, clad+Epsilon, thick-Epsilon), ),
        (( thick+LxPositive-Epsilon, thick-a-Epsilon, thick-Epsilon), ),
        (( thick+LxPositive-Epsilon, thick-a+Epsilon, thick-Epsilon), ),
        )
        CELLS = cells+cells2
    region = PipePart.Set(cells=CELLS, name='BaseSteel_Set')
    PipePart.SectionAssignment(region=region, sectionName='BaseSteel_Section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)



    WeldSteelMaterial=PipeModel.Material('WeldSteel')
    WeldSteelMaterial.DeformationPlasticity(table=((WeldSteelE, WeldSteelPratio, Yield_Stress_Weld, N_Weld, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='WeldSteel_Section', material='WeldSteel', thickness=None) 
    ###
    cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=0.0,zMax=thick*0.99)
    region = PipePart.Set(cells=cells, name='WeldSteel_Set')
    PipePart.SectionAssignment(region=region, sectionName='WeldSteel_Section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)


    CladSteelWeldMaterial=PipeModel.Material('CladSteelWeld')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelWeldE, CladSteelWeldPratio, Yield_Stress_CladSteelWeld, N_CladWeld, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeld_section', material='CladSteelWeld', thickness=None)   
    cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=0.0,zMax=Weld_Fixed_Width+rr6)
    region = PipePart.Set(cells=cells, name='CladSteelWeld_Set')
    PipePart.SectionAssignment(region=region, sectionName='CladSteelWeld_section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)


    CladSteelWeldMaterial=PipeModel.Material('CladSteelBase')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelBaseE, CladSteelBasePratio, Yield_Stress_CladSteelBase, N_CladSteelBase, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelBase_section', material='CladSteelBase', thickness=None) 
    if Weld_Angle == 0:
        cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=Weld_Fixed_Width,zMax=thick+L)
    else:
        cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=thick,zMax=thick+L)
        cells2 = Cells.findAt(
        ((-C + Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick - Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick + Epsilon, clad-Epsilon, thick-Epsilon), ),
        ((thick + LxPositive - Epsilon, clad-Epsilon, thick-Epsilon), ),
        )
        cells = cells1+ cells2

    region = PipePart.Set(cells=cells, name='CladSteelBase_Set')
    PipePart.SectionAssignment(region=region, sectionName='CladSteelBase_section', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)
  
# # ###############################################################################
# ### STEP

if LoadingType == 'Bending':
    PipeModel.StaticStep(name='Bending_Step', previous='Initial',   description='Remote Bending', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
elif LoadingType == 'Int_Pressure':
    PipeModel.StaticStep(name='Int_Pressure_Step', previous='Initial',   description='Internal Pressure', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
elif LoadingType == 'Tensile':
    PipeModel.StaticStep(name='Tensile_Step', previous='Initial',   description='Axial Tensile', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
elif LoadingType == 'Torsion':
    PipeModel.StaticStep(name='Torsion_Step', previous='Initial',   description='Torsion', 
	    timePeriod=100.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)


# # ##############################################################################
# # # ###############################################################################
# # # # ###ASSEMBLY

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

# ########### Surface Set for Internal Pressure

BALA = C + thick + LxPositive

PipeAssambly = PipeModel.rootAssembly
Faces = PipeAssambly.instances['Pipe_Instance'].faces
side1Faces1 = Faces.findAt(
    ((-thick, BALA - Epsilon, -Epsilon), ), 
    ((-thick, BALA - Epsilon, -Weld_Fixed_Width-Epsilon), ), 
    ((-thick, BALA - Epsilon, -thick-Epsilon), ), 
    ((-thick, BALA - Epsilon, -thick-LDiv1-Epsilon), ), 
    ((-thick, BALA - Epsilon, -thick-LDiv2-Epsilon), ), 
    ((-thick, BALA - Epsilon, -thick-LDiv3-Epsilon), ), 

    ((-thick, C + thick + Epsilon, -Epsilon), ), 
    ((-thick, C + thick + Epsilon, -Weld_Fixed_Width-Epsilon), ), 
    ((-thick, C + thick + Epsilon, -thick-Epsilon), ), 
    ((-thick, C + thick + Epsilon, -thick-LDiv1-Epsilon), ), 
    ((-thick, C + thick + Epsilon, -thick-LDiv2-Epsilon), ), 
    ((-thick, C + thick + Epsilon, -thick-LDiv3-Epsilon), ), 

    ((-thick, C + thick - Epsilon, -Epsilon), ), 
    ((-thick, C + thick - Epsilon, -Weld_Fixed_Width-Epsilon), ), 
    ((-thick, C + thick - Epsilon, -thick-Epsilon), ), 
    ((-thick, C + thick - Epsilon, -thick-LDiv1-Epsilon), ), 
    ((-thick, C + thick - Epsilon, -thick-LDiv2-Epsilon), ), 
    ((-thick, C + thick - Epsilon, -thick-LDiv3-Epsilon), ), 

    ((-thick, C + Epsilon, -Epsilon), ), 
    ((-thick, C + Epsilon, -Weld_Fixed_Width-Epsilon), ), 
    ((-thick, C + Epsilon, -thick-Epsilon), ), 
    ((-thick, C + Epsilon, -thick-LDiv1-Epsilon), ), 
    ((-thick, C + Epsilon, -thick-LDiv2-Epsilon), ), 
    ((-thick, C + Epsilon, -thick-LDiv3-Epsilon), ), 

   ((-thick,  Epsilon, -Epsilon), ), 
    ((-thick,  Epsilon, -Weld_Fixed_Width-Epsilon), ), 
    ((-thick,  Epsilon, -thick-Epsilon), ), 
    ((-thick,  Epsilon, -thick-LDiv1-Epsilon), ), 
    ((-thick,  Epsilon, -thick-LDiv2-Epsilon), ), 
    ((-thick,  Epsilon, -thick-LDiv3-Epsilon), ), 
    )
PipeAssambly.Surface(side1Faces=side1Faces1, name='SurfacePressureSet')

side1Faces2 = Faces.findAt(
        ((-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-a-a2-Epsilon, C*0.1, -L), ), 
        ((-thick+Epsilon, C*0.1, -L), ), 

        ((-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-a-a2-Epsilon, C*1.01, -L), ), 
        ((-thick+Epsilon, C*1.01, -L), ), 

        ((-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick-Epsilon, -L), ),

        ((-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick+Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+Epsilon, -L), ), 

        ((-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        
)
PipeAssambly.Surface(side1Faces=side1Faces2, name='Tensile_Faces_Set')



# # # # # # ##### MPC

PipeAssembly = PipeModel.rootAssembly
RefPoints = PipeAssembly.referencePoints
refPoints1=(RefPoints[5], )
region1=PipeAssembly.Set(referencePoints=refPoints1, name='Ref_Point_Set')
Faces = PipeAssembly.instances['Pipe_Instance'].faces
Edges = PipeAssembly.instances['Pipe_Instance'].edges

if LoadingType == 'Bending':
    faces1 = Faces.findAt(
        ((-Epsilon, C*0.1, -L), ), 
        ((-a-Epsilon, C*0.1, -L), ), 
        ((-a-a2-Epsilon, C*0.1, -L), ), 
        ((-thick+Epsilon, C*0.1, -L), ), 

        ((-Epsilon, C*1.01, -L), ), 
        ((-a-Epsilon, C*1.01, -L), ), 
        ((-a-a2-Epsilon, C*1.01, -L), ), 
        ((-thick+Epsilon, C*1.01, -L), ), 

        ((-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick-Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick-Epsilon, -L), ),

        ((-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick+Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+Epsilon, -L), ), 

        ((-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-a-a2-Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        ((-thick+Epsilon, C+thick+LxPositive-Epsilon, -L), ), 
        
    )
    region2=PipeAssembly.Set( faces=faces1, name='MPC_Faces_Set')
    PipeModel.MultipointConstraint(name='MPC_Ref_Point', controlPoint=region1,surface=region2, mpcType=BEAM_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)
  


#######################################
######### BCs
###################################

PipeAssembly = PipeModel.rootAssembly
Faces = PipeAssembly.instances['Pipe_Instance'].faces

BALA = C+thick+LxPositive
faces1 = Faces.findAt(

        # ((-a+Weld_Fixed_Width*1.1, 0.0, -Weld_Fixed_Width*0.1), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -Weld_Fixed_Width*1.1), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -thick*0.99), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -thick*1.1), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -thick-LDiv1*1.1), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -thick-LDiv2*1.1), ), 
        # ((-a+Weld_Fixed_Width*1.1, 0.0, -thick-LDiv3*1.1), ), 

        ((-a+A1*0.01, 0.0, -A1*0.9), ), 
        # ((-a+A2*0.9, 0.0, -A2*0.1), ), 
        ((-a+A2*0.01, 0.0, -A2*0.9), ), 
        # ((-a+A3*0.9, 0.0, -A3*0.1), ), 
        ((-a+A3*0.01, 0.0, -A3*0.9), ), 
        ((-a+A4*0.9, 0.0, -A4*0.1), ), 
        ((-a+A4*0.01, 0.0, -A4*0.9), ), 
        ((-a+A5*0.9, 0.0, -A5*0.1), ), 
        ((-a+A5*0.01, 0.0, -A5*0.9), ),
        ((-a+A6*0.9, 0.0, -A6*0.1), ), 
        ((-a+A6*0.01, 0.0, -A6*0.9), ),
        ((-a+A7*0.9, 0.0, -A7*0.1), ), 
        ((-a+A7*0.01, 0.0, -A7*0.9), ),
        ((-a+A8*0.9, 0.0, -A8*0.1), ), 
        ((-a+A8*0.01, 0.0, -A8*0.9), ),
        ((-a+A9*0.9, 0.0, -A9*0.1), ), 
        ((-a+A9*0.01, 0.0, -A9*0.9), ), 
        ((-a+A10*0.01, 0.0, -A10*0.9), ), 

        # ######
         ((-a-A1*0.9, 0.0, -A1*0.1), ), 
        ((-a-A1*0.01, 0.0, -A1*0.9), ),
        ((-a-A2*0.9, 0.0, -A2*0.1), ), 
        ((-a-A2*0.01, 0.0, -A2*0.9), ),
        ((-a-A3*0.9, 0.0, -A3*0.1), ), 
        ((-a-A3*0.01, 0.0, -A3*0.9), ),
        ((-a-A4*0.9, 0.0, -A4*0.1), ), 
        ((-a-A4*0.01, 0.0, -A4*0.9), ), 
        ((-a-A5*0.9, 0.0, -A5*0.1), ), 
        ((-a-A5*0.01, 0.0, -A5*0.9), ),
        ((-a-A6*0.9, 0.0, -A6*0.1), ), 
        ((-a-A6*0.01, 0.0, -A6*0.9), ),
        ((-a-A7*0.9, 0.0, -A7*0.1), ), 
        ((-a-A7*0.01, 0.0, -A7*0.9), ),
        ((-a-A8*0.9, 0.0, -A8*0.1), ), 
        ((-a-A8*0.01, 0.0, -A8*0.9), ),
        ((-a-A9*0.9, 0.0, -A9*0.1), ), 
        ((-a-A9*0.01, 0.0, -A9*0.9), ), 
        ((-a-A10*0.01, 0.0, -A10*0.9), ), 

        ((-a*0.1, 0.0, -Weld_Fixed_Width*0.1), ), 
        ((-a*0.1, 0.0, -Weld_Fixed_Width*1.1), ), 
        ((-a*0.1, 0.0, -thick*0.99), ), 
        ((-a*0.1, 0.0, -thick*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv1*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv2*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv3*1.1), ),

        ((-a-Weld_Fixed_Width*0.99, 0.0, -Weld_Fixed_Width*0.1), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -Weld_Fixed_Width*1.1), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -thick*0.99), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -thick*1.1), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -thick-LDiv1*1.1), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -thick-LDiv2*1.1), ), 
        ((-a-Weld_Fixed_Width*0.99, 0.0, -thick-LDiv3*1.1), ),


        ((-a*0.1, 0.0, -Weld_Fixed_Width*0.1), ), 
        ((-a*0.1, 0.0, -Weld_Fixed_Width*1.1), ), 
        ((-a*0.1, 0.0, -thick*0.99), ), 
        ((-a*0.1, 0.0, -thick*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv1*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv2*1.1), ), 
        ((-a*0.1, 0.0, -thick-LDiv3*1.1), ),

        ((-a-Weld_Fixed_Width*0.01, 0.0, -Weld_Fixed_Width*0.99), ), 
        ((-a+Weld_Fixed_Width*0.01, 0.0, -Weld_Fixed_Width*0.99), ), 

        ((-thick+clad*1.01, 0.0, -Weld_Fixed_Width*0.1), ), 
        ((-thick+clad*1.01, 0.0, -Weld_Fixed_Width*1.1), ), 
        ((-thick+clad*1.01, 0.0, -thick*0.99), ), 
        ((-thick+clad*1.01, 0.0, -thick*1.1), ), 
        ((-thick+clad*1.01, 0.0, -thick-LDiv1*1.1), ), 
        ((-thick+clad*1.01, 0.0, -thick-LDiv2*1.1), ), 
        ((-thick+clad*1.01, 0.0, -thick-LDiv3*1.1), ),


        
        ((-thick*0.99, 0.0, -Weld_Fixed_Width*0.01), ), 
        ((-thick*0.99, 0.0, -Weld_Fixed_Width*1.01), ), 
        ((-thick*0.99, 0.0, -thick*0.99), ), 
        ((-thick*0.99, 0.0, -thick*1.1), ), 
        ((-thick*0.99, 0.0, -thick-LDiv1*1.1), ), 
        ((-thick*0.99, 0.0, -thick-LDiv2*1.1), ), 
        ((-thick*0.99, 0.0, -thick-LDiv3*1.1), ),


        ###########################################

     
        ((-a*0.1, BALA, -Weld_Fixed_Width*0.1), ), 
        ((-a*0.1, BALA, -Weld_Fixed_Width*1.1), ), 
        ((-a*0.1, BALA, -thick*0.99), ), 
        ((-a*0.1, BALA, -thick*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv1*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv2*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv3*1.1), ),

        ((-a*0.1, BALA, -Weld_Fixed_Width*0.1), ), 
        ((-a*0.1, BALA, -Weld_Fixed_Width*1.1), ), 
        ((-a*0.1, BALA, -thick*0.99), ), 
        ((-a*0.1, BALA, -thick*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv1*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv2*1.1), ), 
        ((-a*0.1, BALA, -thick-LDiv3*1.1), ),

        ((-thick+clad*1.01, BALA, -Weld_Fixed_Width*0.1), ), 
        ((-thick+clad*1.01, BALA, -Weld_Fixed_Width*1.1), ), 
        ((-thick+clad*1.01, BALA, -thick*0.99), ), 
        ((-thick+clad*1.01, BALA, -thick*1.1), ), 
        ((-thick+clad*1.01, BALA, -thick-LDiv1*1.1), ), 
        ((-thick+clad*1.01, BALA, -thick-LDiv2*1.1), ), 
        ((-thick+clad*1.01, BALA, -thick-LDiv3*1.1), ),
        
        ((-thick*0.99, BALA, -Weld_Fixed_Width*0.1), ), 
        ((-thick*0.99, BALA, -Weld_Fixed_Width*1.01), ), 
        ((-thick*0.99, BALA, -thick*0.99), ), 
        ((-thick*0.99, BALA, -thick*1.1), ), 
        ((-thick*0.99, BALA, -thick-LDiv1*1.1), ), 
        ((-thick*0.99, BALA, -thick-LDiv2*1.1), ), 
        ((-thick*0.99, BALA, -thick-LDiv3*1.1), ),
     
)
region = PipeAssembly.Set(faces=faces1, name='Set_Faces_XSYM')

if LoadingType == 'Bending':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Bending_Step', region=region, localCsys=None)
elif LoadingType == 'Int_Pressure':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Int_Pressure_Step', region=region, localCsys=None)
elif LoadingType == 'Tensile':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Tensile_Step', region=region, localCsys=None)
elif LoadingType == 'Torsion':
    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='Torsion_Step', region=region, localCsys=None)


PipeAssembly = PipeModel.rootAssembly
Faces = PipeAssembly.instances['Pipe_Instance'].faces
if KeyHole ==0:
    KK = 0
else:
    KK = KeyHole
faces1 = Faces.findAt(
     ((-thick*0.99, C*0.01, 0.0), ), 
     ((-thick+clad*1.01, C*0.01, 0.0), ), 
     ((-a-Weld_Fixed_Width*0.99, C*0.01, 0.0), ), 

    ((-a-A10*0.99, C*0.01, 0.0), ), 
    ((-a-A9*0.99, C*0.01, 0.0), ), 
    ((-a-A8*0.99, C*0.01, 0.0), ), 
    ((-a-A7*0.99, C*0.01, 0.0), ), 
    ((-a-A6*0.99, C*0.01, 0.0), ), 
    ((-a-A5*0.99, C*0.01, 0.0), ), 
    ((-a-A4*0.99, C*0.01, 0.0), ), 
    ((-a-A3*0.99, C*0.01, 0.0), ), 
    ((-a-A2*0.99, C*0.01, 0.0), ), 
    ((-a-A1*0.99, C*0.01, 0.0), ), 

    ((-thick*0.99, C*1.01, 0.0), ), 
    ((-thick+clad*1.01, C*1.01, 0.0), ), 
    ((-a-Weld_Fixed_Width*0.99, C*1.01, 0.0), ), 

    ((-a-A10*0.99, C*1.01, 0.0), ), 
    ((-a-A9*0.99, C*1.01, 0.0), ), 
    ((-a-A8*0.99, C*1.01, 0.0), ), 
    ((-a-A7*0.99, C*1.01, 0.0), ), 
    ((-a-A6*0.99, C*1.01, 0.0), ), 
    ((-a-A5*0.99, C*1.01, 0.0), ), 
    ((-a-A4*0.99, C*1.01, 0.0), ), 
    ((-a-A3*0.99, C*1.01, 0.0), ), 
    ((-a-A2*0.99, C*1.01, 0.0), ), 
    ((-a-A1*0.99, C*1.01, 0.0), ), 

    ((-thick*0.99, C+a+Weld_Fixed_Width*1.01, 0.0), ), 
    ((-thick+clad*1.01, C+a+Weld_Fixed_Width*1.01, 0.0), ), 
    ((-a-Weld_Fixed_Width*0.99, C+a+Weld_Fixed_Width*1.01, 0.0), ), 
    
    ((-thick*0.99, C+thick*1.01, 0.0), ), 
    ((-thick+clad*1.01, C+thick*1.01, 0.0), ), 
    ((-a-Weld_Fixed_Width*0.99, C+thick*1.01, 0.0), ), 

    ((-thick*0.99, C+thick+LxPositive*0.9, 0.0), ), 
    ((-thick+clad*1.01, C+thick+LxPositive*0.9, 0.0), ), 
    ((-a-Weld_Fixed_Width*0.99, C+thick+LxPositive*0.9, 0.0), ), 
    )
region = PipeAssembly.Set(faces=faces1, name='Set_Faces_ZSYM')

if LoadingType == 'Bending':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Bending_Step', region=region, localCsys=None)
elif LoadingType == 'Int_Pressure':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Int_Pressure_Step', region=region, localCsys=None)
elif LoadingType == 'Tensile':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Tensile_Step', region=region, localCsys=None)
elif LoadingType == 'Torsion':
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='Torsion_Step', region=region, localCsys=None)


# # ##################### Applied LOADs
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

# # # ##### Crack definition

PipeAssembly = PipeModel.rootAssembly
PipeAssembly.makeIndependent(instances=(PipeAssembly.instances['Pipe_Instance'], ))

PipeAssembly = PipeModel.rootAssembly
Edges = PipeAssembly.instances['Pipe_Instance'].edges
edges1 = Edges.findAt(
	((-CR1, C*0.1, -KeyHole), ), 
	((-CR1*cos(pi/8), C+CR1*sin(pi/8), -KeyHole), ),
    ((-CR1*cos(3*pi/8), C+CR1*sin(3*pi/8), -KeyHole), )

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
##############################################

# # # ######### OutPut J-Integral

if LoadingType == 'Bending':
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Bending_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)  
    PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='Bending_Step',timeInterval=1.0, contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)
    # ###History OUTPUT 
    ###OUTPUT Reaction for and displacment
    regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
    PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='Bending_Step', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)  
    ##History output for Rotational Moment and Displacmente 
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Bending_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Bending_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'Int_Pressure':
    PipeAssembly = PipeModel.rootAssembly
    v1 = PipeAssembly.instances['Pipe_Instance'].vertices
    verts1 = v1.findAt(((0.0, 0.0, -KeyHole), ))
    PipeAssembly.Set(vertices=verts1, name='Ref_Point_Set_t')

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
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set_t']
    PipeModel.HistoryOutputRequest(name='rotation', 
        createStepName='Int_Pressure_Step', variables=('UR','RM', ), timeInterval=1.0, 
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    ###Field OUTPUT : stress $ strains
    PipeModel.FieldOutputRequest(name='SS_results', createStepName='Int_Pressure_Step', 
        variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)
elif LoadingType == 'Tensile':
    PipeAssembly = PipeModel.rootAssembly
    v1 = PipeAssembly.instances['Pipe_Instance'].vertices
    verts1 = v1.findAt(((0.0, 0.0, -KeyHole), ))
    PipeAssembly.Set(vertices=verts1, name='Ref_Point_Set_t')

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
    regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set_t']
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


del PipeModel.fieldOutputRequests['F-Output-1']
del PipeModel.historyOutputRequests['H-Output-1']


# # # # ###############################
NCPUs = int(NCPUs)
radvar = str(radius)
vari = str(180/(radius*pi))
PipeModel.keywordBlock.synchVersions(storeNodesAndElements=False)
PipeModel.keywordBlock.insert(35, 
  '\n*NMAP, NSET=ALLN1, TYPE=RECTANGULAR\n%s,0.,0.,999999.,0.,0.\n%s,1.,0.\n1.,1.,1.\n*NMAP, NSET=ALLN1, TYPE=CYLINDRICAL\n0.,0.,0.,0.,0.,1.\n0.,1.,0.\n1.,%s,1.\n*NMAP, NSET=ALLN1, TYPE=RECTANGULAR\n0.,0.,0., 1.,0.,0.\n1.,1.,0.' %(radvar,radvar,vari) )

mdb.Job(name='Ext_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Mat))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick))), model='Pipe_Model_External', description='', type=ANALYSIS, 
     atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
     memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
     explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
     modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
     scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=NCPUs, 
     activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=NCPUs)
PipeModel.setValues(noPartsInputFile=ON)
mdb.jobs['Ext_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Mat))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].writeInput(consistencyChecking=OFF)
#mdb.jobs['Ext_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Mat))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].submit(consistencyChecking=OFF)
#mdb.jobs['Ext_'+LoadingType+'_WA'+str(int(Weld_Angle))+'_N'+str(int(N_Mat))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(theta_pi*100))+'_AT'+str(int(round(10*a/thick)))].waitForCompletion()



