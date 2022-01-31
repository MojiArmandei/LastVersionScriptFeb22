### Moji 

##### New Interface Crack 20/07/2021


#importing main packages
from abaqus import *
from abaqusConstants import *
import regionToolset
import numpy as np 

import sketch
import part
import fields
import mesh
import material
import assembly
import step
import load
import job
import os
import section
import regionToolset

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

##############################
## Checking errors due to Compatability of Commands
backwardCompatibility.setValues(includeDeprecated=True, reportDeprecated=False)
#### This command forces Abaqus to generate Reply file with coordinates and not only index of geometry entity
cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)""")

######################################
# jUST TO RESET ALL PREVIOUS MODELS
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
###############################################
cae_no_parts_input_file=ON 


################################################################################################################
################ FIXED VALUES  ###################################################################################
#################################################################################################################

##### THE WALL THICKNESS SIZE 
thick = 20.6 # thickness in all direction of main body
##### EXTERNAL RADIUS OF PIPE
radius = 206
### Crack Center position
CR_POS = thick * 0.48  #### 0.2 < CR_POS <0.8
##### Postion of Weld Strip from two ends of tube
WELD_POS = 0.5 #### 0.15 < CR_POS <0.85
###############################################################################################################
################  USER INPUT SECTION  ##########################################################################
###############################################################################################################

External = 'YES'

#### 1 )  NUMBER OF CPU USED FOR ANALYSIS (Depends on your computer)
NCPUs = 6

#### 2)  SELECT TYPE OF THE MATERIAL DEFINITION USED IN THE MODEL
MaterialType = 'Plastic-Deformation'  # DEFAULT 5 parametros para 3 materials : Base , Weld , clad , CladWled
# MaterialType = 'Elastic-Plastic' # 3 parametros e uma curva para 3 materials: Base , solda , Clad

#### 3) SELCET THE CRACK DEPTH SIZE 

a = 2.06
# a= 4.12
# a= 6.18
# a= 8.24
# a= 10.3 

#### 4) SELCET THE THETA/PI RATIO FOR C SIZE

Theta_Pi = 0.04
# Theta_Pi = 0.12
# Theta_Pi = 0.2

###############################
########################
# LoadingType = 'Bending'
LoadingType = 'Int_Pressure'
# LoadingType = 'Tensile' ### Tensile applied by Pressure/
# LoadingType = 'Torsion'
# LoadingType = 'ALL'

Applied_Bending = 0.24
Applied_Int_Pressure = 10
Applied_Tensile_Pressure = 5
Applied_Torsion = 0.1

Applied_Bending = 0.24

##### 6 ) SELECT CLAD LAYER THICKNESS
clad = thick * 0.1
#clad = 2.00





#######################################################################################################################
################# FROM HERE ALL VARAIBLES IS CALCULATED AUTOMATICALLY #################################################
########################################################################################################################

#### P.S: If you want to change variable, make sure you know what you are doing
#### If you need any help, please do not hesitate to contact me: mojtaba.armandei@usp.br

D_ext = 2 * radius
L = 5 * D_ext 
C= Theta_Pi *2* D_ext #49.5 #16.5
L1 = L * WELD_POS
L2 = L - L1
LxPositive=((2*pi*radius)/2.0)-(thick+C)

LDiv1_1 = thick
LDiv2_1 = thick + L1 * 0.2
LDiv3_1 = thick + L1 * 0.6 
LDiv1_2 = thick
LDiv2_2 =  L2 * 0.2
LDiv3_2 =  L2 * 0.6 

# clad =2 mm ate 4 mm 
# side of rectangle area arround of the key hole
 
BOX = a * 1.5


A1 = a/20 * 1 
A2 = a/20 * 2
A3 = a/20 * 3
A4 = a/20 * 4
A5 = a/20 * 5
A6 = a/20 * 6
A7 = a/20 * 7
A8 = a/20 * 8
A9 = a/20 * 9

# if a/thick >= 0.49:
#     thick = thick + A1
#### Clad is always internal also is bigger than "a"
a2= a + A1
a3=   a + 2 * A1

nseedCurveCrack = 5 ## curved lines in XY crack (crack face)
nseedcurved = 8 ## curved lines in XZ plane (sweep lines)
nseed_Zero_CR_POS = 4
nseed_CR_POS_thick = 4
nseed_LDiv1_1 = 4
nseed_LDiv2_1 = 4
nseed_LDiv3_1 = 4
nseed_L1 = 4
nseed_LDiv1_2 = 4
nseed_LDiv2_2 = 4
nseed_LDiv3_2 = 4
nseed_L2 = 4

nseeda = 5
nseedC = 10
nseedthickH = 4
nseedthickV = 4

nseedLxPos = 7
nseedBOX = 4
nseedclad = 3


nseedLDiv1 = 4
nseedLDiv2 = 4
nseedLDiv3 = 4
nseedL = 4
nseedMid = 4 ## distance between clad and (thick -a3)



########################### INICIO de ALTERACAO

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

################ FINAL de ALTERACAO



##################################################################################
###############   Crack Block
##################################################################################
#---------
#create the model
mdb.models.changeKey(fromName='Model-1',toName='EX_MOD')
PipeModel=mdb.models['EX_MOD']


#---------
#Create the part
import sketch
import part

#### crack block 1

# a) sketch the pipe plate before rolling using rectangle tool
Orig=[0.0,0.0,0.0]
PipeSketch=PipeModel.ConstrainedSketch(name='Pipe Plate Sketch', sheetSize=100.0)
PipeSketch.setPrimaryObject(option=STANDALONE)

PipeSketch.rectangle(point1=(Orig[0], Orig[1]), point2=(thick, thick)) 

#b) creating the Main part
PipePart = PipeModel.Part(name='CrBlock1', dimensionality=THREE_D,type=DEFORMABLE_BODY)
PipePart = PipeModel.parts['CrBlock1']
PipePart.BaseSolidExtrude(sketch=PipeSketch, depth= thick )
PipeSketch.unsetPrimaryObject()

# # a-1) make a sketch central crack curve 

PipePart = PipeModel.parts['CrBlock1']

Faces=PipePart.faces
Edges=PipePart.edges
Cells=PipePart.cells

Face_Sketch_plane=Faces.findAt(coordinates=(thick*0.3, thick*0.3, 0.0))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, thick*0.3, 0.0))
Origin_Sketch_plane=(0.0, 0.0, 0.0)

transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)



Sketch1.ArcByCenterEnds(center=(0.0, thick), point1=(0.0, thick-a), point2=(-a,  thick), direction=CLOCKWISE)

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane,  faces=Face_Sketch_plane, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


##########################
### Sketching lateral crack circles
#############################################################


PipePart = PipeModel.parts['CrBlock1']

Faces=PipePart.faces
Edges=PipePart.edges
Cells=PipePart.cells


Face_Sketch_plane=Faces.findAt(coordinates=(0., thick*0.3, thick*0.3))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, thick*0.3, 0.0))
Origin_Sketch_plane=(0.0, 0.0, 0.0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Spot(point=(CR_POS, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A1, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A2, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A3, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A4, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A5, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A6, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A7, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A8, thick-a))
Sketch1.CircleByCenterPerimeter(center=(CR_POS, thick-a), point1=(CR_POS-A9, thick-a))


PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane,  faces=Face_Sketch_plane, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

############################################
##### Extruding crack face curve

PipePart = PipeModel.parts['CrBlock1']

Faces=PipePart.faces
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((0.0, thick*0.99, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(a*cos(pi/4), thick - a*sin(pi/4), 0.0)), 
    )
ExtrudeLine= Edges.findAt(coordinates=(0.0, thick, thick*0.01))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#######
### sweeping A1 circles by crack central curve

### A1
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A1*cos(pi/4), CR_POS-A1*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A1*cos(pi/4),CR_POS-A1*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)

### A2
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A2*cos(pi/4), CR_POS-A2*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A2*cos(pi/4),CR_POS-A2*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)


### A3
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A3*cos(pi/4), CR_POS-A3*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A3*cos(pi/4),CR_POS-A3*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)

### A4
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A4*cos(pi/4), CR_POS-A4*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A4*cos(pi/4),CR_POS-A4*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)

### A5
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A5*cos(pi/4), CR_POS-A5*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A5*cos(pi/4),CR_POS-A5*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)

### A6
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A6*cos(pi/4), CR_POS-A6*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A6*cos(pi/4),CR_POS-A6*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)


### A7
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A7*cos(pi/4), CR_POS-A7*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A7*cos(pi/4),CR_POS-A7*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)


### A8
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A8*cos(pi/4), CR_POS-A8*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A8*cos(pi/4),CR_POS-A8*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)

### A9
PipePart = PipeModel.parts['CrBlock1']
Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.01, 0.0, thick*0.01), ),
    ((0.0, thick, thick*0.01), ),    
    )

pickedEdges =(
    Edges.findAt(coordinates=(0, thick - a+A9*cos(pi/4), CR_POS-A9*sin(pi/4))), 
    Edges.findAt( coordinates=( 0, thick - a-A9*cos(pi/4),CR_POS-A9*sin(pi/4)))
    )
ExtrudeLine = Edges.findAt(coordinates=(a * cos(pi/4), thick -a * sin(pi/4), 0.0))

PipePart.PartitionCellBySweepEdge(sweepPath=ExtrudeLine, cells=pickedCells, edges=pickedEdges)



###################################################################################
#### Creating Partition of a vertical line to form crack line

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells
Face_Sketch_plane=Faces.findAt(coordinates=(thick, thick*0.01, thick*0.01))
Edge_Sketch_plane=Edges.findAt(coordinates=(thick, thick*0.01, thick))
Origin_Sketch_plane=(thick, 0.0, 0.0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(CR_POS, 0.0), point2=(CR_POS, -thick))


PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane,  faces=Face_Sketch_plane, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

##### Extruding vertical partion

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=0.0,xMax=thick,yMin=0.0,yMax=thick,zMin=0.0,zMax=thick)

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.9, CR_POS)), 
    )
ExtrudeLine= Edges.findAt(coordinates=(thick*0.99, thick, thick))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)



#####################################################
####################################################
### Creating partion of BOX for improving mesh


PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

Face_Sketch_plane=Faces.findAt(coordinates=(thick*0.01, thick*0.01, 0.0))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, thick*0.01, 0.0))
Origin_Sketch_plane=(0.0, 0.0, 0.0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=61.8, gridSpacing=1.54, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(-BOX, thick), point2=(-BOX, 0.0))
Sketch1.Line(point1=(0.0, thick-BOX), point2=(-thick, thick-BOX))


PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane,  faces=Face_Sketch_plane, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

##### #xtrduing vertical box line

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.9, thick*0.9, thick*0.9), ), 
    ((thick*0.9, thick*0.9, 0.0), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(BOX, thick*0.9, 0.0)), 
    Edges.findAt( coordinates=(BOX, thick*0.1, 0.0))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick, thick, thick*0.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### #xtrduing horizontal box line

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((thick*0.9, thick*0.9, thick*0.9), ), 
    ((thick*0.9, thick*0.9, 0.0), ),

    ((thick*0.01, thick*0.1, thick*0.9), ), 
    ((thick*0.01, thick*0.1, thick*0.), ), 

    )
pickedEdges =(
    Edges.findAt(coordinates=(thick*0.01, thick-BOX, 0.0)), 
    Edges.findAt( coordinates=(thick*0.9, thick-BOX, 0.0))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick, thick, thick*0.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)




##################################################################
#####  Solid Extrude for Lx-Negative (C)
###################################################################

PipePart = PipeModel.parts['CrBlock1']


Face_Sketch_plane=Faces.findAt(coordinates=(0.0, thick*0.01, thick*0.01))
Edge_Sketch_plane=Edges.findAt(coordinates=(0.0, thick*0.01, 0.0))
Origin_Sketch_plane=(0.0, 0.0, 0.0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)
Sketch1.Line(point1=(0.0, thick), point2=(thick, thick))

Sketch1.Line(point1=(thick, thick), point2=(thick, 0.0))
Sketch1.Line(point1=(thick, 0.0), point2=(0.0, 0.0))
Sketch1.Line(point1=(0.0, 0.0), point2=(0.0, thick))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=Sketch1, depth=C, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch1.unsetPrimaryObject()
del mdb.models['EX_MOD'].sketches['__profile__']


#######################################
#### Extruding 


#### A1

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(pi/4), CR_POS+A1*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(3*pi/4), CR_POS+A1*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(3*pi/4), CR_POS-A1*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A1*cos(pi/4), CR_POS-A1*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### A2

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(pi/4), CR_POS+A2*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(3*pi/4), CR_POS+A2*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(3*pi/4), CR_POS-A2*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A2*cos(pi/4), CR_POS-A2*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### A3

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(pi/4), CR_POS+A3*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(3*pi/4), CR_POS+A3*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(3*pi/4), CR_POS-A3*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A3*cos(pi/4), CR_POS-A3*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### A4

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(pi/4), CR_POS+A4*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(3*pi/4), CR_POS+A4*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(3*pi/4), CR_POS-A4*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A4*cos(pi/4), CR_POS-A4*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### A5

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(pi/4), CR_POS+A5*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(3*pi/4), CR_POS+A5*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(3*pi/4), CR_POS-A5*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A5*cos(pi/4), CR_POS-A5*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### A6

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(pi/4), CR_POS+A6*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(3*pi/4), CR_POS+A6*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(3*pi/4), CR_POS-A6*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A6*cos(pi/4), CR_POS-A6*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### A7

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(pi/4), CR_POS+A7*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(3*pi/4), CR_POS+A7*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(3*pi/4), CR_POS-A7*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A7*cos(pi/4), CR_POS-A7*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### A8

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(pi/4), CR_POS+A8*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(3*pi/4), CR_POS+A8*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(3*pi/4), CR_POS-A8*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A8*cos(pi/4), CR_POS-A8*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### A9

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(pi/4), CR_POS+A9*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(3*pi/4), CR_POS+A9*sin(pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(3*pi/4), CR_POS-A9*sin(3*pi/4))), 
    Edges.findAt(coordinates=(0.0, thick-a+A9*cos(pi/4), CR_POS-A9*sin(3*pi/4))), 
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


##### Extruding vertical partion

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick,yMin=0.0,yMax=thick,zMin=0.0,zMax=thick)

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.9, CR_POS)), 
    Edges.findAt(coordinates=(thick, thick*0.01, CR_POS)), 

    )

ExtrudeLine= Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### Extruding horizontal partion

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick,yMin=0.0,yMax=thick,zMin=0.0,zMax=thick)

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-a, thick*0.01)), 

    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A9+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A8+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A7+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A6+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A5+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A4+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A3+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A2+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS-A1/2)), 

    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A2-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A3-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A4-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A5-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A6-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A7-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A8-A1/2)), 
    Edges.findAt(coordinates=(0.0, thick-a, CR_POS+A9-A1/2)), 

    Edges.findAt(coordinates=(0.0, thick-a, thick*0.99)), 

    )

ExtrudeLine= Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


##### Extruding horizontal BOX partion

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells

pickedCells = Cells.findAt(
    ((-thick*0.01, thick*0.01, thick*0.01), ),
    ((-thick*0.01, thick*0.01, thick*0.9), )

    )
pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick-BOX,thick*0.01 )), 
    Edges.findAt(coordinates=(0.0, thick-BOX,thick*0.99 )), 

  
    )
ExtrudeLine = Edges.findAt(coordinates=(-thick*0.01, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine,  cells=pickedCells, edges=pickedEdges, sense=FORWARD)


###################################################################
##### Solid Extrude of LxPositive


PipePart = PipeModel.parts['CrBlock1']


Face_Sketch_plane=Faces.findAt(coordinates=(thick, thick*0.01, thick*0.01))
Edge_Sketch_plane=Edges.findAt(coordinates=(thick, thick*0.01,thick))
Origin_Sketch_plane=(thick, 0.0, 0.0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=71.36, gridSpacing=1.78, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(0.0, 0.0), point2=(thick, 0.0))
Sketch1.Line(point1=(thick, 0.0), point2=(thick, -thick))
Sketch1.Line(point1=(thick, -thick), point2=(0.0, -thick))
Sketch1.Line(point1=(0.0, -thick), point2=(0.0, 0.0))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=Sketch1, depth=LxPositive, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch1.unsetPrimaryObject()
del mdb.models['EX_MOD'].sketches['__profile__']


################################################################
#####
### Extrdue vertical line through Lx postive

PipePart = PipeModel.parts['CrBlock1']

pickedCells = Cells.findAt(
    ((thick*1.1, thick, thick*0.99), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.99, CR_POS)),
    Edges.findAt(coordinates=(thick, thick*0.01, CR_POS))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick*1.1, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

### Extrdue vertical line through Lx postive

PipePart = PipeModel.parts['CrBlock1']

pickedCells = Cells.findAt(
    ((thick*1.1, thick*0.1, thick*0.99), ),
    ((thick*1.1, thick*0.1, thick*0.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(thick, thick -BOX, thick*0.1)),
    Edges.findAt(coordinates=(thick, thick -BOX, thick*0.9))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick*1.1, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


##################################
#### Sketchin LXpostive/2 partition

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells
Faces=PipePart.faces

Face_Sketch_plane=Faces.findAt(coordinates=(thick*1.5, thick*0.01, thick))
Edge_Sketch_plane=Edges.findAt(coordinates=(thick+LxPositive, thick*0.01,thick))
Origin_Sketch_plane=(thick+LxPositive, 0.0, thick)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=575.55, gridSpacing=14.38, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(-(thick + LxPositive/2), thick), point2=((-(thick + LxPositive/2), 0.0)))
pickedFaces = Faces.findAt(
    ((thick*1.1, thick*0.01, thick), ), 
    ((thick*1.1, thick*0.99, thick), )
    )

ExtrudeLine = Edges.findAt(coordinates=(thick+LxPositive, thick*0.01, thick))

PipePart.PartitionFaceBySketch(sketchUpEdge=ExtrudeLine, faces=pickedFaces, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


#### Extrude LxPositive/2 line

PipePart = PipeModel.parts['CrBlock1']

Edges=PipePart.edges
Cells=PipePart.cells
Faces=PipePart.faces

pickedCells = Cells.findAt(
    ((thick*1.1, thick*0.99, thick*0.99), ), 
    ((thick*1.1, thick*0.99, thick*0.01), ), 
    ((thick*1.1, thick*0.01, thick*0.99), ), 
    ((thick*1.1, thick*0.01, thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=( LxPositive/2, thick*0.99, thick)), 
    Edges.findAt(coordinates=(LxPositive/2, thick*0.01, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(thick+LxPositive, thick,  thick*0.99))

PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#############################################################################################
############## Solid Extrude in Z direction L1


PipePart = PipeModel.parts['CrBlock1']


Face_Sketch_plane=Faces.findAt(coordinates=(thick*1.1, thick*0.1, thick))
Edge_Sketch_plane=Edges.findAt(coordinates=(LxPositive+thick, thick*0.1, thick))
Origin_Sketch_plane=(thick+LxPositive, 0.0, thick)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=329.56, gridSpacing=8.23, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)


Sketch1.Line(point1=(0.0, 0.0), point2=(0.0, thick))
Sketch1.Line(point1=(0.0, thick), point2=(-(LxPositive+thick+C), thick))
Sketch1.Line(point1=(-(LxPositive+thick+C), thick), point2=(-(LxPositive+thick+C), 0.0))
Sketch1.Line(point1=(-(LxPositive+thick+C), 0.0), point2=(0.0, 0.0))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=Sketch1, depth=L1, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch1.unsetPrimaryObject()
del mdb.models['EX_MOD'].sketches['__profile__']


########################################################
#### Extrduing LxPositive/2 line across L1

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.99, thick*0.99, thick*1.1), )
    )
pickedEdges =(
    Edges.findAt(coordinates=(LxPositive/2, thick*0.01, thick)), 
    Edges.findAt(coordinates=(LxPositive/2, thick*0.99, thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### Extrduing thick line across L1

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.99, thick*0.99, thick*1.1), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.1, thick)), 
    Edges.findAt(coordinates=(thick, thick*0.9, thick))
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### Extrduing BOX line across L1

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.99, thick*0.99, thick*1.1), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(BOX, thick*0.1, thick)), 
    Edges.findAt(coordinates=(BOX, thick*0.9, thick))
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

##### Extruding Curved curve acroos L1
PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.99, thick*0.99, thick*1.1), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(-C*0.9, thick-a, thick)), 
    Edges.findAt(coordinates=(a*cos(pi/4), thick-a*sin(pi/4), thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### Extrduing X = 0 line across L1

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.99, thick*0.99, thick*1.1), ),
    ((-C*0.99, thick*0.1, thick*1.1), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick*0.1, thick)), 
    Edges.findAt(coordinates=(0.0, thick-a*1.1, thick)), 
    Edges.findAt(coordinates=(0.0, thick*0.9, thick))
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### Extrduing Horizontal BOX line across L1

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0.0,zMax=thick+L1)



pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive*0.9, thick-BOX, thick)), 
    Edges.findAt( coordinates=(thick+LxPositive*0.2, thick-BOX, thick)), 
    Edges.findAt(coordinates=(thick*0.9, thick-BOX, thick)), 
    Edges.findAt( coordinates=(thick*0.01, thick-BOX, thick)),
    Edges.findAt(coordinates=(-C*0.9, thick-BOX,  thick)), 
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

####################################################################
############## Solid Extrude in Z direction L2

PipePart = PipeModel.parts['CrBlock1']


Face_Sketch_plane=Faces.findAt(coordinates=(thick*1.1, thick*0.1, 0))
Edge_Sketch_plane=Edges.findAt(coordinates=(LxPositive+thick, thick*0.1, 0))
Origin_Sketch_plane=(thick+LxPositive, 0.0, 0)


transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=329.56, gridSpacing=8.23, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)


Sketch1.Line(point1=(0.0, 0.0), point2=(-(LxPositive+thick+C), 0.0))
Sketch1.Line(point1=(-(LxPositive+thick+C), 0.0), point2=(-(LxPositive+thick+C), -thick))
Sketch1.Line(point1=(-(LxPositive+thick+C), -thick), point2=(0.0, -thick))
Sketch1.Line(point1=(0.0, -thick), point2=(0.0, 0.0))


PipePart.SolidExtrude(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane,sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=Sketch1, depth=L2, flipExtrudeDirection=OFF, keepInternalBoundaries=ON)
Sketch1.unsetPrimaryObject()
del mdb.models['EX_MOD'].sketches['__profile__']


# #### Extrduing LxPositive/2 line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.9, 0.0, -thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(LxPositive/2, thick*0.99, 0.0)), 
    Edges.findAt(coordinates=(LxPositive/2, thick*0.1, 0.0)), 
 
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


# #### Extrduing thick line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.9, 0.0, -thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(thick, thick*0.99, 0.0)), 
    Edges.findAt(coordinates=(thick, thick*0.1, 0.0)), 
 
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### Extrduing BOX line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.9, 0.0, -thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(BOX, thick*0.99, 0.0)), 
    Edges.findAt(coordinates=(BOX, thick*0.1, 0.0)), 
 
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


#### Extrduing Curved crack line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.9, 0.0, -thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(a*cos(pi/4), thick-a*sin(pi/4), 0.0)), 
    Edges.findAt(coordinates=(-C*0.9, thick-a, 0.0)),  
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

#### Extrduing X=0  line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges

pickedCells = Cells.findAt(
    ((-C*0.9, 0.0, -thick*0.01), ),
    ((-C*0.9, thick*0.99, -thick*0.01), )
    )

pickedEdges =(
    Edges.findAt(coordinates=(0.0, thick*0.99, 0.0)), 
    Edges.findAt(coordinates=(0.0, thick-a*1.1, 0.0)), 
    Edges.findAt(coordinates=(0.0, thick*0.1, 0.0)),  
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)



#### Extrduing Horizontal BOX line across L2

PipePart = PipeModel.parts['CrBlock1']
Cells = PipePart.cells
Edges = PipePart.edges
pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=-thick - L2,zMax=thick)

pickedEdges =(
    Edges.findAt(coordinates=(thick+LxPositive*0.9, thick-BOX, 0)), 
    Edges.findAt( coordinates=(thick+LxPositive*0.2, thick-BOX, 0)), 
    Edges.findAt(coordinates=(thick*0.9, thick-BOX, 0)), 
    Edges.findAt( coordinates=(thick*0.01, thick-BOX, 0)),
    Edges.findAt(coordinates=(-C*0.9, thick-BOX,  0)), 
    )

ExtrudeLine = Edges.findAt(coordinates=(-C, thick, -thick*1.1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

################################################################################
############   Creating Clad line
##################################################################################

PipePart = PipeModel.parts['CrBlock1']

Faces=PipePart.faces
Edges=PipePart.edges
Cells=PipePart.cells

Face_Sketch_plane= Faces.findAt(coordinates=(thick+LxPositive*0.7, thick*0.01, thick+L1))
Edge_Sketch_plane = Edges.findAt(coordinates=(0.0, thick*0.01, thick+L1))
Origin_Sketch_plane =(0.0, 0.0, thick+L1)
# f.findAt(coordinates=(252.486694, 5.493333, 144.2))
# e.findAt(coordinates=(0.0, 6.18, 144.2))
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch = mdb.models['EX_MOD'].ConstrainedSketch(name='__profile__', sheetSize=678.74, gridSpacing=16.96, transform=transformXY)
Sketch.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

Sketch.Line(point1=(C, -clad), point2=(-(thick+LxPositive), -clad))

pickedFaces = Faces.findAt(
    ((thick+LxPositive*0.7, thick*0.01, thick+L1), ), 
    ((thick+LxPositive*0.02, thick*0.01, thick+L1), ), 
    ((thick*0.8, thick*0.01, thick+L1), ), 
    ((BOX * 0.2, thick*0.01,  thick+L1), ), 
    ((-C * 0.4, thick*0.01, thick+L1), )
    )


PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=Sketch)
Sketch.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

### Extruding clad line

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=-thick - L2,zMax=thick+L1)
pickedEdges =(
    Edges.findAt(coordinates=(-C*0.1, clad, thick+L1)), 
    Edges.findAt(coordinates=(BOX*0.1, clad, thick+L1)), 
    Edges.findAt(coordinates=(thick*0.9, clad, thick+L1)), 
    Edges.findAt(coordinates=(thick*1.9, clad, thick+L1)), 
    Edges.findAt(coordinates=(thick+LxPositive*0.9, clad, thick+L1)), 
   
    )


ExtrudeLine = Edges.findAt(coordinates=(thick+LxPositive, 0.0, thick*1.1))

PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)


# # ##################################################
# # ####3 Creating divisions in L1 & L2

PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

Face_Sketch_plane= Faces.findAt(coordinates=(-C, thick-a*0.3,L1 * 0.2))
Edge_Sketch_plane =Edges.findAt(coordinates=(-C, thick*0.1, thick))
Origin_Sketch_plane = (-C, 0.0, thick)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)

Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=825.4, gridSpacing=20.63, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(-LDiv1_1, 0.0), point2=(-LDiv1_1, -thick))
Sketch1.Line(point1=(-LDiv2_1, 0.0), point2=(-LDiv2_1, -thick))
Sketch1.Line(point1=(-LDiv3_1, 0.0), point2=(-LDiv3_1, -thick))

pickedFaces = Faces.findAt(
    ((-C, thick-a*0.3, L1 * 0.2), ),
    ((-C, thick-a*1.1, L1 * 0.2),), 
    ((-C, clad*1.1, L1 * 0.2), ), 
    ((-C, clad*0.3, L1 * 0.2), )
    )

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']

# # ### extruding divisions of L1

### LDiv1_1
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=thick ,zMax=L1+thick)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, LDiv1_1+thick)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, LDiv1_1+thick)), 
    Edges.findAt(coordinates=(-C, clad*1.1, LDiv1_1+thick)),
    Edges.findAt(coordinates=(-C, clad*0.3, LDiv1_1+thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, thick+L1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

## LDiv2_1
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=thick ,zMax=L1+thick)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, LDiv2_1+thick)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, LDiv2_1+thick)), 
    Edges.findAt(coordinates=(-C, clad*1.1, LDiv2_1+thick)),
    Edges.findAt(coordinates=(-C,clad*0.5,LDiv2_1+thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, thick+L1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

### LDiv3_1
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=thick ,zMax=L1+thick)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, LDiv3_1+thick)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, LDiv3_1+thick)), 
    Edges.findAt(coordinates=(-C, clad*1.1, LDiv3_1+thick)),
    Edges.findAt(coordinates=(-C, clad*0.5, LDiv3_1+thick))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, thick+L1))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ##################################################

# ###### sketching divisions in L2

PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

Face_Sketch_plane= Faces.findAt(coordinates=(-C, thick-a*0.3,  - L2 * 0.2))
Edge_Sketch_plane =Edges.findAt(coordinates=(-C, thick*0.1, 0.0))
Origin_Sketch_plane = (-C, 0.0, 0.0)
transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
Sketch1 = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=865.9, gridSpacing=21.64, transform=transformXY)
Sketch1.setPrimaryObject(option=SUPERIMPOSE)
PipePart.projectReferencesOntoSketch(sketch=Sketch1, filter=COPLANAR_EDGES)

Sketch1.Line(point1=(-LDiv1_2, thick), point2=(-LDiv1_2, 0.0))
Sketch1.Line(point1=(-LDiv2_2, thick), point2=(-LDiv2_2, 0.0))
Sketch1.Line(point1=(-LDiv3_2, thick), point2=(-LDiv3_2, 0.0))


pickedFaces = Faces.findAt(
    ((-C,thick-a*0.3, -L2 * 0.2), ),
    ((-C, thick-a*1.1, -L2 * 0.2),), 
    ((-C, clad*1.1, -L2 * 0.2), ), 
    ((-C,clad*0.3, -L2 * 0.2), )
    )

PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=Sketch1)
Sketch1.unsetPrimaryObject()
del PipeModel.sketches['__profile__']


# ### extruding divisions of L2

### LDiv1_2
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=-L2 ,zMax=0)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, -LDiv1_2)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, -(LDiv1_2))), 
    Edges.findAt(coordinates=(-C, clad*1.1, -(LDiv1_2))),
    Edges.findAt(coordinates=(-C, clad*0.3, -(LDiv1_2)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)

# ### LDiv2_2
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=-L2 ,zMax=0)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, -LDiv2_2)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, -(LDiv2_2))), 
    Edges.findAt(coordinates=(-C, clad*1.1, -(LDiv2_2))),
    Edges.findAt(coordinates=(-C, clad*0.3, -(LDiv2_2)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


### LDiv3_2
PipePart = PipeModel.parts['CrBlock1']
Faces= PipePart.faces
Edges= PipePart.edges
Cells = PipePart.cells

pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=-L2 ,zMax=0)

pickedEdges =(
    Edges.findAt(coordinates=(-C, thick-a*0.3, -LDiv3_2)),
    Edges.findAt(coordinates=(-C, thick-a*1.1, -(LDiv3_2))), 
    Edges.findAt(coordinates=(-C, clad*1.1, -(LDiv3_2))),
    Edges.findAt(coordinates=(-C, clad*0.3, -(LDiv3_2)))
    )
ExtrudeLine = Edges.findAt(coordinates=(-C*0.99, thick, 0.0))
PipePart.PartitionCellByExtrudeEdge(line=ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=FORWARD)


# ###### Creating partion LBox arround the crack

if a <= 0.399:


    PipePart = PipeModel.parts['CrBlock1']
    Faces= PipePart.faces
    Edges= PipePart.edges
    Cells = PipePart.cells

    Face_Sketch_plane= Faces.findAt(coordinates=(-C, thick*0.99, CR_POS*0.1))
    Edge_Sketch_plane = Edges.findAt(coordinates=(-C, thick*0.99, thick))
    Origin_Sketch_plane = (-C, thick, 0.0)
    transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
    Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
    Sketch.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

    Sketch.Line(point1=(CR_POS-a, 0), point2=(CR_POS-a, -thick))
    Sketch.Line(point1=(CR_POS+a, 0.0), point2=(CR_POS+a, -thick))


    pickedFaces = Faces.findAt(
        ((-C, thick*0.99, CR_POS*0.1), ), 
        ((-C,thick-a*1.1, CR_POS*0.1), ), 
        ((-C, clad*1.1, CR_POS*0.1), ), 
        ((-C, clad*0.2,CR_POS*0.1), ), 

        ((-C, thick*0.99,thick*0.9), ), 
        ((-C, thick-a*1.1, thick*0.9), ),
        ((-C, clad*1.1, thick*0.9), ),
        ((-C, clad*0.2, thick*0.9), )
        )
    PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=Sketch)
    Sketch.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']


    pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0 ,zMax=thick)
    pickedEdges =(
        Edges.findAt(coordinates=(-C, thick*0.99, CR_POS+a)), 
        Edges.findAt(coordinates=(-C, thick-a*1.1, CR_POS+a)), 
        Edges.findAt(coordinates=(-C, clad*1.1, CR_POS+a)),
        Edges.findAt(coordinates=(-C, clad*0.2, CR_POS+a))
        )
    ExtrudeLine = Edges.findAt(coordinates=(-C*0.2, thick, thick))
    PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

    pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0 ,zMax=thick)
    pickedEdges =(
        Edges.findAt(coordinates=(-C, thick*0.99, CR_POS-a)), 
        Edges.findAt(coordinates=(-C, thick-a*1.1, CR_POS-a)), 
        Edges.findAt(coordinates=(-C, clad*1.1, CR_POS-a)),
        Edges.findAt(coordinates=(-C, clad*0.2, CR_POS-a))
        )
    ExtrudeLine = Edges.findAt(coordinates=(-C*0.2, thick, thick))
    PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

else:
    PipePart = PipeModel.parts['CrBlock1']
    Faces= PipePart.faces
    Edges= PipePart.edges
    Cells = PipePart.cells

    Face_Sketch_plane= Faces.findAt(coordinates=(-C, thick*0.99, CR_POS*0.1))
    Edge_Sketch_plane = Edges.findAt(coordinates=(-C, thick*0.99, thick))
    Origin_Sketch_plane = (-C, thick, 0.0)
    transformXY = PipePart.MakeSketchTransform(sketchPlane=Face_Sketch_plane, sketchUpEdge=Edge_Sketch_plane, sketchPlaneSide=SIDE1, origin=Origin_Sketch_plane)
    Sketch = PipeModel.ConstrainedSketch(name='__profile__', sheetSize=58.26, gridSpacing=1.45, transform=transformXY)
    Sketch.setPrimaryObject(option=SUPERIMPOSE)
    PipePart.projectReferencesOntoSketch(sketch=Sketch, filter=COPLANAR_EDGES)

    Sketch.Line(point1=(CR_POS-a*0.6, 0), point2=(CR_POS-a*0.6, -thick))
    Sketch.Line(point1=(CR_POS+a*0.6, 0.0), point2=(CR_POS+a*0.6, -thick))


    pickedFaces = Faces.findAt(
        ((-C, thick*0.99, CR_POS*0.1), ), 
        ((-C,thick-a*1.1, CR_POS*0.1), ), 
        ((-C, clad*1.1, CR_POS*0.1), ), 
        ((-C, clad*0.2,CR_POS*0.1), ), 

        ((-C, thick*0.99,thick*0.9), ), 
        ((-C, thick-a*1.1, thick*0.9), ),
        ((-C, clad*1.1, thick*0.9), ),
        ((-C, clad*0.2, thick*0.9), )
        )
    PipePart.PartitionFaceBySketch(sketchUpEdge=Edge_Sketch_plane, faces=pickedFaces, sketch=Sketch)
    Sketch.unsetPrimaryObject()
    del PipeModel.sketches['__profile__']


    pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0 ,zMax=thick)
    pickedEdges =(
        Edges.findAt(coordinates=(-C, thick*0.99, CR_POS+a*0.6)), 
        Edges.findAt(coordinates=(-C, thick-a*1.1, CR_POS+a*0.6)), 
        Edges.findAt(coordinates=(-C, clad*1.1, CR_POS+a*0.6)),
        Edges.findAt(coordinates=(-C, clad*0.2, CR_POS+a*0.6))
        )
    ExtrudeLine = Edges.findAt(coordinates=(-C*0.2, thick, thick))
    PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)

    pickedCells=Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=thick,zMin=0 ,zMax=thick)
    pickedEdges =(
        Edges.findAt(coordinates=(-C, thick*0.99, CR_POS-a*0.6)), 
        Edges.findAt(coordinates=(-C, thick-a*1.1, CR_POS-a*0.6)), 
        Edges.findAt(coordinates=(-C, clad*1.1, CR_POS-a*0.6)),
        Edges.findAt(coordinates=(-C, clad*0.2, CR_POS-a*0.6))
        )
    ExtrudeLine = Edges.findAt(coordinates=(-C*0.2, thick, thick))
    PipePart.PartitionCellByExtrudeEdge(line= ExtrudeLine, cells=pickedCells, edges=pickedEdges, sense=REVERSE)



##########################################################


# Create Material


if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='BaseSteel')
    PipeModel.materials['BaseSteel'].Density(table=((BaseSteelDens, ), ))
    PipeModel.materials['BaseSteel'].Elastic(table=((BaseSteelE,  BaseSteelPratio), ))
    PipeModel.materials['BaseSteel'].Plastic(table=((SS_Curve_Base[0][0], SS_Curve_Base[0][1]),  (SS_Curve_Base[1][0], SS_Curve_Base[1][1])))
    PipeModel.HomogeneousSolidSection(name='BaseSteelSection', material='BaseSteel', thickness=None)

if MaterialType == 'Plastic-Deformation':
    BaseSteelMaterial=PipeModel.Material('BaseSteel')
    BaseSteelMaterial.DeformationPlasticity(table=((BaseSteelE, BaseSteelPratio, Yield_Stress_Base, N_Base, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='BaseSteelSection', material='BaseSteel', thickness=None)

cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=-thick-L2,zMax=0)
cells2 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=CR_POS,zMax=thick+L1)
cells = cells1+cells2

region = PipePart.Set(cells=cells, name='SteelBase_Set')
PipePart.SectionAssignment(region=region, sectionName='BaseSteelSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)



if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='WeldSteel')
    PipeModel.materials['WeldSteel'].Density(table=((WeldSteelDens, ), ))
    PipeModel.materials['WeldSteel'].Elastic(table=((WeldSteelE,  WeldSteelPratio), ))
    PipeModel.materials['WeldSteel'].Plastic(table=((SS_Curve_Weld[0][0], SS_Curve_Weld[0][1]),  (SS_Curve_Weld[1][0] , SS_Curve_Weld[1][1] )))
    PipeModel.HomogeneousSolidSection(name='WeldSteelSection', material='WeldSteel', thickness=None)

if MaterialType == 'Plastic-Deformation':
    WeldSteelMaterial=PipeModel.Material('WeldSteel')
    WeldSteelMaterial.DeformationPlasticity(table=((WeldSteelE, WeldSteelPratio, Yield_Stress_Weld, N_Weld, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='WeldSteelSection', material='WeldSteel', thickness=None) 

cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=clad,yMax=thick,zMin=0.0,zMax=CR_POS)
region = PipePart.Set(cells=cells, name='WeldSteel_Set')
PipePart.SectionAssignment(region=region, sectionName='WeldSteelSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)




if MaterialType == 'Elastic-Plastic':
    PipeModel.Material(name='CladSteelWeld')
    PipeModel.materials['CladSteelWeld'].Density(table=((CladSteelWeldDens, ), ))
    PipeModel.materials['CladSteelWeld'].Elastic(table=((CladSteelWeldE,  CladSteelWeldPratio), ))
    PipeModel.materials['CladSteelBase'].Plastic(table=((SS_Curve_CladWeld[0][0], SS_Curve_CladWeld[0][1]),  (SS_Curve_CladWeld[1][0] , SS_Curve_CladWeld[1][1])))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeldSection', material='CladSteelWeld', thickness=None)

if MaterialType == 'Plastic-Deformation':
    CladSteelWeldMaterial=PipeModel.Material('CladSteelWeld')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelWeldE, CladSteelWeldPratio, Yield_Stress_CladSteelWeld, N_clad, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelWeldSection', material='CladSteelWeld', thickness=None) 

cells = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=0.0,zMax=CR_POS)
region = PipePart.Set(cells=cells, name='CladSteelWeld_Set')
PipePart.SectionAssignment(region=region, sectionName='CladSteelWeldSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

if MaterialType == 'Elastic-Plastic':
    PipeModel.Material('CladSteelBase')
    PipeModel.materials['CladSteelBase'].Density(table=((CladSteelBaseDens, ), ))
    PipeModel.materials['CladSteelBase'].Elastic(table=((CladSteelBaseE,CladSteelBasePratio), ))
    PipeModel.materials['CladSteelBase'].Plastic(table=((SS_Curve_Clad[0][0], SS_Curve_Clad[0][1]),  (SS_Curve_Clad[1][0] , SS_Curve_Clad[1][1])))
    PipeModel.HomogeneousSolidSection(name='CladSteelBaseSection', material='CladSteelBase', thickness=None)

if MaterialType == 'Plastic-Deformation':
    CladSteelWeldMaterial=PipeModel.Material('CladSteelBase')
    CladSteelWeldMaterial.DeformationPlasticity(table=((CladSteelBaseE, CladSteelBasePratio, Yield_Stress_CladSteelBase, N_clad, Alpha), ))
    PipeModel.HomogeneousSolidSection(name='CladSteelBaseSection', material='CladSteelBase', thickness=None)   


cells1 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=-thick-L2,zMax=0)
cells2 = Cells.getByBoundingBox(xMin=-C,xMax=thick+LxPositive,yMin=0.0,yMax=clad,zMin=CR_POS,zMax=thick+L1)
cells = cells1+cells2

region = PipePart.Set(cells=cells, name='CladSteelBase_Set')
PipePart.SectionAssignment(region=region, sectionName='CladSteelBaseSection', offset=0.0,  offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

#######################################################################

#####################################################
### Seeding the circles arround 

## partial global meshing
PipePart.seedPart(size=31.0, deviationFactor=0.1, minSizeFactor=0.1)
PipePart.generateMesh()


PipePart = PipeModel.parts['CrBlock1']
Edges = PipePart.edges
pickedEdges = Edges.findAt(
    ### A1
    ((a-A1*cos(pi/4), thick, CR_POS- A1*sin(pi/4)), ),
    ((a-A1*cos(3*pi/4), thick, CR_POS- A1*sin(pi/4)), ),
    ((a-A1*cos(3*pi/4), thick, CR_POS+ A1*sin(pi/4)), ),
    ((a-A1*cos(pi/4), thick, CR_POS+ A1*sin(pi/4)), ),

        ### A2
    ((a-A2*cos(pi/4), thick, CR_POS- A2*sin(pi/4)), ),
    ((a-A2*cos(3*pi/4), thick, CR_POS- A2*sin(pi/4)), ),
    ((a-A2*cos(3*pi/4), thick, CR_POS+ A2*sin(pi/4)), ),
    ((a-A2*cos(pi/4), thick, CR_POS+ A2*sin(pi/4)), ),

     ### A3
    ((a-A3*cos(pi/4), thick, CR_POS- A3*sin(pi/4)), ),
    ((a-A3*cos(3*pi/4), thick, CR_POS- A3*sin(pi/4)), ),
    ((a-A3*cos(3*pi/4), thick, CR_POS+ A3*sin(pi/4)), ),
    ((a-A3*cos(pi/4), thick, CR_POS+ A3*sin(pi/4)), ),

    ### A4
    ((a-A4*cos(pi/4), thick, CR_POS- A4*sin(pi/4)), ),
    ((a-A4*cos(3*pi/4), thick, CR_POS- A4*sin(pi/4)), ),
    ((a-A4*cos(3*pi/4), thick, CR_POS+ A4*sin(pi/4)), ),
    ((a-A4*cos(pi/4), thick, CR_POS+ A4*sin(pi/4)), ),

    ### A5
    ((a-A5*cos(pi/4), thick, CR_POS- A5*sin(pi/4)), ),
    ((a-A5*cos(3*pi/4), thick, CR_POS- A5*sin(pi/4)), ),
    ((a-A5*cos(3*pi/4), thick, CR_POS+ A5*sin(pi/4)), ),
    ((a-A5*cos(pi/4), thick, CR_POS+ A5*sin(pi/4)), ),

    ### A6
    ((a-A6*cos(pi/4), thick, CR_POS- A6*sin(pi/4)), ),
    ((a-A6*cos(3*pi/4), thick, CR_POS- A6*sin(pi/4)), ),
    ((a-A6*cos(3*pi/4), thick, CR_POS+ A6*sin(pi/4)), ),
    ((a-A6*cos(pi/4), thick, CR_POS+ A6*sin(pi/4)), ),

    ### A7
    ((a-A7*cos(pi/4), thick, CR_POS- A7*sin(pi/4)), ),
    ((a-A7*cos(3*pi/4), thick, CR_POS- A7*sin(pi/4)), ),
    ((a-A7*cos(3*pi/4), thick, CR_POS+ A7*sin(pi/4)), ),
    ((a-A7*cos(pi/4), thick, CR_POS+ A7*sin(pi/4)), ),

    ### A8
    ((a-A8*cos(pi/4), thick, CR_POS- A8*sin(pi/4)), ),
    ((a-A8*cos(3*pi/4), thick, CR_POS- A8*sin(pi/4)), ),
    ((a-A8*cos(3*pi/4), thick, CR_POS+ A8*sin(pi/4)), ),
    ((a-A8*cos(pi/4), thick, CR_POS+ A8*sin(pi/4)), ),

    ### A9
    ((a-A9*cos(pi/4), thick, CR_POS- A9*sin(pi/4)), ),
    ((a-A9*cos(3*pi/4), thick, CR_POS- A9*sin(pi/4)), ),
    ((a-A9*cos(3*pi/4), thick, CR_POS+ A9*sin(pi/4)), ),
    ((a-A9*cos(pi/4), thick, CR_POS+ A9*sin(pi/4)), ),

    ######## At x=0
     ### A1
    ((0.0, thick-(a-A1*cos(pi/4)), CR_POS- A1*sin(pi/4)), ),
    ((0.0, thick-(a-A1*cos(3*pi/4)), CR_POS- A1*sin(pi/4)), ),
    ((0.0, thick-(a-A1*cos(3*pi/4)), CR_POS+ A1*sin(pi/4)), ),
    ((0.0, thick-(a-A1*cos(pi/4)), CR_POS+ A1*sin(pi/4)), ), 
     ### A2
    ((0.0, thick-(a-A2*cos(pi/4)), CR_POS- A2*sin(pi/4)), ),
    ((0.0, thick-(a-A2*cos(3*pi/4)), CR_POS- A2*sin(pi/4)), ),
    ((0.0, thick-(a-A2*cos(3*pi/4)), CR_POS+ A2*sin(pi/4)), ),
    ((0.0, thick-(a-A2*cos(pi/4)), CR_POS+ A2*sin(pi/4)), ),

     ### A3
    ((0.0, thick-(a-A3*cos(pi/4)), CR_POS- A3*sin(pi/4)), ),
    ((0.0, thick-(a-A3*cos(3*pi/4)), CR_POS- A3*sin(pi/4)), ),
    ((0.0, thick-(a-A3*cos(3*pi/4)), CR_POS+ A3*sin(pi/4)), ),
    ((0.0, thick-(a-A3*cos(pi/4)), CR_POS+ A3*sin(pi/4)), ),

     ### A4
    ((0.0, thick-(a-A4*cos(pi/4)), CR_POS- A4*sin(pi/4)), ),
    ((0.0, thick-(a-A4*cos(3*pi/4)), CR_POS- A4*sin(pi/4)), ),
    ((0.0, thick-(a-A4*cos(3*pi/4)), CR_POS+ A4*sin(pi/4)), ),
    ((0.0, thick-(a-A4*cos(pi/4)), CR_POS+ A4*sin(pi/4)), ),
     ### A5
    ((0.0, thick-(a-A5*cos(pi/4)), CR_POS- A5*sin(pi/4)), ),
    ((0.0, thick-(a-A5*cos(3*pi/4)), CR_POS- A5*sin(pi/4)), ),
    ((0.0, thick-(a-A5*cos(3*pi/4)), CR_POS+ A5*sin(pi/4)), ),
    ((0.0, thick-(a-A5*cos(pi/4)), CR_POS+ A5*sin(pi/4)), ),

     ### A6
    ((0.0, thick-(a-A6*cos(pi/4)), CR_POS- A6*sin(pi/4)), ),
    ((0.0, thick-(a-A6*cos(3*pi/4)), CR_POS- A6*sin(pi/4)), ),
    ((0.0, thick-(a-A6*cos(3*pi/4)), CR_POS+ A6*sin(pi/4)), ),
    ((0.0, thick-(a-A6*cos(pi/4)), CR_POS+ A6*sin(pi/4)), ),

     ### A7
    ((0.0, thick-(a-A7*cos(pi/4)), CR_POS- A7*sin(pi/4)), ),
    ((0.0, thick-(a-A7*cos(3*pi/4)), CR_POS- A7*sin(pi/4)), ),
    ((0.0, thick-(a-A7*cos(3*pi/4)), CR_POS+ A7*sin(pi/4)), ),
    ((0.0, thick-(a-A7*cos(pi/4)), CR_POS+ A7*sin(pi/4)), ),

    ### A8
    ((0.0, thick-(a-A8*cos(pi/4)), CR_POS- A8*sin(pi/4)), ),
    ((0.0, thick-(a-A8*cos(3*pi/4)), CR_POS- A8*sin(pi/4)), ),
    ((0.0, thick-(a-A8*cos(3*pi/4)), CR_POS+ A8*sin(pi/4)), ),
    ((0.0, thick-(a-A8*cos(pi/4)), CR_POS+ A8*sin(pi/4)), ),

     ### A9
    ((0.0, thick-(a-A9*cos(pi/4)), CR_POS- A9*sin(pi/4)), ),
    ((0.0, thick-(a-A9*cos(3*pi/4)), CR_POS- A9*sin(pi/4)), ),
    ((0.0, thick-(a-A9*cos(3*pi/4)), CR_POS+ A9*sin(pi/4)), ),
    ((0.0, thick-(a-A9*cos(pi/4)), CR_POS+ A9*sin(pi/4)), ),

    ##### At x= -C
     ### A1
    ((-C, thick-(a-A1*cos(pi/4)), CR_POS- A1*sin(pi/4)), ),
    ((-C, thick-(a-A1*cos(3*pi/4)), CR_POS- A1*sin(pi/4)), ),
    ((-C, thick-(a-A1*cos(3*pi/4)), CR_POS+ A1*sin(pi/4)), ),
    ((-C, thick-(a-A1*cos(pi/4)), CR_POS+ A1*sin(pi/4)), ), 
     ### A2
    ((-C, thick-(a-A2*cos(pi/4)), CR_POS- A2*sin(pi/4)), ),
    ((-C, thick-(a-A2*cos(3*pi/4)), CR_POS- A2*sin(pi/4)), ),
    ((-C, thick-(a-A2*cos(3*pi/4)), CR_POS+ A2*sin(pi/4)), ),
    ((-C, thick-(a-A2*cos(pi/4)), CR_POS+ A2*sin(pi/4)), ),

     ### A3
    ((-C, thick-(a-A3*cos(pi/4)), CR_POS- A3*sin(pi/4)), ),
    ((-C, thick-(a-A3*cos(3*pi/4)), CR_POS- A3*sin(pi/4)), ),
    ((-C, thick-(a-A3*cos(3*pi/4)), CR_POS+ A3*sin(pi/4)), ),
    ((-C, thick-(a-A3*cos(pi/4)), CR_POS+ A3*sin(pi/4)), ),

     ### A4
    ((-C, thick-(a-A4*cos(pi/4)), CR_POS- A4*sin(pi/4)), ),
    ((-C, thick-(a-A4*cos(3*pi/4)), CR_POS- A4*sin(pi/4)), ),
    ((-C, thick-(a-A4*cos(3*pi/4)), CR_POS+ A4*sin(pi/4)), ),
    ((-C, thick-(a-A4*cos(pi/4)), CR_POS+ A4*sin(pi/4)), ),
     ### A5
    ((-C, thick-(a-A5*cos(pi/4)), CR_POS- A5*sin(pi/4)), ),
    ((-C, thick-(a-A5*cos(3*pi/4)), CR_POS- A5*sin(pi/4)), ),
    ((-C, thick-(a-A5*cos(3*pi/4)), CR_POS+ A5*sin(pi/4)), ),
    ((-C, thick-(a-A5*cos(pi/4)), CR_POS+ A5*sin(pi/4)), ),

     ### A6
    ((-C, thick-(a-A6*cos(pi/4)), CR_POS- A6*sin(pi/4)), ),
    ((-C, thick-(a-A6*cos(3*pi/4)), CR_POS- A6*sin(pi/4)), ),
    ((-C, thick-(a-A6*cos(3*pi/4)), CR_POS+ A6*sin(pi/4)), ),
    ((-C, thick-(a-A6*cos(pi/4)), CR_POS+ A6*sin(pi/4)), ),

     ### A7
    ((-C, thick-(a-A7*cos(pi/4)), CR_POS- A7*sin(pi/4)), ),
    ((-C, thick-(a-A7*cos(3*pi/4)), CR_POS- A7*sin(pi/4)), ),
    ((-C, thick-(a-A7*cos(3*pi/4)), CR_POS+ A7*sin(pi/4)), ),
    ((-C, thick-(a-A7*cos(pi/4)), CR_POS+ A7*sin(pi/4)), ),

    ### A8
    ((-C, thick-(a-A8*cos(pi/4)), CR_POS- A8*sin(pi/4)), ),
    ((-C, thick-(a-A8*cos(3*pi/4)), CR_POS- A8*sin(pi/4)), ),
    ((-C, thick-(a-A8*cos(3*pi/4)), CR_POS+ A8*sin(pi/4)), ),
    ((-C, thick-(a-A8*cos(pi/4)), CR_POS+ A8*sin(pi/4)), ),

     ### A9
    ((-C, thick-(a-A9*cos(pi/4)), CR_POS- A9*sin(pi/4)), ),
    ((-C, thick-(a-A9*cos(3*pi/4)), CR_POS- A9*sin(pi/4)), ),
    ((-C, thick-(a-A9*cos(3*pi/4)), CR_POS+ A9*sin(pi/4)), ),
    ((-C, thick-(a-A9*cos(pi/4)), CR_POS+ A9*sin(pi/4)), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedCurveCrack, constraint=FINER)


pickedEdges = Edges.findAt(
   ((a*cos(pi/4), thick-(a*cos(pi/4)), CR_POS), ),

    (((a+A1)*cos(pi/4), thick-((a+A1)*cos(pi/4)), CR_POS), ),
    (((a+A2)*cos(pi/4), thick-((a+A2)*cos(pi/4)), CR_POS), ),
    (((a+A3)*cos(pi/4), thick-((a+A3)*cos(pi/4)), CR_POS), ),
    (((a+A4)*cos(pi/4), thick-((a+A4)*cos(pi/4)), CR_POS), ),
    (((a+A5)*cos(pi/4), thick-((a+A5)*cos(pi/4)), CR_POS), ),
    (((a+A6)*cos(pi/4), thick-((a+A6)*cos(pi/4)), CR_POS), ),
    (((a+A7)*cos(pi/4), thick-((a+A7)*cos(pi/4)), CR_POS), ),
    (((a+A8)*cos(pi/4), thick-((a+A8)*cos(pi/4)), CR_POS), ),
    (((a+A9)*cos(pi/4), thick-((a+A9)*cos(pi/4)), CR_POS), ),

     (((a+A1)*cos(pi/4), thick-((a+A1)*cos(pi/4)), CR_POS), ),
    (((a+A2)*cos(pi/4), thick-((a+A2)*cos(pi/4)), CR_POS), ),
    (((a+A3)*cos(pi/4), thick-((a+A3)*cos(pi/4)), CR_POS), ),
    (((a+A4)*cos(pi/4), thick-((a+A4)*cos(pi/4)), CR_POS), ),
    (((a+A5)*cos(pi/4), thick-((a+A5)*cos(pi/4)), CR_POS), ),
    (((a+A6)*cos(pi/4), thick-((a+A6)*cos(pi/4)), CR_POS), ),
    (((a+A7)*cos(pi/4), thick-((a+A7)*cos(pi/4)), CR_POS), ),
    (((a+A8)*cos(pi/4), thick-((a+A8)*cos(pi/4)), CR_POS), ),
    (((a+A9)*cos(pi/4), thick-((a+A9)*cos(pi/4)), CR_POS), ),

    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A1), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A2), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A3), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A4), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A5), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A6), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A7), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A8), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS+A9), ),

    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A1), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A2), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A3), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A4), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A5), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A6), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A7), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A8), ),
    (((a)*cos(pi/4), thick-((a)*cos(pi/4)), CR_POS-A9), ),

   ((a*cos(pi/4), thick-(a*cos(pi/4)), thick), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), 0.0), ),

    ((a*cos(pi/4), thick-(a*cos(pi/4)), thick+LDiv1_1), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), thick+LDiv2_1), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), thick+LDiv3_1), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), thick+L1), ),

    ((a*cos(pi/4), thick-(a*cos(pi/4)), -LDiv1_2), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), -LDiv2_2), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), -LDiv3_2), ),
    ((a*cos(pi/4), thick-(a*cos(pi/4)), -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedcurved, constraint=FINER)


###########################################

pickedEdges = Edges.findAt(
    ((-C, thick, CR_POS*0.1), ),
    ((-C, thick-a, CR_POS*0.1), ),
    ((-C, thick-BOX, CR_POS*0.1), ),
    ((-C, clad, CR_POS*0.1), ),
    ((-C, 0.0, CR_POS*0.1), ),

    ((0.0, thick, CR_POS*0.1), ),
    ((0.0, thick-a, CR_POS*0.1), ),
    ((0.0, thick-BOX, CR_POS*0.1), ),
    ((0.0, clad, CR_POS*0.1), ),
    ((0.0, 0.0, CR_POS*0.1), ),

    ((a, thick, CR_POS*0.1), ),
    ((BOX, thick, CR_POS*0.1), ),


    ((thick, thick, CR_POS*0.1), ),
    ((thick, thick-BOX, CR_POS*0.1), ),
    ((thick, clad, CR_POS*0.1), ),
    ((thick, 0.0, CR_POS*0.1), ),

    ((LxPositive/2, thick, CR_POS*0.1), ),
    ((LxPositive/2, thick-BOX, CR_POS*0.1), ),
    ((LxPositive/2, clad, CR_POS*0.1), ),
    ((LxPositive/2, 0.0, CR_POS*0.1), ),

    ((thick+LxPositive, thick, CR_POS*0.1), ),
    ((thick+LxPositive, thick-BOX, CR_POS*0.1), ),
    ((thick+LxPositive, clad, CR_POS*0.1), ),
    ((thick+LxPositive, 0.0, CR_POS*0.1), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_Zero_CR_POS, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, thick-CR_POS*0.01), ),
    ((-C, thick-a, thick-CR_POS*0.01), ),
    ((-C, thick-BOX, thick-CR_POS*0.01), ),
    ((-C, clad, thick-CR_POS*0.01), ),
    ((-C, 0.0, thick-CR_POS*0.01), ),

    ((0.0, thick, thick-CR_POS*0.01), ),
    ((0.0, thick-a, thick-CR_POS*0.01), ),
    ((0.0, thick-BOX, thick-CR_POS*0.01), ),
    ((0.0, clad, thick-CR_POS*0.01), ),
    ((0.0, 0.0, thick-CR_POS*0.01), ),

    ((a, thick, thick-CR_POS*0.01), ),
    ((BOX, thick, thick-CR_POS*0.01), ),


    ((thick, thick, thick-CR_POS*0.01), ),
    ((thick, thick-BOX, thick-CR_POS*0.01), ),
    ((thick, clad, thick-CR_POS*0.01), ),
    ((thick, 0.0, thick-CR_POS*0.01), ),

    ((LxPositive/2, thick, thick-CR_POS*0.01), ),
    ((LxPositive/2, thick-BOX, thick-CR_POS*0.01), ),
    ((LxPositive/2, clad, thick-CR_POS*0.01), ),
    ((LxPositive/2, 0.0, thick-CR_POS*0.01), ),

    ((thick+LxPositive, thick, thick-CR_POS*0.01), ),
    ((thick+LxPositive, thick-BOX, thick-CR_POS*0.01), ),
    ((thick+LxPositive, clad, thick-CR_POS*0.01), ),
    ((thick+LxPositive, 0.0, thick-CR_POS*0.01), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_CR_POS_thick, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, thick+LDiv1_1*0.9), ),
    ((-C, thick-a, thick+LDiv1_1*0.9), ),
    ((-C, thick-BOX, thick+LDiv1_1*0.9), ),
    ((-C, clad, thick+LDiv1_1*0.9), ),
    ((-C, 0.0, thick+LDiv1_1*0.9), ),

    ((0.0, thick, thick+LDiv1_1*0.9), ),
    ((0.0, thick-a, thick+LDiv1_1*0.9), ),
    ((0.0, thick-BOX, thick+LDiv1_1*0.9), ),
    ((0.0, clad, thick+LDiv1_1*0.9), ),
    ((0.0, 0.0, thick+LDiv1_1*0.9), ),

    ((a, thick, thick+LDiv1_1*0.9), ),
    ((BOX, thick, thick+LDiv1_1*0.9), ),


    ((thick, thick, thick+LDiv1_1*0.9), ),
    ((thick, thick-BOX, thick+LDiv1_1*0.9), ),
    ((thick, clad, thick+LDiv1_1*0.9), ),
    ((thick, 0.0, thick+LDiv1_1*0.9), ),

    ((LxPositive/2, thick, thick+LDiv1_1*0.9), ),
    ((LxPositive/2, thick-BOX, thick+LDiv1_1*0.9), ),
    ((LxPositive/2, clad, thick+LDiv1_1*0.9), ),
    ((LxPositive/2, 0.0, thick+LDiv1_1*0.9), ),

    ((thick+LxPositive, thick, thick+LDiv1_1*0.9), ),
    ((thick+LxPositive, thick-BOX, thick+LDiv1_1*0.9), ),
    ((thick+LxPositive, clad, thick+LDiv1_1*0.9), ),
    ((thick+LxPositive, 0.0, thick+LDiv1_1*0.9), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv1_1, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, thick+LDiv2_1*0.9), ),
    ((-C, thick-a, thick+LDiv2_1*0.9), ),
    ((-C, thick-BOX, thick+LDiv2_1*0.9), ),
    ((-C, clad, thick+LDiv2_1*0.9), ),
    ((-C, 0.0, thick+LDiv2_1*0.9), ),

    ((0.0, thick, thick+LDiv2_1*0.9), ),
    ((0.0, thick-a, thick+LDiv2_1*0.9), ),
    ((0.0, thick-BOX, thick+LDiv2_1*0.9), ),
    ((0.0, clad, thick+LDiv2_1*0.9), ),
    ((0.0, 0.0, thick+LDiv2_1*0.9), ),

    ((a, thick, thick+LDiv2_1*0.9), ),
    ((BOX, thick, thick+LDiv2_1*0.9), ),


    ((thick, thick, thick+LDiv2_1*0.9), ),
    ((thick, thick-BOX, thick+LDiv2_1*0.9), ),
    ((thick, clad, thick+LDiv2_1*0.9), ),
    ((thick, 0.0, thick+LDiv2_1*0.9), ),

    ((LxPositive/2, thick, thick+LDiv2_1*0.9), ),
    ((LxPositive/2, thick-BOX, thick+LDiv2_1*0.9), ),
    ((LxPositive/2, clad, thick+LDiv2_1*0.9), ),
    ((LxPositive/2, 0.0, thick+LDiv2_1*0.9), ),

    ((thick+LxPositive, thick, thick+LDiv2_1*0.9), ),
    ((thick+LxPositive, thick-BOX, thick+LDiv2_1*0.9), ),
    ((thick+LxPositive, clad, thick+LDiv2_1*0.9), ),
    ((thick+LxPositive, 0.0, thick+LDiv2_1*0.9), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv2_1, constraint=FINER)



pickedEdges = Edges.findAt(
    ((-C, thick, thick+LDiv3_1*0.9), ),
    ((-C, thick-a, thick+LDiv3_1*0.9), ),
    ((-C, thick-BOX, thick+LDiv3_1*0.9), ),
    ((-C, clad, thick+LDiv3_1*0.9), ),
    ((-C, 0.0, thick+LDiv3_1*0.9), ),

    ((0.0, thick, thick+LDiv3_1*0.9), ),
    ((0.0, thick-a, thick+LDiv3_1*0.9), ),
    ((0.0, thick-BOX, thick+LDiv3_1*0.9), ),
    ((0.0, clad, thick+LDiv3_1*0.9), ),
    ((0.0, 0.0, thick+LDiv3_1*0.9), ),

    ((a, thick, thick+LDiv3_1*0.9), ),
    ((BOX, thick, thick+LDiv3_1*0.9), ),


    ((thick, thick, thick+LDiv3_1*0.9), ),
    ((thick, thick-BOX, thick+LDiv3_1*0.9), ),
    ((thick, clad, thick+LDiv3_1*0.9), ),
    ((thick, 0.0, thick+LDiv3_1*0.9), ),

    ((LxPositive/2, thick, thick+LDiv3_1*0.9), ),
    ((LxPositive/2, thick-BOX, thick+LDiv3_1*0.9), ),
    ((LxPositive/2, clad, thick+LDiv3_1*0.9), ),
    ((LxPositive/2, 0.0, thick+LDiv3_1*0.9), ),

    ((thick+LxPositive, thick, thick+LDiv3_1*0.9), ),
    ((thick+LxPositive, thick-BOX, thick+LDiv3_1*0.9), ),
    ((thick+LxPositive, clad, thick+LDiv3_1*0.9), ),
    ((thick+LxPositive, 0.0, thick+LDiv3_1*0.9), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv3_1, constraint=FINER)


pickedEdges = Edges.findAt(
    ((-C, thick, thick+LDiv3_1*1.2), ),
    ((-C, thick-a, thick+LDiv3_1*1.2), ),
    ((-C, thick-BOX, thick+LDiv3_1*1.2), ),
    ((-C, clad, thick+LDiv3_1*1.2), ),
    ((-C, 0.0, thick+LDiv3_1*1.2), ),

    ((0.0, thick, thick+LDiv3_1*1.2), ),
    ((0.0, thick-a, thick+LDiv3_1*1.2), ),
    ((0.0, thick-BOX, thick+LDiv3_1*1.2), ),
    ((0.0, clad, thick+LDiv3_1*1.2), ),
    ((0.0, 0.0, thick+LDiv3_1*1.2), ),

    ((a, thick, thick+LDiv3_1*1.2), ),
    ((BOX, thick, thick+LDiv3_1*1.2), ),


    ((thick, thick, thick+LDiv3_1*1.2), ),
    ((thick, thick-BOX, thick+LDiv3_1*1.2), ),
    ((thick, clad, thick+LDiv3_1*1.2), ),
    ((thick, 0.0, thick+LDiv3_1*1.2), ),

    ((LxPositive/2, thick, thick+LDiv3_1*1.2), ),
    ((LxPositive/2, thick-BOX, thick+LDiv3_1*1.2), ),
    ((LxPositive/2, clad, thick+LDiv3_1*1.2), ),
    ((LxPositive/2, 0.0, thick+LDiv3_1*1.2), ),

    ((thick+LxPositive, thick, thick+LDiv3_1*1.2), ),
    ((thick+LxPositive, thick-BOX, thick+LDiv3_1*1.2), ),
    ((thick+LxPositive, clad, thick+LDiv3_1*1.2), ),
    ((thick+LxPositive, 0.0, thick+LDiv3_1*1.2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_L1, constraint=FINER)


pickedEdges = Edges.findAt(
    ((-C, thick, -LDiv1_2*0.99), ),
    ((-C, thick-a, -LDiv1_2*0.99), ),
    ((-C, thick-BOX, -LDiv1_2*0.99), ),
    ((-C, clad, -LDiv1_2*0.99), ),
    ((-C, 0.0, -LDiv1_2*0.99), ),

    ((0.0, thick, -LDiv1_2*0.99), ),
    ((0.0, thick-a, -LDiv1_2*0.99), ),
    ((0.0, thick-BOX, -LDiv1_2*0.99), ),
    ((0.0, clad, -LDiv1_2*0.99), ),
    ((0.0, 0.0, -LDiv1_2*0.99), ),

    ((a, thick, -LDiv1_2*0.99), ),
    ((BOX, thick, -LDiv1_2*0.99), ),


    ((thick, thick, -LDiv1_2*0.99), ),
    ((thick, thick-BOX, -LDiv1_2*0.99), ),
    ((thick, clad, -LDiv1_2*0.99), ),
    ((thick, 0.0, -LDiv1_2*0.99), ),

    ((LxPositive/2, thick, -LDiv1_2*0.99), ),
    ((LxPositive/2, thick-BOX, -LDiv1_2*0.99), ),
    ((LxPositive/2, clad, -LDiv1_2*0.99), ),
    ((LxPositive/2, 0.0, -LDiv1_2*0.99), ),

    ((thick+LxPositive, thick, -LDiv1_2*0.99), ),
    ((thick+LxPositive, thick-BOX, -LDiv1_2*0.99), ),
    ((thick+LxPositive, clad, -LDiv1_2*0.99), ),
    ((thick+LxPositive, 0.0, -LDiv1_2*0.99), ),
)

PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv1_2, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, -LDiv1_2*1.01), ),
    ((-C, thick-a, -LDiv1_2*1.01), ),
    ((-C, thick-BOX, -LDiv1_2*1.01), ),
    ((-C, clad, -LDiv1_2*1.01), ),
    ((-C, 0.0, -LDiv1_2*1.01), ),

    ((0.0, thick, -LDiv1_2*1.01), ),
    ((0.0, thick-a, -LDiv1_2*1.01), ),
    ((0.0, thick-BOX, -LDiv1_2*1.01), ),
    ((0.0, clad, -LDiv1_2*1.01), ),
    ((0.0, 0.0, -LDiv1_2*1.01), ),

    ((a, thick, -LDiv1_2*1.01), ),
    ((BOX, thick, -LDiv1_2*1.01), ),


    ((thick, thick, -LDiv1_2*1.01), ),
    ((thick, thick-BOX, -LDiv1_2*1.01), ),
    ((thick, clad, -LDiv1_2*1.01), ),
    ((thick, 0.0, -LDiv1_2*1.01), ),

    ((LxPositive/2, thick, -LDiv1_2*1.01), ),
    ((LxPositive/2, thick-BOX, -LDiv1_2*1.01), ),
    ((LxPositive/2, clad, -LDiv1_2*1.01), ),
    ((LxPositive/2, 0.0, -LDiv1_2*1.01), ),

    ((thick+LxPositive, thick, -LDiv1_2*1.01), ),
    ((thick+LxPositive, thick-BOX, -LDiv1_2*1.01), ),
    ((thick+LxPositive, clad, -LDiv1_2*1.01), ),
    ((thick+LxPositive, 0.0, -LDiv1_2*1.01), ),
)

PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv2_2, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, -LDiv2_2*1.01), ),
    ((-C, thick-a, -LDiv2_2*1.01), ),
    ((-C, thick-BOX, -LDiv2_2*1.01), ),
    ((-C, clad, -LDiv2_2*1.01), ),
    ((-C, 0.0, -LDiv2_2*1.01), ),

    ((0.0, thick, -LDiv2_2*1.01), ),
    ((0.0, thick-a, -LDiv2_2*1.01), ),
    ((0.0, thick-BOX, -LDiv2_2*1.01), ),
    ((0.0, clad, -LDiv2_2*1.01), ),
    ((0.0, 0.0, -LDiv2_2*1.01), ),

    ((a, thick, -LDiv2_2*1.01), ),
    ((BOX, thick, -LDiv2_2*1.01), ),


    ((thick, thick, -LDiv2_2*1.01), ),
    ((thick, thick-BOX, -LDiv2_2*1.01), ),
    ((thick, clad, -LDiv2_2*1.01), ),
    ((thick, 0.0, -LDiv2_2*1.01), ),

    ((LxPositive/2, thick, -LDiv2_2*1.01), ),
    ((LxPositive/2, thick-BOX, -LDiv2_2*1.01), ),
    ((LxPositive/2, clad, -LDiv2_2*1.01), ),
    ((LxPositive/2, 0.0, -LDiv2_2*1.01), ),

    ((thick+LxPositive, thick, -LDiv2_2*1.01), ),
    ((thick+LxPositive, thick-BOX, -LDiv2_2*1.01), ),
    ((thick+LxPositive, clad, -LDiv2_2*1.01), ),
    ((thick+LxPositive, 0.0, -LDiv2_2*1.01), ),
)

PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_LDiv3_2, constraint=FINER)

pickedEdges = Edges.findAt(
    ((-C, thick, -LDiv3_2*1.01), ),
    ((-C, thick-a, -LDiv3_2*1.01), ),
    ((-C, thick-BOX, -LDiv3_2*1.01), ),
    ((-C, clad, -LDiv3_2*1.01), ),
    ((-C, 0.0, -LDiv3_2*1.01), ),

    ((0.0, thick, -LDiv3_2*1.01), ),
    ((0.0, thick-a, -LDiv3_2*1.01), ),
    ((0.0, thick-BOX, -LDiv3_2*1.01), ),
    ((0.0, clad, -LDiv3_2*1.01), ),
    ((0.0, 0.0, -LDiv3_2*1.01), ),

    ((a, thick, -LDiv3_2*1.01), ),
    ((BOX, thick, -LDiv3_2*1.01), ),


    ((thick, thick, -LDiv3_2*1.01), ),
    ((thick, thick-BOX, -LDiv3_2*1.01), ),
    ((thick, clad, -LDiv3_2*1.01), ),
    ((thick, 0.0, -LDiv3_2*1.01), ),

    ((LxPositive/2, thick, -LDiv3_2*1.01), ),
    ((LxPositive/2, thick-BOX, -LDiv3_2*1.01), ),
    ((LxPositive/2, clad, -LDiv3_2*1.01), ),
    ((LxPositive/2, 0.0, -LDiv3_2*1.01), ),

    ((thick+LxPositive, thick, -LDiv3_2*1.01), ),
    ((thick+LxPositive, thick-BOX, -LDiv3_2*1.01), ),
    ((thick+LxPositive, clad, -LDiv3_2*1.01), ),
    ((thick+LxPositive, 0.0, -LDiv3_2*1.01), ),
)

PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseed_L2, constraint=FINER)


pickedEdges = Edges.findAt(

    ((-C*0.9, thick, 0.0), ),
    ((-C*0.9, thick-a, 0.0), ),
    ((-C*0.9, thick-BOX, 0.0), ),
    ((-C*0.9, clad, 0.0), ),
    ((-C*0.9, 0.0, 0.0), ),

    ((-C*0.9, thick, CR_POS), ),
    ((-C*0.9, thick-a+A9, CR_POS), ),
    ((-C*0.9, thick-a+A8, CR_POS), ),
    ((-C*0.9, thick-a+A7, CR_POS), ),
    ((-C*0.9, thick-a+A6, CR_POS), ),
    ((-C*0.9, thick-a+A5, CR_POS), ),
    ((-C*0.9, thick-a+A4, CR_POS), ),
    ((-C*0.9, thick-a+A3, CR_POS), ),
    ((-C*0.9, thick-a+A2, CR_POS), ),
    ((-C*0.9, thick-a+A1, CR_POS), ),

    ((-C*0.9, thick-a, CR_POS), ),
    ((-C*0.9, thick-a, CR_POS-A9), ),
    ((-C*0.9, thick-a, CR_POS-A8), ),
    ((-C*0.9, thick-a, CR_POS-A7), ),
    ((-C*0.9, thick-a, CR_POS-A6), ),
    ((-C*0.9, thick-a, CR_POS-A5), ),
    ((-C*0.9, thick-a, CR_POS-A4), ),
    ((-C*0.9, thick-a, CR_POS-A3), ),
    ((-C*0.9, thick-a, CR_POS-A2), ),
    ((-C*0.9, thick-a, CR_POS-A1), ),

    ((-C*0.9, thick-a, CR_POS+A1), ),
    ((-C*0.9, thick-a, CR_POS+A2), ),
    ((-C*0.9, thick-a, CR_POS+A3), ),
    ((-C*0.9, thick-a, CR_POS+A4), ),
    ((-C*0.9, thick-a, CR_POS+A5), ),
    ((-C*0.9, thick-a, CR_POS+A6), ),
    ((-C*0.9, thick-a, CR_POS+A7), ),
    ((-C*0.9, thick-a, CR_POS+A8), ),
    ((-C*0.9, thick-a, CR_POS+A9), ),

    ((-C*0.9, thick-a-A1, CR_POS), ),
    ((-C*0.9, thick-a-A2, CR_POS), ),
    ((-C*0.9, thick-a-A3, CR_POS), ),
    ((-C*0.9, thick-a-A4, CR_POS), ),
    ((-C*0.9, thick-a-A5, CR_POS), ),
    ((-C*0.9, thick-a-A6, CR_POS), ),
    ((-C*0.9, thick-a-A7, CR_POS), ),
    ((-C*0.9, thick-a-A8, CR_POS), ),
    ((-C*0.9, thick-a-A9, CR_POS), ),

    ((-C*0.9, thick-BOX, CR_POS), ),
    ((-C*0.9, clad, CR_POS), ),
    ((-C*0.9, 0.0, CR_POS), ),

    ((-C*0.9, thick, thick), ),
    ((-C*0.9, thick-a, thick), ),
    ((-C*0.9, thick-BOX, thick), ),
    ((-C*0.9, clad, thick), ),
    ((-C*0.9, 0.0, thick), ),

    ((-C*0.9, thick, thick+LDiv1_1), ),
    ((-C*0.9, thick-a, thick+LDiv1_1), ),
    ((-C*0.9, thick-BOX, thick+LDiv1_1), ),
    ((-C*0.9, clad, thick+LDiv1_1), ),
    ((-C*0.9, 0.0, thick+LDiv1_1), ),

    ((-C*0.9, thick, thick+LDiv2_1), ),
    ((-C*0.9, thick-a, thick+LDiv2_1), ),
    ((-C*0.9, thick-BOX, thick+LDiv2_1), ),
    ((-C*0.9, clad, thick+LDiv2_1), ),
    ((-C*0.9, 0.0, thick+LDiv2_1), ),

    ((-C*0.9, thick, thick+LDiv3_1), ),
    ((-C*0.9, thick-a, thick+LDiv3_1), ),
    ((-C*0.9, thick-BOX, thick+LDiv3_1), ),
    ((-C*0.9, clad, thick+LDiv3_1), ),
    ((-C*0.9, 0.0, thick+LDiv3_1), ),

    ((-C*0.9, thick, thick+L1), ),
    ((-C*0.9, thick-a, thick+L1), ),
    ((-C*0.9, thick-BOX, thick+L1), ),
    ((-C*0.9, clad, thick+L1), ),
    ((-C*0.9, 0.0, thick+L1), ),

    ((-C*0.9, thick, -LDiv1_2), ),
    ((-C*0.9, thick-a, -LDiv1_2), ),
    ((-C*0.9, thick-BOX, -LDiv1_2), ),
    ((-C*0.9, clad, -LDiv1_2), ),
    ((-C*0.9, 0.0, -LDiv1_2), ),

    ((-C*0.9, thick, -LDiv2_2), ),
    ((-C*0.9, thick-a, -LDiv2_2), ),
    ((-C*0.9, thick-BOX, -LDiv2_2), ),
    ((-C*0.9, clad, -LDiv2_2), ),
    ((-C*0.9, 0.0, -LDiv2_2), ),

    ((-C*0.9, thick, -LDiv3_2), ),
    ((-C*0.9, thick-a, -LDiv3_2), ),
    ((-C*0.9, thick-BOX, -LDiv3_2), ),
    ((-C*0.9, clad, -LDiv3_2), ),
    ((-C*0.9, 0.0, -LDiv3_2), ),

    ((-C*0.9, thick, -L2), ),
    ((-C*0.9, thick-a, -L2), ),
    ((-C*0.9, thick-BOX, -L2), ),
    ((-C*0.9, clad, -L2), ),
    ((-C*0.9, 0.0, -L2), ),







)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedC, constraint=FINER)

pickedEdges = Edges.findAt(
    ((a*0.01, thick, 0.0), ),
    ((a*0.01, thick-BOX, 0.0), ),
    ((a*0.01, clad, 0.0), ),
    ((a*0.01, 0.0, 0.0), ),

    ((a*0.01, thick, CR_POS), ),
    ((a*0.01, thick-BOX, CR_POS), ),
    ((a*0.01, clad, CR_POS), ),
    ((a*0.01, 0.0, CR_POS), ),

    ((a*0.01, thick, thick), ),
    ((a*0.01, thick-BOX, thick), ),
    ((a*0.01, clad, thick), ),
    ((a*0.01, 0.0, thick), ),

    ((a*0.01, thick, thick+LDiv1_1), ),
    ((a*0.01, thick-BOX, thick+LDiv1_1), ),
    ((a*0.01, clad, thick+LDiv1_1), ),
    ((a*0.01, 0.0, thick+LDiv1_1), ),

    ((a*0.01, thick, thick+LDiv2_1), ),
    ((a*0.01, thick-BOX, thick+LDiv2_1), ),
    ((a*0.01, clad, thick+LDiv2_1), ),
    ((a*0.01, 0.0, thick+LDiv2_1), ),

    ((a*0.01, thick, thick+LDiv3_1), ),
    ((a*0.01, thick-BOX, thick+LDiv3_1), ),
    ((a*0.01, clad, thick+LDiv3_1), ),
    ((a*0.01, 0.0, thick+LDiv3_1), ),

    ((a*0.01, thick, thick+L1), ),
    ((a*0.01, thick-BOX, thick+L1), ),
    ((a*0.01, clad, thick+L1), ),
    ((a*0.01, 0.0, thick+L1), ),


    ((a*0.01, thick, -LDiv1_2), ),
    ((a*0.01, thick-BOX, -LDiv1_2), ),
    ((a*0.01, clad, -LDiv1_2), ),
    ((a*0.01, 0.0, -LDiv1_2), ),

    ((a*0.01, thick, -LDiv2_2), ),
    ((a*0.01, thick-BOX, -LDiv2_2), ),
    ((a*0.01, clad, -LDiv2_2), ),
    ((a*0.01, 0.0, -LDiv2_2), ),

    ((a*0.01, thick, -LDiv3_2), ),
    ((a*0.01, thick-BOX, -LDiv3_2), ),
    ((a*0.01, clad, -LDiv3_2), ),
    ((a*0.01, 0.0, -LDiv3_2), ),

    ((a*0.01, thick, -L2), ),
    ((a*0.01, thick-BOX, -L2), ),
    ((a*0.01, clad, -L2), ),
    ((a*0.01, 0.0, -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseeda, constraint=FINER)

pickedEdges = Edges.findAt(

    ((thick*0.99, thick, -L2), ),
    ((thick*0.99, thick-BOX, -L2), ),
    ((thick*0.99, clad, -L2), ),
    ((thick*0.99, 0.0, -L2), ),

    ((thick*0.99, thick, -LDiv3_2), ),
    ((thick*0.99, thick-BOX, -LDiv3_2), ),
    ((thick*0.99, clad, -LDiv3_2), ),
    ((thick*0.99, 0.0, -LDiv3_2), ),

    ((thick*0.99, thick, -LDiv2_2), ),
    ((thick*0.99, thick-BOX, -LDiv2_2), ),
    ((thick*0.99, clad, -LDiv2_2), ),
    ((thick*0.99, 0.0, -LDiv2_2), ),

    ((thick*0.99, thick, -LDiv1_2), ),
    ((thick*0.99, thick-BOX, -LDiv1_2), ),
    ((thick*0.99, clad, -LDiv1_2), ),
    ((thick*0.99, 0.0, -LDiv1_2), ),

    ((thick*0.99, thick, 0.0), ),
    ((thick*0.99, thick-BOX, 0.0), ),
    ((thick*0.99, clad, 0.0), ),
    ((thick*0.99, 0.0, 0.0), ),

    ((thick*0.99, thick, CR_POS), ),
    ((thick*0.99, thick-BOX, CR_POS), ),
    ((thick*0.99, clad, CR_POS), ),
    ((thick*0.99, 0.0, CR_POS), ),

    ((thick*0.99, thick, thick), ),
    ((thick*0.99, thick-BOX, thick), ),
    ((thick*0.99, clad, thick), ),
    ((thick*0.99, 0.0, thick), ),

    ((thick*0.99, thick, thick+LDiv1_1), ),
    ((thick*0.99, thick-BOX, thick+LDiv1_1), ),
    ((thick*0.99, clad, thick+LDiv1_1), ),
    ((thick*0.99, 0.0, thick+LDiv1_1), ),

    
    ((thick*0.99, thick, thick+LDiv2_1), ),
    ((thick*0.99, thick-BOX, thick+LDiv2_1), ),
    ((thick*0.99, clad, thick+LDiv2_1), ),
    ((thick*0.99, 0.0, thick+LDiv2_1), ),

    
    ((thick*0.99, thick, thick+LDiv3_1), ),
    ((thick*0.99, thick-BOX, thick+LDiv3_1), ),
    ((thick*0.99, clad, thick+LDiv3_1), ),
    ((thick*0.99, 0.0, thick+LDiv3_1), ),

    
    ((thick*0.99, thick, thick+L1), ),
    ((thick*0.99, thick-BOX, thick+L1), ),
    ((thick*0.99, clad, thick+L1), ),
    ((thick*0.99, 0.0, thick+L1), ),

)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedthickH, constraint=FINER)

pickedEdges = Edges.findAt(

    ((thick*1.01, thick, thick+L1), ),
    ((thick*1.01, thick-BOX, thick+L1), ),
    ((thick*1.01, clad, thick+L1), ),
    ((thick*1.01, 0.0, thick+L1), ),


    ((thick*1.01, thick, thick+LDiv3_1), ),
    ((thick*1.01, thick-BOX, thick+LDiv3_1), ),
    ((thick*1.01, clad, thick+LDiv3_1), ),
    ((thick*1.01, 0.0, thick+LDiv3_1), ),

    ((thick*1.01, thick, thick+LDiv2_1), ),
    ((thick*1.01, thick-BOX, thick+LDiv2_1), ),
    ((thick*1.01, clad, thick+LDiv2_1), ),
    ((thick*1.01, 0.0, thick+LDiv2_1), ),

    ((thick*1.01, thick, thick+LDiv1_1), ),
    ((thick*1.01, thick-BOX, thick+LDiv1_1), ),
    ((thick*1.01, clad, thick+LDiv1_1), ),
    ((thick*1.01, 0.0, thick+LDiv1_1), ),

    ((thick*1.01, thick, thick), ),
    ((thick*1.01, thick-BOX, thick), ),
    ((thick*1.01, clad, thick), ),
    ((thick*1.01, 0.0, thick), ),

    ((thick*1.01, thick, CR_POS), ),
    ((thick*1.01, thick-BOX, CR_POS), ),
    ((thick*1.01, clad, CR_POS), ),
    ((thick*1.01, 0.0, CR_POS), ),

    ((thick*1.01, thick, 0.0), ),
    ((thick*1.01, thick-BOX, 0.0), ),
    ((thick*1.01, clad, 0.0), ),
    ((thick*1.01, 0.0, 0.0), ),

    ((thick*1.01, thick, -LDiv1_2), ),
    ((thick*1.01, thick-BOX, -LDiv1_2), ),
    ((thick*1.01, clad, -LDiv1_2), ),
    ((thick*1.01, 0.0, -LDiv1_2), ),

    ((thick*1.01, thick, -LDiv2_2), ),
    ((thick*1.01, thick-BOX, -LDiv2_2), ),
    ((thick*1.01, clad, -LDiv2_2), ),
    ((thick*1.01, 0.0, -LDiv2_2), ),

    ((thick*1.01, thick, -LDiv3_2), ),
    ((thick*1.01, thick-BOX, -LDiv3_2), ),
    ((thick*1.01, clad, -LDiv3_2), ),
    ((thick*1.01, 0.0, -LDiv3_2), ),

    ((thick*1.01, thick, -L2), ),
    ((thick*1.01, thick-BOX, -L2), ),
    ((thick*1.01, clad, -L2), ),
    ((thick*1.01, 0.0, -L2), ),

########

    ((thick+LxPositive*0.98, thick, thick+L1), ),
    ((thick+LxPositive*0.98, thick-BOX, thick+L1), ),
    ((thick+LxPositive*0.98, clad, thick+L1), ),
    ((thick+LxPositive*0.98, 0.0, thick+L1), ),


    ((thick+LxPositive*0.98, thick, thick+LDiv3_1), ),
    ((thick+LxPositive*0.98, thick-BOX, thick+LDiv3_1), ),
    ((thick+LxPositive*0.98, clad, thick+LDiv3_1), ),
    ((thick+LxPositive*0.98, 0.0, thick+LDiv3_1), ),

    ((thick+LxPositive*0.98, thick, thick+LDiv2_1), ),
    ((thick+LxPositive*0.98, thick-BOX, thick+LDiv2_1), ),
    ((thick+LxPositive*0.98, clad, thick+LDiv2_1), ),
    ((thick+LxPositive*0.98, 0.0, thick+LDiv2_1), ),

    ((thick+LxPositive*0.98, thick, thick+LDiv1_1), ),
    ((thick+LxPositive*0.98, thick-BOX, thick+LDiv1_1), ),
    ((thick+LxPositive*0.98, clad, thick+LDiv1_1), ),
    ((thick+LxPositive*0.98, 0.0, thick+LDiv1_1), ),

    ((thick+LxPositive*0.98, thick, thick), ),
    ((thick+LxPositive*0.98, thick-BOX, thick), ),
    ((thick+LxPositive*0.98, clad, thick), ),
    ((thick+LxPositive*0.98, 0.0, thick), ),

    ((thick+LxPositive*0.98, thick, CR_POS), ),
    ((thick+LxPositive*0.98, thick-BOX, CR_POS), ),
    ((thick+LxPositive*0.98, clad, CR_POS), ),
    ((thick+LxPositive*0.98, 0.0, CR_POS), ),

    ((thick+LxPositive*0.98, thick, 0.0), ),
    ((thick+LxPositive*0.98, thick-BOX, 0.0), ),
    ((thick+LxPositive*0.98, clad, 0.0), ),
    ((thick+LxPositive*0.98, 0.0, 0.0), ),

    ((thick+LxPositive*0.98, thick, -LDiv1_2), ),
    ((thick+LxPositive*0.98, thick-BOX, -LDiv1_2), ),
    ((thick+LxPositive*0.98, clad, -LDiv1_2), ),
    ((thick+LxPositive*0.98, 0.0, -LDiv1_2), ),

    ((thick+LxPositive*0.98, thick, -LDiv2_2), ),
    ((thick+LxPositive*0.98, thick-BOX, -LDiv2_2), ),
    ((thick+LxPositive*0.98, clad, -LDiv2_2), ),
    ((thick+LxPositive*0.98, 0.0, -LDiv2_2), ),

    ((thick+LxPositive*0.98, thick, -LDiv3_2), ),
    ((thick+LxPositive*0.98, thick-BOX, -LDiv3_2), ),
    ((thick+LxPositive*0.98, clad, -LDiv3_2), ),
    ((thick+LxPositive*0.98, 0.0, -LDiv3_2), ),

    ((thick+LxPositive*0.98, thick, -L2), ),
    ((thick+LxPositive*0.98, thick-BOX, -L2), ),
    ((thick+LxPositive*0.98, clad, -L2), ),
    ((thick+LxPositive*0.98, 0.0, -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedLxPos, constraint=FINER)

pickedEdges = Edges.findAt(


    ((-C, thick-BOX*1.01, thick+L1), ),
    ((0, thick-BOX*1.01, thick+L1), ),
    ((BOX, thick-BOX*1.01, thick+L1), ),
    ((thick, thick-BOX*1.01, thick+L1), ),
    ((LxPositive/2, thick-BOX*1.01, thick+L1), ),
    ((thick+LxPositive, thick-BOX*1.01, thick+L1), ),

    ((-C, thick-BOX*1.01, thick+LDiv3_1), ),
    ((0, thick-BOX*1.01, thick+LDiv3_1), ),
    ((BOX, thick-BOX*1.01, thick+LDiv3_1), ),
    ((thick, thick-BOX*1.01, thick+LDiv3_1), ),
    ((LxPositive/2, thick-BOX*1.01, thick+LDiv3_1), ),
    ((thick+LxPositive, thick-BOX*1.01, thick+LDiv3_1), ),

    ((-C, thick-BOX*1.01, thick+LDiv2_1), ),
    ((0, thick-BOX*1.01, thick+LDiv2_1), ),
    ((BOX, thick-BOX*1.01, thick+LDiv2_1), ),
    ((thick, thick-BOX*1.01, thick+LDiv2_1), ),
    ((LxPositive/2, thick-BOX*1.01, thick+LDiv2_1), ),
    ((thick+LxPositive, thick-BOX*1.01, thick+LDiv2_1), ),

    ((-C, thick-BOX*1.01, thick+LDiv1_1), ),
    ((0, thick-BOX*1.01, thick+LDiv1_1), ),
    ((BOX, thick-BOX*1.01, thick+LDiv1_1), ),
    ((thick, thick-BOX*1.01, thick+LDiv1_1), ),
    ((LxPositive/2, thick-BOX*1.01, thick+LDiv1_1), ),
    ((thick+LxPositive, thick-BOX*1.01, thick+LDiv1_1), ),

    ((-C, thick-BOX*1.01, thick), ),
    ((0, thick-BOX*1.01, thick), ),
    ((BOX, thick-BOX*1.01, thick), ),
    ((thick, thick-BOX*1.01, thick), ),
    ((LxPositive/2, thick-BOX*1.01, thick), ),
    ((thick+LxPositive, thick-BOX*1.01, thick), ),

    ((-C, thick-BOX*1.01, CR_POS), ),
    ((0, thick-BOX*1.01, CR_POS), ),
    ((BOX, thick-BOX*1.01, CR_POS), ),
    ((thick, thick-BOX*1.01, CR_POS), ),
    ((LxPositive/2, thick-BOX*1.01, CR_POS), ),
    ((thick+LxPositive, thick-BOX*1.01, CR_POS), ),

    ((-C, thick-BOX*1.01, 0.), ),
    ((0, thick-BOX*1.01, 0.), ),
    ((BOX, thick-BOX*1.01, 0.), ),
    ((thick, thick-BOX*1.01, 0.), ),
    ((LxPositive/2, thick-BOX*1.01, 0.), ),
    ((thick+LxPositive, thick-BOX*1.01, 0.), ),

    ((-C, thick-BOX*1.01, -LDiv1_2), ),
    ((0, thick-BOX*1.01, -LDiv1_2), ),
    ((BOX, thick-BOX*1.01, -LDiv1_2), ),
    ((thick, thick-BOX*1.01, -LDiv1_2), ),
    ((LxPositive/2, thick-BOX*1.01, -LDiv1_2), ),
    ((thick+LxPositive, thick-BOX*1.01, -LDiv1_2), ),

    ((-C, thick-BOX*1.01, -LDiv2_2), ),
    ((0, thick-BOX*1.01, -LDiv2_2), ),
    ((BOX, thick-BOX*1.01, -LDiv2_2), ),
    ((thick, thick-BOX*1.01, -LDiv2_2), ),
    ((LxPositive/2, thick-BOX*1.01, -LDiv2_2), ),
    ((thick+LxPositive, thick-BOX*1.01, -LDiv2_2), ),

    ((-C, thick-BOX*1.01, -LDiv3_2), ),
    ((0, thick-BOX*1.01, -LDiv3_2), ),
    ((BOX, thick-BOX*1.01, -LDiv3_2), ),
    ((thick, thick-BOX*1.01, -LDiv3_2), ),
    ((LxPositive/2, thick-BOX*1.01, -LDiv3_2), ),
    ((thick+LxPositive, thick-BOX*1.01, -LDiv3_2), ),

    ((-C, thick-BOX*1.01, -L2), ),
    ((0, thick-BOX*1.01, -L2), ),
    ((BOX, thick-BOX*1.01, -L2), ),
    ((thick, thick-BOX*1.01, -L2), ),
    ((LxPositive/2, thick-BOX*1.01, -L2), ),
    ((thick+LxPositive, thick-BOX*1.01, -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedthickV, constraint=FINER)


###%%%

pickedEdges = Edges.findAt(


    ((-C, thick-a*0.2, thick+L1), ),
    ((0, thick-a*0.2, thick+L1), ),
    ((BOX, thick-a*0.2, thick+L1), ),
    ((thick, thick-a*0.2, thick+L1), ),
    ((LxPositive/2, thick-a*0.2, thick+L1), ),
    ((thick+LxPositive, thick-a*0.2, thick+L1), ),

    ((-C, thick-a*0.2, thick+LDiv3_1), ),
    ((0, thick-a*0.2, thick+LDiv3_1), ),
    ((BOX, thick-a*0.2, thick+LDiv3_1), ),
    ((thick, thick-a*0.2, thick+LDiv3_1), ),
    ((LxPositive/2, thick-a*0.2, thick+LDiv3_1), ),
    ((thick+LxPositive, thick-a*0.2, thick+LDiv3_1), ),

    ((-C, thick-a*0.2, thick+LDiv2_1), ),
    ((0, thick-a*0.2, thick+LDiv2_1), ),
    ((BOX, thick-a*0.2, thick+LDiv2_1), ),
    ((thick, thick-a*0.2, thick+LDiv2_1), ),
    ((LxPositive/2, thick-a*0.2, thick+LDiv2_1), ),
    ((thick+LxPositive, thick-a*0.2, thick+LDiv2_1), ),

    ((-C, thick-a*0.2, thick+LDiv1_1), ),
    ((0, thick-a*0.2, thick+LDiv1_1), ),
    ((BOX, thick-a*0.2, thick+LDiv1_1), ),
    ((thick, thick-a*0.2, thick+LDiv1_1), ),
    ((LxPositive/2, thick-a*0.2, thick+LDiv1_1), ),
    ((thick+LxPositive, thick-a*0.2, thick+LDiv1_1), ),

    ((-C, thick-a*0.2, thick), ),
    ((0, thick-a*0.2, thick), ),
    ((BOX, thick-a*0.2, thick), ),
    ((thick, thick-a*0.2, thick), ),
    ((LxPositive/2, thick-a*0.2, thick), ),
    ((thick+LxPositive, thick-a*0.2, thick), ),

    ((-C, thick-a*0.2, CR_POS), ),
    ((0, thick-a*0.2, CR_POS), ),
    ((BOX, thick-a*0.2, CR_POS), ),
    ((thick, thick-a*0.2, CR_POS), ),
    ((LxPositive/2, thick-a*0.2, CR_POS), ),
    ((thick+LxPositive, thick-a*0.2, CR_POS), ),

    ((-C, thick-a*0.2, 0.), ),
    ((0, thick-a*0.2, 0.), ),
    ((BOX, thick-a*0.2, 0.), ),
    ((thick, thick-a*0.2, 0.), ),
    ((LxPositive/2, thick-a*0.2, 0.), ),
    ((thick+LxPositive, thick-a*0.2, 0.), ),

    ((-C, thick-a*0.2, -LDiv1_2), ),
    ((0, thick-a*0.2, -LDiv1_2), ),
    ((BOX, thick-a*0.2, -LDiv1_2), ),
    ((thick, thick-a*0.2, -LDiv1_2), ),
    ((LxPositive/2, thick-a*0.2, -LDiv1_2), ),
    ((thick+LxPositive, thick-a*0.2, -LDiv1_2), ),

    ((-C, thick-a*0.2, -LDiv2_2), ),
    ((0, thick-a*0.2, -LDiv2_2), ),
    ((BOX, thick-a*0.2, -LDiv2_2), ),
    ((thick, thick-a*0.2, -LDiv2_2), ),
    ((LxPositive/2, thick-a*0.2, -LDiv2_2), ),
    ((thick+LxPositive, thick-a*0.2, -LDiv2_2), ),

    ((-C, thick-a*0.2, -LDiv3_2), ),
    ((0, thick-a*0.2, -LDiv3_2), ),
    ((BOX, thick-a*0.2, -LDiv3_2), ),
    ((thick, thick-a*0.2, -LDiv3_2), ),
    ((LxPositive/2, thick-a*0.2, -LDiv3_2), ),
    ((thick+LxPositive, thick-a*0.2, -LDiv3_2), ),

    ((-C, thick-a*0.2, -L2), ),
    ((0, thick-a*0.2, -L2), ),
    ((BOX, thick-a*0.2, -L2), ),
    ((thick, thick-a*0.2, -L2), ),
    ((LxPositive/2, thick-a*0.2, -L2), ),
    ((thick+LxPositive, thick-a*0.2, -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedBOX, constraint=FINER)
#####$$$$

pickedEdges = Edges.findAt(


    ((-C, clad*0.3, thick+L1), ),
    ((0, clad*0.3, thick+L1), ),
    ((BOX, clad*0.3, thick+L1), ),
    ((thick, clad*0.3, thick+L1), ),
    ((LxPositive/2, clad*0.3, thick+L1), ),
    ((thick+LxPositive, clad*0.3, thick+L1), ),

    ((-C, clad*0.3, thick+LDiv3_1), ),
    ((0, clad*0.3, thick+LDiv3_1), ),
    ((BOX, clad*0.3, thick+LDiv3_1), ),
    ((thick, clad*0.3, thick+LDiv3_1), ),
    ((LxPositive/2, clad*0.3, thick+LDiv3_1), ),
    ((thick+LxPositive, clad*0.3, thick+LDiv3_1), ),

    ((-C, clad*0.3, thick+LDiv2_1), ),
    ((0, clad*0.3, thick+LDiv2_1), ),
    ((BOX, clad*0.3, thick+LDiv2_1), ),
    ((thick, clad*0.3, thick+LDiv2_1), ),
    ((LxPositive/2, clad*0.3, thick+LDiv2_1), ),
    ((thick+LxPositive, clad*0.3, thick+LDiv2_1), ),

    ((-C, clad*0.3, thick+LDiv1_1), ),
    ((0, clad*0.3, thick+LDiv1_1), ),
    ((BOX, clad*0.3, thick+LDiv1_1), ),
    ((thick, clad*0.3, thick+LDiv1_1), ),
    ((LxPositive/2, clad*0.3, thick+LDiv1_1), ),
    ((thick+LxPositive, clad*0.3, thick+LDiv1_1), ),

    ((-C, clad*0.3, thick), ),
    ((0, clad*0.3, thick), ),
    ((BOX, clad*0.3, thick), ),
    ((thick, clad*0.3, thick), ),
    ((LxPositive/2, clad*0.3, thick), ),
    ((thick+LxPositive, clad*0.3, thick), ),

    ((-C, clad*0.3, CR_POS), ),
    ((0, clad*0.3, CR_POS), ),
    ((BOX, clad*0.3, CR_POS), ),
    ((thick, clad*0.3, CR_POS), ),
    ((LxPositive/2, clad*0.3, CR_POS), ),
    ((thick+LxPositive, clad*0.3, CR_POS), ),

    ((-C, clad*0.3, 0.), ),
    ((0, clad*0.3, 0.), ),
    ((BOX, clad*0.3, 0.), ),
    ((thick, clad*0.3, 0.), ),
    ((LxPositive/2, clad*0.3, 0.), ),
    ((thick+LxPositive, clad*0.3, 0.), ),

    ((-C, clad*0.3, -LDiv1_2), ),
    ((0, clad*0.3, -LDiv1_2), ),
    ((BOX, clad*0.3, -LDiv1_2), ),
    ((thick, clad*0.3, -LDiv1_2), ),
    ((LxPositive/2, clad*0.3, -LDiv1_2), ),
    ((thick+LxPositive, clad*0.3, -LDiv1_2), ),

    ((-C, clad*0.3, -LDiv2_2), ),
    ((0, clad*0.3, -LDiv2_2), ),
    ((BOX, clad*0.3, -LDiv2_2), ),
    ((thick, clad*0.3, -LDiv2_2), ),
    ((LxPositive/2, clad*0.3, -LDiv2_2), ),
    ((thick+LxPositive, clad*0.3, -LDiv2_2), ),

    ((-C, clad*0.3, -LDiv3_2), ),
    ((0, clad*0.3, -LDiv3_2), ),
    ((BOX, clad*0.3, -LDiv3_2), ),
    ((thick, clad*0.3, -LDiv3_2), ),
    ((LxPositive/2, clad*0.3, -LDiv3_2), ),
    ((thick+LxPositive, clad*0.3, -LDiv3_2), ),

    ((-C, clad*0.3, -L2), ),
    ((0, clad*0.3, -L2), ),
    ((BOX, clad*0.3, -L2), ),
    ((thick, clad*0.3, -L2), ),
    ((LxPositive/2, clad*0.3, -L2), ),
    ((thick+LxPositive, clad*0.3, -L2), ),


)
PipePart.seedEdgeByNumber(edges=pickedEdges, number=nseedclad, constraint=FINER)




elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,  secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)

cells = Cells.getByBoundingBox(xMin=-20000, xMax=20000,yMin=-20000,yMax=20000,zMin=-20000,zMax=20000)
pickedRegions =(cells,)
PipePart.setElementType(regions=pickedRegions, elemTypes=(elemType1,elemType2,))
PipePart.generateMesh()

####################################################################
###############################################################################
### STEP

PipeModel.StaticStep(name='BendingStep', previous='Initial',   description='Bending', 
	timePeriod=30.0, maxNumInc=1000, initialInc=1.0,   minInc=0.00001, maxInc=1.0)
##############################################################################
###############################################################################
# ###ASSEMBLY
PipeAssambly = PipeModel.rootAssembly
PipeAssambly.DatumCsysByDefault(CARTESIAN)
PipePart = PipeModel.parts['CrBlock1']
PipeAssambly.Instance(name='Pipe_Instance', part=PipePart, dependent=ON)

n1 = PipeAssambly.instances['Pipe_Instance'].nodes
nodes1 = n1[0:10000000000]
PipeAssambly.Set(nodes=nodes1, name='ALLN')


PipePart = PipeModel.parts['CrBlock1']
Edges = PipePart.edges

Pipessembly = PipeModel.rootAssembly
v11 = Pipessembly.instances['Pipe_Instance'].vertices
Pipessembly.ReferencePoint(point=v11.findAt(coordinates=(-C, 0.0, -L2)))


# Movement of Instance- Assembly

PipeAssambly.translate(instanceList=('Pipe_Instance', ), vector=(C, 0.0, 0))
PipeAssembly = PipeModel.rootAssembly
PipeAssembly.rotate(instanceList=('Pipe_Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(   0.0, 0.0, 1.2), angle=90.0)
PipeModel.setValues(noPartsInputFile=ON)

# ##### MPC

PipeAssembly = PipeModel.rootAssembly
RefPoints = PipeAssembly.referencePoints
refPoints1=(RefPoints[5], )
region1=PipeAssembly.Set(referencePoints=refPoints1, name='Ref_Point_Set')
Faces = PipeAssembly.instances['Pipe_Instance'].faces
Edges = PipeAssembly.instances['Pipe_Instance'].edges

faces1 = Faces.findAt(
	((-clad*0.2, C*0.1, -L2), ), 
	((-clad*1.2, C*0.1, -L2), ),
	((-thick+a*1.01, C*0.1, -L2), ), 
 	((-thick+a*0.9, C*0.1, -L2), ), 

	((-clad*0.2, C+a*0.1, -L2), ), 
	((-clad*1.2, C+a*0.1, -L2), ),
	((-thick+a*1.01, C+a*0.1, -L2), ), 
 	((-thick+a*0.9, C+a*0.1, -L2), ), 

	((-clad*0.2, C+thick*0.9, -L2), ), 
	((-clad*1.2, C+thick*0.9, -L2), ),
	((-thick+a*1.01, C+thick*0.9, -L2), ), 
 	((-thick+a*0.9, C+thick*0.9, -L2), ), 

	((-clad*0.2, C+thick*1.1, -L2), ), 
	((-clad*1.2, C+thick*1.1, -L2), ),
	((-thick+a*1.01, C+thick*1.1, -L2), ), 
 	((-thick+a*0.9, C+thick*1.1, -L2), ), 

	((-clad*0.2, C+thick+LxPositive*0.9, -L2), ), 
	((-clad*1.2, C+thick+LxPositive*0.9, -L2), ),
	((-thick+a*1.01, C+thick+LxPositive*0.9, -L2), ), 
 	((-thick+a*0.9, C+thick+LxPositive*0.9, -L2), ), 

	)

region2=PipeAssembly.Set( faces=faces1, name='MPC_Faces_Set')
PipeModel.MultipointConstraint(name='MPC_Ref_Point', controlPoint=region1,surface=region2, mpcType=BEAM_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)

##### Crack definition

PipeAssembly = PipeModel.rootAssembly
Edges = PipeAssembly.instances['Pipe_Instance'].edges
edges1 = Edges.findAt(
	((-(thick - a), C*0.1, CR_POS), ), 
	((-thick + a*cos(pi/4), C+a*cos(pi/4), CR_POS), )
	)
crackFront = regionToolset.Region(edges=edges1)
crackTip = regionToolset.Region(edges=edges1)
PipeAssembly.engineeringFeatures.ContourIntegral(name='Crack_Interface', symmetric=ON, crackFront=crackFront, crackTip=crackTip, 
    extensionDirectionMethod=CRACK_NORMAL, crackNormal=((0.0, 0.0, 0.0), (0.0, 0.0, -1.0)), midNodePosition=0.25, collapsedElementAtTip=NONE, crackTipSense=FORWARD)
###################################################

if External == 'YES':
    PipeAssembly = PipeModel.rootAssembly
    PipeAssembly.rotate(instanceList=('Pipe_Instance', ), axisPoint=(0.0, 0.0, 0.0),  axisDirection=(0.0, 8.24, 0.0), angle=180.0)
    Faces = PipeAssembly.instances['Pipe_Instance'].faces
    Faces = PipeAssembly.instances['Pipe_Instance'].faces
    faces1 = Faces.findAt(
        ((clad*0.2, 0.0, -thick-L1*0.99), ), 
        ((clad*1.2, 0.0, -thick -L1*0.99), ), 
        ((thick-a*1.01, 0.0, -thick - L1*0.99), ), 
        ((thick-a*0.99, 0.0, -thick - L1*0.99), ), 

        ((clad*0.2, 0.0, -thick-LDiv3_1*0.99), ), 
        ((clad*1.2, 0.0, -thick-LDiv3_1*0.99), ), 
        ((thick-a*1.01, 0.0, -thick-LDiv3_1*0.99), ), 
        ((thick-a*0.99, 0.0, -thick-LDiv3_1*0.99), ), 

        ((clad*0.2, 0.0, -thick-LDiv2_1*0.9), ), 
        ((clad*1.2, 0.0, -thick-LDiv2_1*0.9), ), 
        ((thick-a*1.01, 0.0, -thick-LDiv2_1*0.9), ), 
        ((thick-a*0.99, 0.0, -thick-LDiv2_1*0.9), ), 

        ((clad*0.2, 0.0, -thick-LDiv1_1*0.9), ), 
        ((clad*1.2, 0.0, -thick-LDiv1_1*0.9), ), 
        ((thick-a*1.01, 0.0, -thick-LDiv1_1*0.9), ), 
        ((thick-a*0.99, 0.0, -thick-LDiv1_1*0.9), ), 

        ((clad*0.2, 0.0, -thick*0.9), ), 
        ((clad*1.2, 0.0, -thick*0.9), ), 
        ((thick-a*1.01, 0.0, -thick*0.9), ), 
        ((thick-a*0.99, 0.0, -thick*0.9), ), 

        ((clad*0.2, 0.0, -CR_POS-a*0.55), ), 
        ((clad*1.2, 0.0, -CR_POS-a*0.55), ), 
        ((thick-a*1.01, 0.0, -CR_POS-a*0.55), ), 
         ((thick-a*0.99, 0.0, -CR_POS-a*0.55), ), 

        ((clad*0.2, 0.0, -CR_POS+a*0.55), ), 
        ((clad*1.2, 0.0, -CR_POS+a*0.55), ), 
        ((thick-a*1.01, 0.0, -CR_POS+a*0.55), ), 
        ((thick-a*0.99, 0.0, -CR_POS+a*0.55), ), 

        ((clad*0.2, 0.0, -a*0.02), ), 
        ((clad*1.2, 0.0, -a*0.02), ), 
        ((thick-a*1.01, 0.0, -a*0.02), ), 
        ((thick-a*0.99, 0.0, -a*0.02), ), 


        ((clad*0.2, 0.0, LDiv1_2*0.9), ), 
        ((clad*1.2, 0.0, LDiv1_2*0.9), ), 
        ((thick-a*1.01, 0.0, LDiv1_2*0.9), ), 
        ((thick-a*0.99, 0.0, LDiv1_2*0.9), ), 

        ((clad*0.2, 0.0, LDiv2_2*0.9), ), 
        ((clad*1.2, 0.0, LDiv2_2*0.9), ), 
        ((thick-a*1.01, 0.0, LDiv2_2*0.9), ), 
        ((thick-a*0.99, 0.0, LDiv2_2*0.9), ), 

        ((clad*0.2, 0.0, LDiv3_2*0.9), ), 
        ((clad*1.2, 0.0, LDiv3_2*0.9), ), 
        ((thick-a*1.01, 0.0, LDiv3_2*0.9), ), 
        ((thick-a*0.99, 0.0, LDiv3_2*0.9), ), 

        ((clad*0.2, 0.0, L2*0.9), ), 
        ((clad*1.2, 0.0, L2*0.9), ), 
        ((thick-a*1.01, 0.0, L2*0.9), ), 
        ((thick-a*0.99, 0.0, L2*0.9), ), 

###################$4444
        ((clad*0.2, C+thick+LxPositive, -thick-L1*0.99), ), 
        ((clad*1.2, C+thick+LxPositive, -thick -L1*0.99), ), 
        ((thick-a*1.01, C+thick+LxPositive, -thick - L1*0.99), ), 
        ((thick-a*0.99, C+thick+LxPositive, -thick - L1*0.99), ), 

        ((clad*0.2, C+thick+LxPositive, -thick-LDiv3_1*0.99), ), 
        ((clad*1.2, C+thick+LxPositive, -thick-LDiv3_1*0.99), ), 
        ((thick-a*1.01, C+thick+LxPositive, -thick-LDiv3_1*0.99), ), 
        ((thick-a*0.99, C+thick+LxPositive, -thick-LDiv3_1*0.99), ), 

        ((clad*0.2, C+thick+LxPositive, -thick-LDiv2_1*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, -thick-LDiv2_1*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, -thick-LDiv2_1*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, -thick-LDiv2_1*0.9), ), 

        ((clad*0.2, C+thick+LxPositive, -thick-LDiv1_1*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, -thick-LDiv1_1*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, -thick-LDiv1_1*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, -thick-LDiv1_1*0.9), ), 

        ((clad*0.2, C+thick+LxPositive, -thick*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, -thick*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, -thick*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, -thick*0.9), ), 

        ((clad*0.2, C+thick+LxPositive, -CR_POS-a*0.55), ), 
        ((clad*1.2, C+thick+LxPositive, -CR_POS-a*0.55), ), 
        ((thick-a*1.01, C+thick+LxPositive, -CR_POS-a*0.55), ), 
        ((thick-a*0.99, C+thick+LxPositive, -CR_POS-a*0.55), ), 

        ((clad*0.2, C+thick+LxPositive, -CR_POS+a*0.55), ), 
        ((clad*1.2, C+thick+LxPositive, -CR_POS+a*0.55), ), 
        ((thick-a*1.01, C+thick+LxPositive, -CR_POS+a*0.55), ), 
        ((thick-a*0.99, C+thick+LxPositive, -CR_POS+a*0.55), ), 

        ((clad*0.2, C+thick+LxPositive, -a*0.02), ), 
        ((clad*1.2, C+thick+LxPositive, -a*0.02), ), 
        ((thick-a*1.01, C+thick+LxPositive, -a*0.02), ), 
        ((thick-a*0.99, C+thick+LxPositive, -a*0.02), ), 


        ((clad*0.2, C+thick+LxPositive, LDiv1_2*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, LDiv1_2*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, LDiv1_2*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, LDiv1_2*0.9), ), 

        ((clad*0.2, C+thick+LxPositive, LDiv2_2*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, LDiv2_2*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, LDiv2_2*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, LDiv2_2*0.9), ), 

        ((clad*0.2, C+thick+LxPositive, LDiv3_2*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, LDiv3_2*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, LDiv3_2*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, LDiv3_2*0.9), ), 


        ((clad*0.2, C+thick+LxPositive, L2*0.9), ), 
        ((clad*1.2, C+thick+LxPositive, L2*0.9), ), 
        ((thick-a*1.01, C+thick+LxPositive, L2*0.9), ), 
        ((thick-a*0.99, C+thick+LxPositive, L2*0.9), ), 


    )
    region = PipeAssembly.Set(faces=faces1, name='Set_Faces_XSYM')

    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='BendingStep', region=region, localCsys=None)

    PipeAssembly = PipeModel.rootAssembly
    Faces = PipeAssembly.instances['Pipe_Instance'].faces
    faces1 = Faces.findAt(
        ((clad*0.2, C*0.1, -CR_POS), ), 
        ((clad*1.1, C*0.1, -CR_POS), ), 
         ((thick-a-A1*0.9, C*0.1, -CR_POS), ), 
         ((thick-a-A1*1.1, C*0.1, -CR_POS), ),
         ((thick-a-A3*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A4*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A5*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A6*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A7*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A8*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A9*0.9, C*0.1, -CR_POS), ),
         ((thick-a-A9*1.1, C*0.1, -CR_POS), ),


         ((clad*0.2, C*1.1, -CR_POS), ), 
         ((clad*1.1, C*1.1, -CR_POS), ), 
         ((thick-a-A1*0.9, C*1.01, -CR_POS), ), 
         ((thick-a-A1*1.1, C*1.01, -CR_POS), ),
         ((thick-a-A3*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A4*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A5*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A6*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A7*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A8*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A9*0.9, C*1.01, -CR_POS), ),
         ((thick-a-A9*1.1, C*1.01, -CR_POS), ),

         
        ((clad*0.2, C+thick*0.9, -CR_POS), ), 
        ((clad*1.1, C+thick*0.9, -CR_POS), ), 
        ((thick*0.9, C+thick*0.9, -CR_POS), ), 

        ((clad*0.2, C+thick*1.1, -CR_POS), ), 
        ((clad*1.1, C+thick*1.1, -CR_POS), ), 
        ((thick*0.9, C+thick*1.1, -CR_POS), ), 

        ((clad*0.2, C+thick+LxPositive*0.9, -CR_POS), ), 
        ((clad*1.1, C+thick+LxPositive*0.9, -CR_POS), ), 
        ((thick*0.9, C+thick+LxPositive*0.9, -CR_POS), ), 

      
        )
    region = PipeAssembly.Set(faces=faces1, name='Set_Faces_ZSYM')
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='BendingStep', region=region, localCsys=None)

else:

    PipeAssembly = PipeModel.rootAssembly

    Faces = PipeAssembly.instances['Pipe_Instance'].faces
    faces1 = Faces.findAt(
        ((-clad*0.2, 0.0, -L2*0.9), ), 
        ((-clad*1.2, 0.0, -L2*0.9), ), 
        ((-thick+a*1.01, 0.0, -L2*0.9), ), 
        ((-thick+a*0.99, 0.0, -L2*0.9), ), 

        ((-clad*0.2, 0.0, -LDiv3_2*0.9), ), 
        ((-clad*1.2, 0.0, -LDiv3_2*0.9), ), 
        ((-thick+a*1.01, 0.0, -LDiv3_2*0.9), ), 
        ((-thick+a*0.99, 0.0, -LDiv3_2*0.9), ), 

        ((-clad*0.2, 0.0, -LDiv2_2*0.9), ), 
        ((-clad*1.2, 0.0, -LDiv2_2*0.9), ), 
        ((-thick+a*1.01, 0.0, -LDiv2_2*0.9), ), 
        ((-thick+a*0.99, 0.0, -LDiv2_2*0.9), ), 

        ((-clad*0.2, 0.0, -LDiv1_2*0.9), ), 
        ((-clad*1.2, 0.0, -LDiv1_2*0.9), ), 
        ((-thick+a*1.01, 0.0, -LDiv1_2*0.9), ), 
        ((-thick+a*0.99, 0.0, -LDiv1_2*0.9), ), 

        ((-clad*0.2, 0.0, thick*0.1), ), 
        ((-clad*1.2, 0.0, thick*0.1), ), 
        ((-thick+a*1.01, 0.0, thick*0.1), ), 
        ((-thick+a*0.99, 0.0, thick*0.1), ), 

        ((-clad*0.2, 0.0, CR_POS-a*0.55), ), 
        ((-clad*1.2, 0.0, CR_POS-a*0.55), ), 
        ((-thick+a*1.01, 0.0, CR_POS-a*0.55), ), 
        ((-thick+a*0.99, 0.0, CR_POS-a*0.55), ), 

        ((-clad*0.2, 0.0, CR_POS+a*0.55), ), 
        ((-clad*1.2, 0.0, CR_POS+a*0.55), ), 
        ((-thick+a*1.01, 0.0, CR_POS+a*0.55), ), 
        ((-thick+a*0.99, 0.0, CR_POS+a*0.55), ), 

        ((-clad*0.2, 0.0, thick*0.9), ), 
        ((-clad*1.2, 0.0, thick*0.9), ), 
        ((-thick+a*1.01, 0.0, thick*0.9), ), 
        ((-thick+a*0.99, 0.0, thick*0.9), ), 

        ((-clad*0.2, 0.0, thick+LDiv1_1*0.9), ), 
        ((-clad*1.2, 0.0, thick+LDiv1_1*0.9), ), 
        ((-thick+a*1.01, 0.0, thick+LDiv1_1*0.9), ), 
        ((-thick+a*0.99, 0.0, thick+LDiv1_1*0.9), ), 

        ((-clad*0.2, 0.0, thick+LDiv2_1*0.9), ), 
        ((-clad*1.2, 0.0, thick+LDiv2_1*0.9), ), 
        ((-thick+a*1.01, 0.0, thick+LDiv2_1*0.9), ), 
        ((-thick+a*0.99, 0.0, thick+LDiv2_1*0.9), ), 

        ((-clad*0.2, 0.0, thick+LDiv3_1*0.9), ), 
        ((-clad*1.2, 0.0, thick+LDiv3_1*0.9), ), 
        ((-thick+a*1.01, 0.0, thick+LDiv3_1*0.9), ), 
        ((-thick+a*0.99, 0.0, thick+LDiv3_1*0.9), ), 


        ((-clad*0.2, 0.0, thick+L1*0.9), ), 
        ((-clad*1.2, 0.0, thick+L1*0.9), ), 
        ((-thick+a*1.01, 0.0, thick+L1*0.9), ), 
        ((-thick+a*0.99, 0.0, thick+L1*0.9), ), 

    ###############$4444
        ((-clad*0.2, C+thick+LxPositive, -L2*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, -L2*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, -L2*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, -L2*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, -LDiv3_2*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, -LDiv3_2*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, -LDiv3_2*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, -LDiv3_2*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, -LDiv2_2*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, -LDiv2_2*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, -LDiv2_2*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, -LDiv2_2*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, -LDiv1_2*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, -LDiv1_2*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, -LDiv1_2*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, -LDiv1_2*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, thick*0.1), ), 
        ((-clad*1.2, C+thick+LxPositive, thick*0.1), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick*0.1), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick*0.1), ), 

        ((-clad*0.2, C+thick+LxPositive, CR_POS-a*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, CR_POS-a*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, CR_POS-a*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, CR_POS-a*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, CR_POS+a*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, CR_POS+a*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, CR_POS+a*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, CR_POS+a*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, thick*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, thick*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, thick+LDiv1_1*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, thick+LDiv1_1*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick+LDiv1_1*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick+LDiv1_1*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, thick+LDiv2_1*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, thick+LDiv2_1*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick+LDiv2_1*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick+LDiv2_1*0.9), ), 

        ((-clad*0.2, C+thick+LxPositive, thick+LDiv3_1*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, thick+LDiv3_1*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick+LDiv3_1*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick+LDiv3_1*0.9), ), 


        ((-clad*0.2, C+thick+LxPositive, thick+L1*0.9), ), 
        ((-clad*1.2, C+thick+LxPositive, thick+L1*0.9), ), 
        ((-thick+a*1.01, C+thick+LxPositive, thick+L1*0.9), ), 
        ((-thick+a*0.99, C+thick+LxPositive, thick+L1*0.9), ), 
    

        )
    region = PipeAssembly.Set(faces=faces1, name='Set_Faces_XSYM')

    PipeModel.XsymmBC(name='X_SYM_BC', createStepName='BendingStep', region=region, localCsys=None)



    PipeAssembly = PipeModel.rootAssembly
    Faces = PipeAssembly.instances['Pipe_Instance'].faces
    faces1 = Faces.findAt(
        ((-clad*0.2, C*0.1, CR_POS), ), 
        ((-clad*1.1, C*0.1, CR_POS), ), 
        ((-thick+a+A1*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A2*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A3*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A4*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A5*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A6*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A7*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A8*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A9*0.9, C*0.1, CR_POS), ), 
        ((-thick+a+A9*1.01, C*0.1, CR_POS), ), 

        ((-clad*0.2, C*1.01, CR_POS), ), 
        ((-clad*1.1, C*1.01, CR_POS), ), 
        ((-thick+a+A1*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A2*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A3*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A4*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A5*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A6*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A7*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A8*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A9*0.9, C*1.01, CR_POS), ), 
        ((-thick+a+A9*1.01, C*1.01, CR_POS), ), 

        ((-clad*0.2, C+thick*0.9, CR_POS), ), 
        ((-clad*1.1, C+thick*0.9, CR_POS), ), 
        ((-thick*0.9, C+thick*0.9, CR_POS), ),

        
        ((-clad*0.2, C+thick*1.01, CR_POS), ), 
        ((-clad*1.1, C+thick*1.01, CR_POS), ), 
        ((-thick*0.9, C+thick*1.01, CR_POS), ),

         ((-clad*0.2, C+thick+LxPositive*0.91, CR_POS), ), 
        ((-clad*1.1, C+thick+LxPositive*0.91, CR_POS), ), 
        ((-thick*0.9, C+thick+LxPositive*0.91, CR_POS), ), 


        )
    region = PipeAssembly.Set(faces=faces1, name='Set_Faces_ZSYM')
    PipeModel.ZsymmBC(name='Z_SYM_BC', createStepName='BendingStep', region=region, localCsys=None)


PipeAssembly = PipeModel.rootAssembly
region = PipeAssembly.sets['Ref_Point_Set']
PipeModel.DisplacementBC(name='Y_Fixed_RefP', createStepName='BendingStep', 
    region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

PipeAssembly = PipeModel.rootAssembly
region = PipeAssembly.sets['Ref_Point_Set']
PipeModel.DisplacementBC(name='Bending', createStepName='BendingStep', 
    region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=-Applied_Bending, ur2=UNSET, 
    ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
    fieldName='', localCsys=None)

######### OutPut J-Integral
PipeModel.HistoryOutputRequest(name='J_integral', 
    createStepName='BendingStep', contourIntegral='Crack_Interface', 
    sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=30)

###OUTPUT RESULTADOS
#OUTPUT DOS DESLOCAMENTOS VERTICAIS DOS NODES DA FACE DA TRINCA
regionDef=PipeModel.rootAssembly.sets['Set_Faces_ZSYM']
PipeModel.HistoryOutputRequest(name='face_displ', 
    createStepName='BendingStep', variables=('U2','RF3', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

#OUTPUT DA ROTACAO
regionDef=PipeModel.rootAssembly.sets['Ref_Point_Set']
PipeModel.HistoryOutputRequest(name='rotation', 
    createStepName='BendingStep', variables=('UR','RM', ), timeInterval=1.0, 
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

#OUTPUT DAS TENSOES E DEFORMACOES
PipeModel.FieldOutputRequest(name='SS_results', createStepName='BendingStep', 
    variables=('MISES', 'E','PEEQ','U'),timeInterval=1.0)

del PipeModel.fieldOutputRequests['F-Output-1']
del PipeModel.historyOutputRequests['H-Output-1']


##########################################


radvar = str(radius)
vari = str(180/(radius*pi))
PipeModel.keywordBlock.synchVersions(storeNodesAndElements=False)
PipeModel.keywordBlock.insert(35, 
   '\n*NMAP, NSET=ALLN, TYPE=RECTANGULAR\n%s,0.,0.,999999.,0.,0.\n%s,1.,0.\n1.,1.,1.\n*NMAP, NSET=ALLN, TYPE=CYLINDRICAL\n0.,0.,0.,0.,0.,1.\n0.,1.,0.\n1.,%s,1.\n*NMAP, NSET=ALLN, TYPE=RECTANGULAR\n0.,0.,0., 1.,0.,0.\n1.,1.,0.' %(radvar,radvar,vari) )

if External == 'YES':
    ModelType = 'EXTERNAL'
else:
    ModelType = 'INTERNAL'

mdb.Job(name='Interface'+ModelType+'_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(Theta_Pi*100))+'_AT'+str(int(round(10*a/thick))),
      model='EX_MOD', description='', type=ANALYSIS, 
      atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
      memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
      explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
      modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
      scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=NCPUs, 
      activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=NCPUs)
PipeModel.setValues(noPartsInputFile=ON)
mdb.jobs['Interface'+ModelType+ '_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(Theta_Pi*100))+'_AT'+str(int(round(10*a/thick)))].writeInput(consistencyChecking=OFF)

#mdb.jobs['Interface'+ModelType+ '_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(Theta_Pi*100))+'_AT'+str(int(round(10*a/thick)))].submit(consistencyChecking=OFF)
#mdb.jobs['Interface'+ModelType+ '_N'+str(int(N_Base))+'_DT'+str(int(round(radius*2/thick)))+'_TP'+str(int(Theta_Pi*100))+'_AT'+str(int(round(10*a/thick)))].waitForCompletion()



