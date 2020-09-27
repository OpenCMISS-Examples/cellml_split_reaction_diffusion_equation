#!/usr/bin/env python

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron



#------------------------------------------------------------------------------------------------
# INPUTS
#------------------------------------------------------------------------------------------------

# Geometric Inputs
numberGlobalXElements = 10
length = 100.0

# Reaction Diffusion Inputs
CellMLFile = "constant_rate.xml"         # CellML Model File Name  
Diffusivity = 0.5                                            # Diffusivity
store_coeff = 1.0                                   # Defining Storage Coefficient
Ca_init = 0.0
cellMLModelsField_Init = 1

# Simulation Runtime Inputs
startT = 0.0
endT   = 0.5
Tstep  = 0.01
ODE_TIME_STEP = 0.0000001
#outputfreq = 10


#------------------------------------------------------------------------------------------------
# USER NUMBERS
#------------------------------------------------------------------------------------------------

# Defining User Numbers
coordinateSystemUserNumber          = 10
regionUserNumber                    = 20
basisUserNumber                     = 30
generatedMeshUserNumber             = 40
meshUserNumber                      = 50
decompositionUserNumber             = 60
geometricFieldUserNumber            = 70
EquationsSetFieldUserNumber         = 80
dependentFieldUserNumber            = 90
MaterialsFieldUserNumber            = 100
EquationsSetUserNumber              = 110
problemUserNumber                = 120
SourceFieldUserNumber               = 130
cellMLUserNumber                    = 140
cellMLModelsFieldUserNumber         = 150
cellMLStateFieldUserNumber          = 160
cellMLIntermediateFieldUserNumber   = 170
cellMLParametersFieldUserNumber     = 180


#------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#------------------------------------------------------------------------------------------------

# iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()


#------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#------------------------------------------------------------------------------------------------

# One Dimensional Coordinate System
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 1
coordinateSystem.CreateFinish()


#------------------------------------------------------------------------------------------------
# REGION
#------------------------------------------------------------------------------------------------

# Start Region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "cellml_split_reaction_diffusion_equation"
region.coordinateSystem = coordinateSystem
region.CreateFinish()


#------------------------------------------------------------------------------------------------
# BASIS
#------------------------------------------------------------------------------------------------

# Simplex Basis Reaction Diffusion
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 1
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
basis.CreateFinish()


#------------------------------------------------------------------------------------------------
# MESH
#------------------------------------------------------------------------------------------------

# Initialise Mesh
mesh = iron.Mesh()
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [length]
generatedMesh.numberOfElements = [numberGlobalXElements]    
generatedMesh.CreateFinish(meshUserNumber,mesh)


#------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#------------------------------------------------------------------------------------------------

# Parallelization
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()


#------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#------------------------------------------------------------------------------------------------

# Geometric Field
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
# geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Coordinate")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
# geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
# geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

generatedMesh.GeometricParametersCalculate(geometricField)


#-------------------------------------------------------------------------------------------
# EQUATION SETS 
#-------------------------------------------------------------------------------------------

# Create standard Laplace equations set
EquationsSetField         = iron.Field()
EquationsSet              = iron.EquationsSet()
EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
EquationsSet.CreateStart(EquationsSetUserNumber,region,
                                   geometricField,EquationsSetSpecification,
                                   EquationsSetFieldUserNumber,EquationsSetField)
EquationsSet.CreateFinish()


#-------------------------------------------------------------------------------------------
# DEPENDENT FIELD 
#-------------------------------------------------------------------------------------------

# Start and label  Field Dependent Field
dependentField = iron.Field()
EquationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"U")
dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"DELUDELN")
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)

# Finish Field Dependent Field
EquationsSet.DependentCreateFinish()

# Initialise Field Dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,Ca_init)


#-------------------------------------------------------------------------------------------
# MATERIAL FIELD - 
#-------------------------------------------------------------------------------------------

MaterialsField = iron.Field()
EquationsSet.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsField)
MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
#geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
#MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
EquationsSet.MaterialsCreateFinish()


# Diffusion Coefficient in X
MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Diffusivity)

# Storage Coefficient
MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,store_coeff)


#-------------------------------------------------------------------------------------------
# SOURCE FIELD
#-------------------------------------------------------------------------------------------

# Setting up the source field for reaction diffusion equation set.
# For split problem subtype, source field is not used
SourceField = iron.Field()
EquationsSet.SourceCreateStart(SourceFieldUserNumber, SourceField)
SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Source")
geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
SourceField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
EquationsSet.SourceCreateFinish()

# Initialise to 0.0
SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


# Update Source Field
SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

# Update Dependent Field
dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)

dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,1,2,1,0.0)
                



#-------------------------------------------------------------------------------------------
# CELLML FIELD
#-------------------------------------------------------------------------------------------

#Initialise cellML
cellML = iron.CellML()
cellML.CreateStart(cellMLUserNumber,region)

#Importing the cellML model
constantModelIndex = cellML.ModelImport(CellMLFile)

#Parameters that will not change --> parameters field
cellML.VariableSetAsKnown(constantModelIndex,"dude/param")
cellML.VariableSetAsWanted(constantModelIndex,"dude/intmd")

#Finish cellML
cellML.CreateFinish()


#-------------------------------------------------------------------------------------------
# CELLML OPENCMISS FIELD MAPS
#-------------------------------------------------------------------------------------------

#Initialise Field Maps
cellML.FieldMapsCreateStart()


#----------------------------#
# CELLML Dependent Variables #
#----------------------------#

#
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "dude/ca",iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(constantModelIndex,"dude/ca",
                              iron.FieldParameterSetTypes.VALUES,dependentField,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)


#Finish Field Maps
cellML.FieldMapsCreateFinish()


#-------------------------------------------------------------------------------------------
# CELLML MODEL FIELD 
#-------------------------------------------------------------------------------------------

#Initialise the Model Field
cellMLModelsField = iron.Field()
cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber,cellMLModelsField)
cellML.ModelsFieldCreateFinish()

cellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,cellMLModelsField_Init)


#--------------------------------------------------------------------------------------------
# CELLML STATE FIELD
#--------------------------------------------------------------------------------------------

#Initialise the State Field
cellMLStateField = iron.Field()
cellML.StateFieldCreateStart(cellMLStateFieldUserNumber,cellMLStateField)
cellML.StateFieldCreateFinish()

#-------------------------------------------------------------------------------------------
# CELLML INTERMEDIATE FIELD 
#-------------------------------------------------------------------------------------------

#Initialise the Parameters Field
cellMLIntermediateField = iron.Field()
cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber,cellMLIntermediateField)
cellML.IntermediateFieldCreateFinish()


#-------------------------------------------------------------------------------------------
# CELLML PARAMETERS FIELD 
#-------------------------------------------------------------------------------------------

#Initialise the Parameters Field
cellMLParametersField = iron.Field()
cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber,cellMLParametersField)
cellML.ParametersFieldCreateFinish()


#------------------------------------------------------------------------------------------------
# EQUATIONS
#------------------------------------------------------------------------------------------------

# Equations Set
Equations = iron.Equations()
EquationsSet.EquationsCreateStart(Equations)
Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Equations.outputType = iron.EquationsOutputTypes.NONE
#Equations.outputType = iron.EquationsOutputTypes.MATRIX

#Finish
EquationsSet.EquationsCreateFinish()


#------------------------------------------------------------------------------------------------
# PROBLEM 
#------------------------------------------------------------------------------------------------
    
# Create Problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                           iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                           iron.ProblemSubtypes.CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()


# Create control loops
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()

problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
# start time, stop time, time increment
controlLoop.TimesSet(startT,endT,Tstep)


#Control Loop Outputs
#controlLoop.TimeOutputSet(outputfreq)
#controlLoop.LoadOutputSet(1)
#controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
#controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)

problem.ControlLoopCreateFinish()


#------------------------------------------------------------------------------------------------
# SOLVER 
#------------------------------------------------------------------------------------------------

#
#    1st Solver --> DAE
#         |
#         v
#    2nd Solver --> Dynamic 
#         |
#         v      
#    3rd Solver --> DAE
#

#Create problem solver for Strang splitting
problem.SolversCreateStart()


#Create first solver --> DAE Solver
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
solver.DAETimeStepSet(ODE_TIME_STEP)
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Create second solver --> Dynamic solver for parabolic equation
solver = iron.Solver()
linearsolver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solver)

#Set theta - backward vs forward time step parameter
solver.DynamicThetaSet([1.0])

#Set output type
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Obtain dynamic linear solver from the solver
solver.DynamicLinearSolverGet(linearsolver)

#Set Library
#solver.LibraryTypeSet(iron.SolverLibraries.LAPACK)
#solver.LibraryTypeSet(iron.SolverLibraries.CMISS)
#solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
linearsolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)

#solver.LinearDirectTypeSet(iron.DirectLinearSolverTypes.LU)

#linearsolver.LinearIterativeMaximumIterationsSet(1000)
#linearsolver.linearIterativeAbsoluteTolerance = 1.0E-12
#linearsolver.linearIterativeRelativeTolerance = 1.0E-12


#Create third solver --> Another DAE Solver
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,solver)
solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
solver.DAETimeStepSet(ODE_TIME_STEP)
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Finish the problem
problem.SolversCreateFinish()


#------------------------------------------------------------------------------------------------
# SOLVER CELLML EQUATIONS 
#------------------------------------------------------------------------------------------------

#Start solver CellML Equations
problem.CellMLEquationsCreateStart()

#Create first solver
#cellML equations
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
cellMLEquations = iron.CellMLEquations()
solver.CellMLEquationsGet(cellMLEquations)
#Add into CellML Environment
cellMLIndex = cellMLEquations.CellMLAdd(cellML)

#Create third solver
#cellML equations
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,solver)
cellMLEquations = iron.CellMLEquations()
solver.CellMLEquationsGet(cellMLEquations)
#Add into CellML Environment
cellMLIndex = cellMLEquations.CellMLAdd(cellML)

#Finish the solver CellML Equations
problem.CellMLEquationsCreateFinish()


#------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS - REACTION DIFFUSION
#------------------------------------------------------------------------------------------------

#Start solver equations
problem.SolverEquationsCreateStart()

#Create second solver
#Solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solver)
solver.SolverEquationsGet(solverEquations)


solverEquations.sparsityTypeSet = iron.SolverEquationsSparsityTypes.SPARSE

equationsSetIndex = solverEquations.EquationsSetAdd(EquationsSet)

problem.SolverEquationsCreateFinish()


#------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS - REACTION DIFFUSION
#------------------------------------------------------------------------------------------------

# Set up Boundary Conditions & Perturbation Initial Conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

condition = iron.BoundaryConditionsTypes.FIXED
value = 1.5
BCNodes = [1,11]

for node in range(1,12):
    nodeDomain = decomposition.NodeDomainGet(node,1)
    if nodeDomain==computationalNodeNumber:    
        if node in BCNodes:
            #boundary Conditions
            boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,
                                        1,1,node,1,condition,value) 
            cellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES,
                                                                   1,1,node,1,0)

solverEquations.BoundaryConditionsCreateFinish()



#--------------------------------------------------------------------------------------------
# OUTPUT 
#--------------------------------------------------------------------------------------------

problem.Solve()


#--------------------------------------------------------------------------------------------
# OUTPUT 
#--------------------------------------------------------------------------------------------

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.ElementsExport("cellml_split_reaction_diffusion_equation","FORTRAN")
fields.NodesExport("cellml_split_reaction_diffusion_equation","FORTRAN")
fields.Finalise()

# Finalise OpenCMISS-Iron
iron.Finalise()