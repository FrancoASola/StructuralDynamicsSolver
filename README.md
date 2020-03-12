# StructuralDynamicsSolver
Structural Dynamics Solver for a structure with 3 degrees of freedoms and coupled damping. <br>

# About
This was a term project for Structural Dynamics, which was an expansion of a mid-term project (Sindgle DOF solver). <br>
The solver uses Newton-Raphson iterations to solve a differential equation at each time step. This project could easily be expanded to include more degrees of freedom and could even be coded to input an arbitrary amount of DOFs (currently 3 DOFS are hard coded).<br>
<br>
This does not use Modal Decomposition, since damping is coupled. I have coded a version that uses Modal Decomposition to solve the differential equation, but I can't seem to find the file. Naturally, Modal Decomposition speeds up the process since you can ignore a portion of modes that do not contribute much to the behavior. However, this would only work for decoupled systems.


