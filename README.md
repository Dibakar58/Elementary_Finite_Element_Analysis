# Elementary_Finite_Element_Analysis

Hello there! If you are learning FEM from ground zero this is a great place for you to start. This is a repository of MATLAB codes that contains some popular FEM problems one usually encounters when learning FEA for the first time. This repository includes beginner-friendly and readable codes (with adequate comments) for some elementary FEM problems. As one would normally encounter in any FEM course, I have included FEM for the 1D Truss,2D Truss,3D Truss , 2D Plates (plane stress condition), Plate with a Hole and some 3D Geometries and finally the results are compared with ANSYS Structural (APDL). Since this is absolutely elementary, I have considered all problems as Linear . Though, it might be a great idea to include a fast iterative solvers and new FEM schemes . It shall soon be added, stay tuned! For details on which code to refer for a specific problem, read on further below:

# Finite Difference Method- 1D Truss
![](Images/Truss_1D.png)

# Finite Difference Method- 2D Truss
## Example 1


![](Images/2Dtruss1.png)
<br>
![](Images/2Dtruss_compare.png)

## Example 2


![](Images/res.png)
<br>
![](Images/2Dtruss2.png)

# Finite Element Method- 2D Plate
## Example 1
### 2D Mesh
![](Images/p1m.png)
### Results
![](Images/p1u.png)
![](Images/p1v.png)
## Example 2
### 2D Mesh
![](Images/mes.png)

### Results
![](Images/um.png)
![](Images/vm.png)

# # Finite Element Method- 3D cantilever 
Using FEM on a 3-D structured mesh . The results look something like 
![](Images/Disp_FEM.png)
![](Images/Disp_ANSYS.png)

![](Images/U_FEM.png)
![](Images/U_ANSYS.png)

![](Images/V_FEM.png)
![](Images/V_ANSYS.png)
### Von Mises Stress
![](Images/S_FEM.png)
![](Images/S_ANSYS.png)
