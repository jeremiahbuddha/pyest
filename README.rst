.. role:: m(raw)
    :format: latex

===============================================================================
 ASEN-5070 Term Project
===============================================================================
:Name: Jonathon Smith
:Date: April 2012

Overview
========
The OD project will consist of deriving the theory and algorithms and developing 
a computer program to establish the trajectory of an Earth-orbiting satellite 
given range and range rate observations from three ground stations.

Requirements
============

1) Derive the equations of motion for a satellite moving under the influence of 
   drag and :m:`$J_2$`
2) Derive the :m:`$A$` and :m:`$\tilde{H}$` matrices for the linearized system
3) Set up the batch algorithm to solve for where is the deviation of X from a 
   reference or nominal solution. Account for a priori values of and P(to).
4) Set up the sequential algorithm to solve for where is the deviation of X 
   from a reference or nominal solution. Account for a priori values of and P(to).

A minimum requirement to obtain a B in the course is to develop a computer 
program and successfully execute the batch processor algorithm. Additional 
credit will be given for implementing the sequential algorithm. You should 
write a final report containing the mathematical development, flowcharts and 
data processing results of your program.

State Vector
============

The state vector being estimated, along with the stm, is given below.
0. X
1. Y 
2. Z
3. dX
4. dY
5. dZ
6. mu
7. J_2
8. C_d
9. X_1
10. Y_1
11. Z_1
12. X_2
13. Y_2
14. Z_2
15. X_3
16. Y_3
17. Z_3
18. phi00
19. phi01
20. phi02
21. phi03
22. phi04
23. phi05
24. phi06
25. phi07
26. phi08
27. phi09
28. phi010
29. phi011
30. phi012
31. phi013
32. phi014
33. phi015
34. phi016
35. phi017
36. phi10
37. phi11
38. phi12
39. phi13
40. phi14



