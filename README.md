## derivative-estimation
implementation of derivative estimation by quotient difference - c++

This project implements the derivative estimation method in random design described in the following paper :    	
"Derivative Estimation in Random Design", from Yu Liu, Kris De Brabanter

## INSTALLATION
You will need the Eigen library, available here : http://eigen.tuxfamily.org/dox/eigen-doc.tgz 		
One only needs to place the eigen-3.x.x folder (containing the Eigen folder) inside the include folder. 	
The project is built using CMake (more information https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2019)_
I personally used Visual Studio software on a Windows 10 distribution.

## USE
The derivative.cpp file containing the main function will allow you to create random data, apply the method for derivative estimation
and write the results in csv files.
You could also easily use your own data by following the instructions and explanations contained in derivative.cpp, annex.h and annex.cpp

## References
I entirely used the method from "Derivative Estimation in Random Design", from Yu Liu, Kris De Brabanter.
The library I used were Eigen 3.3.7 for Linear algebra calculations, and KDE by Tim Nugent (c) 2014 for kernel density estimation.

## Contact
Feel free to contact me if you have any questions at
andrew.caunes@gmail.com or andrew.caunes@etu.emse.fr
