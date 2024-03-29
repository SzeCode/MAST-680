# MAST-680
A repository to show my work done for the course MAST 680, Data Driven Methods in Dynamical Systems, from Concordia University in Montreal taken in Winter 2023. The course was taught by Jason Bramburger.

The repository contains the files for:

Assignment 1: Background Subtraction in Video.
- Description: An assignment that incorporates dynamic mode decomposition (DMD) with low-rank and sparse reconstruction to find the moving foreground in a video.

Assignment 2: Learning Dynamics from Video.
- Description: An assignment that uses sparse identification of non-linear dynamics (SINDY) to discover the governing equations of a mass-spring system in a video.

Assignment 3: Forecasting Chaos with Neural Networks.
- Description: An assignment that uses neural networks to forcast the dynamics of the Lorenz 	equations.

Final Project:  Dynamic Mode Decomposition (DMD) in Salient Region Detection.
- Description: A project which explores the use of dynamic mode decomposition (DMD) to help identify the salient regions in a static image. This is an attempt to recreate the results of the work of Sikha et al (https://doi.org/10.1016/j.jocs.2017.07.007)


Instructions on how to run the files.

Assignment 1: 
- Run the .m file directly in Octave or Matlab. The outputs are displayed when you run the program.

Assignment 2: 
- Run all three of the centroid tracking files Assignment2_CamX.m, then run the main file PCA_SINDy_VersFinal.m. The outputs are displayed when you run the program.

Assignment 3: 
- Run the .sln (Microsoft visual studios 2022) file in Main_forecast for predictions of the Lorenz system. The outputs are displayed when you run the program. Alternatively, run the .m file using Octave or Matlab in the outputs folder.
- Run the .sln (Microsoft visual studios 2022) file in Lobe_Transition for predictions of the time to transition from one lobe to another in the Lorenz system. Alternatively, run the .m file using Octave or Matlab in the outputs folder.
- Run the .sln (Microsoft visual studios 2022) file in LoadModel to load saved models generated from the Neural Networks.

Final Project: 
- Run the MAIN_DMD_SALIENT_EXTRACT.m file directly in Octave or Matlab. Specify manualy the image name to change image.


Feedbacks from the professor is also included in each folder. They will be implemented in the near future.
