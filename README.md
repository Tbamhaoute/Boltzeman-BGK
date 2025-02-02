- Ce code présente l’application de la méthode des différences finies en une dimension pour résoudre une équation cinétique simplifiée de type BGK décrivant un gaz fictif avec et sans terme de collision, et une comparaison avec les équations d'euler en fonction du nombre de Knudsen.

- Pour la compilation du code éxecuter:
   -  make
   -  gfortran -fopenmp -o bgk bgk.f90
   -  ./bgk
