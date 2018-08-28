## GIS
Grafy i sieci 17Z.

### Temat
Problem komiwojazera w obecnosci ulic jednokierunkowych - implementacja algorytmu rozwiązującego problem komiwojażera oraz zbadanie, czy można modelować ulice jednokierunkowe (identyfikacja warunków, w których algorytm wybiera drogę pod prąd mimo istnienia 'legalnej' drogi).

### Źródła
Oryginalna implementacja [tabu search](https://github.com/CaoManhDat/TSP-TabuSearch)

### Dane
Struktura pliku danych:
```
5           // liczba wierzchołków / miast
0 6 10 1 9  // a_{i, j}  - koszt przejścia z miasta i do j (macierz kosztów C)
1 0 7 10 4
10 7 0 5 1
2 10 1 0 4
9 1 3 4 0
```
