CFLAGS=-g -Wall
#CFLAGS=-O2 -Wall

CC = gcc -Wall -Wextra
CXX = g++ -Wall -Wextra

all: mesure_similarite

run: mesure_similarite
	./mesure_similarite 28973 4672 2 10 1600

init: 
	rm -rf Dossier
	rm -f dernier
	rm -f prochain
	rm -f tmp
	rm -f fichiers/liste_molecules_mces
	rm -f fichiers/liste_molecules_mces.save
	echo "0 0" > prochain
	head -500 fichiers/molecules.data > fichiers/liste_molecules_mces
	cp fichiers/liste_molecules_mces fichiers/liste_molecules_mces.save
	mkdir Dossier
	touch Dossier/similarite.result
	touch Dossier/temps.result
	
val: mesure_similarite
	valgrind --leak-check=full --show-leak-kinds=all ./mesure_similarite 28973 4672  2 20 600	

mesure_similarite: mesure_similarite.o lecture_molecule_sdf.o fonctions_mces.o feuilles.o helpers/graph.o helpers/cliquerecursif.o
	gcc ${CFLAGS} -o $@ $^

mesure_similarite_scip : mesure_similarite.o lecture_molecule_sdf.o fonctions_mces.o feuilles.o helpers/graph.o helpers/proglin_helper_scip.o helpers/sciplib.a helpers/cliquescip.o 
	$(CXX) -I helpers/scip -o $@ $^ -lpopt -lgmp -lm -lz -lreadline -lncurses 	

mesure_similarite.o: mesure_similarite.c
	gcc ${CFLAGS} -c mesure_similarite.c

fonctions_mces.o: fonctions_mces.c fonctions_mces.h 
	gcc ${CFLAGS} -c fonctions_mces.c
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h
	gcc ${CFLAGS} -c lecture_molecule_sdf.c

helpers/proglin_helper_scip.o : helpers/proglin_helper_scip.c helpers/proglin_helper.h helpers/sciplib.a
	$(CC) -I helpers/scip -o $@ -c $<

feuilles.o: feuilles.c fonctions_mces.h
	gcc ${CFLAGS} -c feuilles.c

	
clean: 
	rm -f mesure_similarite
	rm -f *.o
	rm -f helpers/*.o

