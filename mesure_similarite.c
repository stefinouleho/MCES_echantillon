#include "fonctions_mces.h"

double last_chrono;
int total_molecules = 500;
float date = 20.00;
int taille_limite = 100000;

int * lecture_liste_molecules()
{
	int *liste = malloc ( total_molecules * sizeof(int));
	FILE *F = fopen("fichiers/liste_molecules_mces","r");
	
	if( F == NULL)
	{
		fprintf(stdout, " Cannot open the file fichiers/liste_molecules_mces");
		exit(10);
	}
	int i; 
	for ( i = 0; i < total_molecules; i++)
	{
		fscanf(F, "%d",&liste[i]);
	}
	
	fclose(F);
	
	return liste;
			
}

void sauvegarde_compteur(int i , int j)
{
//sauvegarde des indices au cas ou le ptogramme serait interompu 
			
	FILE *G = fopen("tmp","w");
	if( G == NULL)
	{
		fprintf(stdout, " Cannot open the file tmp");
		exit(11);
	}
			
	if( j == i - 1)
	{
		fprintf(G,"%d	%d", i + 1 , 0);
	}
			
	else
	{
		fprintf(G,"%d	%d", i , j + 1);
	}
			
	fclose(G);
	
	
	FILE *F = fopen("prochain","r");
	if( F == NULL)
	{
		fprintf(stdout, " Cannot open the file prochain");
		exit(12);
	}
	int numero, indice;
	
	fscanf(F,"%d", &numero);
	fscanf(F,"%d", &indice);
	
	fclose(F);
	
	FILE *H = fopen("dernier","w");
	if( H == NULL)
	{
		fprintf(stdout, " Cannot write in the file dernier");
		exit(13);
	}
	
	fprintf(H,"%d	%d", numero, indice );
	
	fclose(H);
	
	G = fopen("tmp","r");
	if( G == NULL)
	{
		fprintf(stdout, " Cannot open the file tmp");
		exit(14);
	}
	fscanf(G,"%d", &numero);
	fscanf(G,"%d", &indice);
	
	fclose(G);
	
	F = fopen("prochain","w");
	if( F == NULL)
	{
		fprintf(stdout, " Cannot open the file prochain");
		exit(15);
	}
	fprintf(F,"%d	%d", numero, indice );
	
	fclose(F);
			
	// fin de la sauvegarde 
}

int main(int argc, char *argv[])
{
	
	
	if( argc != 6){
		fprintf(stdout,"Missing arguments (num chebi1 , num chebi 2 ( 0  for all) ,type, date , lenght limit )\n");
		exit(20);
	}

	//lecture des molecules dans le fichier chebi_lite.sdf
	struct molecule *M = lecture_fichier_chebi();
	
	
	FILE *F = fopen("prochain","r");
	if( F == NULL)
	{
		fprintf(stdout, " Cannot open the file prochain");
		exit(10);
	}
	int numero, indice;
	
	fscanf(F,"%d", &numero);
	fscanf(F,"%d", &indice);
	fclose(F);
	
	int *liste_molecules = lecture_liste_molecules();
	
	int i,j;
	
	float r,start,stop;
	
	int pos1,pos2;
	for ( i = numero ; i < total_molecules; i++)
	{
		
		fprintf(stdout,"\r%5d / %d",i,total_molecules);
		fflush(stdout); 
		for ( j = indice ; j < i ; j++)
		{
			
			
			pos1 = position_M(liste_molecules[i],M);
			pos2 = position_M(liste_molecules[j],M);
			M[pos1]= elimination_feuilles(M[pos1]);
			M[pos2]= elimination_feuilles(M[pos2]);
			start = chrono();
			r = mesure_similarite( liste_molecules[i], liste_molecules[j],M, date, taille_limite);
			stop = chrono();
			
			char fichier[256];
			//sprintf(fichier,"Dossier/%d_%d.result",i,j);
			sprintf(fichier,"Dossier/similarite.result");
			
			F = fopen(fichier,"a");
			if( F == NULL)
			{
				fprintf(stdout, "Cannot open the file %s\n",fichier);
				exit(19);
			}
			
			sprintf(fichier,"Dossier/temps.result");
			FILE *G = fopen(fichier,"a");
			if( G == NULL)
			{
				fprintf(stdout, " Cannot open the file %s\n",fichier);
				exit(19);
			}
			fprintf(F,"%.2f	",r);
			fprintf(G,"%.2f	",stop -start);
			
			if( j == i - 1)
			{
				fprintf(F,"\n");
				fprintf(G,"\n");
			}
				
			fclose(F);
			fclose(G);
			
			sauvegarde_compteur(i,j);
		}
		indice = 0;
		
	} 
	
	
	
	//liberation de la memoire
	free(liste_molecules);
	int nb_mol;
	printf("\n3. Libération de la mémoire : %.3lf s\n",chrono());

	for(nb_mol=0 ; nb_mol < NB_MOLECULES ; nb_mol++) 
		liberer_molecule(M[nb_mol]);
	free(M);

	exit(0);
}
