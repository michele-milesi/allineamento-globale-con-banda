#include <stdlib.h>
#include <stdio.h>
#include "blosum62.h"

int main(){
	//dichiarazione variabili
	int** M;	//matrice di programmazione dinamica
	short** P;	//matrice dei predecessori, per ricostruire l'allineamento
				//0 -> up-left, 1 -> left, 2 -> up
	char* s1;	//stringa in input	
	char* s2;	//stringa in input
	char* alignment1;	//allineamneto di s1 con s2
	char* alignment2;	//allineamneto di s2 con s1
	int l1, l2;	//lunghezze delle stringhe in input
	int base, extra = 1, k, band_width;
	int i, j, index;	
	int best_alignment_score = 0, upper_bound, score;
	int first_time = 1;
	
	FILE* input_file = fopen("input.txt", "r");
	FILE* output_file;
	
	//input
	fscanf(input_file, "%d", &l1);
	s1 = (char*) malloc(sizeof(char) * (l1 + 1));
	fscanf(input_file, "%s", s1);
	
	fscanf(input_file, "%d", &l2);
	s2 = (char*) malloc(sizeof(char) * (l2 + 1));
	fscanf(input_file, "%s", s2);
	
	//calcolo parametri per la larghezza della banda
	base = abs(l1 - l2);
	band_width = base + 2 * extra;
	
	//calcolo valore allineamento nel caso in cui s1 e s2 coincidano
	for(i=0; i<l1; i++)
		best_alignment_score += d[index_of(s1[i])][index_of(s1[i])];
		
	
	M = (int**) malloc(sizeof(int*) * (l1 + 1));
	P = (short**) malloc(sizeof(short*) * (l1 + 1));
	
	do{
		for(i = 0; i <= l1; i++) {
			if (first_time != 1) {
				free((void*) M[i]);
				free((void*) P[i]);
			}
			M[i] = (int*) malloc(sizeof(int) * band_width);
			P[i] = (short*) malloc(sizeof(short) * band_width);
		}
		first_time = 0;
		
		//calcolo allineamento
		k = band_width / 2;	
		M[0][k] = 0;
		P[0][k] = -1;
		for(i = 0; i <= l1; i++) {
			for(j = 0; j < band_width; j++) {
				int s1_index = i - 1;	//carattere di s2 che si sta analizzando
				int s2_index = j + i - k - 1;	//carattere di s2 che si sta analizzando
					
				//condizione al contorno (left)
				if(s1_index == -1 && s2_index >= 0 && s2_index < l2) {
					M[i][j] = M[i][j-1] + d[index_of('*')][index_of(s2[s2_index])];
					P[i][j] = 1;
				}	
				//condizione al contorno (up)
				if(s2_index == -1 && s1_index >= 0) {
					M[i][j] = M[i-1][j+1] + d[index_of(s1[s1_index])][index_of('*')];
					P[i][j] = 2;
				}
				
				//passo ricorsivo
				if(s1_index >= 0 && s2_index >= 0 && s2_index < l2) {	
					//up-left
					M[i][j] = M[i-1][j] + d[index_of(s1[s1_index])][index_of(s2[s2_index])];
					P[i][j] = 0;
					//left
					if(j > 0 && M[i][j] < M[i][j-1] + d[index_of('*')][index_of(s2[s2_index])]) {
						M[i][j] = M[i][j-1] + d[index_of('*')][index_of(s2[s2_index])];	
						P[i][j] = 1;
					}
					//up
					if(j < band_width - 1 && M[i][j] < M[i-1][j+1] + d[index_of(s1[s1_index])][index_of('*')]) {
						M[i][j] = M[i-1][j+1] + d[index_of(s1[s1_index])][index_of('*')];
						P[i][j] = 2;
					}
				}
			}
		}
		score = M[l1][l2 - l1 + k];	//valore calcolato dell'allineamento
		upper_bound = best_alignment_score -(4 * (band_width + 1));	//upper_bound nel caso migliore (s1 = s2) e k + 1 indel
		extra *= 2;	
		band_width = base + 2 * extra;
		
	//si continua a raddoppiare la banda fino a quando o la grandezza della matrice è maggiore o uguale a 
	//quella dell'allineamento globale "classico"
	//oppure quando ci si accorge che raddoppiare la banda non può più portare ad un miglioramento del valore trovato
	}while(score <= upper_bound && k <= (l1>l2?l1:l2));
	
	//ricostruzione di un allineamento ottimo
	alignment1 = (char*) malloc(sizeof(char) * (l1 + l2 + 1));	//l'allineamento non può essere più lungo di |s1| + |s2|
	alignment2 = (char*) malloc(sizeof(char) * (l1 + l2 + 1));
	i = l1;	
	j = l2 - l1 + k;	//j parte dall'elemento della matrice P che corrisponde all'ultimo carattere di s2
	index = 0;	//tiene conto dei caratteri inseriti nell'allineamento
	while(P[i][j] == 0 || P[i][j] == 1 || P[i][j] == 2)
	{
		//up
		if(P[i][j] == 2){
			alignment1[index] = s1[i-1];
			alignment2[index] = '-';
			i--, j++;
		}
		//up-left
		else if(P[i][j] == 0){
			alignment1[index] = s1[i-1];
			alignment2[index] = s2[j + i - k - 1];			
			i--;
		}
		//left
		else if(P[i][j] == 1){
			alignment1[index] = '-';
			alignment2[index] = s2[j + i - k - 1];
			j--;
		}
		index++;
	}
	alignment1[index] = '\0';
	alignment2[index] = '\0';
	alignment1 = (char*) strrev(alignment1);
	alignment2 = (char*) strrev(alignment2);
	
	//scrittura risultati su file
	output_file = fopen("result.txt", "w");
	fprintf(output_file, "Alignment Score: %d\n\n", score);
	fprintf(output_file, "Alignment:\n");
	fprintf(output_file, "%s\n", alignment1);
	fprintf(output_file, "%s", alignment2);
	
	//chiusura file input e output
	fclose(output_file);
	fclose(input_file);
	
	//libera spazio nello heap
	for(i = 0; i <= l1; i++) {
		free((void*) M[i]);
		free((void*) P[i]);
	}
	free((void*) M);
	free((void*) P);
	free((void*) s1);
	free((void*) s2);
	free((void*) alignment1);
	free((void*) alignment2);
	
	return 0;
}
