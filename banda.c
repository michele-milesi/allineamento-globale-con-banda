/*
	Michele Milesi 844682
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "blosum62.h"

int main(int argc, char **argv){
	//dichiarazione variabili
	int64_t** M;	//matrice di programmazione dinamica
	short** P;	//matrice dei predecessori, per ricostruire l'allineamento
				//0 -> up-left
				//1 -> left
				// 2 -> up
				
	char* s1, *s2;	//stringhe in input	
	char* alignment1;	//allineamneto di s1 con s2
	char* alignment2;	//allineamneto di s2 con s1
	int l1, l2;	//lunghezze delle stringhe in input
	int base, extra = 1, band_width;
	int k;	
	int i, j, index;	
	int64_t best_alignment_score = 0, upper_bound, score;
	short first_time = 1;	//variabile che serve per capire se bisogna rilasciare memoria
							//occupata dalla matrice di programmazione dinamica
							//prima di allocare nuovo spazio per una nuova iterazione
	
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
	base = abs(l1 - l2) + 1;
	band_width = base + 2 * extra;
	
	//scambio s1 e s2 nel caso in cui |s1| > |s1|
	if(l1 > l2) {
		char* tmp = s1;
		s1 = s2;
		s2 = tmp;
		int t = l1;
		l1 = l2;
		l2 = t;
	}
	
	//calcolo valore allineamento nel caso in cui s1 e s2 coincidano
	for(i=0; i<l1; i++)
		best_alignment_score += d[index_of(s1[i])][index_of(s1[i])];
		
		
	M = (int64_t**) malloc(sizeof(int64_t*) * (l1 + 1));
	P = (short**) malloc(sizeof(short*) * (l1 + 1));
	
	do{
		//allocazione spazio necessario
		for(i = 0; i <= l1; i++) {
			if (first_time != 1) {
				free(M[i]);
				free(P[i]);
			}
			M[i] = (int64_t*) malloc(sizeof(int64_t) * band_width);
			P[i] = (short*) malloc(sizeof(short) * band_width);
		}
		first_time = 0;	//imposto che d'ora in poi, se è necessario raddopiare la banda
						//prima di allocare spazio, bisogna rilasciare quello già occupato
		
		//calcolo allineamento
		k = extra;	//la prima cella da considerare della matrice di programmazione dinamica è in posizione [0, extra]
		M[0][k] = 0;	//imposto caso base: la cella M[0, extra] della matrice coincide con
						//la cella in posizione [0, 0] della matrice di programmazione dinamica nell'algoritmo
						//di Smith-Waterman
						
		P[0][k] = -1;	//mi segno dove mi devo fermare nel ricostruire l'allineamento
		
		
		for(i = 0; i <= l1; i++) {
			for(j = 0; j < band_width; j++) {
				int s1_index = i - 1;	//carattere di s1 che si sta analizzando
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
					//up -> j non deve sforare i bordi della matrice
					if(j < band_width - 1 && M[i][j] < M[i-1][j+1] + d[index_of(s1[s1_index])][index_of('*')]) {
						M[i][j] = M[i-1][j+1] + d[index_of(s1[s1_index])][index_of('*')];
						P[i][j] = 2;
					}
				}
			}
		}
		score = M[l1][l2 - l1 + k];	//valore calcolato dell'allineamento
		upper_bound = best_alignment_score - (4 * (k + 1));	
		extra *= 2;	
		band_width = base + 2 * extra;
		
	//si continua a raddoppiare la banda fino a quando o la grandezza della matrice è maggiore a 
	//quella dell'algoritmo dell'allineamento globale "classico" (Smith-Waterman)
	//oppure quando ci si accorge che raddoppiare la banda non può più portare ad un miglioramento del valore trovato
	}while(score <= upper_bound && k <= l2);
	
	
	//ricostruzione di un allineamento ottimo
	//l'allineamento non può essere più lungo di |s1| + |s2|
	//+ 1 per contenere il carattere di fine stringa '/0'
	alignment1 = (char*) malloc(sizeof(char) * (l1 + l2 + 1));
	alignment2 = (char*) malloc(sizeof(char) * (l1 + l2 + 1));
	i = l1;	
	j = l2 - l1 + k;	//j parte dall'elemento della matrice P che corrisponde all'ultimo carattere di s2
	index = 0;	//tiene conto dei caratteri inseriti nell'allineamento
	while(P[i][j] == 0 || P[i][j] == 1 || P[i][j] == 2)
	{
		//up
		if(P[i][j] == 2){
			alignment1[index] = s1[i - 1];
			alignment2[index] = '-';
			i--, j++;
		}
		//up-left
		else if(P[i][j] == 0){
			alignment1[index] = s1[i - 1];
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
	
	//scrittura risultati su file
	output_file = fopen("result.txt", "w");
	fprintf(output_file, "Alignment Score: %ld\n\n", score);
	fprintf(output_file, "Alignment:\n");
	fprintf(output_file, "%s\n", strrev(alignment1));
	fprintf(output_file, "%s", strrev(alignment2));
	
	//chiusura file input e output
	fclose(output_file);
	fclose(input_file);
	
	//libera spazio nello heap
	for(i = 0; i <= l1; i++) {
		free(M[i]);
		free(P[i]);
	}
	free(M);
	free(P);
	free(s1);
	free(s2);
	free(alignment1);
	free(alignment2);
	
	return 0;
}
