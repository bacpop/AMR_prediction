CC = g++
CFLAGS = -Wall -g -lz -lsdsl -ldivsufsort -ldivsufsort64
INC = -I ~/include
LIB = -L ~/lib

amr_prediction: amr_pred.o
	$(CC) $(INC) $(LIB) $(CFLAGS) -o amr_prediction amr_pred.o 



