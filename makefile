CC = g++
CFLAGS = -Wall -g -lz -lsdsl -ldivsufsort -ldivsufsort64
INC = -I ~/include
LIB = -L ~/lib

amr_prediction: amr_pred.o
	$(CC) $(INC) $(LIB) $(CFLAGS) -o amr_prediction amr_pred.o 


# works in the command line:
# g++ -I ~/include -L ~/lib -g amr_pred.cpp -lz -o amr_pred -lsdsl -ldivsufsort -ldivsufsort64
