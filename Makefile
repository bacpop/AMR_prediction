CXX = g++
CXXFLAGS = -Wall -g -O3 -Wextra
CPPFLAGS = -I ${HOME}/include
LDFLAGS = -L ${HOME}/lib  
LDLIBS = -lz -lsdsl -ldivsufsort -ldivsufsort64

amr_prediction: amr_pred.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@ 

clean:
	$(RM) *.o ~* amr_prediction


# works in the command line:
# g++ -I ~/include -L ~/lib -g amr_pred.cpp -lz -o amr_pred -lsdsl -ldivsufsort -ldivsufsort64
