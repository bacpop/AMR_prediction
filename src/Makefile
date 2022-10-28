CXX = g++
CXXFLAGS = -Wall -g -O0 
CPPFLAGS = -I ${HOME}/include -I . -I ..
LDFLAGS = -L ${HOME}/lib  
LDLIBS = -lz #-lsdsl -ldivsufsort -ldivsufsort64

amr_prediction: amr_pred.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@ 

clean:
	$(RM) *.o ~* amr_prediction

amr_pred_for_CLI: amr_pred_for_CLI.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@ 

# web specific options
web: CXX = em++
web: CPPFLAGS += -DWEB
web: CXXFLAGS = -O3 -s ASSERTIONS=1  \
				-DJSON_NOEXCEPTION \
				-s DISABLE_EXCEPTION_CATCHING=1 \
				-flto \
				--bind -s STRICT=1 \
				-s ALLOW_MEMORY_GROWTH=1 \
				-s USE_ZLIB=1 \
				-s "EXPORTED_FUNCTIONS=['_malloc']" \
				-s 'EXPORTED_RUNTIME_METHODS=["FS"]' \
				-s EXPORT_NAME=AMRprediction \
				-s MODULARIZE=1 \
				-Wall -std=c++14 \
				--preload-file ./files

web: LDFLAGS = -lnodefs.js -lworkerfs.js

WEB_OUT=../build/web_amr_prediction
WEB_OBJS=${WEB_OUT}.js ${WEB_OUT}.html ${WEB_OUT}.wasm

web: amr_pred.cpp
	$(LINK.cpp) $^ $(LDLIBS) -o ${WEB_OUT}.js
	sed -i.old '1s;^;\/* eslint-disable *\/;' ${WEB_OUT}.js

check:
	./amr_prediction