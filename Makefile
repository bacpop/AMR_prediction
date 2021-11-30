CXX = g++
CXXFLAGS = -Wall -g -O0 
CPPFLAGS = -I ${HOME}/include -I ${HOME}/cpp_projects/amr_pred
LDFLAGS = -L ${HOME}/lib  
LDLIBS = -lz #-lsdsl -ldivsufsort -ldivsufsort64

amr_prediction: amr_pred.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@ 

clean:
	$(RM) *.o ~* amr_prediction


# web specific options
web: CXX = em++

web: CXXFLAGS = -O3 -s ASSERTIONS=1  \
				--bind -s STRICT=1 \
				-s ALLOW_MEMORY_GROWTH=1 \
				-s USE_ZLIB=1 \
				-s "EXPORTED_FUNCTIONS=['_malloc']" \
				-s 'EXTRA_EXPORTED_RUNTIME_METHODS=["FS"]' \
				-s EXPORT_NAME=HelloWorld \
				-s MODULARIZE=1 \
				-Wall -std=c++14 \
				--preload-file files

web: LDFLAGS = -lnodefs.js -lworkerfs.js

WEB_OUT=web/web_amr_prediction
WEB_OBJS=${WEB_OUT}.js ${WEB_OUT}.html ${WEB_OUT}.wasm

web: amr_pred.o
	$(LINK.cpp) $^ $(LDLIBS) -o ${WEB_OUT}.js
	sed -i.old '1s;^;\/* eslint-disable *\/;' ${WEB_OUT}.js