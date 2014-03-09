CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX -I eigen -I CImg
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-L /opt/X11 -L /opt/X11/lib\
    	-lGL -lGLU -lm -lstdc++ -lX11

else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I glut-3.7.6-bin eigen
	LDFLAGS = -lglut -lGLU
endif
	


RM = /bin/rm -f 
all: main 
main: example_02.o 
	$(CC) $(CFLAGS) -o rt example_02.o $(LDFLAGS) 
example_02.o: example_02.cpp
	$(CC) $(CFLAGS) -c example_02.cpp -o example_02.o
clean: 
	$(RM) *.o rt
 


