CFlags = -ansi -g -Wall -Wextra -Werror -pedantic-errors -lm 
spkmeans: spkmeans.o wam.o ddg.o gl.o jacobi.o  
	gcc -o spkmeans spkmeans.o wam.o ddg.o gl.o jacobi.o $(CFlags)

spkmeans.o: spkmeans.c jacobi.o gl.o ddg.o wam.o spkmeans.h
	gcc -c spkmeans.c $(CFlags)

jacobi.o: jacobi.c spkmeans.h
	gcc -c jacobi.c $(CFlags)

gl.o: gl.c spkmeans.h
	gcc -c gl.c $(CFlags)

ddg.o: ddg.c spkmeans.h
	gcc -c ddg.c $(CFlags)      
        
wam.o: wam.c spkmeans.h
	gcc -c wam.c $(CFlags)


